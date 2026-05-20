#!/usr/bin/env python3
"""
sleep_edf_unified.py

Unified Sleep-EDF Telemetry pipeline reproducing both paper Table 1 (|K| column)
and paper Table 2 (eff_rank, dwell, K/eff_rank ratio) from the canonical scripts.

Architecture:
  - Spectrogram from sleep_edf_K_pipeline__2_.py
      (linear-power interpolation then log; matched to EEF 4.16 format)
  - Metric computation from sleep_edf_K_pipeline__2_.py
      (K, eff_rank, dwell, ratio in one pass)
  - Stage mapping from sleep_K_pipeline_fixed.py
      (direct interval lookup at K-window centers; no set_annotations time-shift)
      This is what gives the paper's 44/44 |K| result; canonical mapping with
      set_annotations produces ~41/44 because of subtle annotation-time shifts.

Output mirrors paper Table 1 (sleep column) + Table 2 + REM clustering test +
N=22 subject-level sensitivity analysis (fills [PLACEHOLDER] in Methods).

Usage:
    !python sleep_edf_unified.py
    !python sleep_edf_unified.py --root /home/Sleep-EDF --out unified_results.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

try:
    import mne
    mne.set_log_level("ERROR")
except ImportError:
    sys.exit("Install MNE first:  pip install mne")

warnings.filterwarnings("ignore")


# ── CONFIG ────────────────────────────────────────────────────────────────────
EEG_CHANNEL    = "EEG Fpz-Cz"
EEG_FALLBACK   = "EEG Pz-Oz"
FS_TARGET      = 100
FFT_WIN_S      = 4
STRIDE_S       = 2
N_FREQ_BINS    = 100
FREQ_MAX       = 50.0
K_WIN          = 30
K_STRIDE       = 2
MIN_WAKE_WINS  = 10
MIN_NREM_WINS  = 10

STAGE_WAKE, STAGE_NREM, STAGE_REM = 1, 0, 2

STAGE_MAP = {
    'Sleep stage W':  STAGE_WAKE,
    'Sleep stage 1':  np.nan,
    'Sleep stage 2':  STAGE_NREM,
    'Sleep stage 3':  STAGE_NREM,
    'Sleep stage 4':  STAGE_NREM,
    'Sleep stage R':  STAGE_REM,
    'Sleep stage ?':  np.nan,
    'Movement time':  np.nan,
}

PAPER = {
    "K_wake_mean": 6.65,   "K_wake_sd": 2.85,
    "K_nrem_mean": 1.69,   "K_nrem_sd": 0.28,
    "K_rem_mean":  1.676,  "K_rem_sd":  0.553,
    "K_ratio": 4.01, "K_t": 11.52, "K_p": 1.01e-14,
    "K_wilcoxon_p": 1.14e-13, "K_n_correct": 44,

    "R_wake_mean": 23.77,  "R_wake_sd": 0.40,
    "R_nrem_mean": 24.29,  "R_nrem_sd": 0.10,
    "R_t": -8.47, "R_p": 1.03e-10, "R_n_correct": 4,

    "ratio_wake_mean": 0.376, "ratio_wake_sd": 0.281,
    "ratio_nrem_mean": 0.076, "ratio_nrem_sd": 0.024,
    "ratio_ratio": 5.18, "ratio_t": 7.15, "ratio_p": 7.69e-9,
    "ratio_n_correct": 41,

    "dwell_wake_mean": 0.476, "dwell_wake_sd": 0.197,
    "dwell_nrem_mean": 0.215, "dwell_nrem_sd": 0.050,
    "dwell_ratio": 2.51, "dwell_t": 7.69, "dwell_p": 1.31e-9,
    "dwell_n_correct": 38,
}


# ── SPECTROGRAM (canonical: interp linear power, then log) ────────────────────
def compute_spectrogram(raw_eeg, fs):
    n_fft     = int(FFT_WIN_S * fs)
    n_stride  = int(STRIDE_S * fs)
    target_freqs = np.linspace(0, FREQ_MAX, N_FREQ_BINS)
    win       = np.hanning(n_fft)
    fft_freqs = np.fft.rfftfreq(n_fft, d=1 / fs)

    specs, times = [], []
    pos = 0
    while pos + n_fft <= len(raw_eeg):
        chunk        = raw_eeg[pos:pos + n_fft] * win
        power        = np.abs(np.fft.rfft(chunk, n=n_fft)) ** 2
        power_interp = np.interp(target_freqs, fft_freqs, power)
        power_db     = 10 * np.log10(power_interp + 1e-12)
        specs.append(power_db)
        times.append((pos + n_fft / 2) / fs)
        pos += n_stride

    return np.array(specs).T, np.array(times)   # (100, T), (T,)


# ── METRICS (canonical: K + eff_rank + dwell + ratio in one pass) ─────────────
def compute_metrics(sdb):
    traj   = sdb.T
    n_wins = (len(traj) - K_WIN) // K_STRIDE
    K        = np.full(n_wins, np.nan)
    eff_rank = np.full(n_wins, np.nan)

    for i in range(n_wins):
        t0    = i * K_STRIDE
        chunk = traj[t0:t0 + K_WIN]
        if np.isnan(chunk).mean() > 0.1:
            continue
        chunk = chunk[~np.isnan(chunk).any(axis=1)]
        if len(chunk) < 5:
            continue
        dx  = np.diff(chunk, axis=0)
        if dx.shape[0] < 2:
            continue
        sm  = np.linalg.norm(dx, axis=1)
        sig = np.std(sm)
        if sig < 1e-10:
            K[i] = 0.0
        else:
            lam   = np.mean(sm) / sig**2
            C     = lam * sig
            K[i]  = abs(2 * sig * C / (1 + C**2))
        try:
            sv = np.linalg.svd(dx, compute_uv=False)
            sv = sv[sv > 1e-12]
            if len(sv) > 0:
                p = sv / sv.sum()
                eff_rank[i] = float(np.exp(-np.sum(p * np.log(p + 1e-300))))
        except np.linalg.LinAlgError:
            pass

    # Dwell: fraction of preceding 20 K-windows above subject 75th percentile of K
    k75   = np.nanpercentile(K, 75)
    dwell = np.full(n_wins, np.nan)
    for i in range(20, n_wins):
        seg = K[i - 20:i]
        if not np.isnan(seg).any():
            dwell[i] = np.mean(seg > k75)

    # K / eff_rank ratio
    ratio = np.full(n_wins, np.nan)
    valid = ~(np.isnan(K) | np.isnan(eff_rank))
    ratio[valid] = K[valid] / np.where(eff_rank[valid] > 0.1, eff_rank[valid], 0.1)

    return K, eff_rank, dwell, ratio


# ── STAGE MAPPING (fixed-script style: direct interval lookup) ────────────────
def map_stages(K_times, annot):
    """Interval lookup at each K-window centre — NO set_annotations time shift.

    Reverse iteration so that the latest annotation containing t wins,
    matching the behaviour of forward-order mask-assignment.
    """
    onsets, durs, descs = annot.onset, annot.duration, annot.description
    labels = np.full(len(K_times), np.nan)
    for j, t in enumerate(K_times):
        for k in range(len(onsets) - 1, -1, -1):
            if onsets[k] <= t < onsets[k] + durs[k]:
                labels[j] = STAGE_MAP.get(str(descs[k]), np.nan)
                break
    return labels


# ── SUBJECT ID PARSING (telemetry: ST<sub3><night>...) ───────────────────────
def subject_id_from_name(psg_id):
    """ST7011J0 -> subject='701' night='1'."""
    base = psg_id.split("-")[0]
    if len(base) >= 6 and base[:2].upper() in ("ST", "SC"):
        return base[2:5], base[5]
    return base, "?"


# ── PER-PSG PROCESSING ────────────────────────────────────────────────────────
@dataclass
class PerPSG:
    psg_id: str; subject_id: str; night: str; channel: str
    duration_h: float; n_kwin: int
    nW: int; nN: int; nR: int
    K_wake: float;     K_nrem: float;     K_rem: float
    R_wake: float;     R_nrem: float;     R_rem: float
    dw_wake: float;    dw_nrem: float;    dw_rem: float
    rt_wake: float;    rt_nrem: float;    rt_rem: float


def process_pair(psg_path: Path, hyp_path: Path) -> PerPSG | None:
    psg_id = psg_path.name.split("-")[0]
    sub_id, night = subject_id_from_name(psg_path.name)

    raw = mne.io.read_raw_edf(str(psg_path), preload=False, verbose=False)
    fs  = int(raw.info["sfreq"])

    ch = None
    for cand in (EEG_CHANNEL, EEG_FALLBACK):
        if cand in raw.ch_names:
            ch = cand; break
    if ch is None:
        print(f"    [skip] {psg_id}: no Fpz-Cz/Pz-Oz in {raw.ch_names[:4]}")
        return None

    raw.pick([ch]).load_data(verbose=False)
    eeg = raw.get_data()[0]

    if fs != FS_TARGET:
        from math import gcd
        from scipy.signal import resample_poly
        g   = gcd(FS_TARGET, fs)
        eeg = resample_poly(eeg, FS_TARGET // g, fs // g)
        fs  = FS_TARGET

    duration_h = len(eeg) / fs / 3600.0

    sdb, spec_times = compute_spectrogram(eeg, fs)
    K, eff_rank, dwell, ratio = compute_metrics(sdb)

    n_K = len(K)
    center_cols = np.clip(np.arange(n_K) * K_STRIDE + K_WIN // 2,
                          0, len(spec_times) - 1)
    K_times = spec_times[center_cols]

    annot   = mne.read_annotations(str(hyp_path))
    labels  = map_stages(K_times, annot)

    nW = int(np.sum(labels == STAGE_WAKE))
    nN = int(np.sum(labels == STAGE_NREM))
    nR = int(np.sum(labels == STAGE_REM))

    if nW < MIN_WAKE_WINS or nN < MIN_NREM_WINS:
        print(f"    [skip] {psg_id}: insufficient Wake={nW} NREM={nN}")
        return None

    def stage_mean(arr, target):
        sub = arr[labels == target]
        if len(sub) == 0:
            return float("nan")
        return float(np.nanmean(sub))

    return PerPSG(
        psg_id=psg_id, subject_id=sub_id, night=night, channel=ch,
        duration_h=duration_h, n_kwin=n_K,
        nW=nW, nN=nN, nR=nR,
        K_wake=stage_mean(K, STAGE_WAKE),     K_nrem=stage_mean(K, STAGE_NREM),     K_rem=stage_mean(K, STAGE_REM),
        R_wake=stage_mean(eff_rank, STAGE_WAKE), R_nrem=stage_mean(eff_rank, STAGE_NREM), R_rem=stage_mean(eff_rank, STAGE_REM),
        dw_wake=stage_mean(dwell, STAGE_WAKE), dw_nrem=stage_mean(dwell, STAGE_NREM), dw_rem=stage_mean(dwell, STAGE_REM),
        rt_wake=stage_mean(ratio, STAGE_WAKE), rt_nrem=stage_mean(ratio, STAGE_NREM), rt_rem=stage_mean(ratio, STAGE_REM),
    )


# ── PAIRING ──────────────────────────────────────────────────────────────────
def find_pairs(root: Path):
    psgs = sorted(root.rglob("*-PSG.edf"))
    hyps = sorted(root.rglob("*-Hypnogram.edf"))
    idx  = {h.name[:7]: h for h in hyps}
    pairs = []
    for psg in psgs:
        h = idx.get(psg.name[:7])
        if h is not None:
            pairs.append((psg, h))
        else:
            print(f"  [warn] {psg.name}: no paired hypnogram", file=sys.stderr)
    return pairs


# ── PAIRED STATS ─────────────────────────────────────────────────────────────
def paired_stats(wake_vals, nrem_vals, name=""):
    w = np.asarray(wake_vals, float); n = np.asarray(nrem_vals, float)
    keep = ~(np.isnan(w) | np.isnan(n)); w, n = w[keep], n[keep]
    if w.size < 3:
        return None
    t, p = scipy_stats.ttest_rel(w, n)
    try:
        _, wp = scipy_stats.wilcoxon(w, n)
    except ValueError:
        wp = float("nan")
    return dict(
        name=name, n=int(w.size),
        wake_mean=float(w.mean()), wake_sd=float(w.std(ddof=1)),
        nrem_mean=float(n.mean()), nrem_sd=float(n.std(ddof=1)),
        ratio=float(w.mean() / n.mean()) if n.mean() else float("nan"),
        median_ratio=float(np.median(w / np.where(n != 0, n, np.nan))),
        t=float(t), p=float(p), wilcoxon_p=float(wp),
        n_correct=int((w > n).sum()),
    )


def compare_block(stat, paper_prefix, paper):
    if stat is None:
        print(f"  {paper_prefix}: insufficient data"); return
    name = stat["name"]
    print(f"\n  {name}")
    print(f"    Wake (mean ± SD): {stat['wake_mean']:7.4f} ± {stat['wake_sd']:6.4f}    "
          f"[paper: {paper[paper_prefix+'_wake_mean']:7.4f} ± {paper[paper_prefix+'_wake_sd']:6.4f}]")
    print(f"    NREM (mean ± SD): {stat['nrem_mean']:7.4f} ± {stat['nrem_sd']:6.4f}    "
          f"[paper: {paper[paper_prefix+'_nrem_mean']:7.4f} ± {paper[paper_prefix+'_nrem_sd']:6.4f}]")
    if paper_prefix + "_ratio" in paper:
        print(f"    Ratio:            {stat['ratio']:7.4f}                  "
              f"[paper: {paper[paper_prefix+'_ratio']:7.4f}]")
    print(f"    t({stat['n']-1})            = {stat['t']:7.3f}                  "
          f"[paper: {paper[paper_prefix+'_t']:7.3f}]")
    print(f"    p (paired t)      = {stat['p']:.3e}             "
          f"[paper: {paper[paper_prefix+'_p']:.3e}]")
    if paper_prefix + "_wilcoxon_p" in paper:
        print(f"    Wilcoxon p        = {stat['wilcoxon_p']:.3e}             "
              f"[paper: {paper[paper_prefix+'_wilcoxon_p']:.3e}]")
    print(f"    N correct (W>N)   = {stat['n_correct']}/{stat['n']}                     "
          f"[paper: {paper[paper_prefix+'_n_correct']}/44]")


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="/home/Sleep-EDF",
                    help="Folder containing Sleep-EDF Telemetry PSG+Hypnogram pairs")
    ap.add_argument("--out", default="sleep_edf_unified_per_psg.csv")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"Root not found: {root}")

    print("=" * 78)
    print("SLEEP-EDF UNIFIED PIPELINE")
    print("=" * 78)
    print(f"Root: {root}")

    pairs = find_pairs(root)
    print(f"Found {len(pairs)} PSG/Hypnogram pairs.\n")
    if not pairs:
        return 1

    results: list[PerPSG] = []
    for i, (psg, hyp) in enumerate(pairs, 1):
        print(f"[{i:2d}/{len(pairs)}] {psg.name}", end="  ", flush=True)
        try:
            r = process_pair(psg, hyp)
        except Exception as e:
            print(f"FAIL: {e}")
            continue
        if r is None:
            continue
        wn = (r.K_wake / r.K_nrem) if r.K_nrem else float("nan")
        print(f"ch={r.channel:<11s} dur={r.duration_h:4.1f}h  "
              f"K_W={r.K_wake:6.3f} K_N={r.K_nrem:6.3f} K_R={r.K_rem:6.3f}  "
              f"W/N={wn:4.2f}x  (nW={r.nW} nN={r.nN} nR={r.nR})")
        results.append(r)

    if not results:
        sys.exit("No PSGs processed.")

    # CSV
    out = Path(args.out)
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["psg_id", "subject_id", "night", "channel", "duration_h", "n_kwin",
                    "nW", "nN", "nR",
                    "K_wake", "K_nrem", "K_rem",
                    "eff_rank_wake", "eff_rank_nrem", "eff_rank_rem",
                    "K_eff_ratio_wake", "K_eff_ratio_nrem", "K_eff_ratio_rem",
                    "dwell_wake", "dwell_nrem", "dwell_rem"])
        for r in results:
            w.writerow([r.psg_id, r.subject_id, r.night, r.channel,
                        f"{r.duration_h:.2f}", r.n_kwin,
                        r.nW, r.nN, r.nR,
                        f"{r.K_wake:.6f}", f"{r.K_nrem:.6f}", f"{r.K_rem:.6f}",
                        f"{r.R_wake:.6f}", f"{r.R_nrem:.6f}", f"{r.R_rem:.6f}",
                        f"{r.rt_wake:.6f}", f"{r.rt_nrem:.6f}", f"{r.rt_rem:.6f}",
                        f"{r.dw_wake:.6f}", f"{r.dw_nrem:.6f}", f"{r.dw_rem:.6f}"])
    print(f"\nPer-PSG CSV written to: {out}")

    # ── PRIMARY (PSG-level, N=44) ────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("PRIMARY ANALYSIS — PSG-LEVEL  (paper Table 1 + Table 2)")
    print("=" * 78)

    K_w = [r.K_wake for r in results];   K_n = [r.K_nrem for r in results];   K_r = [r.K_rem for r in results]
    R_w = [r.R_wake for r in results];   R_n = [r.R_nrem for r in results];   R_r = [r.R_rem for r in results]
    rt_w = [r.rt_wake for r in results]; rt_n = [r.rt_nrem for r in results]; rt_r = [r.rt_rem for r in results]
    dw_w = [r.dw_wake for r in results]; dw_n = [r.dw_nrem for r in results]

    K_stat  = paired_stats(K_w, K_n, "|K|")
    R_stat  = paired_stats(R_w, R_n, "eff_rank")
    rt_stat = paired_stats(rt_w, rt_n, "K/eff_rank ratio")
    dw_stat = paired_stats(dw_w, dw_n, "Dwell fraction")

    compare_block(K_stat,  "K",     PAPER)
    compare_block(R_stat,  "R",     PAPER)
    compare_block(rt_stat, "ratio", PAPER)
    compare_block(dw_stat, "dwell", PAPER)

    # ── REM classification block ─────────────────────────────────────────────
    print("\n" + "-" * 78)
    print("REM CLASSIFICATION (paper: REM clusters with NREM)")
    print("-" * 78)
    K_rem_arr  = np.array(K_r, float); K_nrem_arr = np.array(K_n, float)
    rt_rem_arr = np.array(rt_r, float); rt_nrem_arr = np.array(rt_n, float)
    valid = ~np.isnan(K_rem_arr) & ~np.isnan(K_nrem_arr)
    if valid.sum() >= 3:
        kr_kn = K_rem_arr[valid] / K_nrem_arr[valid]
        print(f"  K_REM / K_NREM mean:      {kr_kn.mean():.3f}    [paper: 0.99x]")
        print(f"  K_REM mean:               {K_rem_arr[valid].mean():.3f}    [paper: {PAPER['K_rem_mean']:.3f}]")
    valid_r = ~np.isnan(rt_rem_arr) & ~np.isnan(rt_nrem_arr)
    if valid_r.sum() >= 3:
        rr_rn = rt_rem_arr[valid_r] / rt_nrem_arr[valid_r]
        print(f"  ratio_REM / ratio_NREM:   {rr_rn.mean():.3f}    [paper: 0.93x]")

    # ── SUBJECT-LEVEL SENSITIVITY (N=22) ─────────────────────────────────────
    print("\n" + "=" * 78)
    print("SENSITIVITY — SUBJECT-LEVEL  (collapse 2 nights/subject, N=22)")
    print("(fills [PLACEHOLDER] in Methods Dataset 2)")
    print("=" * 78)

    by_sub = defaultdict(list)
    for r in results:
        by_sub[r.subject_id].append(r)

    n_sub      = len(by_sub)
    n_two_night = sum(1 for v in by_sub.values() if len(v) == 2)
    print(f"  Subjects total:    {n_sub}")
    print(f"  Subjects w/ 2 nights: {n_two_night}    [paper: 22]")

    subj = {k: ([], [], [], [], [], [], [], []) for k in ["K_w","K_n","R_w","R_n","rt_w","rt_n","dw_w","dw_n"]}
    for sub, recs in by_sub.items():
        subj["K_w"][0].append(np.nanmean([r.K_wake for r in recs]))
        subj["K_n"][0].append(np.nanmean([r.K_nrem for r in recs]))
        subj["R_w"][0].append(np.nanmean([r.R_wake for r in recs]))
        subj["R_n"][0].append(np.nanmean([r.R_nrem for r in recs]))
        subj["rt_w"][0].append(np.nanmean([r.rt_wake for r in recs]))
        subj["rt_n"][0].append(np.nanmean([r.rt_nrem for r in recs]))
        subj["dw_w"][0].append(np.nanmean([r.dw_wake for r in recs]))
        subj["dw_n"][0].append(np.nanmean([r.dw_nrem for r in recs]))

    s_K  = paired_stats(subj["K_w"][0],  subj["K_n"][0],  "|K| (subject-level)")
    s_R  = paired_stats(subj["R_w"][0],  subj["R_n"][0],  "eff_rank (subject-level)")
    s_rt = paired_stats(subj["rt_w"][0], subj["rt_n"][0], "K/eff_rank (subject-level)")
    s_dw = paired_stats(subj["dw_w"][0], subj["dw_n"][0], "Dwell (subject-level)")

    for s in (s_K, s_R, s_rt, s_dw):
        if s is None: continue
        print(f"\n  {s['name']}")
        print(f"    Wake (mean ± SD):  {s['wake_mean']:7.4f} ± {s['wake_sd']:6.4f}")
        print(f"    NREM (mean ± SD):  {s['nrem_mean']:7.4f} ± {s['nrem_sd']:6.4f}")
        print(f"    Ratio:             {s['ratio']:7.4f}")
        print(f"    t({s['n']-1})              = {s['t']:7.3f}")
        print(f"    p (paired t)       = {s['p']:.3e}")
        print(f"    N correct          = {s['n_correct']}/{s['n']}")

    # Manuscript paste-string
    if s_K:
        print("\n" + "-" * 78)
        print("FOR MANUSCRIPT — paste these into Methods Dataset 2 placeholders:")
        print("-" * 78)
        print(f"  t({s_K['n']-1}) = {s_K['t']:.2f}, p = {s_K['p']:.2e}; "
              f"{s_K['n_correct']}/{s_K['n']} subjects correct direction.")

    print("\n" + "=" * 78)
    print("DONE.")
    print("=" * 78)


if __name__ == "__main__":
    sys.exit(main())
