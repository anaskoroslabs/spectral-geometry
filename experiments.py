"""
experiments.py — Empirical checks reported in the revised paper (Sec. 8).

E1  Machine-check of the worked examples and the Appendix C counterexample.
E2  Theorem 1 validation on random VAMP-valid proof objects (alpha <= B at
    every node; zero violations expected) and mutation testing (a single
    grade perturbed above its verification bound must be flagged).
E3  Runtime scaling. (a) Width scaling at fixed depth: audit time linear in
    DAG size. (b) Diamond-ladder depth scaling: tree-unfolding is
    exponential in d, the DAG-native audit is not; budget bit-length grows
    linearly in d (the O((n+e)(d+b)) bit-cost statement).
E4  Certifiability horizon: certification rate vs. proof depth under a
    random mix of plain and high-tension steps (audited-rule regime).
"""

from __future__ import annotations
import json
import random
import time
from fractions import Fraction

from vamp_audit import (Node, ProofObject, vamp_audit, F,
                        worked_example_rejected, worked_example_certified,
                        counterexample_dag)

random.seed(396)
RESULTS = {}


# ----------------------------------------------------------------- E1 ----
def e1_machine_checks():
    out = {}
    pi, T, theta, Rh, eta = worked_example_rejected()
    r = vamp_audit(pi, T, theta, Rh, eta)
    out["rejected_variant"] = {
        "kappa": r.kappa,
        "beta_P6": str(r.beta["P6"]), "beta_P7": str(r.beta["P7"]),
        "theta_star": str(r.theta_star),
    }
    assert r.kappa == "rejected" and r.beta["P7"] == F("0.70")

    pi, T, theta, Rh, eta = worked_example_certified()
    r = vamp_audit(pi, T, theta, Rh, eta)
    out["certified_variant"] = {
        "kappa": r.kappa, "beta_P6p": str(r.beta["P6"]),
    }
    assert r.kappa == "certified" and r.beta["P6"] == F("0.40")

    pi, T, naive = counterexample_dag(F(1, 10), F(1, 5))
    r = vamp_audit(pi, T, {"x": F(2)}, {"x"}, F(1))
    out["counterexample"] = {
        "alpha_v": str(pi.nodes["v"].alpha),
        "B_tree_v": str(r.beta["v"]),
        "naive_B_DAG_v": str(naive),
        "naive_unsound": bool(pi.nodes["v"].alpha > naive),
        "tree_budget_dominates": bool(pi.nodes["v"].alpha <= r.beta["v"]),
    }
    assert pi.nodes["v"].alpha > naive
    assert pi.nodes["v"].alpha <= r.beta["v"]
    RESULTS["E1"] = out


# ----------------------------------------------------------------- E2 ----
def random_valid_object(depth: int, tension: Fraction, saturate=False):
    """Random VAMP-valid proof DAG, audited-rule regime."""
    rnd = lambda hi: F(random.randint(0, int(hi * 100)), 100)
    nodes, frontier, T = [], [], {}
    for i in range(random.randint(2, 4)):
        n = Node(f"p{i}", "premise", rnd(0.10))
        nodes.append(n); frontier.append(n)
    for d in range(depth):
        cls = f"c{d}"
        T[cls] = tension if random.random() < 0.5 else F(0)
        m = min(random.choice([1, 1, 2]), len(frontier))
        asserts = random.sample(frontier, m)
        w = Node(f"w{d}", "premise", rnd(0.05), is_rule=True)
        nodes.append(w)
        bound = sum((a.alpha for a in asserts), F(0)) + w.alpha + T[cls]
        if saturate:
            alpha = min(bound, F(1))
        else:
            alpha = min(F(random.randint(0, int(bound * 100)), 100), F(1))
        v = Node(f"v{d}", "vamp", alpha,
                 premises=[a.nid for a in asserts], imp=w.nid, cls=cls,
                 is_action=(d == depth - 1))
        nodes.append(v); frontier.append(v)
    return ProofObject(nodes), T


def e2_theorem_validation(n_objects=10000, mutants=2000):
    viol, checked = 0, 0
    for _ in range(n_objects):
        pi, T = random_valid_object(random.randint(1, 8), F(3, 10))
        r = vamp_audit(pi, T, {"x": F(10)}, {"x"}, F(2))  # thresholds vacuous
        assert r.kappa == "certified", r.failure
        for nid, n in pi.nodes.items():
            checked += 1
            if n.alpha > r.beta[nid]:
                viol += 1
    detected = 0
    for _ in range(mutants):
        pi, T = random_valid_object(random.randint(1, 6), F(3, 10))
        derived = [n for n in pi.nodes.values() if n.kind == "vamp"]
        m = random.choice(derived)
        asserts = [pi.nodes[p].alpha for p in m.premises]
        bound = sum(asserts, F(0)) + pi.nodes[m.imp].alpha + T[m.cls]
        m.alpha = bound + F(1, 100)          # push past the bound
        r = vamp_audit(pi, T, {"x": F(10)}, {"x"}, F(2))
        if r.kappa == "rejected":
            detected += 1
    RESULTS["E2"] = {
        "objects": n_objects, "nodes_checked": checked,
        "budget_violations": viol,
        "mutants": mutants, "mutants_detected": detected,
        "detection_rate": detected / mutants,
    }


# ----------------------------------------------------------------- E3 ----
def width_dag(n_blocks: int):
    """Fixed depth 2, growing width: n_blocks independent 3-node blocks
    merged pairwise; DAG size ~ 5*n_blocks."""
    nodes = []
    tops = []
    for i in range(n_blocks):
        p = Node(f"p{i}", "premise", F(1, 20))
        w = Node(f"w{i}", "premise", F(0), is_rule=True)
        v = Node(f"v{i}", "vamp", F(1, 20), premises=[p.nid], imp=w.nid,
                 cls="p")
        nodes += [p, w, v]; tops.append(v)
    for i in range(0, n_blocks - 1, 2):
        w = Node(f"m{i}", "premise", F(0), is_rule=True)
        v = Node(f"z{i}", "vamp", F(1, 10),
                 premises=[tops[i].nid, tops[i + 1].nid], imp=w.nid, cls="p")
        nodes += [w, v]
    return ProofObject(nodes), {"p": F(0)}


def ladder_dag(d: int, t=F(3, 10)):
    """Diamond ladder of depth d: x0 -> (a_i, b_i) -> x_{i+1}. Tree
    unfolding has ~3*2^d nodes; the DAG has 3d+1 (+rules)."""
    nodes = [Node("x0", "premise", F(1, 100))]
    prev = "x0"
    for i in range(d):
        wa = Node(f"wa{i}", "premise", F(0), is_rule=True)
        wb = Node(f"wb{i}", "premise", F(0), is_rule=True)
        wm = Node(f"wm{i}", "premise", F(0), is_rule=True)
        a = Node(f"a{i}", "vamp", F(1, 100), premises=[prev], imp=wa.nid,
                 cls="H")
        b = Node(f"b{i}", "vamp", F(1, 100), premises=[prev], imp=wb.nid,
                 cls="H")
        x = Node(f"x{i+1}", "vamp", F(1, 100), premises=[a.nid, b.nid],
                 imp=wm.nid, cls="H")
        nodes += [wa, wb, wm, a, b, x]
        prev = x.nid
    return ProofObject(nodes), {"H": t}


def e3_runtime():
    out = {"width": [], "ladder": []}
    for nb in [200, 2000, 20000]:
        pi, T = width_dag(nb)
        n = len(pi.nodes)
        t0 = time.perf_counter()
        r = vamp_audit(pi, T, {"x": F(10)}, {"x"}, F(2))
        dt = time.perf_counter() - t0
        out["width"].append({"dag_nodes": n, "seconds": round(dt, 4),
                             "us_per_node": round(1e6 * dt / n, 2)})
    for d in [100, 1000, 4000]:
        pi, T = ladder_dag(d)
        n = len(pi.nodes)
        t0 = time.perf_counter()
        r = vamp_audit(pi, T, {"x": F(10 ** 30)}, {"x"}, F(2))
        dt = time.perf_counter() - t0
        top = r.beta[f"x{d}"]
        bits = top.numerator.bit_length()
        out["ladder"].append({
            "depth": d, "dag_nodes": n,
            "unfolded_nodes_approx": f"~2^{d + 2}",
            "seconds": round(dt, 4),
            "budget_bits": bits,
        })
    RESULTS["E3"] = out


# ----------------------------------------------------------------- E4 ----
def e4_horizon(trials=2000):
    theta_star = F(6, 10)
    rows = []
    for depth in range(1, 9):
        cert = 0
        for _ in range(trials):
            pi, T = random_valid_object(depth, F(3, 10), saturate=False)
            r = vamp_audit(pi, T, {"esc": theta_star}, {"esc"}, F(2))
            cert += (r.kappa == "certified")
        rows.append({"depth": depth, "cert_rate": round(cert / trials, 3)})
    # High-tension-step count view (chains, certified rules, saturating):
    ks = []
    for k in range(0, 5):
        # chain with k tension steps (T=0.3) and 4-k plain steps
        nodes = [Node("p", "premise", F(1, 20))]
        prev, T = "p", {}
        order = ["H"] * k + ["P"] * (4 - k)
        random.shuffle(order)
        for i, c in enumerate(order):
            cls = f"c{i}"
            T[cls] = F(3, 10) if c == "H" else F(0)
            w = Node(f"w{i}", "premise", F(0), is_rule=True)
            prev_alpha = nodes[-1].alpha if i else nodes[0].alpha
            v = Node(f"v{i}", "vamp",
                     min(F(1), prev_alpha + T[cls]),
                     premises=[prev], imp=w.nid, cls=cls,
                     is_action=(i == 3))
            nodes += [w, v]; prev = v.nid
        pi = ProofObject(nodes)
        r = vamp_audit(pi, T, {"esc": F(6, 10)}, {"esc"}, F(2))
        ks.append({"high_tension_steps": k,
                   "B_at_action": str(r.beta[prev]),
                   "kappa": r.kappa})
    RESULTS["E4"] = {"random_dags_T0.3_theta0.6": rows,
                     "chain_by_tension_count": ks}


if __name__ == "__main__":
    e1_machine_checks()
    e2_theorem_validation()
    e3_runtime()
    e4_horizon()
    with open("results.json", "w") as f:
        json.dump(RESULTS, f, indent=2)
    print(json.dumps(RESULTS, indent=2))
