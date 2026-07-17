# VAMP Repair Program — AIES-26 submission 396 → revision R1

Inputs: main paper (v5), supplementary materials, three reviews (18 itemized
comments R1–R18; Reviewer A long-form, reject-leaning; Reviewer B, 4
strengths / 5 weaknesses). Deliverables: `paper_r1.tex/.pdf`,
`supplement_r1.tex/.pdf`, `vamp_audit.py`, `experiments.py`, `results.json`.

## Load-bearing design decisions

**D1 — m-ary VAMP.** VAMP generalized to arity m ≥ 1: from assertions
A_1..A_m and implication (⋀A_j → B), derive B. Activation: every assertion
grade < η. Verification: α(B) ≤ Σα(A_j) + α(⋀A_j→B) + T(H(s)). m=1 recovers
the original rule. Gödel-conjunction non-expansion (supp Lemma 2) keeps the
bundled antecedent's grade ≤ max of the parts, so activation semantics are
undisturbed. Kills R18 and legitimizes the worked example and every step in
Appendix C.

**D2 — Explicit rule regimes.** The implication premise of a VAMP step is a
proof node. Two declared regimes: *audited-rule* (implications carry their
own α, entering verification and budget) and *certified-rule* (designer
declares α = 0 for the rule base; declaration is itself part of the audit
input). Corollary 2 now states both bounds: (2d+1)·T(H) general,
(d+1)·T(H) certified-rule. Kills R15 / Reviewer A's Corollary-2 objection /
overall-assessment item 1.

**D3 — Rule-aware budget (Definition 4').** Budget recursion = verification
condition with grades replaced by budgets. VAMP node: Σ budgets of all
premises (assertions + implication) + T. VGP node: min of the pair's budgets
+ T (sound by negation preservation; tighter than summing; supp Lemma 3
aligned). Eliminates the formula-vs-practice mismatches the reviewers kept
finding (R15, R18, supp A.3 vs Def 4).

**D4 — DAG-native linear audit (Proposition 5 restated — an upgrade, not a
patch).** The tree-unfolded budget satisfies its own recursion, so a single
memoized topological pass over the DAG computes it exactly; no unfolding is
ever materialized (supp Lemma, DAG-evaluation). Audit runtime is O(n+e) in
the *DAG*, full stop. What is exponential in depth is the budget *value*
(certifiability), not the computation. The original paper's "linear in the
unfolded object, which may be exponentially larger" claim was needlessly
pessimistic and is deleted. Kills R16; reframes R13; empirically confirmed
(E3: depth-4000 ladder, unfolding ≈ 2^4002 nodes, audited in 0.16 s).

**D5 — Certifiability horizon, owned.** Closed forms stated: chains certify
iff α(v_0) + Σ T(H(s_i)) < ϑ*; balanced binary trees iff
2^d·a + (2^d−1)·T < ϑ*. The "discrimination = counting tension exposure"
observation (Reviewer A's "thin basis") is now stated by the paper itself as
the declared mechanism, with the unsoundness counterexample explaining why
nothing tighter is available under VAMP. Neutral wording replaces the
unsupported "real proof objects rarely realize this worst case" (R13).

**D6 — Action typing + fail-closed deployment.** Action expressions R(A)
live in L_act; action nodes are terminal (sinks; audited structurally);
canonical default valuation v(R(A)) := v(A) provided, or designer-supplied
and checked by schema-boundedness. Dead zone [η, ϑ*): any attempted step
fails activation → audit rejects → default action suppressed → declared
fallback F (escalate-to-human / refuse / defer). Kills R10, R11; answers
Reviewer A's "are action rules implications or not" via a
"Schemas, not conditionals" remark tying Prop 3 to schemas.

**D7 — Reference implementation + empirical checks (new Section 8).**
Exact-rational audit implementation; worked examples machine-checked;
Theorem 1 validated on 10^4 random VAMP-valid objects (119,876 nodes, 0
violations); 2,000/2,000 single-grade mutants detected; 80,000-node audit
in 0.53 s (~6 µs/node); horizon table reproducing Reviewer A's
once-certified/twice-rejected arithmetic as designed behavior. Converts
Reviewer B #1 from "no implementation" to "reference implementation +
synthetic evaluation" (honestly scoped: no live-system integration).

**D8 — Related work section (new §2) + scope honesty.** Annotated
paraconsistent logic programming, bilattice logics, graded/weighted
argumentation, possibilistic logic, resource-bounded/substructural logics,
runtime verification, shielding and correct-by-construction runtime
enforcement (both reviewer-named citations included), calibration, CoT
faithfulness. New subsections: "What Theorem 1 does and does not establish"
(Reviewer A ask #1, verbatim honored: qualified the negative-to-positive
claim); "Deployment architectures" (Reviewer B #3, faithfulness assumption
explicit with citation).

## Concern → fix ledger

| # (sev) | Concern | Fix | Where |
|---|---|---|---|
| R1 (Min) | Abstract/Intro identical opening | Abstract opening rewritten | Abstract |
| R2 (Min) | Jargon in contribution 1 | Parenthetical glosses added | §1 |
| R3 (Min) | Formula/value same symbol | v(A) = (τ_A, φ_A) notation throughout | §3 |
| R4 (Mod) | Heterogeneous-source calibration | Calibration-assumed paragraph + refs + named future work | §3, §9 |
| R5 (Min) | "geometric" | "structural" | §4 |
| R6 (Min) | "flows through inference" | "accumulates through inference" everywhere | Abs., §4, §9 |
| R7 (Min) | Plain-text pseudo-formula | Formal display: α(B) ≤ max_{A∈Δ} α(A) | §4 |
| R8 (Maj) | Ch.4 prose contradicts additive formula | Prose rewritten: "exceed the **sum** of the premises' grades by at most T" | §5 |
| R9 (Mod) | "α(premises)" ambiguous | Unified claim restated in plain text; per-rule forms displayed | §5 |
| R10 (Maj) | α of R(A) undefined; type discipline | D6: L_act, default v(R(A)):=v(A), terminality | §5, §7 |
| R11 (Maj) | Dead-zone operational behavior | D6: fail-closed deployment definition; Algorithm routes to fallback | §5, §7 |
| R12 (Min) | R_high = ∅ vacuous coordination | R_high ≠ ∅ added to Def. 2/3 + Algorithm step 0 | §5, §7 |
| R13 (Maj) | Exponential blowup on CoT DAGs; "rare" unsupported | D4 (no computational blowup) + D5 (horizon owned, neutral wording, counterexample = why conservatism is forced) | §6–§8, supp C |
| R14 (Min) | "once per use" ambiguous | "counted once per path / duplicated across each distinct path" | §6 |
| R15 (Maj) | Corollary 2 ignores implication budget | D2: both regimes, both bounds; example declares regime; arithmetic redone and machine-checked | §6, §7 |
| R16 (Mod) | Prop 5 statement vs caveat | D4: Prop 5 restated linear-on-DAG; caveat now about value, in the statement; bit-complexity footnote | §7, supp C |
| R17 (Min) | Reals→rationals unjustified | Decidability sentence at audit inputs | §7 |
| R18 (Maj) | s_3 three premises unlicensed | D1: m-ary VAMP; example steps re-typed s_1..s_4 as VAMP_1/VAMP_3; supp C steps re-licensed | §5, §7, supp C |
| Rev A | Theorem near-tautological; qualify claims | "Does and does not establish" subsection; claim reworded | §6 |
| Rev A | Action rules: implications or not | "Schemas, not conditionals" remark | §4 |
| Rev A | Situate vs annotated paraconsistent / RV / resource logics | Related Work §2 | §2 |
| Rev A | Connect valuations or reframe | Scope: certificate relative to declared valuations; calibration open | §3, §9 |
| Rev A | Sum-not-max conservatism | D5 + supp C: tighter is unsound; conservatism forced | §6, supp C |
| Rev B #1 | No implementation | D7 | §8 |
| Rev B #2 | Grades' source deferred | Same as Rev A valuations | §3, §9 |
| Rev B #3 | Who constructs the proof object | Deployment architectures + faithfulness citation | §7 |
| Rev B #4 | Counterexample not in package | Supp C included and machine-checked; merged summary in §6 | supp C |
| Rev B #5 | Thin related work | §2 (reviewer-named refs included) | §2 |

## Residual risks (a revision cannot fix these)

1. **Significance.** Theorem 1 remains a soundness-of-accounting induction.
   The revision owns this instead of overclaiming; a reviewer who wants
   mathematical depth will still not find it. Mitigation only.
2. **Live-system gap.** Reference implementation + synthetic evaluation is
   not an integration with a real reasoner. Faithful proof-object
   construction from LLM traces is named as an open problem, not solved.
3. **Calibration.** Valuation sourcing remains exogenous by design; stated,
   cited, and scoped — not resolved.
4. **Bib placeholders.** Könighofer et al. 2022 author list abbreviated
   ("et al.") — complete before submission. Verify Dunne et al. 2011 and
   Baroni/Rago/Toni 2018 details against the originals.

## Venue notes

- Two-column `article` approximation of AAAI format; re-flow into the
  official style file at submission. Current length will exceed a 7-page
  limit: compress §2 (related work to ~0.5 col), §6 horizon discussion, and
  the certified variant of the worked example first.
- AAAI-27 window (abs 7/21, paper 7/28) collides with the corrigibility
  paper. Alternatives: SafeAI/AAAI workshop, KR, JELIA, or a journal (JAIR)
  where the logic-first framing is native. Reviewer A's "logic paper" remark
  is a venue signal, not only a criticism.
