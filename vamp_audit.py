r"""
vamp_audit.py — Reference implementation of VAMP-AUDIT (revised calculus).

Implements the repaired LPA-alpha audit:
  * m-ary VAMP steps: assertions A_1..A_m + implication premise (/\A_j -> B).
  * VGP steps: contradictory pair (A, ~A) -> R(A), min-budget convention.
  * Rule-aware tree-unfolded budgets computed DAG-natively in one
    topological pass (no explicit unfolding is materialized).
  * Fail-closed decision: certified | rejected (+ fallback routing).
  * Exact rational arithmetic (fractions.Fraction).

Node kinds:
  premise : leaf; carries alpha. May be flagged is_rule (an implication
            premise). In the certified-rule regime, rule premises are
            declared with alpha = 0.
  vamp    : derived by an m-ary VAMP step. premises = [a_1..a_m],
            imp = id of the implication premise node, cls = collapse class.
  vgp     : derived by a VGP step. premises = [u, u_neg], cls = collapse
            class, schema = name of consequence schema R.

Action-licensing nodes (is_action=True) must be sinks (terminality check).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional

F = Fraction


@dataclass
class Node:
    nid: str
    kind: str                      # 'premise' | 'vamp' | 'vgp'
    alpha: Fraction
    premises: list = field(default_factory=list)   # assertion premise ids
    imp: Optional[str] = None      # implication premise id (vamp only)
    cls: Optional[str] = None      # collapse class id (derived nodes)
    schema: Optional[str] = None   # consequence schema name (vgp only)
    is_action: bool = False
    is_rule: bool = False          # implication/rule premise node


@dataclass
class AuditResult:
    kappa: str                     # 'certified' | 'rejected'
    sigma: dict                    # step id -> 'pass'|'fail:<reason>'
    beta: dict                     # node id -> Fraction budget
    rho: dict                      # node id -> 'below'|'at-or-above'
    theta_star: Fraction
    fallback_engaged: bool
    failure: Optional[str] = None  # first failure description


class ProofObject:
    def __init__(self, nodes: list[Node]):
        self.nodes = {n.nid: n for n in nodes}
        if len(self.nodes) != len(nodes):
            raise ValueError("duplicate node ids")

    def all_premise_ids(self, n: Node) -> list[str]:
        ids = list(n.premises)
        if n.imp is not None:
            ids.append(n.imp)
        return ids

    def topological_order(self) -> list[str]:
        indeg = {nid: 0 for nid in self.nodes}
        children = {nid: [] for nid in self.nodes}
        for n in self.nodes.values():
            for p in self.all_premise_ids(n):
                if p not in self.nodes:
                    raise ValueError(f"{n.nid}: missing premise {p}")
                indeg[n.nid] += 1
                children[p].append(n.nid)
        queue = [nid for nid, d in indeg.items() if d == 0]
        order = []
        while queue:
            v = queue.pop()
            order.append(v)
            for c in children[v]:
                indeg[c] -= 1
                if indeg[c] == 0:
                    queue.append(c)
        if len(order) != len(self.nodes):
            raise ValueError("proof object is cyclic")
        return order

    def sinks(self) -> set[str]:
        used = set()
        for n in self.nodes.values():
            used.update(self.all_premise_ids(n))
        return set(self.nodes) - used


def vamp_audit(pi: ProofObject,
               T: dict[str, Fraction],
               theta: dict[str, Fraction],
               R_high: set[str],
               eta: Fraction) -> AuditResult:
    """Algorithm 1 (revised): single topological pass over the proof DAG."""
    if not R_high:
        raise ValueError("well-formedness: R_high must be nonempty")
    theta_star = min(theta[r] for r in R_high)

    sigma: dict[str, str] = {}
    beta: dict[str, Fraction] = {}
    rho: dict[str, str] = {}

    def reject(reason: str) -> AuditResult:
        return AuditResult("rejected", sigma, beta, rho, theta_star,
                           fallback_engaged=True, failure=reason)

    # Structural well-formedness: action nodes are terminal (sinks).
    sink_set = pi.sinks()
    for n in pi.nodes.values():
        if n.is_action and n.nid not in sink_set:
            return reject(f"type discipline: action node {n.nid} is not terminal")

    for nid in pi.topological_order():
        n = pi.nodes[nid]
        if n.kind == "premise":
            beta[nid] = n.alpha
            continue

        if n.kind == "vamp":
            step = f"s({nid})"
            asserts = [pi.nodes[p] for p in n.premises]
            imp = pi.nodes[n.imp] if n.imp is not None else None
            if imp is None:
                sigma[step] = "fail:missing-implication-premise"
                return reject(f"{step}: VAMP step lacks implication premise")
            # Activation: every assertion premise strictly below eta.
            for a in asserts:
                if not (a.alpha < eta):
                    sigma[step] = "fail:activation"
                    return reject(f"{step}: activation alpha({a.nid})="
                                  f"{a.alpha} >= eta={eta}")
            # Verification: alpha(B) <= sum alpha(A_j) + alpha(imp) + T.
            bound = sum((a.alpha for a in asserts), F(0)) + imp.alpha + T[n.cls]
            if not (n.alpha <= bound):
                sigma[step] = "fail:verification"
                return reject(f"{step}: alpha={n.alpha} > bound={bound}")
            sigma[step] = "pass"
            beta[nid] = sum((beta[a.nid] for a in asserts), F(0)) \
                + beta[imp.nid] + T[n.cls]
            continue

        if n.kind == "vgp":
            step = f"s({nid})"
            if len(n.premises) != 2:
                sigma[step] = "fail:arity"
                return reject(f"{step}: VGP requires the pair (A, ~A)")
            u, un = (pi.nodes[p] for p in n.premises)
            # Negation preservation (Lemma 1): alpha(A) == alpha(~A).
            if u.alpha != un.alpha:
                sigma[step] = "fail:negation-preservation"
                return reject(f"{step}: alpha({u.nid}) != alpha({un.nid})")
            if n.schema not in theta:
                sigma[step] = "fail:unknown-schema"
                return reject(f"{step}: schema {n.schema} has no threshold")
            # Activation: alpha(A) >= theta(R).
            if not (u.alpha >= theta[n.schema]):
                sigma[step] = "fail:activation"
                return reject(f"{step}: alpha={u.alpha} < theta({n.schema})")
            # Schema boundedness: alpha(R(A)) <= alpha(A) + T.
            if not (n.alpha <= u.alpha + T[n.cls]):
                sigma[step] = "fail:schema-boundedness"
                return reject(f"{step}: alpha={n.alpha} > "
                              f"{u.alpha + T[n.cls]}")
            sigma[step] = "pass"
            # Min-budget convention (Lemma 3'): the pair carries one grade.
            beta[nid] = min(beta[u.nid], beta[un.nid]) + T[n.cls]
            continue

        return reject(f"unknown node kind {n.kind} at {nid}")

    # High-stakes scan over action-licensing nodes.
    kappa = "certified"
    for nid, n in pi.nodes.items():
        rho[nid] = "below" if beta[nid] < theta_star else "at-or-above"
        if n.is_action and rho[nid] == "at-or-above":
            kappa = "rejected"
    return AuditResult(kappa, sigma, beta, rho, theta_star,
                       fallback_engaged=(kappa == "rejected"))


# ---------------------------------------------------------------------------
# Worked example (Section 6): cross-jurisdictional content moderation.
# Certified-rule regime: all rule-base implications carry alpha = 0.
# ---------------------------------------------------------------------------

def worked_example_rejected() -> tuple[ProofObject, dict, dict, set, Fraction]:
    a = lambda s: F(s)
    nodes = [
        Node("P1", "premise", a("0.05")),
        Node("P2", "premise", a("0.05")),
        Node("P3", "premise", a("0")),
        Node("w1", "premise", a("0"), is_rule=True),   # P1 -> P4
        Node("w2", "premise", a("0"), is_rule=True),   # P2 -> P5
        Node("w3", "premise", a("0"), is_rule=True),   # (P4^P5^P3) -> P6
        Node("w4", "premise", a("0"), is_rule=True),   # P6 -> P7
        Node("P4", "vamp", a("0.05"), premises=["P1"], imp="w1", cls="plain"),
        Node("P5", "vamp", a("0.05"), premises=["P2"], imp="w2", cls="plain"),
        Node("P6", "vamp", a("0.30"), premises=["P4", "P5", "P3"],
             imp="w3", cls="J1J2"),
        Node("P7", "vamp", a("0.55"), premises=["P6"], imp="w4",
             cls="J1J2", is_action=True),
    ]
    T = {"plain": a("0"), "J1J2": a("0.30")}
    theta = {"escalate": a("0.60")}
    return ProofObject(nodes), T, theta, {"escalate"}, a("0.5")


def worked_example_certified() -> tuple[ProofObject, dict, dict, set, Fraction]:
    a = lambda s: F(s)
    nodes = [
        Node("P1", "premise", a("0.05")),
        Node("P2", "premise", a("0.05")),
        Node("P3", "premise", a("0")),
        Node("w1", "premise", a("0"), is_rule=True),
        Node("w2", "premise", a("0"), is_rule=True),
        Node("w3", "premise", a("0"), is_rule=True),
        Node("P4", "vamp", a("0.05"), premises=["P1"], imp="w1", cls="plain"),
        Node("P5", "vamp", a("0.05"), premises=["P2"], imp="w2", cls="plain"),
        Node("P6", "vamp", a("0.20"), premises=["P4", "P5", "P3"],
             imp="w3", cls="J1J2", is_action=True),
    ]
    T = {"plain": a("0"), "J1J2": a("0.30")}
    theta = {"escalate": a("0.60")}
    return ProofObject(nodes), T, theta, {"escalate"}, a("0.5")


# ---------------------------------------------------------------------------
# Appendix C counterexample: naive DAG-aware budgets are unsound.
# q used twice through u1, u2; VAMP_2 merge at v. eps=0.1, t=0.2.
# ---------------------------------------------------------------------------

def counterexample_dag(eps: Fraction = F(1, 10), t: Fraction = F(1, 5)):
    nodes = [
        Node("q",  "premise", eps),
        Node("wq1", "premise", F(0), is_rule=True),
        Node("wq2", "premise", F(0), is_rule=True),
        Node("wv",  "premise", F(0), is_rule=True),
        Node("u1", "vamp", eps + t, premises=["q"], imp="wq1", cls="H"),
        Node("u2", "vamp", eps + t, premises=["q"], imp="wq2", cls="H"),
        Node("v",  "vamp", 2 * eps + 3 * t, premises=["u1", "u2"],
             imp="wv", cls="H"),
    ]
    T = {"H": t}
    pi = ProofObject(nodes)
    # Naive DAG-aware budget: each distinct ancestor charged once.
    naive = eps + 3 * t
    return pi, T, naive


if __name__ == "__main__":
    pi, T, theta, Rh, eta = worked_example_rejected()
    r = vamp_audit(pi, T, theta, Rh, eta)
    print("rejected-variant:", r.kappa, "| beta(P7) =", r.beta["P7"])
    pi, T, theta, Rh, eta = worked_example_certified()
    r = vamp_audit(pi, T, theta, Rh, eta)
    print("certified-variant:", r.kappa, "| beta(P6') =", r.beta["P6"])
    pi, T, naive = counterexample_dag()
    r = vamp_audit(pi, T, {"escalate": F(2)}, {"escalate"}, F(1))
    print("counterexample: alpha(v) =", pi.nodes["v"].alpha,
          "| B_tree(v) =", r.beta["v"], "| naive B_DAG(v) =", naive,
          "| unsound:", pi.nodes["v"].alpha > naive)
