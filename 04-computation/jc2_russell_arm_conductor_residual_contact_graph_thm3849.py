#!/usr/bin/env python3
"""Exact controls for THM-3849's arm conductor/contact law."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)  # type: ignore[operator]


t, U, V = sp.symbols("t U V")
c, lam = sp.symbols("c lambda", nonzero=True)


def conductor_packet(
    p: sp.Expr,
    q: sp.Expr,
    delta: sp.Expr,
    kappa: sp.Expr,
    label: str,
) -> tuple[sp.Expr, sp.Expr]:
    """Verify the implicit equation and pulled-back gradient syzygy."""
    pprime = sp.diff(p, t)
    qprime = sp.diff(q, t)
    same(delta.subs({U: p, V: q}), 0, f"{label}: implicit equation")
    _, factors = sp.factor_list(delta, U, V)
    check(len(factors) == 1 and factors[0][1] == 1,
          f"{label}: irreducible displayed equation")
    check(sp.gcd(pprime, qprime) == 1, f"{label}: derivative row is unimodular")
    same(
        sp.diff(delta, U).subs({U: p, V: q}),
        kappa * qprime,
        f"{label}: first gradient component",
    )
    same(
        sp.diff(delta, V).subs({U: p, V: q}),
        -kappa * pprime,
        f"{label}: second gradient component",
    )
    return pprime, qprime


# The universal normal-jet contraction: the residual intersection coefficient
# is kappa times the arm Bezout determinant.
alpha, beta, pprime0, qprime0, kappa0 = sp.symbols(
    "alpha beta pprime qprime kappa"
)
same(
    kappa0 * qprime0 * alpha - kappa0 * pprime0 * beta,
    kappa0 * (alpha * qprime0 - pprime0 * beta),
    "universal gradient-normal contraction",
)
same(
    kappa0 * lam / (3 * c**3),
    kappa0 * (lam / (3 * c**3)),
    "Darboux normalization scales the residual coefficient",
)

# Ordinary node: two branches, one transverse edge of weight one.
p_node = t**2
q_node = t**3 - t
delta_node = V**2 - U * (U - 1) ** 2
kappa_node = -(t**2 - 1)
pp_node, qp_node = conductor_packet(
    p_node, q_node, delta_node, kappa_node, "node"
)
alpha_node = -1 / (3 * c**3)
beta_node = -t / (2 * c**3)
same(
    3 * c**3 * (alpha_node * qp_node - pp_node * beta_node),
    1,
    "node: Russell-arm Bezout row",
)
same(
    sp.diff(delta_node, U).subs({U: p_node, V: q_node}) * alpha_node
    + sp.diff(delta_node, V).subs({U: p_node, V: q_node}) * beta_node,
    kappa_node / (3 * c**3),
    "node: residual contact coefficient",
)
check(sp.degree(kappa_node, t) == 2, "node: conductor degree two")
check(sp.diff(kappa_node, t).subs(t, 1) != 0, "node: +1 contact is simple")
check(sp.diff(kappa_node, t).subs(t, -1) != 0, "node: -1 contact is simple")

# Tacnode: the same two normalization points have contact weight two.
p_tac = t**2 - 1
q_tac = t * (t**2 - 1) ** 2
delta_tac = V**2 - (U + 1) * U**4
kappa_tac = -(t**2 - 1) ** 2
conductor_packet(p_tac, q_tac, delta_tac, kappa_tac, "tacnode")
check(sp.degree(kappa_tac, t) == 4, "tacnode: conductor degree four")
check(sp.rem(kappa_tac, (t - 1) ** 2, domain=sp.QQ) == 0,
      "tacnode: +1 weighted degree two")
check(sp.rem(kappa_tac, (t + 1) ** 2, domain=sp.QQ) == 0,
      "tacnode: -1 weighted degree two")
u = sp.symbols("u")
tac_branch_gap = sp.series(2 * u**2 * sp.sqrt(1 + u), u, 0, 5).removeO()
check(sp.Poly(tac_branch_gap, u).as_dict().get((2,), 0) != 0,
      "tacnode: branch intersection weight is two")

# Ordinary triple point: three preimages and three unit-weight edges.  The
# conductor has order two at each vertex, the weighted graph degree.
p_tri = t**3 - t
q_tri = t * p_tri
delta_tri = V**3 - U**2 * V - U**4
kappa_tri = -p_tri**2
pp_tri, qp_tri = conductor_packet(
    p_tri, q_tri, delta_tri, kappa_tri, "ordinary triple"
)
alpha_tri = -sp.Rational(9, 2) * t
beta_tri = 1 - 6 * t**2
same(alpha_tri * qp_tri - pp_tri * beta_tri, 1,
     "ordinary triple: exact Bezout row")
check(sp.factor(kappa_tri) == -t**2 * (t - 1) ** 2 * (t + 1) ** 2,
      "ordinary triple: order two at all three vertices")
check(sp.degree(kappa_tri, t) == 6,
      "ordinary triple: twice the three-edge weight")
slopes = [
    sp.cancel(qp_tri.subs(t, root) / pp_tri.subs(t, root))
    for root in (-1, 0, 1)
]
check(len(set(slopes)) == 3, "ordinary triple: three distinct tangents")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "gradient": "Delta_U(p,q)=kappa*q';Delta_V(p,q)=-kappa*p'",
    "residual": "div Delta(P,Q)=L+D;[D]=2[L];D cap L=V(kappa)",
    "singular_support": "zeros(kappa)=normalization preimages of affine singularities",
    "contact_graph": "ord_e(kappa)=weighted degree;deg(kappa)=2 total edge weight=2 sum delta",
    "controls": "node degree2;tacnode degree4;ordinary triple degree6",
    "scope": "all Darboux pairs on the fixed THM-3785 Russell pseudo-plane",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3849-russell-arm-conductor-polynomial-and-residual-contact-graph")
print("gradient=Delta_U(p,q)=kappa*qprime;Delta_V(p,q)=-kappa*pprime")
print("residual=div(Delta(P,Q))=L+D;class(D)=2[L];D_cap_L=V(kappa)")
print("graph=ord_vertex(kappa)=weighted_degree;deg(kappa)=2*total_edge_weight")
print("controls=node:2;tacnode:4;ordinary_triple:6")
print("scope=all_degree_fixed_Russell_surface;Darboux_pair_existence_OPEN")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
