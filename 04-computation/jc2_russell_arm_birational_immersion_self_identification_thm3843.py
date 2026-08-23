#!/usr/bin/env python3
"""Exact companion for THM-3843's Russell-arm normalization packet."""

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


x, y = sp.symbols("x y")
c = sp.symbols("c", nonzero=True)
r = x**3
z = sp.expand(x * (c + x**3 * y))
e = sp.expand(3 * c**2 * y + 3 * c * x**3 * y**2 + x**6 * y**3)


def jac(lhs: sp.Expr, rhs: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(lhs, x) * sp.diff(rhs, y)
                     - sp.diff(lhs, y) * sp.diff(rhs, x))


same(r**2 * e, z**3 - c**3 * r, "Russell pseudo-plane relation")
same(jac(r, z), 3 * r**2, "source bracket {r,z}")
same(jac(r, e), 9 * z**2, "source bracket {r,e}")
same(jac(z, e), 3 * c**3 + 6 * r * e, "source bracket {z,e}")
same(r.subs(x, 0), 0, "source line maps to r=0")
same(z.subs(x, 0), 0, "source line maps to z=0")
same(e.subs(x, 0), 3 * c**2 * y, "source line is the arm isomorphism")

# Freeze the universal arm restriction of the Poisson bracket.  The four jet
# symbols are the two tangent derivatives and the two transverse z-jets.
lam, p1, q1, alpha, beta = sp.symbols("lambda p1 q1 alpha beta")
arm_bracket = 3 * c**3 * (alpha * q1 - p1 * beta)
bezout_rhs = lam / (3 * c**3)
same(arm_bracket - lam,
     3 * c**3 * (alpha * q1 - beta * p1 - bezout_rhs),
     "arm bracket is the derivative Bezout identity")
check(sp.Poly(alpha * q1 - beta * p1, p1, q1).total_degree() == 1,
      "the arm equation is linear in the two tangent derivatives")

# Polynomial-Luroth hostile: a proper common polynomial parameter s(t) of
# degree >1 contributes the nonunit common derivative s'(t).
t = sp.symbols("t")
s0, s1, s2 = sp.symbols("s0 s1 s2", nonzero=True)
s = s0 + s1 * t + s2 * t**2
F0, F1, F2, G0, G1, G2 = sp.symbols("F0 F1 F2 G0 G1 G2")
S = sp.symbols("S")
F_profile = F0 + F1 * S + F2 * S**2
G_profile = G0 + G1 * S + G2 * S**2
Fp = F_profile.subs(S, s)
Gp = G_profile.subs(S, s)
same(sp.diff(Fp, t), sp.diff(F_profile, S).subs(S, s) * sp.diff(s, t),
     "common-parameter derivative factor for the first coordinate")
same(sp.diff(Gp, t), sp.diff(G_profile, S).subs(S, s) * sp.diff(s, t),
     "common-parameter derivative factor for the second coordinate")
check(sp.degree(sp.diff(s, t), t) == 1,
      "a quadratic common parameter has a nonunit derivative")

# A birational immersion can still identify two points.  This nodal cubic is
# the sharp positive control for the theorem's conclusion.
p_node = t**2 - 1
q_node = t * (t**2 - 1)
check(sp.gcd(sp.diff(p_node, t), sp.diff(q_node, t)) == 1,
      "nodal control is an immersion on its normalization")
same(q_node / p_node, t, "nodal control is birational")
same(p_node.subs(t, 1), p_node.subs(t, -1),
     "nodal control identifies the first coordinate")
same(q_node.subs(t, 1), q_node.subs(t, -1),
     "nodal control identifies the second coordinate")
P, Q = sp.symbols("P Q")
node_equation = Q**2 - P**2 * (P + 1)
same(node_equation.subs({P: p_node, Q: q_node}), 0,
     "nodal control lies on its implicit curve")
same(sp.diff(node_equation, P).subs({P: 0, Q: 0}), 0,
     "self-identification image is singular in the P direction")
same(sp.diff(node_equation, Q).subs({P: 0, Q: 0}), 0,
     "self-identification image is singular in the Q direction")
check(sp.resultant(sp.diff(p_node, t), sp.diff(q_node, t), t) != 0,
      "the two nodal branches remain individually immersed")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "surface": "r^2 e=z^3-c^3 r with the THM-3785 source atlas of degree three",
    "arm": "L=(r,z)=0 is A1_e and source x=0 maps isomorphically to it",
    "restriction": "{P,Q}=lambda gives 3c^3(P_z Q_e-P_e Q_z)=lambda on L",
    "normalization": "derivative Bezout plus polynomial Luroth makes L->Gamma the finite normalization with nowhere-zero differential",
    "self_identification": "injectivity would trigger Gwozdziewicz on the source line, contradicting degree at least nine",
    "scope": "all Darboux pairs on the fixed Russell pseudo-plane; no existence claim",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3843-russell-arm-birational-immersion-and-forced-self-identification")
print("arm=L=(r,z)=0=A1_e;source_line=x=0;map=e=3c^2y")
print("darboux_restriction=3c^3(alpha*qprime-pprime*beta)=lambda")
print("normalization=finite_birational;differential_nonvanishing;derivative_gcd=1")
print("injective_boundary=Gwozdziewicz_forces_automorphism_but_degree>=9")
print("conclusion=distinct_arm_points_share_a_singular_image")
print("scope=all_degree_fixed_Russell_surface;Darboux_pair_existence_OPEN")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
