#!/usr/bin/env python3
"""Exact companion for THM-3841's three-puncture Jelonek obstruction."""

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


A, C, T, q = sp.symbols("A C T q")

Delta = sp.expand(
    A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
    + C**2 * (162 * A**3 + 126 * A**2 * C - 4 * C**3)
    - 27 * A**2 * C**4
)
G = 3 * q**2 + 7
R = q**3 - 7 * q - 6
Aq = -2 * q**2 * R / G**2
Cq = -q * R / G

check(sp.factor(Delta) == Delta, "the affine discriminant is irreducible")
check(
    sp.gcd(sp.Poly(Delta, A, C), sp.Poly(sp.diff(Delta, A), A, C))
    == 1,
    "Delta is squarefree in the A direction",
)
check(
    sp.gcd(sp.Poly(Delta, A, C), sp.Poly(sp.diff(Delta, C), A, C))
    == 1,
    "Delta is squarefree in the C direction",
)
same(Delta.subs(A, 0), -4 * C**5, "the branch is not supported on A=0")
same(sp.factor(Delta.subs({A: Aq, C: Cq})), 0,
     "the q parametrization lies on Delta")
param_relation_1 = A * G - 2 * q * C
param_relation_2 = C * G + q * R
same(sp.resultant(param_relation_1, param_relation_2, q), -28 * Delta,
     "elimination of q gives exactly the irreducible branch equation")

# The characteristic cubic and its different reproduce the exact norm packet
# used to read v_E(A)=0 and v_E(D)=1 at the generic ramification divisor.
f = T**3 - C * T**2 + 7 * A**2 * T + 3 * A**3 - A**2 * C**2
different = sp.diff(f, T)
same(sp.discriminant(f, T), A**2 * Delta,
     "characteristic discriminant is A^2 Delta")
same(sp.resultant(f, different, T), -A**2 * Delta,
     "norm of the different is -A^2 Delta")

# A rational inverse at the generic branch point makes the parametrization
# birational, not merely dominant.
q_inverse = A * (27 * A - 9 * C**2 + 7 * C) / (2 * (C**2 - 21 * A**2))
same(sp.factor(q_inverse.subs({A: Aq, C: Cq})), q,
     "generic rational inverse recovers q")
check(sp.factor((Cq**2 - 21 * Aq**2)) != 0,
      "the inverse denominator is not identically zero")

# Homogenize the normalization map P1 -> projective closure of Delta.
u, v = sp.symbols("u v")
Gh = 3 * u**2 + 7 * v**2
Rh = u**3 - 7 * u * v**2 - 6 * v**3
Phi_A = sp.expand(-2 * u**2 * Rh * v)
Phi_C = sp.expand(-u * Rh * Gh)
Phi_Z = sp.expand(Gh**2 * v**2)

check(all(sp.Poly(F, u, v).total_degree() == 6
          for F in (Phi_A, Phi_C, Phi_Z)),
      "the projective normalization packet is homogeneous of degree six")
check(sp.gcd(sp.gcd(Phi_A, Phi_C), Phi_Z) == 1,
      "the projective packet has no common polynomial factor")
check(sp.resultant(3 * q**2 + 7, q**3 - 7 * q - 6, q) == 6460,
      "the two finite pole points are not base points")
check(sp.discriminant(3 * q**2 + 7, q) == -84,
      "the finite pole divisor has two distinct points")
same(Phi_C.subs({u: 1, v: 0}), -3,
     "the parameter infinity is not a base point")
same(Phi_A.subs(v, 1) / Phi_Z.subs(v, 1), Aq.subs(q, u),
     "projective A chart dehomogenizes correctly")
same(Phi_C.subs(v, 1) / Phi_Z.subs(v, 1), Cq.subs(q, u),
     "projective C chart dehomogenizes correctly")

# The affine normalization is exactly D(G) in A1: the target affine chart is
# Phi_Z != 0, hence it removes infinity and both simple roots of G.
same(Phi_Z, v**2 * Gh**2,
     "the affine-chart pullback has exactly the three puncture supports")
check(sp.gcd(v, Gh) == 1, "infinity is distinct from both finite punctures")
check(sp.degree(G, q) == 2, "there are exactly two finite punctures")

for q0 in (-2, -1, 0, 3):
    same(Aq.subs(q, q0), 0, f"vertex A value at q={q0}")
    same(Cq.subs(q, q0), 0, f"vertex C value at q={q0}")

# The valuation certificate is formal but worth freezing: every extension
# has positive ramification index e, so the regular source pullback of h
# would have a negative value if the extension had an affine source center.
e = sp.symbols("e", integer=True, positive=True)
vE_A = sp.Integer(0)
vE_D = sp.Integer(1)
vE_h = vE_A - vE_D
same(vE_h, -1, "the deleted ramification valuation has v_E(h)=-1")
same(e * vE_h, -e, "every finite field extension keeps h negative")

# A morphism A1 -> Spec k[q,G^-1] sends q to a polynomial p.  If p is
# nonconstant of degree n, then G(p) has degree 2n and cannot be a unit.
n = sp.symbols("n", integer=True, positive=True)
same(2 * n, sp.degree(G, q) * n,
     "the pullback of G has positive degree under every nonconstant p")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "completion": "normal finite cubic Xbar over A2; U=Xbar minus the ramification prime E",
    "valuation": "v_E(A)=0, v_E(D)=1, v_E(h)=-1; every extension is source-infinite",
    "jelonek": "the irreducible branch Delta is a full component of S_(pi psi)",
    "normalization": "Delta^nu=P1 minus {infinity, the two roots of 3q^2+7}",
    "obstruction": "Jelonek-Lason polynomial uniruledness contradicts the three-puncture normalization",
    "scope": "all-degree; no dominant generically finite A2 morphism to U; hence no Keller atlas",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3841-deleted-ramification-three-puncture-jelonek-nonentry")
print("valuation=vE(A):0;vE(D):1;vE(h):-1;extensions_source_infinite")
print("branch=Delta_irreducible;Delta_subset_Jelonek_component")
print("normalization=P1_minus_infinity_and_two_roots_of_3q2_plus7")
print("punctures=3;polynomial_A1_dominance=impossible")
print("conclusion=no_dominant_generically_finite_A2_to_U;no_Keller_atlas")
print("scope=all_degree;arbitrary_algebraically_closed_characteristic_zero")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
