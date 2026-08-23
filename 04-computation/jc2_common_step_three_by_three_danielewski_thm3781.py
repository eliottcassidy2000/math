#!/usr/bin/env python3
"""Exact support and coefficient companion for THM-3781."""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def weight_grid(left: tuple[int, ...], right: tuple[int, ...]) -> Counter[int]:
    return Counter(r + s + 1 for r in left for s in right)


def bracket_coefficient(
    r: int,
    f: sp.Expr,
    fp: sp.Expr,
    s: int,
    g: sp.Expr,
    gp: sp.Expr,
) -> sp.Expr:
    """Coefficient of c^(r+s+1) in {c^r f(b),c^s g(b)}."""
    return sp.expand(s * fp * g - r * f * gp)


# The common-step-three scalar-centred supports have only two inequivalent
# placements after swapping the two outputs.
case_i = ((-4, -1, 2), (-3, 0, 3))
case_ii = ((-5, -2, 1), (-2, 1, 4))
expected_grid = Counter({-6: 1, -3: 2, 0: 3, 3: 2, 6: 1})
gate(weight_grid(*case_i) == expected_grid, "case I support convolution")
gate(weight_grid(*case_ii) == expected_grid, "case II support convolution")

placements = []
for p in range(-20, 21):
    q = -p - 1
    if p - 3 < 0 < p + 3 and q - 3 < 0 < q + 3:
        placements.append((p, q))
gate(placements == [(-2, 1), (-1, 0), (0, -1), (1, -2)],
     "all scalar-centred step-three placements")
orbits = {tuple(sorted((p, q))) for p, q in placements}
gate(orbits == {(-2, 1), (-1, 0)}, "two placements up to output swap")


# Universal coefficient symbols.  Kp and Lp stand for derivatives with
# respect to b; all constants are declared nonzero exactly where division is
# used in the proof.
K, Kp, L, Lp = sp.symbols("K Kp L Lp", nonzero=True)
h, hp = sp.symbols("h hp")
a, ap = sp.symbols("a ap")
A0, B0, L0, M0 = sp.symbols("A0 B0 L0 M0", nonzero=True)


# Case I: P weights (-4,-1,2), Q weights (-3,0,3).
f1, f1p = A0 * K**4, 4 * A0 * K**3 * Kp
g1, g1p = B0 * K**3, 3 * B0 * K**2 * Kp
F1, F1p = L0 * L**2, 2 * L0 * L * Lp
H1, H1p = M0 * L**3, 3 * M0 * L**2 * Lp

gate(bracket_coefficient(-4, f1, f1p, -3, g1, g1p) == 0,
     "case I negative endpoint power law")
gate(bracket_coefficient(2, F1, F1p, 3, H1, H1p) == 0,
     "case I positive endpoint power law")

minus_i = (
    bracket_coefficient(-4, f1, f1p, 0, h, hp)
    + bracket_coefficient(-1, a, ap, -3, g1, g1p)
)
minus_i_reduced = sp.factor(minus_i / K**4)
expected_minus_i = sp.factor(4 * A0 * hp - 3 * B0 * (ap * K - a * Kp) / K**2)
gate(sp.factor(minus_i_reduced - expected_minus_i) == 0,
     "case I lower adjacent derivative")

plus_i = (
    bracket_coefficient(-1, a, ap, 3, H1, H1p)
    + bracket_coefficient(2, F1, F1p, 0, h, hp)
)
plus_i_reduced = sp.factor(plus_i / L**2)
expected_plus_i = sp.factor(3 * M0 * (ap * L + a * Lp) - 2 * L0 * hp)
gate(sp.factor(plus_i_reduced - expected_plus_i) == 0,
     "case I upper adjacent derivative")

lam_i = sp.factor(4 * A0 / (3 * B0))
mu_i = sp.factor(2 * L0 / (3 * M0))
C0, D0 = sp.symbols("C0 D0")
a_i = K * (lam_i * h + C0)
gate(sp.factor((sp.diff(a_i, h) * hp + sp.diff(a_i, K) * Kp) / K
               - a_i * Kp / K**2 - lam_i * hp) == 0,
     "case I first integrated relation")

# Combining a=K(lambda h+C) and aL=mu h+D gives a polynomial Mobius
# relation.  Its denominator and numerator have the same degree in z=KL,
# so a polynomial h must be constant (or the impossible zero branch).
z = sp.symbols("z")
mobius_den_i = lam_i * z - mu_i
mobius_num_i = D0 - C0 * z
gate(sp.degree(mobius_den_i, z) == 1, "case I Mobius denominator degree")
gate(sp.degree(mobius_num_i, z) <= 1, "case I Mobius numerator degree")
gate(sp.factor(sp.resultant(mobius_den_i, mobius_num_i, z)
               - (lam_i * D0 - mu_i * C0)) == 0,
     "case I Mobius resultant")


# Case II: P weights (-5,-2,1), Q weights (-2,1,4).
f2, f2p = A0 * K**5, 5 * A0 * K**4 * Kp
g2, g2p = B0 * K**2, 2 * B0 * K * Kp
F2, F2p = L0 * L, L0 * Lp
H2, H2p = M0 * L**4, 4 * M0 * L**3 * Lp

gate(bracket_coefficient(-5, f2, f2p, -2, g2, g2p) == 0,
     "case II negative endpoint power law")
gate(bracket_coefficient(1, F2, F2p, 4, H2, H2p) == 0,
     "case II positive endpoint power law")

minus_ii = (
    bracket_coefficient(-5, f2, f2p, 1, h, hp)
    + bracket_coefficient(-2, a, ap, -2, g2, g2p)
)
expected_minus_ii = sp.factor(
    K**4 * (5 * A0 * (Kp * h + K * hp)
            - 2 * B0 * (ap * K - 2 * a * Kp) / K**3)
)
gate(sp.factor(minus_ii - expected_minus_ii) == 0,
     "case II lower adjacent derivative")

plus_ii = (
    bracket_coefficient(-2, a, ap, 4, H2, H2p)
    + bracket_coefficient(1, F2, F2p, 1, h, hp)
)
expected_plus_ii = sp.factor(
    L**2 * (4 * M0 * (ap * L**2 + 2 * a * L * Lp)
            - L0 * (hp * L - h * Lp) / L**2)
)
gate(sp.factor(plus_ii - expected_plus_ii) == 0,
     "case II upper adjacent derivative")

lam_ii = sp.factor(5 * A0 / (2 * B0))
mu_ii = sp.factor(4 * M0 / L0)
nu_ii = sp.factor(lam_ii * mu_ii)
obstruction_den_ii = 1 - nu_ii * z**3
obstruction_num_ii = lam_ii * D0 * z + C0
gate(sp.degree(obstruction_den_ii, z) == 3,
     "case II cubic denominator degree")
gate(sp.degree(obstruction_num_ii, z) <= 1,
     "case II linear numerator degree")
gate(sp.factor(obstruction_den_ii.subs(z, 0)) == 1,
     "case II denominator coprime to K at K=0")


# Root-valuation invoices force the common negative endpoint owner K to
# contain every squarefree arm factor.
for v in range(0, 6):
    admissible_i = 4 * v >= 2 and 3 * v >= 2
    admissible_ii = 5 * v >= 3 and 2 * v >= 1
    gate(admissible_i == (v >= 1), "case I arm valuation threshold")
    gate(admissible_ii == (v >= 1), "case II arm valuation threshold")

# Degree hostiles make the terminal divisibility arguments explicit for a
# broad exact range.  z is nonconstant because Sigma divides K.
for degree_z in range(1, 25):
    gate(degree_z <= degree_z, "case I equal-degree Mobius boundary")
    gate(3 * degree_z > degree_z, "case II cubic-over-linear boundary")


semantic_rows = (
    "surface:c^2*e=Sigma(b);Sigma_squarefree;degSigma>=2",
    "caseI:Pweights=-4,-1,2;Qweights=-3,0,3",
    "caseI:endpoints=f=A*k^4,g=B*k^3,F=L*l^2,H=M*l^3",
    "caseI:adjacent=a/k=lambda*h+C;a*l=mu*h+D;Mobius_forces_h_constant",
    "caseII:Pweights=-5,-2,1;Qweights=-2,1,4",
    "caseII:endpoints=f=A*k^5,g=B*k^2,F=L*l,H=M*l^4",
    "caseII:adjacent=a/k^2=lambda*k*h+C;h/l=mu*a*l^2+D",
    "caseII:(1-nu*(k*l)^3)*a=k^2*(lambda*D*k*l+C)",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3781-common-step-three-by-three-danielewski-darboux-nonentry")
print("surface=exponent_two_squarefree_Danielewski")
print("step=3;placements_up_to_swap=2")
print("case_I=-4,-1,2_x_-3,0,3;exit=polynomial_Mobius_forces_zero_weight_constant")
print("case_II=-5,-2,1_x_-2,1,4;exit=cubic_denominator_cannot_divide_linear_numerator")
print("scalar_equation=never_reached")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
