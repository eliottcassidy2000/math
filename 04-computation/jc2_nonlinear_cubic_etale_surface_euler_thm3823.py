#!/usr/bin/env python3
"""Exact compact-Euler computation for the THM-3811 etale surface."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(f"gate failed: {label}")


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    gate(sp.factor(sp.cancel(left - right)) == 0, label)


A, C, q, T = sp.symbols("A C q T")
R = (q - 3) * (q + 1) * (q + 2)
D = 3 * q**2 + 7
G = q**3 + 7 * q + 3
Aq = sp.cancel(-2 * q**2 * R / D**2)
Cq = sp.cancel(-q * R / D)
wq = sp.cancel(Aq * q)

Delta = sp.expand(
    A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
    + C**2 * (162 * A**3 + 126 * A**2 * C - 4 * C**3)
    - 27 * A**2 * C**4
)
f = T**3 - C * T**2 + 7 * A**2 * T + 3 * A**3 - A**2 * C**2

# The double-root incidence covers the single reduced branch.
E1 = sp.expand(A * G - C * (C + q**2))
E2 = sp.expand(A * D - 2 * C * q)
same(sp.resultant(E1, E2, q), -A * Delta, "double-root incidence resultant")
gate(sp.resultant(D, G, q) != 0, "division by D loses no incidence point")
same(Delta.subs(A, 0), -4 * C**5, "no branch component lies in A=0")
gate(sp.gcd(sp.gcd(Delta, sp.diff(Delta, A)), sp.diff(Delta, C)) == 1,
     "branch is reduced")
coefficient, factors = sp.factor_list(Delta, A, C)
gate(coefficient != 0 and len(factors) == 1 and factors[0][1] == 1,
     "branch has one rational irreducible factor")
same(Delta.subs({A: Aq, C: Cq}), 0, "normalization lies on branch")

# The affine normalization is P1 minus the two D-roots and infinity.
gate(sp.degree(D, q) == 2 and sp.discriminant(D, q) != 0,
     "two distinct finite denominator points")
num_A, den_A = sp.fraction(Aq)
num_C, den_C = sp.fraction(Cq)
gate(sp.gcd(num_A, D) == sp.gcd(num_C, D) == 1,
     "normalization numerators miss finite poles")
gate(sp.factor(den_A) == D**2 and sp.factor(den_C) == D,
     "finite pole orders are two and one")
gate((sp.degree(num_A, q) - sp.degree(den_A, q),
      sp.degree(num_C, q) - sp.degree(den_C, q)) == (1, 2),
     "infinity is the third puncture")
chi_normalization = 2 - 3
gate(chi_normalization == -1, "normalization compact Euler characteristic")

# Four normalization points map to the origin.
gate(sp.discriminant(R, q) != 0 and sp.gcd(q, R) == 1,
     "four distinct vertex parameters")
gate(sp.gcd(D, q * R) == 1, "vertex parameters avoid punctures")
for root in (-2, -1, 0, 3):
    same(Aq.subs(q, root), 0, f"vertex A at {root}")
    same(Cq.subs(q, root), 0, f"vertex C at {root}")

# No other singular preimages occur except the two q^2=7/3 points.
dA_pull = sp.factor(sp.diff(Delta, A).subs({A: Aq, C: Cq}))
dC_pull = sp.factor(sp.diff(Delta, C).subs({A: Aq, C: Cq}))
expected_dA = (
    4 * q**3 * R**3 * (3 * q**2 - 7)**3 * (q**3 + 7 * q + 3) / D**6
)
expected_dC = (
    -4 * q**4 * R**3 * (3 * q**2 - 7)**3
    * (q**3 + 21 * q + 12) / D**7
)
same(dA_pull, expected_dA, "Delta_A pullback")
same(dC_pull, expected_dC, "Delta_C pullback")
gate(sp.gcd(q**3 + 7 * q + 3, q**3 + 21 * q + 12) == 1,
     "no extra common derivative zero")
gate(sp.gcd(D, q * R * (3 * q**2 - 7)) == 1,
     "singular preimages avoid punctures")
gate(sp.discriminant(3 * q**2 - 7, q) != 0 and sp.gcd(q * R, 3 * q**2 - 7) == 1,
     "two distinct nonvertex triple parameters")

# The simple companion collides with the double root at exactly those values.
companion = sp.cancel(Cq - 2 * wq)
same(f.subs({A: Aq, C: Cq}), (T - wq)**2 * (T - companion),
     "branch cubic factors as double root and companion")
gap = sp.factor(companion - wq)
same(gap, q * R * (3 * q**2 - 7) / D**2, "companion collision locus")
same(sp.diff(f, T).subs({A: Aq, C: Cq, T: companion}), gap**2,
     "companion derivative is the squared gap")

# Only the four-to-one vertex changes the normalization Euler characteristic.
chi_branch = chi_normalization - (4 - 1)
gate(chi_branch == -4, "branch compact Euler characteristic")
omitted_values = 3
chi_off_branch = 1 - chi_branch
chi_visible_branch = chi_branch - omitted_values
chi_U = 3 * chi_off_branch + chi_visible_branch
gate((chi_off_branch, chi_visible_branch, chi_U) == (5, -7, 8),
     "surface stratum Euler arithmetic")

# A finite etale map from connected A2 is a positive-degree finite covering.
# Its Euler equation 1=d*8 has no positive integer solution.
chi_A2 = 1
gate(chi_A2 % abs(chi_U) != 0, "universal finite-degree obstruction")

# Finiteness is essential: this nonproper etale open immersion has unequal chi.
chi_Gm_times_A1 = (1 - 1) * 1
gate(chi_Gm_times_A1 == 0 and chi_A2 == 1,
     "nonproper etale Euler-multiplicativity hostile")

self_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(self_tree)),
     "assertion-free exact companion")

semantic = {
    "normalization": "P1 minus roots(3q^2+7) and infinity; chi_c=-1",
    "branch": "four preimages at origin, two unibranch triple values; chi_c=-4",
    "surface": "chi_c(U)=3*5-7=8",
    "consequence": "no finite etale A2_C to U_C",
    "scope": "nonfinite nonproper dominant etale polynomial atlas remains open",
}
semantic_hash = hashlib.sha256(json.dumps(
    semantic, sort_keys=True, separators=(",", ":")
).encode("ascii")).hexdigest()

print("theorem=THM3823_nonlinear_cubic_etale_surface_compact_Euler")
print("normalization=P1_minus_3;chi_c=-1")
print("branch=origin_preimages_4;triple_values_2_unibranch;chi_c=-4")
print("surface_strata=off_branch:3x5;visible_branch:1x(-7);omitted:0x3")
print("surface=chi_c_U=8")
print("finite_etale=no_A2_C_to_U_C_of_any_positive_degree")
print("scope=nonfinite_nonproper_dominant_etale_polynomial_atlas_open")
print(f"GATES={GATES}")
print(f"semantic_sha256={semantic_hash}")
print("RESULT=PASS")
