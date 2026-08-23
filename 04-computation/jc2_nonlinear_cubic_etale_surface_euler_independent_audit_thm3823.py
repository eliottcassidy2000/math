#!/usr/bin/env python3
"""Independent hostile audit of the THM-3823 compact-Euler obstruction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    ROOT / "04-computation/jc2_nonlinear_cubic_etale_surface_euler_thm3823.py":
        "11922872c37d824cd22da0875bca48044afe7267f92ad2c67be3fd5d4088c815",
    ROOT / "05-knowledge/results/jc2_nonlinear_cubic_etale_surface_euler_thm3823.out":
        "1007441d08053e669ca299ca6022c6b2462682ef3fa54e1ebe08e3850fb64b30",
    ROOT / "01-canon/theorems/THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet.md":
        "e367b22565fa5bca53a48d798250f72e56cdb11c4c260381ab45b2f612d2590c",
}

GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(f"audit gate failed: {label}")


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    gate(sp.factor(sp.cancel(left - right)) == 0, label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


for input_path, expected_hash in INPUTS.items():
    gate(input_path.is_file(), f"input exists: {input_path.name}")
    gate(sha256(input_path) == expected_hash, f"input hash: {input_path.name}")


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

# Independent branch coverage and rational inverse.
E1 = sp.expand(A * G - C * (C + q**2))
E2 = sp.expand(A * D - 2 * C * q)
same(sp.resultant(E1, E2, q), -A * Delta, "incidence resultant")
gate(sp.resultant(D, G, q) != 0, "D division loses no incidence point")
same(Delta.subs(A, 0), -4 * C**5, "no A=0 branch component")
gate(sp.gcd(sp.gcd(Delta, sp.diff(Delta, A)), sp.diff(Delta, C)) == 1,
     "branch reduced")
coefficient, factors = sp.factor_list(Delta, A, C)
gate(coefficient != 0 and len(factors) == 1 and factors[0][1] == 1,
     "branch irreducible over Q")
same(Delta.subs({A: Aq, C: Cq}), 0, "parametrization on branch")
double_root = A**2 * (27 * A - 9 * C**2 + 7 * C) / (2 * (C**2 - 21 * A**2))
same(double_root.subs({A: Aq, C: Cq}), wq,
     "generic inverse recovers Aq")

# Exactly three normalization punctures.
gate(sp.degree(D, q) == 2 and sp.discriminant(D, q) != 0,
     "two finite punctures")
num_A, den_A = sp.fraction(Aq)
num_C, den_C = sp.fraction(Cq)
gate(sp.gcd(num_A, D) == sp.gcd(num_C, D) == 1,
     "finite poles do not cancel")
gate(sp.factor(den_A) == D**2 and sp.factor(den_C) == D,
     "finite pole orders")
gate((sp.degree(num_A, q) - sp.degree(den_A, q),
      sp.degree(num_C, q) - sp.degree(den_C, q)) == (1, 2),
     "infinity pole orders")
chi_normalization = 2 - 3
gate(chi_normalization == -1, "normalization chi_c")

# The vertex has four genuine branches with distinct tangent directions.
gate(sp.discriminant(R, q) != 0 and sp.gcd(q, R) == 1,
     "four distinct vertex parameters")
gate(sp.gcd(D, q * R) == 1, "vertex parameters affine")
for root in (-2, -1, 0, 3):
    same(Aq.subs(q, root), 0, f"vertex A at {root}")
    same(Cq.subs(q, root), 0, f"vertex C at {root}")
ratio_CA = sp.cancel(Cq / Aq)
gate(tuple(sp.limit(ratio_CA, q, root) for root in (-2, -1, 0, 3)) ==
     (-sp.Rational(19, 4), -5, sp.oo, sp.Rational(17, 3)),
     "four tangent directions")

# Complete affine singular-preimage census.
dA_pull = sp.factor(sp.diff(Delta, A).subs({A: Aq, C: Cq}))
dC_pull = sp.factor(sp.diff(Delta, C).subs({A: Aq, C: Cq}))
same(dA_pull,
     4 * q**3 * R**3 * (3 * q**2 - 7)**3 * (q**3 + 7 * q + 3) / D**6,
     "Delta_A pullback")
same(dC_pull,
     -4 * q**4 * R**3 * (3 * q**2 - 7)**3
     * (q**3 + 21 * q + 12) / D**7,
     "Delta_C pullback")
gate(sp.gcd(q**3 + 7 * q + 3, q**3 + 21 * q + 12) == 1,
     "terminal cubics coprime")
num_dA = sp.fraction(sp.cancel(dA_pull))[0]
num_dC = sp.fraction(sp.cancel(dC_pull))[0]
common_zero = sp.Poly(sp.gcd(num_dA, num_dC), q).monic().as_expr()
expected_zero = sp.Poly(q**3 * R**3 * (3 * q**2 - 7)**3, q).monic().as_expr()
same(common_zero, expected_zero, "complete common derivative divisor")
gate(sp.gcd(D, q * R * (3 * q**2 - 7)) == 1,
     "singular preimages avoid punctures")

triple_points = []
for sign in (-1, 1):
    sigma = sign * sp.sqrt(21)
    triple_q = sigma / 3
    triple_A = sp.Rational(1, 7) + sigma / 27
    triple_C = sp.Rational(7, 9) + sigma / 7
    same(Aq.subs(q, triple_q), triple_A, f"triple A sign {sign}")
    same(Cq.subs(q, triple_q), triple_C, f"triple C sign {sign}")
    gate(triple_A != 0 and triple_C != 0, f"triple point nonzero sign {sign}")
    same((Cq / (3 * Aq)).subs(q, triple_q), triple_q,
         f"triple point has unique normalization preimage sign {sign}")
    triple_points.append((triple_A, triple_C))
gate(any(sp.simplify(x - y) != 0 for x, y in zip(*triple_points)),
     "triple images distinct")

# Arithmetic controls for normalization identification at good reductions.
def delta_mod(a_value: int, c_value: int, prime: int) -> int:
    return int(Delta.subs({A: a_value, C: c_value})) % prime


good_counts = []
for prime in (11, 13, 23, 29, 31):
    actual = sum(delta_mod(a_value, c_value, prime) == 0
                 for a_value in range(prime) for c_value in range(prime))
    d_roots = sum((3 * value * value + 7) % prime == 0 for value in range(prime))
    predicted = prime - d_roots - 3
    gate(actual == predicted, f"normalization point count mod {prime}")
    good_counts.append((prime, actual, d_roots))

bad_counts = []
for prime in (17, 19):
    actual = sum(delta_mod(a_value, c_value, prime) == 0
                 for a_value in range(prime) for c_value in range(prime))
    d_roots = sum((3 * value * value + 7) % prime == 0 for value in range(prime))
    predicted = prime - d_roots - 3
    gate(actual != predicted, f"bad-reduction hostile mod {prime}")
    bad_counts.append((prime, actual, predicted))

# Branch fibres: companion survives except at the three singular values.
companion = sp.cancel(Cq - 2 * wq)
same(f.subs({A: Aq, C: Cq}), (T - wq)**2 * (T - companion),
     "double-root plus companion factorization")
gap = sp.factor(companion - wq)
same(gap, q * R * (3 * q**2 - 7) / D**2, "companion collision divisor")
same(sp.diff(f, T).subs({A: Aq, C: Cq, T: companion}), gap**2,
     "companion derivative square")

# Euler integration and universal finite-degree obstruction.
chi_branch = chi_normalization - (4 - 1)
gate(chi_branch == -4, "branch chi_c")
chi_A2 = 1
chi_off_branch = chi_A2 - chi_branch
chi_visible_branch = chi_branch - 3
chi_U = 3 * chi_off_branch + chi_visible_branch
gate((chi_off_branch, chi_visible_branch, chi_U) == (5, -7, 8),
     "surface Euler integration")
gate(chi_A2 % abs(chi_U) != 0,
     "all finite positive covering degrees excluded")

# Nonproper etale maps do not obey finite-cover Euler multiplicativity.
chi_Gm_times_A1 = (1 - 1) * 1
gate(chi_Gm_times_A1 == 0 and chi_A2 == 1,
     "D(x) open-immersion hostile")

self_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(self_tree)),
     "independent audit assertion-free")

semantic = {
    "normalization": {"punctures": 3, "chi_c": -1},
    "branch": {"origin_preimages": 4, "triple_values": 2, "chi_c": -4},
    "surface": {"strata": ((5, 3), (-7, 1), (3, 0)), "chi_c": 8},
    "finite_field_good": good_counts,
    "finite_field_bad_hostiles": bad_counts,
    "consequence": "no finite etale A2_C to U_C",
    "nonclaim": "nonfinite nonproper dominant etale polynomial atlas remains open",
}
semantic_hash = hashlib.sha256(json.dumps(
    semantic, sort_keys=True, separators=(",", ":")
).encode("ascii")).hexdigest()

print("audit=THM3823_independent_compact_Euler_obstruction")
print("normalization=P1_minus_3;chi_c=-1")
print("branch=irreducible;origin_preimages_4;triple_values_2_unibranch;chi_c=-4")
print(f"finite_field_controls=good:{good_counts};bad_hostiles:{bad_counts}")
print("surface_strata=off_branch:3x5;visible_branch:1x(-7);omitted:0x3")
print("surface=chi_c_U=8")
print("finite_etale=excluded_all_positive_degrees;equation_1=8d_impossible")
print("nonfinite_control=D(x)_in_A2;chi_c_0_to_1;atlas_lane_open")
print(f"GATES={GATES}")
print(f"semantic_sha256={semantic_hash}")
print("RESULT=PASS")
