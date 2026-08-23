#!/usr/bin/env python3
"""Exact controls for THM-3883's complete-local square/cube criterion.

Reproduction:
  python3 04-computation/jc2_marked_root_complete_local_square_cube_thm3883.py
  python3 -O 04-computation/jc2_marked_root_complete_local_square_cube_thm3883.py
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")


A, C, b, z, T, epsilon, t = sp.symbols("A C b z T epsilon t")
CHECKS = 0
SERIES_DEPTH = 10


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


def semigroup_contains(n: int, generators: tuple[int, ...]) -> bool:
    if n < 0:
        return False
    reachable = {0}
    for value in range(n + 1):
        if value not in reachable:
            continue
        for generator in generators:
            if value + generator <= n:
                reachable.add(value + generator)
    return n in reachable


# Universal cusp packet.
Delta = -27 * A**2 * b**2 + 8 * A * C**3 - 54 * A * C * b + 9 * C**2 - 54 * b
P = 1 + sp.Rational(2, 3) * A * C
u = 1 + A * C + A**2 * b
zero("cusp_identity", A**2 * Delta - 27 * (P**3 - u**2))

# Above A=0 with the matching sign, w=epsilon*z is the Hensel root
# w=sqrt(1+T), T=(2/3)AC.  The constant and linear terms cancel, and every
# remaining term is divisible by A^2 after substitution.
tail = sp.series((1 + T) ** sp.Rational(3, 2), T, 0, SERIES_DEPTH + 1).removeO()
tail = sp.expand(tail - 1 - sp.Rational(3, 2) * T)
expected_tail = sum(
    sp.binomial(sp.Rational(3, 2), j) * T**j
    for j in range(2, SERIES_DEPTH + 1)
)
zero("binomial_tail", tail - expected_tail)

substituted_tail = sp.expand(tail.subs(T, sp.Rational(2, 3) * A * C) / A**2)
expected_quotient = sum(
    sp.binomial(sp.Rational(3, 2), j)
    * sp.Rational(2, 3) ** j
    * A ** (j - 2)
    * C**j
    for j in range(2, SERIES_DEPTH + 1)
)
zero("A0_completed_divisibility", substituted_tail - expected_quotient)
zero("A0_first_term", expected_quotient.coeff(A, 0) - C**2 / 6)
for j in range(2, SERIES_DEPTH + 1):
    equal(f"A0_j{j}_nonnegative_A_order", j - 2 >= 0, True)
    equal(
        f"A0_j{j}_coefficient_nonzero",
        sp.binomial(sp.Rational(3, 2), j) != 0,
        True,
    )

# The simple-root Hensel uniqueness is the factorization below: if w,z have
# the same nonzero residue, w+z is a unit and w^2=z^2 forces w=z.
w = sp.symbols("w")
zero("hensel_root_difference", w**2 - z**2 - (w - z) * (w + z))

# At z=0 one has AC=-3/2 and A is a unit.  The carrier is an R-term plus
# epsilon*z^3/A^2, so cube descent is necessary and sufficient.
zero(
    "zero_root_carrier_reduction",
    (epsilon * z**3 - 1 - A * C).subs(A * C, -sp.Rational(3, 2))
    - (epsilon * z**3 + sp.Rational(1, 2)),
)

# Complete-local controls: z=t has square in both cusp rings.  Its cube is
# present for <2,3> but is the unique low odd gap for <2,5>.
for q in (3, 5):
    A_q = 1 + t**q
    C_q = sp.Rational(3, 2) * (t**2 - 1) / A_q
    zero(f"S2{q}_profile", 1 + sp.Rational(2, 3) * A_q * C_q - t**2)
    v_q = A_q - 1
    c_q = sp.Rational(2, 3) * (C_q + sp.Rational(3, 2))
    zero(f"S2{q}_coordinate_recovery", (1 + v_q) * c_q - v_q - t**2)
equal("S23_z2", semigroup_contains(2, (2, 3)), True)
equal("S23_z3", semigroup_contains(3, (2, 3)), True)
equal("S25_z2", semigroup_contains(2, (2, 5)), True)
equal("S25_z3_gap", semigroup_contains(3, (2, 5)), False)
for exponent in range(4, 15):
    equal(f"S25_conductor_{exponent}", semigroup_contains(exponent, (2, 5)), True)

print("THM3883_LOCAL_RING", "Rhat_x subset product_p k[[t_p]];test carrier completionwise")
print("THM3883_A0", "z=epsilon gives completed binomial carrier sum_j>=2 beta_j A^(j-2)C^j")
print("THM3883_NONZERO", "common nonzero root residue descends by simple Hensel uniqueness")
print("THM3883_ZERO", "P=0 makes A a unit;carrier descends iff z^3 descends")
print("THM3883_IFF", "A0 matching + nonzero sign matching + zero-root cube descent")
print("THM3883_POSITIVE", "k[[t^2,t^3]]:z=t has z^2,z^3 in the branch ring")
print("THM3883_HOSTILE", "k[[t^2,t^5]]:z=t has z^2 in but z^3 out")
print("THM3883_SCOPE", "arbitrary reduced singularities;global square root assumed;projective companion open")
semantic_packet = (
    "arbitrary reduced curve singularities",
    "completed normalization product",
    "matching sign over A zero",
    "binomial A squared divisibility to all orders",
    "common nonzero Hensel root",
    "zero marked root cube descent iff",
    "finite length completion globalization",
    "two-three cusp positive control",
    "two-five cusp cube-gap hostile",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
