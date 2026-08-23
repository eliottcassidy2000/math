#!/usr/bin/env python3
"""Exact controls for THM-3876's monomial branch-ring non-descent.

The proof is all-m.  This companion checks the universal two-address
collision identity, the marked-profile formula, the m=1 and m=2 boundaries,
and primitive-root hostiles for m=3,...,20.

Reproduction:
  python3 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
  python3 -O 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


t, zeta, U, A, C, b = sp.symbols("t zeta U A C b")
CHECKS = 0
MAX_M = 20


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


# Universal cusp identity and marked-root calculation with U=t^(m+1).
Delta = -27 * A**2 * b**2 + 8 * A * C**3 - 54 * A * C * b + 9 * C**2 - 54 * b
P = 1 + sp.Rational(2, 3) * A * C
u = 1 + A * C + A**2 * b
zero("universal_cusp_identity", A**2 * Delta - 27 * (P**3 - u**2))

P_square = sp.expand(1 + 4 * U * (1 + U))
zero("marked_P_square", P_square - (1 + 2 * U) ** 2)
zero(
    "forced_marked_profile_numerator",
    (1 + 2 * U) ** 3 - 1 - 6 * U * (1 + U) - 2 * U**2 * (3 + 4 * U),
)

# For a primitive mth root zeta and U=-1/(zeta+1), the two normalization
# addresses t and zeta*t have the same (A,C).  Their B-values differ by the
# displayed nonzero universal factor.
collision_U = -1 / (zeta + 1)
coordinate_packet = sp.factor(
    (zeta - 1) + (zeta**2 - 1) * collision_U
)
zero("universal_C_collision", coordinate_packet)

B_difference_over_t2 = sp.factor(
    2 * (3 * (zeta**2 - 1) + 4 * (zeta**3 - 1) * collision_U)
)
zero(
    "universal_B_difference",
    B_difference_over_t2 + 2 * (zeta - 1) ** 3 / (zeta + 1),
)

# Exact endpoint controls.
A1 = t
C1 = 6 * t * (1 + t**2)
B1 = 2 * t**2 * (3 + 4 * t**2)
zero("m1_forward_C", C1 - 6 * A1 * (1 + A1**2))
zero("m1_forward_B", B1 - 2 * A1**2 * (3 + 4 * A1**2))

A2 = t**2
C2 = 6 * t * (1 + t**3)
T2 = sp.expand(C2 / 6 - A2**2)
B2 = 2 * t**2 * (3 + 4 * t**3)
equal("m2_triangular_inverse", T2, t)
zero("m2_marked_profile", B2 - (6 * A2 + 8 * A2**2 * T2))
nonzero("m2_minus_address_separates", C2.subs(t, -t) - C2)

# A primitive mth root exists for every m>=3 and is neither 1 nor -1.
# Exact cyclotomic gcds freeze those exclusions through m=20; the theorem
# uses the same order argument uniformly for all m.
for m in range(3, MAX_M + 1):
    phi = sp.Poly(sp.cyclotomic_poly(m, zeta), zeta)
    equal(f"m{m}_primitive_degree_positive", phi.degree() > 0, True)
    equal(
        f"m{m}_primitive_not_one",
        sp.gcd(phi, sp.Poly(zeta - 1, zeta)).degree(),
        0,
    )
    equal(
        f"m{m}_primitive_not_minus_one",
        sp.gcd(phi, sp.Poly(zeta + 1, zeta)).degree(),
        0,
    )
    nonzero(
        f"m{m}_B_collision_numerator_mod_phi",
        sp.rem(sp.Poly((zeta - 1) ** 3, zeta), phi).as_expr(),
    )

print("THM3876_TOWER", "A=t^m,C=6t(1+t^(m+1)),B=2t^2(3+4t^(m+1))")
print("THM3876_NORMALIZATION", "k[t^m,t+t^(m+2)] has normalization k[t]")
print("THM3876_COLLISION", "zeta^m=1,U=t^(m+1)=-1/(zeta+1)")
print("THM3876_SAME_TARGET", "(A,C)(zeta*t)=(A,C)(t)")
print("THM3876_DIFFERENT_MARK", "B(zeta*t)-B(t)=-2t^2(zeta-1)^3/(zeta+1)!=0")
print("THM3876_M1", "forward graph;THM3866")
print("THM3876_M2", "triangular parabola;THM3873")
print("THM3876_MGE3", "forced B does not descend;no polynomial carrier profile")
print("THM3876_REPLAY_UNIVERSE", "universal identity;primitive roots m=3..20")
print("THM3876_SCOPE", "monomial marked-root tower only;general singular A1 normalization open")
semantic_packet = (
    "monomial marked-root normalization tower",
    "forced marked value on normalization",
    "finite birational but noninjective branch map",
    "explicit primitive-root double address",
    "same target coordinates and unequal marked values",
    "branch-ring descent obstruction",
    "m1 forward graph boundary",
    "m2 first non-graph triangular boundary",
    "all m at least three excluded before transverse search",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
