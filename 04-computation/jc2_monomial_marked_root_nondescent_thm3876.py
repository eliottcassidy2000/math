#!/usr/bin/env python3
"""Exact controls for THM-3876's two-exponent monomial non-descent.

For positive original exponents (m,n), divide by g=gcd(m,n) and write
(M,N)=(m/g,n/g).  The proof is all-(M,N).  This companion checks the
universal marked-root formula, the primitive-root collision, both descent
boundaries M=1,2, and coprime cyclotomic hostiles on a finite exact grid.

Reproduction:
  python3 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
  python3 -O 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
"""

from __future__ import annotations

import hashlib
import math

import sympy as sp


r, zeta, eta, U, A, C, b = sp.symbols("r zeta eta U A C b")
CHECKS = 0
MAX_M = 16
MAX_N = 12


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


# Universal cusp identity and marked-root calculation with U=r^(M+N).
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

# For eta=zeta^N primitive of order M and U=-1/(eta+1), the addresses
# r and zeta*r have the same (A,C), whereas the forced B-values differ.
collision_U = -1 / (eta + 1)
zero(
    "universal_C_collision",
    (eta - 1) + (eta**2 - 1) * collision_U,
)
B_difference_over_r2N = sp.factor(
    2 * (3 * (eta**2 - 1) + 4 * (eta**3 - 1) * collision_U)
)
zero(
    "universal_B_difference",
    B_difference_over_r2N + 2 * (eta - 1) ** 3 / (eta + 1),
)

# Exact M=1 and M=2 endpoint controls for a range of reduced N.  At M=2,
# coprimality forces N odd and S=C/6-A^(N+1)=r^N.
for N in range(1, MAX_N + 1):
    A1 = r
    C1 = 6 * r**N * (1 + r ** (N + 1))
    B1 = 2 * r ** (2 * N) * (3 + 4 * r ** (N + 1))
    zero(f"M1_N{N}_C_descent", C1 - 6 * A1**N * (1 + A1 ** (N + 1)))
    zero(f"M1_N{N}_B_descent", B1 - 2 * A1 ** (2 * N) * (3 + 4 * A1 ** (N + 1)))

    if N % 2 == 1:
        A2 = r**2
        C2 = 6 * r**N * (1 + r ** (N + 2))
        S2 = sp.expand(C2 / 6 - A2 ** (N + 1))
        B2 = 2 * r ** (2 * N) * (3 + 4 * r ** (N + 2))
        zero(f"M2_N{N}_recover_s", S2 - r**N)
        zero(f"M2_N{N}_B_descent", B2 - 2 * S2**2 * (3 + 4 * A2 * S2))

# Exact primitive-root hostiles.  For every coprime (M,N) in the frozen
# grid, eta=zeta^N has exact order M, hence is neither 1 nor -1.  The theorem
# uses this elementary order argument uniformly beyond the grid.
COPRIME_PAIRS = 0
for M in range(3, MAX_M + 1):
    phi = sp.Poly(sp.cyclotomic_poly(M, zeta), zeta)
    equal(f"M{M}_primitive_degree_positive", phi.degree() > 0, True)
    for N in range(1, MAX_N + 1):
        if math.gcd(M, N) != 1:
            continue
        COPRIME_PAIRS += 1
        eta_mn = zeta**N
        equal(
            f"M{M}_N{N}_eta_not_one",
            sp.gcd(phi, sp.Poly(eta_mn - 1, zeta)).degree(),
            0,
        )
        equal(
            f"M{M}_N{N}_eta_not_minus_one",
            sp.gcd(phi, sp.Poly(eta_mn + 1, zeta)).degree(),
            0,
        )
        nonzero(
            f"M{M}_N{N}_B_collision_numerator_mod_phi",
            sp.rem(sp.Poly((eta_mn - 1) ** 3, zeta), phi).as_expr(),
        )

print("THM3876_REDUCTION", "g=gcd(m,n);r=t^g;M=m/g;N=n/g;gcd(M,N)=1")
print("THM3876_TOWER", "A=r^M,C=6r^N(1+r^(M+N)),B=2r^(2N)(3+4r^(M+N))")
print("THM3876_NORMALIZATION", "k[r^M,r^N+r^(M+2N)] has normalization k[r]")
print("THM3876_COLLISION", "eta=zeta^N,U=r^(M+N)=-1/(eta+1)")
print("THM3876_SAME_TARGET", "(A,C)(zeta*r)=(A,C)(r)")
print("THM3876_DIFFERENT_MARK", "B(zeta*r)-B(r)=-2r^(2N)(eta-1)^3/(eta+1)!=0")
print("THM3876_M1", "A=r;forced B descends")
print("THM3876_M2", "N odd;S=C/6-A^(N+1)=r^N;forced B descends")
print("THM3876_MGE3", "forced B does not descend;no polynomial carrier profile")
print("THM3876_IFF", "forced marked value descends iff reduced exponent M<=2")
print(
    "THM3876_REPLAY_UNIVERSE",
    f"universal identities;M<=16,N<=12;coprime_Mge3_pairs={COPRIME_PAIRS}",
)
print("THM3876_SCOPE", "two-exponent monomial marked-root tower only;general singular A1 normalization open")
semantic_packet = (
    "two-exponent monomial marked-root normalization tower",
    "gcd reduction to coprime exponents M N",
    "forced marked value on normalization",
    "finite birational but noninjective branch map",
    "explicit primitive-root double address",
    "same target coordinates and unequal marked values",
    "branch-ring descent obstruction",
    "M1 forward-coordinate boundary",
    "M2 triangular recovery boundary for odd N",
    "descent if and only if reduced A exponent at most two",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
