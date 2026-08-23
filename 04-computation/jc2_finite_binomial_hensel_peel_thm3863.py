#!/usr/bin/env python3
"""Exact finite controls for the all-depth THM-3863 proof candidate.

The proof is symbolic in N.  This companion independently constructs and
checks the first six nontrivial truncations.

Reproduction:
  python3 04-computation/jc2_finite_binomial_hensel_peel_thm3863.py
  python3 -O 04-computation/jc2_finite_binomial_hensel_peel_thm3863.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, Z, sigma = sp.symbols("A C Z sigma")
CHECKS = 0
MAX_N = 6


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


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def beta(j: int) -> sp.Rational:
    return sp.factor(
        sp.binomial(sp.Rational(3, 2), j + 2)
        * sp.Rational(2, 3) ** (j + 2)
    )


def branch(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


F0 = C - 6 * sigma
F = F0 - 6 * A * sigma**2
b0 = 6 * sigma**2 + 8 * A * sigma**3

# Coefficient recursion for Q_*=(b_*-b0)/F.  The exact branch evaluation
# makes every Euclidean remainder zero.
q: list[sp.Expr] = []
for j in range(MAX_N + 1):
    bj = beta(j) * C ** (j + 2)
    if j == 0:
        bj -= 6 * sigma**2
    if j == 1:
        bj -= 8 * sigma**3
    numerator = sp.expand(bj + (6 * sigma**2 * q[j - 1] if j else 0))
    quotient, remainder = sp.div(numerator, F0, C)
    zero(f"q_{j}_division_remainder", remainder)
    qj = sp.factor(quotient)
    q.append(qj)
    equal(f"q_{j}_degree", sp.degree(qj, C), j + 1)
    equal(f"q_{j}_leading", sp.LC(sp.Poly(qj, C)), beta(j))
    nonzero(f"beta_{j}_nonzero", beta(j))

# Check the first six complete truncations, their exact A-adic factor, and
# the projective C-infinity contact Z^N over A=0.
for N in range(1, MAX_N + 1):
    QN = sp.expand(sum(A**j * q[j] for j in range(N)))
    bN = sp.expand(b0 + F * QN)
    DeltaN = branch(bN)
    RN = sp.cancel(DeltaN / (A**N * F))
    zero(f"N{N}_residual_polynomial", sp.denom(RN) - 1)
    RN = sp.expand(RN)

    equal(f"N{N}_factorization", DeltaN, A**N * F * RN)
    equal(f"N{N}_profile_degree", sp.degree(bN, C), N + 1)
    equal(
        f"N{N}_profile_leading",
        sp.LC(sp.Poly(bN, C)),
        beta(N - 1) * A ** (N - 1),
    )
    equal(f"N{N}_residual_degree", sp.degree(RN, C), 2 * N + 1)
    equal(
        f"N{N}_residual_leading",
        sp.LC(sp.Poly(RN, C)),
        -27 * beta(N - 1) ** 2 * A**N,
    )
    equal(f"N{N}_special_residual", RN.subs(A, 0), 54 * q[N])
    equal(f"N{N}_special_degree", sp.degree(RN.subs(A, 0), C), N + 1)
    equal(
        f"N{N}_special_leading",
        sp.LC(sp.Poly(RN.subs(A, 0), C)),
        54 * beta(N),
    )

    D = 2 * N + 1
    RN_projective_C = sp.expand(
        sum(
            sp.Poly(RN, C).coeff_monomial(C**j) * C**j * Z ** (D - j)
            for j in range(D + 1)
        )
    )
    qN_projective = sp.expand(
        sum(
            sp.Poly(q[N], C).coeff_monomial(C**j)
            * C**j
            * Z ** (N + 1 - j)
            for j in range(N + 2)
        )
    )
    equal(
        f"N{N}_projective_special_fibre",
        RN_projective_C.subs(A, 0),
        54 * qN_projective * Z**N,
    )
    nonzero(
        f"N{N}_infinity_contact_exact",
        sp.Poly(qN_projective.subs(C, 1), Z).coeff_monomial(1),
    )

# The first nontrivial two-layer peel in the THM-3852 normalization.
N = 2
QN2 = sp.expand(q[0] + A * q[1])
bN2 = sp.expand((b0 + F * QN2).subs(sigma, -sp.Rational(1, 6)))
F2 = sp.factor(F.subs(sigma, -sp.Rational(1, 6)))
R2 = sp.cancel(branch(bN2) / (A**2 * F2))
equal(
    "normalized_depth_two_profile",
    bN2,
    (
        2 * A**2 * C**2
        - 2 * A**2 * C
        - A**2
        - 12 * A * C**3
        + 108 * C**2
    )
    / 648,
)
equal(
    "normalized_depth_two_special_residual",
    R2.subs(A, 0),
    (3 * C**3 - 3 * C**2 + C + 1) / 12,
)

print("THM3863_FORMAL_LIFT", "b_*=[(1+2AC/3)^(3/2)-1-AC]/A^2")
print("THM3863_GRAPH", "F=C-6sigma(1+A sigma), b0=2sigma^2(3+4A sigma)")
print("THM3863_QUOTIENT", "Q_*=(b_*-b0)/F=sum q_j(C)A^j")
print("THM3863_Q_DEGREE", "deg q_j=j+1, lc q_j=binom(3/2,j+2)(2/3)^(j+2)")
print("THM3863_TRUNCATION", "Delta_N=A^N F R_N")
print("THM3863_DEGREE_DROP", "deg_C R_N=2N+1, deg_C R_N(0,C)=N+1")
print("THM3863_CONTACT", "projective special fibre has exact Z^N contact at C=infinity")
print("THM3863_REPLAY_RANGE", f"N=1..{MAX_N}")
print("THM3863_SCOPE", "canonical finite truncations only; arbitrary C-dependent quotient open")
semantic_packet = (
    "unique binomial Hensel lift",
    "monic graph division recursion",
    "all coefficients polynomial and degree-sharp",
    "finite truncation A-adic branch multiplicity",
    "residual C-degree drop N",
    "projective infinity contact at finite base",
    "componentwise two-puncture consequence",
    "THM3859 Newton edge and THM3860 vertical pole synthesis",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
