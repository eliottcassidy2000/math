#!/usr/bin/env python3
"""Exact controls for the THM-3866 full polynomial-graph closure.

The theorem is all-degree.  This companion checks generic s(A)-jets through
order four and all first-mismatch degree regimes for N=0,...,4 in a fixed
hostile graph normalization.

Reproduction:
  python3 04-computation/jc2_all_polynomial_graph_first_mismatch_thm3866.py
  python3 -O 04-computation/jc2_all_polynomial_graph_first_mismatch_thm3866.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C = sp.symbols("A C")
tau, eta = sp.symbols("tau eta")
s_jet = sp.symbols("s0:6")
CHECKS = 0
MAX_N = 4


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


# Generic polynomial s(A)-jet: monic graph division is polynomial and keeps
# the universal binomial leading coefficient at every order.
S_generic = sum(s_jet[j] * A**j for j in range(MAX_N + 2))
g_generic = sp.expand(6 * S_generic * (1 + A * S_generic))
b0_generic = sp.expand(2 * S_generic**2 * (3 + 4 * A * S_generic))
q_generic: list[sp.Expr] = []
for j in range(MAX_N + 1):
    rhs = beta(j) * C ** (j + 2) - b0_generic.coeff(A, j)
    rhs += sum(
        g_generic.coeff(A, i) * q_generic[j - i]
        for i in range(1, j + 1)
    )
    quotient, remainder = sp.div(sp.expand(rhs), C - 6 * s_jet[0], C)
    zero(f"generic_q_{j}_remainder", remainder)
    qj = sp.factor(quotient)
    q_generic.append(qj)
    equal(f"generic_q_{j}_degree", sp.degree(qj, C), j + 1)
    equal(f"generic_q_{j}_leading", sp.LC(sp.Poly(qj, C)), beta(j))

# The exact perturbation formula driving the first-mismatch proof.
b_symbol, delta_symbol = sp.symbols("b_symbol delta_symbol")
u_symbol = 1 + A * C + A**2 * b_symbol
equal(
    "universal_profile_perturbation",
    branch(b_symbol + delta_symbol) - branch(b_symbol),
    -54 * u_symbol * delta_symbol - 27 * A**2 * delta_symbol**2,
)

Blead = sp.symbols("Blead")
equal(
    "resonant_top_coefficient_square",
    -27 * A**MAX_N * (Blead**2 + 2 * A * Blead * tau + A**2 * tau**2),
    -27 * A**MAX_N * (Blead + A * tau) ** 2,
)

# Hostile exact grid in the THM-3852 normalization sigma=-1/6.  This checks
# d<N, d=N, and d>N (including a larger special-fibre degree) for N=0..4.
sigma = -sp.Rational(1, 6)
F = sp.expand(C - 6 * sigma * (1 + A * sigma))
b0 = sp.expand(2 * sigma**2 * (3 + 4 * A * sigma))
q: list[sp.Expr] = []
for j in range(MAX_N + 1):
    rhs = beta(j) * C ** (j + 2)
    if j == 0:
        rhs -= 6 * sigma**2
    if j == 1:
        rhs -= 8 * sigma**3
    if j:
        rhs += 6 * sigma**2 * q[j - 1]
    quotient, remainder = sp.div(sp.expand(rhs), C - 6 * sigma, C)
    zero(f"hostile_q_{j}_remainder", remainder)
    q.append(sp.expand(quotient))

for N in range(MAX_N + 1):
    QN = sp.expand(sum(A**j * q[j] for j in range(N)))
    bN = sp.expand(b0 + F * QN)
    RN = sp.cancel(branch(bN) / (A**N * F))
    zero(f"N{N}_base_residual_polynomial", sp.denom(RN) - 1)
    RN = sp.expand(RN)
    uN = sp.expand(1 + A * C + A**2 * bN)

    regimes = sorted({max(0, N - 1), N, N + 1, N + 2})
    for d in regimes:
        T = tau * C**d + eta
        Sres = sp.expand(RN - 54 * uN * T - 27 * A ** (N + 2) * F * T**2)
        b = sp.expand(bN + A**N * F * T)
        equal(
            f"N{N}_d{d}_first_mismatch_factorization",
            branch(b),
            A**N * F * Sres,
        )
        equal(f"N{N}_d{d}_special_residual", Sres.subs(A, 0), 54 * (q[N] - T))

        expected_degree = 2 if N == 0 and d == 0 else 2 * max(N, d) + 1
        equal(f"N{N}_d{d}_generic_degree", sp.degree(Sres, C), expected_degree)
        equal(
            f"N{N}_d{d}_special_degree",
            sp.degree(Sres.subs(A, 0), C),
            max(N + 1, d),
        )
        nonzero(
            f"N{N}_d{d}_strict_degree_drop",
            expected_degree - max(N + 1, d),
        )

        leading = sp.LC(sp.Poly(Sres, C))
        if N == 0 and d == 0:
            expected_leading = 8 * A
        elif d < N:
            expected_leading = -27 * beta(N - 1) ** 2 * A**N
        elif d > N:
            expected_leading = -27 * A ** (N + 2) * tau**2
        else:
            expected_leading = -27 * A**N * (beta(N - 1) + A * tau) ** 2
        equal(f"N{N}_d{d}_leading_coefficient", leading, expected_leading)

print("THM3866_GRAPH_CLASSIFICATION", "F=C-6s(A)(1+As(A)), b=b0+QF")
print("THM3866_FORMAL_QUOTIENT", "Q_*=sum q_j(C)A^j, deg q_j=j+1, lc beta_j")
print("THM3866_FIRST_MISMATCH", "Q=Q_<N+A^N T, T(0,C)!=q_N")
print("THM3866_RESIDUAL", "S=R_N-54u_N T-27A^(N+2)F T^2")
print("THM3866_SPECIAL", "S(0,C)=54(q_N-T(0,C))!=0")
print("THM3866_DEGREES", "d<N:2N+1; d>N:2d+1; d=N:2N+1")
print("THM3866_RESONANCE", "lc=-27A^N(beta_(N-1)+A lc(T))^2 !=0")
print("THM3866_REPLAY_UNIVERSE", "generic s-jets j<=4; N=0..4; d=N-1,N,N+1,N+2")
print("THM3866_SCOPE", "all polynomial graph branches; non-graph A1 and A=0 open")
semantic_packet = (
    "all polynomial graph branches",
    "arbitrary polynomial C-dependent transverse quotient",
    "unique first A-adic mismatch",
    "exact nonlinear perturbation identity",
    "three C-degree regimes",
    "resonant leading coefficient perfect square",
    "finite-base projective contact",
    "reduced component two-place obstruction",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
