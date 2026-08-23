#!/usr/bin/env python3
"""Exact controls for THM-3870's vertical-axis first-mismatch theorem.

The proof is all-degree.  This companion verifies the universal vertical
classification and nonlinear response, replays finite truncations through
N=5 in every degree regime around d=N+2, and checks the reverse-graph affine
boundary used in the aggregate two-axis synthesis.

Reproduction:
  python3 04-computation/jc2_vertical_axis_first_mismatch_thm3870.py
  python3 -O 04-computation/jc2_vertical_axis_first_mismatch_thm3870.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C = sp.symbols("A C")
tau, eta, lam = sp.symbols("tau eta lam", nonzero=True)
b_symbol, delta_symbol, Q_symbol = sp.symbols(
    "b_symbol delta_symbol Q_symbol"
)
CHECKS = 0
MAX_N = 5


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


# A divides Delta exactly when the A=0 profile is C^2/6.  The displayed
# quotient is also a hostile control on every sign in the vertical chart.
vertical_profile = C**2 / 6 + A * Q_symbol
vertical_residual = sp.expand(
    -C**3
    - 54 * Q_symbol
    - sp.Rational(3, 4) * A * C**4
    - 54 * A * C * Q_symbol
    - 9 * A**2 * C**2 * Q_symbol
    - 27 * A**3 * Q_symbol**2
)
equal(
    "vertical_profile_factorization",
    branch(vertical_profile),
    A * vertical_residual,
)
equal(
    "vertical_special_classification",
    branch(b_symbol).subs(A, 0),
    9 * C**2 - 54 * b_symbol,
)

# Exact quadratic response in the profile variable.
u_symbol = 1 + A * C + A**2 * b_symbol
equal(
    "universal_profile_response",
    branch(b_symbol + delta_symbol) - branch(b_symbol),
    -54 * u_symbol * delta_symbol - 27 * A**2 * delta_symbol**2,
)

# Every generalized binomial coefficient used by the replay is nonzero.
for j in range(MAX_N + 2):
    nonzero(f"beta_{j}_nonzero", beta(j))

# Truncation and first-mismatch grid.  The all-degree proof is the degree
# comparison encoded here; the finite grid is an adversarial replay only.
for N in range(MAX_N + 1):
    bN = sp.expand(sum(beta(j) * A**j * C ** (j + 2) for j in range(N + 1)))
    uN = sp.expand(1 + A * C + A**2 * bN)
    RN = sp.cancel(branch(bN) / A ** (N + 1))
    equal(f"N{N}_base_residual_denominator", sp.denom(RN), 1)
    RN = sp.expand(RN)
    equal(
        f"N{N}_exact_truncation_factorization",
        branch(bN),
        A ** (N + 1) * RN,
    )
    equal(
        f"N{N}_base_special",
        RN.subs(A, 0),
        54 * beta(N + 1) * C ** (N + 3),
    )
    equal(f"N{N}_base_degree", sp.degree(RN, C), 2 * N + 4)
    equal(
        f"N{N}_base_leading",
        sp.LC(sp.Poly(RN, C)),
        -27 * beta(N) ** 2 * A ** (N + 1),
    )
    equal(f"N{N}_u_degree", sp.degree(uN, C), N + 2)
    equal(
        f"N{N}_u_leading",
        sp.LC(sp.Poly(uN, C)),
        beta(N) * A ** (N + 2),
    )

    regimes = sorted({0, max(0, N + 1), N + 2, N + 3, N + 4})
    for d in regimes:
        T = tau * C**d + eta
        S = sp.expand(RN - 54 * uN * T - 27 * A ** (N + 3) * T**2)
        b = sp.expand(bN + A ** (N + 1) * T)
        equal(
            f"N{N}_d{d}_factorization",
            branch(b),
            A ** (N + 1) * S,
        )
        equal(
            f"N{N}_d{d}_special",
            S.subs(A, 0),
            54 * (beta(N + 1) * C ** (N + 3) - T),
        )

        expected_degree = 2 * max(N + 2, d)
        equal(
            f"N{N}_d{d}_generic_degree",
            sp.degree(S, C),
            expected_degree,
        )
        nonzero(
            f"N{N}_d{d}_strict_degree_drop",
            expected_degree - max(N + 3, d),
        )

        leading = sp.LC(sp.Poly(S, C))
        if d < N + 2:
            expected_leading = -27 * beta(N) ** 2 * A ** (N + 1)
        elif d > N + 2:
            expected_leading = -27 * A ** (N + 3) * tau**2
        else:
            expected_leading = (
                -27 * A ** (N + 1) * (beta(N) + A * tau) ** 2
            )
        equal(f"N{N}_d{d}_leading", leading, expected_leading)

# Reverse polynomial graphs A=h(C) reduce either to A=0 or to an invertible
# affine graph.  The nonvertical solution has z+1=lambda*C.
z = -1 + lam * C
h = sp.expand(sp.Rational(3, 2) * (z**2 - 1) / C)
B = sp.expand(sp.Rational(2, 9) * (2 * z + 1) / lam**2)
equal("reverse_affine_h", h, -3 * lam + sp.Rational(3, 2) * lam**2 * C)
equal("reverse_divisibility_identity", 9 * (z + 1) ** 2 * B, 2 * C**2 * (2 * z + 1))
equal("reverse_graph_branch", branch(B).subs(A, h), 0)
equal(
    "reverse_affine_inverse",
    sp.Rational(2, 3) * (h + 3 * lam) / lam**2,
    C,
)

print("THM3870_VERTICAL_IFF", "A|Delta_b iff b=C^2/6+A Q(A,C)")
print("THM3870_FORMAL_PROFILE", "Q_*=sum beta_(j+1) A^j C^(j+3)")
print("THM3870_FIRST_MISMATCH", "Q=Q_<N+A^N T, T(0,C)!=beta_(N+1)C^(N+3)")
print("THM3870_RESIDUAL", "S=R_N-54u_N T-27A^(N+3)T^2")
print("THM3870_DEGREES", "d<N+2:2N+4; d>N+2:2d; d=N+2:2N+4")
print("THM3870_RESONANCE", "lc=-27A^(N+1)(beta_N+A lc(T))^2 !=0")
print("THM3870_REVERSE_GRAPH", "A=h(C) is A=0 or invertible affine")
print("THM3870_REPLAY_UNIVERSE", "N=0..5; d=0,N+1,N+2,N+3,N+4")
print("THM3870_SCOPE", "vertical axis + reverse graphs; non-graph A1 remains open")
semantic_packet = (
    "vertical branch exact classification",
    "unique binomial formal comparator",
    "first A-adic mismatch",
    "three C-degree regimes",
    "N+2 resonance perfect square",
    "finite-base projective contact",
    "multiplicity-safe reduced companion",
    "reverse graph affine reduction",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
