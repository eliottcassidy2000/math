#!/usr/bin/env python3
"""Exact controls for THM-3873's first non-graph A1 closure.

The proof is all-degree.  This companion checks the triangular normalization,
the formal quotient through order six, the exceptional N=0 residual, and all
three degree regimes around d=N-1 for N=1,...,6.

Reproduction:
  python3 04-computation/jc2_first_nongraph_parabola_thm3873.py
  python3 -O 04-computation/jc2_first_nongraph_parabola_thm3873.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, T, Z, r = sp.symbols("A T Z r")
tau, eta = sp.symbols("tau eta", nonzero=True)
b_symbol, delta_symbol = sp.symbols("b_symbol delta_symbol")
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


C = 6 * (T + A**2)
F = A - T**2
b0 = 6 * A + 8 * A**2 * T


def branch(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


# Triangular normalization and explicit inverse.
A_r = r**2
C_r = 6 * (r + r**4)
equal("normalization_inverse_T", (C_r / 6 - A_r**2), r)
equal("normalization_equation", F.subs({A: A_r, T: r}), 0)
equal(
    "fixed_coordinate_equation",
    36 * A_r - (C_r - 6 * A_r**2) ** 2,
    0,
)
equal(
    "marked_profile_on_normalization",
    b0.subs({A: A_r, T: r}),
    2 * r**2 * (3 + 4 * r**3),
)
nonzero("not_forward_graph_odd_term", C_r.subs(r, 1) - C_r.subs(r, -1))
nonzero("not_reverse_graph_degree_gap", sp.degree(C_r, r) - sp.degree(A_r, r))
equal("marked_branch_factor", sp.rem(branch(b0), F, T), 0)
equal("projective_F_avoids_contact", (A * Z**2 - T**2).subs({A: 0, T: 1, Z: 0}), -1)

# Exact nonlinear response.
u_symbol = 1 + A * C + A**2 * b_symbol
equal(
    "universal_profile_response",
    branch(b_symbol + delta_symbol) - branch(b_symbol),
    -54 * u_symbol * delta_symbol - 27 * A**2 * delta_symbol**2,
)

# Formal quotient Q_*=(b_*-b0)/F.  Division by the monic quadratic in T is
# exact; coefficient recursion freezes polynomiality, degree, and top term.
bstar = sp.expand(
    sum(beta(j) * A**j * C ** (j + 2) for j in range(MAX_N + 1))
)
numerator = sp.expand(bstar - b0)
q: list[sp.Expr] = []
for j in range(MAX_N + 1):
    n_j = numerator.coeff(A, j)
    previous = q[j - 1] if j else 0
    q_j = sp.cancel((previous - n_j) / T**2)
    equal(f"q_{j}_denominator", sp.denom(q_j), 1)
    q_j = sp.expand(q_j)
    q.append(q_j)
    equal(f"q_{j}_division_row", previous - T**2 * q_j, n_j)
    equal(f"q_{j}_degree", sp.degree(q_j, T), j)
    equal(
        f"q_{j}_leading",
        sp.LC(sp.Poly(q_j, T)),
        -beta(j) * 6 ** (j + 2),
    )

equal("q0", q[0], -6)
equal("q1", q[1], 4 * T)
equal("q2", q[2], -6 * T**2)
equal("q3", q[3], 12 * T**3 + 6)
equal("q4", q[4], -28 * T**4 - 12 * T)

# Finite hostile replay of the all-degree first-mismatch proof.
for N in range(MAX_N + 1):
    QN = sp.expand(sum(A**j * q[j] for j in range(N)))
    bN = sp.expand(b0 + F * QN)
    uN = sp.expand(1 + A * C + A**2 * bN)
    RN = sp.cancel(branch(bN) / (A**N * F))
    equal(f"N{N}_base_denominator", sp.denom(RN), 1)
    RN = sp.expand(RN)
    equal(f"N{N}_base_factorization", branch(bN), A**N * F * RN)
    equal(f"N{N}_base_special", RN.subs(A, 0), 54 * q[N])

    if N == 0:
        equal("N0_R_degree", sp.degree(RN, T), 1)
        equal("N0_R_leading", sp.LC(sp.Poly(RN, T)), -1728 * A)
        regimes = [0, 1, 2, 3]
    else:
        gamma = -beta(N - 1) * 6 ** (N + 1)
        equal(f"N{N}_R_degree", sp.degree(RN, T), 2 * N)
        equal(
            f"N{N}_R_leading",
            sp.LC(sp.Poly(RN, T)),
            27 * gamma**2 * A**N,
        )
        equal(f"N{N}_u_degree", sp.degree(uN, T), N + 1)
        equal(
            f"N{N}_u_leading",
            sp.LC(sp.Poly(uN, T)),
            -gamma * A ** (N + 1),
        )
        regimes = sorted({0, max(0, N - 2), N - 1, N, N + 1})

    for d in regimes:
        H = tau if d == 0 else tau * T**d + eta
        S = sp.expand(RN - 54 * uN * H - 27 * A ** (N + 2) * F * H**2)
        b = sp.expand(bN + A**N * F * H)
        equal(f"N{N}_d{d}_factorization", branch(b), A**N * F * S)
        equal(f"N{N}_d{d}_special", S.subs(A, 0), 54 * (q[N] - H))

        if N == 0:
            expected_degree = 2 * d + 2
            special_bound = d
            expected_leading = 27 * A**2 * tau**2
        elif d < N - 1:
            expected_degree = 2 * N
            special_bound = max(N, d)
            expected_leading = 27 * gamma**2 * A**N
        elif d > N - 1:
            expected_degree = 2 * d + 2
            special_bound = max(N, d)
            expected_leading = 27 * A ** (N + 2) * tau**2
        else:
            expected_degree = 2 * N
            special_bound = N
            expected_leading = 27 * A**N * (gamma + A * tau) ** 2

        equal(f"N{N}_d{d}_degree", sp.degree(S, T), expected_degree)
        nonzero(
            f"N{N}_d{d}_strict_drop",
            expected_degree - special_bound,
        )
        equal(
            f"N{N}_d{d}_leading",
            sp.LC(sp.Poly(S, T)),
            expected_leading,
        )

print("THM3873_NORMALIZATION", "A=r^2, C=6(r+r^4), T=C/6-A^2=r")
print("THM3873_BRANCH", "F=A-T^2, b0=6A+8A^2T, b=b0+FQ(A,T)")
print("THM3873_FORMAL", "Q_*=sum q_j(T)A^j, deg q_j=j, lc=-beta_j6^(j+2)")
print("THM3873_FIRST_MISMATCH", "Q=Q_<N+A^N H, H(0,T)!=q_N(T)")
print("THM3873_RESIDUAL", "S=R_N-54u_NH-27A^(N+2)F H^2")
print("THM3873_DEGREES", "N>=1 threshold d=N-1; N=0 is exceptional")
print("THM3873_RESONANCE", "lc=27A^N(gamma_(N-1)+A lc(H))^2 !=0")
print("THM3873_REPLAY_UNIVERSE", "q_0..q_6; N=0..6; d near N-1")
print("THM3873_SCOPE", "first fixed-coordinate non-graph A1; higher triangular depth open")
semantic_packet = (
    "first non-graph marked-root A1",
    "triangular normalization and polynomial inverse",
    "monic quadratic formal division",
    "all-depth polynomial coefficient law",
    "exceptional N=0 residual",
    "N-1 resonance perfect square",
    "multiplicity-safe projective companion",
    "higher triangular depth remains open",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
