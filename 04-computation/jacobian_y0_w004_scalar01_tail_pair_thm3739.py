#!/usr/bin/env python3
"""Exact companion for THM-3739 (W004 scalar-01 tail pair)."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition: object, label: str) -> None:
    if condition is not True and condition != sp.S.true:
        raise ArithmeticError(label)


# W004 fibre word, absolute placements, and inherited boundary gates.
n = sp.symbols("n", integer=True, positive=True)
A_support = (0, n, 3 * n)
B_support = (0, n, 2 * n, 4 * n)
fibres: dict[sp.Expr, list[str]] = defaultdict(list)
for i, x in enumerate(A_support):
    for j, y in enumerate(B_support):
        fibres[sp.expand(x + y)].append(f"{i}{j}")
fibre_word = tuple(
    "+".join(fibres[k]) for k in sorted(fibres, key=lambda z: int(z.coeff(n)))
)
expected_word = ("00", "01+10", "02+11", "12+20", "03+21", "13+22", "23")
require(fibre_word == expected_word, "W004 fibre word")

p_a = (-2, n - 2, 3 * n - 2)
q_a = (1 - n, 1, n + 1, 3 * n + 1)
p_b = (1 - n, 1, 2 * n + 1)
q_b = (-2, n - 2, 2 * n - 2, 4 * n - 2)
require((p_a[0], q_a[1]) == (-2, 1), "family A arm 01")
require((p_b[1], q_b[0]) == (1, -2), "family B arm 10")
for pw, qw, label in ((p_a, q_a, "A"), (p_b, q_b, "B")):
    require(sp.expand(pw[0] + qw[1]) == -1, f"family {label} scalar 01")
    require(sp.expand(pw[1] + qw[0]) == -1, f"family {label} scalar 10")

require((p_a[0], q_a[0].subs(n, 1)) == (-2, 0), "family A n=1 singleton")
require((p_b[0].subs(n, 1), q_b[0]) == (0, -2), "family B n=1 singleton")
parity_exponents = tuple(2 // gcd(2, nv - 1) for nv in range(2, 10))
require(parity_exponents == (2, 1, 2, 1, 2, 1, 2, 1),
        "alternating singleton common-base exponent")
for nv in range(2, 302):
    exponent = 2 // gcd(2, nv - 1)
    require((exponent == 2) == (nv % 2 == 0), f"even parity gate n={nv}")


# Exact W-bracket engine.
a, c, d, t, lam, kap = sp.symbols("a c d t lambda kappa", nonzero=True)
H, K, R, Hp, Kp, Rp = sp.symbols("H K R Hp Kp Rp", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.diff(expr, H) * Hp + sp.diff(expr, K) * Kp + sp.diff(expr, R) * Rp


def bracket(rw: int, sw: int, f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sw * deriv(f) * g - rw * f * deriv(g))


# Complete gcd and identity controls across many odd scales.  The formula
# checks are symbolic identities in H,K,R and their independent jets.
checked_scales: list[int] = []
for nv in range(3, 82, 2):
    m = (nv - 1) // 2
    require(gcd(2, nv - 1) == 2, f"negative endpoint gcd n={nv}")
    require(gcd(3 * nv - 2, 3 * nv + 1) == 1, f"family A high gcd n={nv}")
    require(gcd(2 * nv + 1, 4 * nv - 2) == 1, f"family B high gcd n={nv}")

    # Family A.
    mu = t * sp.Rational(3 * nv + 1, 3 * nv - 2) / d
    U = H * K**2
    f0 = a * H
    g0 = c * H**m
    f2 = d * K ** (3 * nv - 2)
    g3 = t * K ** (3 * nv + 1)
    L = lam * K + a * mu * H * K**3
    M = kap * K ** (nv + 1) + mu * R * K**3

    row_03 = bracket(-2, 3 * nv + 1, f0, g3) + bracket(
        3 * nv - 2, 1, f2, L
    )
    row_13 = bracket(nv - 2, 3 * nv + 1, R, g3) + bracket(
        3 * nv - 2, nv + 1, f2, M
    )
    require(sp.factor(row_03) == 0, f"family A 03+21 transport n={nv}")
    require(sp.factor(row_13) == 0, f"family A 13+22 transport n={nv}")

    phi = R / K ** (nv - 2)
    A = lam + 3 * a * mu * U
    Psi = (nv + 1) * kap + 3 * mu * phi
    C_A = c * d * (3 * nv - 2) * m
    row_02 = bracket(-2, nv + 1, f0, M) + bracket(nv - 2, 1, R, L)
    row_12 = bracket(nv - 2, nv + 1, R, M) + bracket(
        3 * nv - 2, 1 - nv, f2, g0
    )
    reduced_02 = K ** (nv - 1) * (A * deriv(phi) + a * Psi * deriv(U))
    reduced_12 = K ** (2 * nv - 1) * (
        Psi * deriv(phi) - C_A * U ** (m - 1) * deriv(U)
    )
    require(sp.factor(row_02 - reduced_02) == 0,
            f"family A 02+11 reduction n={nv}")
    require(sp.factor(row_12 - reduced_12) == 0,
            f"family A 12+20 reduction n={nv}")

    eta = sp.symbols(f"etaA{nv}")
    phi_solution = (eta - a * (nv + 1) * kap * U) / A
    Psi_solution = sp.factor(Psi.subs(R, phi_solution * K ** (nv - 2)))
    N_A = (nv + 1) * kap * lam + 3 * mu * eta
    require(sp.factor(A * Psi_solution - N_A) == 0,
            f"family A Mobius numerator n={nv}")
    z = sp.symbols(f"zA{nv}")
    A_z = lam + 3 * a * mu * z
    require(sp.Poly(A_z**3 * z ** (m - 1), z).degree() == m + 2,
            f"family A final degree n={nv}")

    # Family B.
    mu = t * sp.Rational(4 * nv - 2, 2 * nv + 1) / d
    V = H * K**2
    U = V**m
    f0 = a * H**m
    g0 = c * H
    f2 = d * K ** (2 * nv + 1)
    g3 = t * K ** (4 * nv - 2)
    L = lam * K ** (nv - 2) + a * mu * H**m * K ** (2 * nv - 3)
    M = kap * K ** (2 * nv - 2) + mu * R * K ** (2 * nv - 3)

    row_03 = bracket(1 - nv, 4 * nv - 2, f0, g3) + bracket(
        2 * nv + 1, nv - 2, f2, L
    )
    row_13 = bracket(1, 4 * nv - 2, R, g3) + bracket(
        2 * nv + 1, 2 * nv - 2, f2, M
    )
    require(sp.factor(row_03) == 0, f"family B 03+21 transport n={nv}")
    require(sp.factor(row_13) == 0, f"family B 13+22 transport n={nv}")

    phi = R / K
    A = (nv - 2) * lam + a * mu * (2 * nv - 3) * U
    Psi = 2 * (nv - 1) * kap + mu * (2 * nv - 3) * phi
    C_B = c * d * (2 * nv + 1)
    row_02 = bracket(1 - nv, 2 * nv - 2, f0, M) + bracket(
        1, nv - 2, R, L
    )
    row_12 = bracket(1, 2 * nv - 2, R, M) + bracket(
        2 * nv + 1, -2, f2, g0
    )
    reduced_02 = K ** (nv - 1) * (A * deriv(phi) + a * Psi * deriv(U))
    reduced_12 = K ** (2 * nv - 1) * (
        Psi * deriv(phi) - C_B * deriv(U) / (m * V ** (m - 1))
    )
    require(sp.factor(row_02 - reduced_02) == 0,
            f"family B 02+11 reduction n={nv}")
    require(sp.factor(row_12 - reduced_12) == 0,
            f"family B 12+20 reduction n={nv}")

    eta = sp.symbols(f"etaB{nv}")
    phi_solution = (eta - 2 * a * (nv - 1) * kap * U) / A
    Psi_solution = sp.factor(Psi.subs(R, phi_solution * K))
    N_B = 2 * (nv - 1) * (nv - 2) * kap * lam + mu * (2 * nv - 3) * eta
    require(sp.factor(A * Psi_solution - N_B) == 0,
            f"family B Mobius numerator n={nv}")
    z = sp.symbols(f"zB{nv}")
    A_z = (nv - 2) * lam + a * mu * (2 * nv - 3) * z**m
    require(sp.Poly(A_z**3, z).degree() == 3 * m,
            f"family B left degree n={nv}")
    require(sp.Poly(z ** (m - 1), z).degree() == m - 1,
            f"family B right degree n={nv}")
    require(3 * m > m - 1, f"family B strict degree gap n={nv}")
    checked_scales.append(nv)


# At n=3 the two absolute placements literally merge.
require(tuple(x.subs(n, 3) if hasattr(x, "subs") else x for x in p_a)
        == tuple(x.subs(n, 3) if hasattr(x, "subs") else x for x in p_b),
        "n=3 P-support merger")
require(tuple(x.subs(n, 3) if hasattr(x, "subs") else x for x in q_a)
        == tuple(x.subs(n, 3) if hasattr(x, "subs") else x for x in q_b),
        "n=3 Q-support merger")


print("THM-3739 exact W004 scalar-01 tail-pair audit")
print("fibre word = " + "|".join(fibre_word))
print("family A = P(-2,n-2,3n-2), Q(1-n,1,n+1,3n+1)")
print("family B = P(1-n,1,2n+1), Q(-2,n-2,2n-2,4n-2)")
print("inherited gates = n=1 singleton failure; every even n parity-square death")
print("odd exact window = 3..81; 40 scales; both orientations")
print("family A mechanism = Mobius first integral then A^3 U^(m-1)=constant")
print("family B mechanism = Mobius first integral then degree 3m versus m-1")
print("n=3 = actual-support merger verified in both orientations")
print("tail consequence = only separately reserved anchor-20 family remains uncaptured")
print("scope = both named W004 scalar fibre 01+10 tail placements")
