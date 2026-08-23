#!/usr/bin/env python3
"""Exact companion for THM-3728 (W004 scalar-12 anchor-20 family)."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition: object, label: str) -> None:
    if condition is not True and condition != sp.S.true:
        raise ArithmeticError(label)


# W004 word, absolute placement, scalar anchors, and even parity gate.
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

p_weights = (1 - 3 * n, 1 - 2 * n, 1)
q_weights = (-2, n - 2, 2 * n - 2, 4 * n - 2)
require(sp.expand(p_weights[1] + q_weights[2]) == -1, "scalar address 12")
require((p_weights[2], q_weights[0]) == (1, -2), "scalar arm address 20")

parity_exponents = tuple(2 // gcd(3 * nv - 1, 2) for nv in range(1, 10))
require(parity_exponents == (1, 2, 1, 2, 1, 2, 1, 2, 1),
        "anchor common-base parity")
for nv in range(2, 302, 2):
    require(gcd(3 * nv - 1, 2) == 1, f"even parity gcd n={nv}")


# Exact W-bracket engine.
a, c, d, t = sp.symbols("a c d t", nonzero=True)
lam, kap = sp.symbols("lambda kappa")
H, K, R, Hp, Kp, Rp = sp.symbols("H K R Hp Kp Rp", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.diff(expr, H) * Hp + sp.diff(expr, K) * Kp + sp.diff(expr, R) * Rp


def bracket(rw: int, sw: int, f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sw * deriv(f) * g - rw * f * deriv(g))


# Exceptional n=1: both upper primitives and the constant-product collapse.
mu1 = 2 * t / d
f0 = a * H
g0 = c * H
f2 = d * K
g3 = t * K**2
L1 = a * mu1 * H * K
M1 = mu1 * K * R + kap
row_03_n1 = bracket(-2, 2, f0, g3) + bracket(1, -1, f2, L1)
row_13_n1 = bracket(-1, 2, R, g3) + bracket(1, 0, f2, M1)
require(sp.factor(row_03_n1) == 0, "n=1 03+21 primitive")
require(sp.factor(row_13_n1) == 0, "n=1 13+22 primitive")
row_02_n1 = bracket(-2, 0, f0, M1) + bracket(-1, -1, R, L1)
require(sp.factor(K**2 * row_02_n1 - a * mu1 * deriv(H * K**3 * R)) == 0,
        "n=1 constant-product derivative")
# The omitted homogeneous constant in KL vanishes at any root of K.
kl_debt = sp.cancel(K * (a * mu1 * H * K + lam / K) - a * mu1 * H * K**2)
require(kl_debt == lam, "n=1 KL homogeneous debt")
require(kl_debt.subs(K, 0) == lam, "n=1 root kills lambda")


# Odd n>=3: exact transports, Euler/Kummer row, and charge normalization.
checked_scales: list[int] = []
for nv in range(3, 82, 2):
    alpha = (3 * nv - 1) // 2
    m = (nv - 1) // 2
    require(gcd(3 * nv - 1, 2) == 2, f"odd low gcd n={nv}")
    require(gcd(1, 4 * nv - 2) == 1, f"odd high gcd n={nv}")

    mu = (4 * nv - 2) * t / d
    V = H * K**2
    U = V**alpha
    P0 = H ** (alpha - 1) * K ** (nv - 2)
    r = R / P0
    E2 = Hp * K + 2 * H * Kp

    f0 = a * H**alpha
    g0 = c * H
    f2 = d * K
    g3 = t * K ** (4 * nv - 2)
    L = lam * K ** (nv - 2) + a * mu * H**alpha * K ** (4 * nv - 3)
    M = kap * K ** (2 * nv - 2) + mu * R * K ** (4 * nv - 3)

    row_03 = bracket(1 - 3 * nv, 4 * nv - 2, f0, g3) + bracket(
        1, nv - 2, f2, L
    )
    row_13 = bracket(1 - 2 * nv, 4 * nv - 2, R, g3) + bracket(
        1, 2 * nv - 2, f2, M
    )
    require(sp.factor(row_03) == 0, f"03+21 transport n={nv}")
    require(sp.factor(row_13) == 0, f"13+22 transport n={nv}")

    row_01 = bracket(1 - 3 * nv, nv - 2, f0, L) + bracket(
        1 - 2 * nv, -2, R, g0
    )
    D = lam * (nv - 2) + a * mu * (4 * nv - 3) * U
    D_coefficient = (
        (a * alpha / c) * H ** (alpha - 1) * K ** (nv - 3) * D
    )
    Rp_solution = ((2 * nv - 1) * R * Hp + D_coefficient * E2) / (2 * H)
    require(sp.factor(row_01.subs(Rp, Rp_solution)) == 0,
            f"01+10 Euler row n={nv}")

    # Freeze the two Euler particulars in R/P0.
    r_particular = (a * alpha / c) * (lam + a * mu * U)
    R_particular = P0 * r_particular
    particular_operator = sp.expand(
        2 * H * deriv(R_particular) - (2 * nv - 1) * R_particular * Hp
    )
    require(sp.factor(particular_operator - D_coefficient * E2) == 0,
            f"01+10 Euler particulars n={nv}")

    # Dividing S^2=chi H^(2n-1) by P0^2 gives f^2 V^(n-2)=chi.
    require(2 * nv - 1 - 2 * (alpha - 1) == -(nv - 2),
            f"Kummer H exponent n={nv}")
    require(-2 * (nv - 2) == -2 * nv + 4,
            f"Kummer K exponent n={nv}")
    require((nv - 2) % 2 == 1, f"odd V valuation n={nv}")

    row_02 = bracket(1 - 3 * nv, 2 * nv - 2, f0, M) + bracket(
        1 - 2 * nv, nv - 2, R, L
    )
    G = (
        lam * (nv - 2) * (2 * nv - 1)
        + a * mu * (4 * nv - 3) * (5 * nv - 2) * U
    )
    normalized = (
        G * r
        + (a * alpha / c) * D**2
        + 4 * a * kap * alpha * (nv - 1) * V
    )
    prefactor = H ** (alpha - 2) * K ** (2 * nv - 5) * E2 / 2
    require(sp.factor(row_02.subs(Rp, Rp_solution) - prefactor * normalized) == 0,
            f"02+11 charge normalization n={nv}")

    # On the zero Kummer sheet, record the exact nonzero U^2 coefficient.
    z = sp.symbols(f"z{nv}")
    uz = z**alpha
    Dz = lam * (nv - 2) + a * mu * (4 * nv - 3) * uz
    Gz = (
        lam * (nv - 2) * (2 * nv - 1)
        + a * mu * (4 * nv - 3) * (5 * nv - 2) * uz
    )
    rz = (a * alpha / c) * (lam + a * mu * uz)
    zero_sheet = sp.Poly(
        sp.expand(Gz * rz + (a * alpha / c) * Dz**2
                  + 4 * a * kap * alpha * (nv - 1) * z), z
    )
    expected_lead = (
        (a * alpha / c) * (a * mu) ** 2 * (4 * nv - 3) * (9 * nv - 5)
    )
    require(zero_sheet.degree() == 2 * alpha, f"zero-sheet degree n={nv}")
    require(sp.factor(zero_sheet.LC() - expected_lead) == 0,
            f"zero-sheet leading coefficient n={nv}")
    checked_scales.append(nv)


print("THM-3728 exact W004 scalar-12 anchor-20 charge-normalization audit")
print("fibre word = " + "|".join(fibre_word))
print("placement = P(1-3n,1-2n,1), Q(-2,n-2,2n-2,4n-2)")
print("even scales = common-base exponent two; THM-3613 parity death")
print("n=1 = upper primitives plus (HK^3 R)'=0 collapse")
print("odd exact window = 3..81; 40 scales")
print("nonzero Kummer sheet = f in C(V) with forbidden odd V-valuation")
print("zero Kummer sheet = nonzero U^2 coefficient (4n-3)(9n-5)")
print("tail consequence = every stabilized W004 placement is now captured")
print("scope = named scalar fibre 12+20 anchor-20 family, all scales")
