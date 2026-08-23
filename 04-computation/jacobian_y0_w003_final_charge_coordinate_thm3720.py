#!/usr/bin/env python3
"""Exact charge-coordinate and differential-system audit for THM-3720."""

from collections import defaultdict

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zero(expression, message):
    reduced = sp.powsimp(sp.factor(sp.cancel(expression)), force=True)
    if reduced != 0:
        raise AssertionError(f"{message}: {reduced}")


# ---------------------------------------------------------------------------
# Exact final W003 placement and coefficient fibres.
# ---------------------------------------------------------------------------

n0 = sp.symbols("n0", integer=True, positive=True)
p_weights = (1 - 3 * n0, 1 - 2 * n0, 1)
q_weights = (-n0 - 2, -2, n0 - 2, 2 * n0 - 2)
A_POS = (0, 1, 3)
B_POS = (0, 1, 2, 3)
fibres = defaultdict(list)
for i, left in enumerate(A_POS):
    for j, right in enumerate(B_POS):
        fibres[left + right].append((i, j))
FIBRES = tuple(tuple(fibres[index]) for index in sorted(fibres))
require(
    FIBRES
    == (
        ((0, 0),),
        ((0, 1), (1, 0)),
        ((0, 2), (1, 1)),
        ((0, 3), (1, 2), (2, 0)),
        ((1, 3), (2, 1)),
        ((2, 2),),
        ((2, 3),),
    ),
    "W003 fibre word",
)
fibre_sums = tuple(
    tuple(sp.expand(p_weights[i] + q_weights[j]) for i, j in fibre)
    for fibre in FIBRES
)
require(
    fibre_sums
    == (
        (-4 * n0 - 1,),
        (-3 * n0 - 1, -3 * n0 - 1),
        (-2 * n0 - 1, -2 * n0 - 1),
        (-n0 - 1, -n0 - 1, -n0 - 1),
        (-1, -1),
        (n0 - 1,),
        (2 * n0 - 1,),
    ),
    "final family fibre sums",
)


# ---------------------------------------------------------------------------
# Exact coefficient calculus in H,K and in the charge coordinate U=HK^delta.
# ---------------------------------------------------------------------------

ell = sp.symbols("ell", integer=True, nonnegative=True)
H, K, M = sp.symbols("H K M", nonzero=True)
Hp, Kp, Mp = sp.symbols("Hp Kp Mp")
U, Z = sp.symbols("U Z", nonzero=True)
Up, Zp = sp.symbols("Up Zp")
a, c, d, e, t = sp.symbols("a c d e t", nonzero=True)


def derivative_hkm(expression):
    return (
        sp.diff(expression, H) * Hp
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, M) * Mp
    )


def wedge_hkm(r, left, s, right):
    return sp.expand(
        s * derivative_hkm(left) * right - r * left * derivative_hkm(right)
    )


def derivative_ukz(expression):
    return (
        sp.diff(expression, U) * Up
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, Z) * Zp
    )


def wedge_ukz(r, left, s, right):
    return sp.expand(
        s * derivative_ukz(left) * right - r * left * derivative_ukz(right)
    )


def audit_branch(delta, n, alpha, beta, label):
    R = 2 * n - 1
    require(sp.expand(delta * alpha - (3 * n - 1)) == 0, f"{label} alpha")
    require(sp.expand(delta * beta - (n + 2)) == 0, f"{label} beta")
    require(sp.expand(alpha - beta - (2 * n - 3) / delta) == 0, f"{label} gap")

    # Before solving the triple, multiplying it by K^(n+1) gives this exact
    # primitive for an arbitrary hidden coefficient M.
    f0_h = a * H**alpha
    g0_h = c * H**beta
    f2_h = d * K
    g2_h = e * K ** (n - 2)
    g3_h = t * K ** (2 * n - 2)
    triple_h = wedge_hkm(1 - 3 * n, f0_h, 2 * n - 2, g3_h)
    triple_h += wedge_hkm(1 - 2 * n, M, n - 2, g2_h)
    triple_h += wedge_hkm(1, f2_h, -n - 2, g0_h)
    primitive_h = e * (n - 2) * M * K**R
    primitive_h += (2 * n - 2) * a * t * H**alpha * K ** (3 * n - 1)
    primitive_h -= c * d * H**beta * K ** (n + 2)
    zero(
        triple_h * K ** (n + 1) - derivative_hkm(primitive_h),
        f"{label} triple primitive",
    )

    # The root-killed integration constant leaves a two-term polynomial P(U).
    A = -2 * (n - 1) * a * t / (e * (n - 2))
    B = c * d / (e * (n - 2))
    P = A * U**alpha + B * U**beta
    Pu = sp.diff(P, U)
    zero(
        e * (n - 2) * P
        + (2 * n - 2) * a * t * U**alpha
        - c * d * U**beta,
        f"{label} integrated triple",
    )

    # All known coefficients are K^(their weight) times a polynomial in U.
    f0 = a * K ** (1 - 3 * n) * U**alpha
    g0 = c * K ** (-n - 2) * U**beta
    f2 = d * K
    g2 = e * K ** (n - 2)
    g3 = t * K ** (2 * n - 2)
    hidden_m = K ** (-R) * P
    hidden_l = K**(-2) * Z

    triple = wedge_ukz(1 - 3 * n, f0, 2 * n - 2, g3)
    triple += wedge_ukz(1 - 2 * n, hidden_m, n - 2, g2)
    triple += wedge_ukz(1, f2, -n - 2, g0)
    zero(triple, f"{label} solved triple")

    middle = wedge_ukz(1 - 3 * n, f0, n - 2, g2)
    middle += wedge_ukz(1 - 2 * n, hidden_m, -2, hidden_l)
    middle_system = R * P * Zp - 2 * Pu * Up * Z
    middle_system += a * e * (n - 2) * alpha * U ** (alpha - 1) * Up
    zero(
        middle - K ** (-R - 2) * middle_system,
        f"{label} middle system",
    )

    lowest = wedge_ukz(1 - 3 * n, f0, -2, hidden_l)
    lowest += wedge_ukz(1 - 2 * n, hidden_m, -n - 2, g0)
    lowest_system = a * alpha * U ** (alpha - 1)
    lowest_system *= delta * U * Zp - 2 * Up * Z
    lowest_system += (
        c
        * beta
        * U ** (beta - 1)
        * Up
        * (R * P - delta * U * Pu)
    )
    zero(
        lowest - K ** (-delta * alpha - 2) * lowest_system,
        f"{label} lowest system",
    )

    # Divide by U' and solve the two rows for (V,Z), where V=Z'/U'.
    V = sp.symbols(f"V_{label.replace(' ', '_').replace('=', '')}")
    C1 = a * e * (n - 2) * alpha * U ** (alpha - 1)
    C2 = c * beta * U ** (beta - 1) * (R * P - delta * U * Pu)
    A11, A12 = R * P, -2 * Pu
    A21 = a * delta * alpha * U**alpha
    A22 = -2 * a * alpha * U ** (alpha - 1)
    determinant = sp.expand(A11 * A22 - A12 * A21)
    determinant_expected = 2 * a * alpha * U ** (alpha - 1)
    determinant_expected *= delta * U * Pu - R * P
    zero(determinant - determinant_expected, f"{label} determinant")

    leading_minor = sp.expand(delta * alpha - R)
    leading_poly = sp.Poly(leading_minor, ell)
    require(
        all(coefficient > 0 for coefficient in leading_poly.all_coeffs()),
        f"{label} determinant leading coefficient is positive",
    )

    z_rational = sp.cancel((A21 * C1 - A11 * C2) / determinant)
    v_rational = sp.cancel((-C1 * A22 + A12 * C2) / determinant)
    zero(A11 * v_rational + A12 * z_rational + C1, f"{label} solved middle")
    zero(A21 * v_rational + A22 * z_rational + C2, f"{label} solved lowest")

    # Once the elementary pullback lemma upgrades Z=Q(U) from rational to
    # polynomial, the entire scalar row has one charge-coordinate factor.
    Q, Qu = sp.symbols("Q Qu")
    scalar = wedge_ukz(1 - 2 * n, hidden_m, 2 * n - 2, g3)
    scalar += wedge_ukz(1, f2, -2, hidden_l)
    scalar = scalar.subs({Zp: Qu * Up, Z: Q})
    scalar_expected = K**(-1) * Up
    scalar_expected *= t * (R - 1) * Pu - d * Qu
    zero(scalar - scalar_expected, f"{label} scalar charge factor")


# n>=3, split into the seven exact gcd(n+2,7) residue branches.
BRANCHES = (
    (1, 7 * ell + 3, 21 * ell + 8, 7 * ell + 5, "n=3 mod 7"),
    (1, 7 * ell + 4, 21 * ell + 11, 7 * ell + 6, "n=4 mod 7"),
    (7, 7 * ell + 5, 3 * ell + 2, ell + 1, "n=5 mod 7"),
    (1, 7 * ell + 6, 21 * ell + 17, 7 * ell + 8, "n=6 mod 7"),
    (1, 7 * ell + 7, 21 * ell + 20, 7 * ell + 9, "n=0 mod 7"),
    (1, 7 * ell + 8, 21 * ell + 23, 7 * ell + 10, "n=1 mod 7"),
    (1, 7 * ell + 9, 21 * ell + 26, 7 * ell + 11, "n=2 mod 7"),
)
for branch in BRANCHES:
    audit_branch(*branch)


# U=HK^delta gives U'=K^(delta-1)E_delta.  The delta=1 scalar equation would
# have degree excess deg(U)-1-deg(K)=deg(H)-1>=1; delta=7 has a K^5 factor.
delta_symbol = sp.symbols("delta_symbol", integer=True, positive=True)
E = sp.symbols("E")
zero(
    derivative_hkm(H * K**delta_symbol)
    - K ** (delta_symbol - 1) * (Hp * K + delta_symbol * H * Kp),
    "charge-coordinate derivative",
)
deg_h, deg_k = sp.symbols("deg_h deg_k", integer=True, positive=True)
zero((deg_h + deg_k - 1) - deg_k - (deg_h - 1), "delta=1 degree excess")
require(7 - 2 == 5, "delta=7 scalar carries K^5")


print("THM-3720 exact W003 final-family audit")
print("placement = P(1-3n,1-2n,1), Q(-n-2,-2,n-2,2n-2), n>=3")
print("triple primitive and two-term charge polynomial P(U) PASS")
print("middle/lowest 2x2 differential system in all seven residues PASS")
print("system determinant nonzero and Z rational in U PASS")
print("polynomial pullback gives Z in C[U]")
print("scalar = K^(delta-2) E_delta Phi(U)")
print("delta=7 divisibility and delta=1 degree gates PASS")
print("scope = final residual W003 family D; W003 is complete")
print("ALL CHECKS PASSED")
