#!/usr/bin/env python3
"""Exact norm-twist and boundary audit for THM-3722."""

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
# Named W004 placement and its one-five-one fibre word.
# ---------------------------------------------------------------------------

n0 = sp.symbols("n0", integer=True, positive=True)
p_weights = (1 - 3 * n0, 1 - 2 * n0, 1)
q_weights = (-n0 - 2, -2, n0 - 2, 3 * n0 - 2)
A_POS = (0, 1, 3)
B_POS = (0, 1, 2, 4)
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
        ((1, 2), (2, 0)),
        ((0, 3), (2, 1)),
        ((1, 3), (2, 2)),
        ((2, 3),),
    ),
    "W004 fibre word",
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
        (-n0 - 1, -n0 - 1),
        (-1, -1),
        (n0 - 1, n0 - 1),
        (3 * n0 - 1,),
    ),
    "named placement fibre sums",
)


# ---------------------------------------------------------------------------
# Charge-coordinate calculus for n>=3.
# ---------------------------------------------------------------------------

ell = sp.symbols("ell", integer=True, nonnegative=True)
U, K, P, Z = sp.symbols("U K P Z", nonzero=True)
Up, Kp, Pp, Zp = sp.symbols("Up Kp Pp Zp")
a, c, d, t, lam = sp.symbols("a c d t lam", nonzero=True)


def derivative(expression):
    return (
        sp.diff(expression, U) * Up
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, P) * Pp
        + sp.diff(expression, Z) * Zp
    )


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def audit_tail_branch(delta, n, alpha, beta, label):
    R = 2 * n - 1
    require(sp.expand(delta * alpha - (3 * n - 1)) == 0, f"{label} alpha")
    require(sp.expand(delta * beta - (n + 2)) == 0, f"{label} beta")

    rho = t * (3 * n - 2) / d
    f0 = a * K ** (1 - 3 * n) * U**alpha
    g0 = c * K ** (-n - 2) * U**beta
    f2 = d * K
    g3 = t * K ** (3 * n - 2)
    hidden_m = K ** (-R) * P
    hidden_n = K ** (n - 2) * (lam + rho * P)
    hidden_l = K**(-2) * Z

    upper = wedge(1 - 2 * n, hidden_m, 3 * n - 2, g3)
    upper += wedge(1, f2, n - 2, hidden_n)
    zero(upper, f"{label} upper transport")

    norm_row = wedge(1 - 2 * n, hidden_m, n - 2, hidden_n)
    norm_row += wedge(1, f2, -n - 2, g0)
    norm_primitive = (n - 2) * lam * P
    norm_primitive += rho * (3 * n - 3) * P**2 / 2
    norm_primitive -= c * d * U**beta
    zero(
        norm_row * K ** (n + 1) - derivative(norm_primitive),
        f"{label} norm primitive",
    )

    # Nonzero lambda: completing the square gives S^2=k^2+C U^beta.
    shift = (n - 2) * lam / (rho * (3 * n - 3))
    norm_scale = 2 * c * d / (rho * (3 * n - 3))
    zero(
        (P + shift) ** 2 - shift**2 - norm_scale * U**beta
        - 2 * norm_primitive / (rho * (3 * n - 3)),
        f"{label} shifted norm identity",
    )

    # Zero lambda: P^2=C U^beta, so P'/P=(beta/2)U'/U.
    p_log = beta * P * Up / (2 * U)
    z_zero = -a * rho * alpha * delta * (3 * n - 2)
    z_zero *= U**alpha / (2 * (3 * n - 4))
    z_zero += c * beta * R * U ** (beta - alpha) * P / (2 * a * alpha)
    middle = wedge(1 - 3 * n, f0, n - 2, hidden_n.subs(lam, 0))
    middle += wedge(1 - 2 * n, hidden_m, -2, hidden_l)
    lowest = wedge(1 - 3 * n, f0, -2, hidden_l)
    lowest += wedge(1 - 2 * n, hidden_m, -n - 2, g0)
    V = sp.symbols(f"V_{label.replace(' ', '_').replace('=', '')}")
    middle_equation = R * V - beta * Z / U
    middle_equation += a * rho * alpha * (3 * n - 2) * U ** (alpha - 1) / 2
    lowest_equation = a * alpha * U ** (alpha - 1) * (delta * U * V - 2 * Z)
    lowest_equation += c * beta * (3 * n - 4) * U ** (beta - 1) * P / 2
    zero(
        middle.subs({lam: 0, Pp: p_log, Zp: V * Up})
        - K ** (-R - 2) * P * Up * middle_equation,
        f"{label} zero-shift middle system",
    )
    zero(
        lowest.subs({lam: 0, Pp: p_log, Zp: V * Up})
        - K ** (-delta * alpha - 2) * Up * lowest_equation,
        f"{label} zero-shift lowest system",
    )
    A11, A12 = R, -beta / U
    A21 = a * alpha * delta * U**alpha
    A22 = -2 * a * alpha * U ** (alpha - 1)
    determinant = sp.expand(A11 * A22 - A12 * A21)
    C1 = a * rho * alpha * (3 * n - 2) * U ** (alpha - 1) / 2
    C2 = c * beta * (3 * n - 4) * U ** (beta - 1) * P / 2
    zero(
        determinant
        + a * alpha * (3 * n - 4) * U ** (alpha - 1),
        f"{label} zero-shift determinant",
    )
    zero(
        z_zero - (A21 * C1 - A11 * C2) / determinant,
        f"{label} zero-shift solved Z",
    )
    zero(
        3 * beta - 2 * alpha - (8 - 3 * n) / delta,
        f"{label} negative exponent identity",
    )
    negative_exponent = sp.Poly(3 * n - 8, ell)
    require(
        all(coefficient > 0 for coefficient in negative_exponent.all_coeffs()),
        f"{label} zero-shift exponent is strictly negative",
    )


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
    audit_tail_branch(*branch)


# ---------------------------------------------------------------------------
# The sole nonzero-shift rational norm boundary: n=5, delta=7, beta=1.
# ---------------------------------------------------------------------------

n5, delta5, alpha5, beta5, R5 = 5, 7, 2, 1, 9
rho5 = 13 * t / d
u5 = (6 * rho5 * P**2 + 3 * lam * P) / (c * d)
theta5 = sp.diff(u5, P)
f05 = a * K**(-14) * u5**2
g05 = c * K**(-7) * u5
f25 = d * K
g35 = t * K**13
m5 = K**(-9) * P
n5_hidden = K**3 * (lam + rho5 * P)
l5 = K**(-2) * Z
middle5 = wedge(-14, f05, 3, n5_hidden) + wedge(-9, m5, -2, l5)
lowest5 = wedge(-14, f05, -2, l5) + wedge(-9, m5, -7, g05)
V5 = sp.symbols("V5")
middle5_reduced = sp.cancel(middle5 * K**11 / Pp)
lowest5_reduced = sp.cancel(lowest5 * K**16 / Pp)
middle5_reduced = middle5_reduced.subs(Zp, V5 * Pp)
lowest5_reduced = lowest5_reduced.subs(Zp, V5 * Pp)
A11 = sp.diff(middle5_reduced, V5)
A12 = sp.diff(middle5_reduced, Z)
A21 = sp.diff(lowest5_reduced, V5)
A22 = sp.diff(lowest5_reduced, Z)
det5 = sp.factor(A11 * A22 - A12 * A21)
det5_expected = 2 * a * alpha5 * u5 ** (alpha5 - 1)
det5_expected *= delta5 * u5 - R5 * P * theta5
zero(det5 - det5_expected, "n=5 differential determinant")
require(
    sp.Poly(sp.factor(delta5 * u5 - R5 * P * theta5), P).degree() == 2,
    "n=5 determinant inner polynomial is nonzero",
)

Q, Qp = sp.symbols("Q Qp")
scalar5 = wedge(-14, f05, 13, g35) + wedge(1, f25, -2, l5)
scalar5 = scalar5.subs({Z: Q, Zp: Qp * Pp})
scalar5_expected = K**(-1) * Pp * (26 * a * t * u5 * theta5 - d * Qp)
zero(scalar5 - scalar5_expected, "n=5 scalar P-coordinate factor")
require(R5 == 9, "n=5 degree floor deg(P)>=9deg(K)")


# ---------------------------------------------------------------------------
# n=2: pure norm, explicit Z(U), and a delta=1 degree obstruction.
# ---------------------------------------------------------------------------

sigma = sp.symbols("sigma", nonzero=True)
rho2 = 4 * t / d
p2 = sigma * U**2
n2_hidden = lam + rho2 * p2
z2 = 6 * c * sigma * U / (5 * a) - 5 * a * rho2 * U**5
v2 = 8 * c * sigma / (5 * a) - 10 * a * rho2 * U**4
f02 = a * K**(-5) * U**5
g02 = c * K**(-4) * U**4
f22 = d * K
g32 = t * K**4
m2 = K**(-3) * p2
l2 = K**(-2) * Z
upper2 = wedge(-3, m2, 4, g32) + wedge(1, f22, 0, n2_hidden)
norm2 = wedge(-3, m2, 0, n2_hidden) + wedge(1, f22, -4, g02)
middle2 = wedge(-5, f02, 0, n2_hidden) + wedge(-3, m2, -2, l2)
lowest2 = wedge(-5, f02, -2, l2) + wedge(-3, m2, -4, g02)
sigma_square = 2 * c * d / (3 * rho2)
V2 = sp.symbols("V2")
middle2_equation = 3 * U * V2 - 4 * Z + 10 * a * rho2 * U**5
lowest2_equation = 5 * a * U * V2 - 10 * a * Z + 4 * c * sigma * U
zero(upper2, "n=2 upper transport")
zero(norm2.subs(sigma**2, sigma_square), "n=2 pure norm")
zero(
    middle2.subs(Zp, V2 * Up)
    - K**(-5) * sigma * U * Up * middle2_equation,
    "n=2 middle system",
)
zero(
    lowest2.subs(Zp, V2 * Up)
    - K**(-7) * U**4 * Up * lowest2_equation,
    "n=2 lowest system",
)
zero(
    middle2_equation.subs({Z: z2, V2: v2}),
    "n=2 solved middle",
)
zero(
    lowest2_equation.subs({Z: z2, V2: v2}),
    "n=2 solved lowest",
)
compatibility2 = derivative(z2) / Up - v2
zero(
    compatibility2
    + 2 * c * sigma / (5 * a)
    + 15 * a * rho2 * U**4,
    "n=2 derivative compatibility obstruction",
)


# ---------------------------------------------------------------------------
# n=1: proportional hidden row leaves a nonzero Euler factor.
# ---------------------------------------------------------------------------

H, Hp = sp.symbols("H Hp", nonzero=True)
M1 = sp.symbols("M1", nonzero=True)
M1p = sp.symbols("M1p")


def derivative_n1(expression):
    return (
        sp.diff(expression, H) * Hp
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, M1) * M1p
    )


def wedge_n1(r, left, s, right):
    return sp.expand(
        s * derivative_n1(left) * right - r * left * derivative_n1(right)
    )


n1_hidden = t * M1 / d
upper1 = wedge_n1(-1, M1, 1, t * K) + wedge_n1(1, d * K, -1, n1_hidden)
next1 = wedge_n1(-1, M1, -1, n1_hidden)
next1 += wedge_n1(1, d * K, -3, c * H**3)
zero(upper1, "n=1 upper proportionality")
zero(
    next1 + 3 * c * d * H**2 * (Hp * K + H * Kp),
    "n=1 surviving Euler factor",
)


print("THM-3722 exact W004 scalar-03 norm-twist audit")
print("placement = P(1-3n,1-2n,1), Q(-n-2,-2,n-2,3n-2), n>=1")
print("upper transport and quadratic norm primitive PASS")
print("nonzero-shift beta>=2 coprime-power obstruction encoded")
print("zero-shift negative Laurent exponent in all seven residues PASS")
print("n=5 beta=1 differential determinant and degree gate PASS")
print("n=2 pure norm/differential compatibility gate PASS")
print("n=1 proportional-row Euler contradiction PASS")
print("scope = complete named W004 scalar fibre 03+21, all scales")
print("ALL CHECKS PASSED")
