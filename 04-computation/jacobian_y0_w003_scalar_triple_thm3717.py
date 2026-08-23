#!/usr/bin/env python3
"""Exact square-law and Euler audit for THM-3717."""

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
# The exact W003 placement and its seven coefficient fibres.
# ---------------------------------------------------------------------------

n0 = sp.symbols("n0", integer=True, positive=True)
p_weights = (-n0 - 2, -2, 2 * n0 - 2)
q_weights = (1 - 2 * n0, 1 - n0, 1, n0 + 1)
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
        (-3 * n0 - 1,),
        (-2 * n0 - 1, -2 * n0 - 1),
        (-n0 - 1, -n0 - 1),
        (-1, -1, -1),
        (n0 - 1, n0 - 1),
        (2 * n0 - 1,),
        (3 * n0 - 1,),
    ),
    "family C fibre sums",
)


# ---------------------------------------------------------------------------
# Formal coefficient calculus.  W_(r,s)(F,G)=s F'G-rFG'.
# ---------------------------------------------------------------------------

ell = sp.symbols("ell", integer=True, positive=True)
H, K, Y = sp.symbols("H K Y", nonzero=True)
Hp, Kp, Yp = sp.symbols("Hp Kp Yp")
a, c, d, e, t = sp.symbols("a c d e t", nonzero=True)


def derivative(expression):
    return (
        sp.diff(expression, H) * Hp
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, Y) * Yp
    )


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def audit_generic(delta, n, alpha, beta, label):
    """Audit n>=4 after the upper double has killed its integration constant."""

    f0 = a * H**alpha
    g0 = c * H**beta
    f2 = d * K ** (2 * n - 2)
    g2 = e * K
    g3 = t * K ** (n + 1)
    rho = (n + 1) * t / (2 * (n - 1) * d)
    M = K ** (n - 3) * Y
    L = rho * Y
    T = K ** (n - 1) * Y

    upper = wedge(-2, M, n + 1, g3)
    upper += wedge(2 * n - 2, f2, 1 - n, L)
    zero(upper, f"{label} upper double")

    middle = wedge(-n - 2, f0, 1, g2)
    middle += wedge(-2, M, 1 - n, L)
    middle_primitive = a * e * H**alpha * K ** (n + 2)
    middle_primitive += rho * (3 - n) * T**2 / 2
    zero(
        middle * K ** (n + 1) - derivative(middle_primitive),
        f"{label} square primitive",
    )

    square_yp = Y * (
        alpha * Hp / (2 * H) + (4 - n) * Kp / (2 * K)
    )
    square_log = alpha * (Hp / H + delta * Kp / K) / 2
    zero(
        derivative(T).subs(Yp, square_yp) - T * square_log,
        f"{label} square logarithmic derivative",
    )

    lower = wedge(-n - 2, f0, 1 - n, L)
    lower += wedge(-2, M, 1 - 2 * n, g0)
    lower = lower.subs(Yp, square_yp)
    euler = Hp * K + delta * H * Kp
    expected = euler * T / 2
    expected *= (
        a
        * rho
        * alpha
        * (4 - n)
        * H ** (alpha - 1)
        * K ** (-n)
        - c
        * beta
        * (n - 2)
        * H ** (beta - 1)
        * K ** (-3)
    )
    zero(lower - expected, f"{label} lower Euler factor")
    zero(beta - alpha - (n - 3) / delta, f"{label} exponent gap")


# The nonexceptional residue branches.  The delta=5 branch is n=3 mod 5;
# n=3 itself is handled below, so ell>=1 here starts at n=8.
BRANCHES = (
    (1, 5 * ell, 5 * ell + 2, 10 * ell - 1, "n=0 mod 5"),
    (1, 5 * ell + 1, 5 * ell + 3, 10 * ell + 1, "n=1 mod 5"),
    (1, 5 * ell + 2, 5 * ell + 4, 10 * ell + 3, "n=2 mod 5"),
    (5, 5 * ell + 3, ell + 1, 2 * ell + 1, "n=3 mod 5"),
    (1, 5 * ell + 4, 5 * ell + 6, 10 * ell + 7, "n=4 mod 5"),
)
for branch in BRANCHES:
    audit_generic(*branch)

# n=4 is the coefficient-degenerate boundary of the lower Euler equation.
audit_generic(1, sp.Integer(4), sp.Integer(6), sp.Integer(7), "n=4")


# ---------------------------------------------------------------------------
# n=3: the hidden middle bracket vanishes, leaving a nonzero Euler factor.
# ---------------------------------------------------------------------------

n3 = sp.Integer(3)
rho3 = t / d
f03 = a * H
M3 = Y
L3 = rho3 * Y
g23 = e * K
g33 = t * K**4
f23 = d * K**4
n3_upper = wedge(-2, M3, 4, g33) + wedge(4, f23, -2, L3)
n3_middle_hidden = wedge(-2, M3, -2, L3)
n3_middle_outer = wedge(-5, f03, 1, g23)
zero(n3_upper, "n=3 upper double")
zero(n3_middle_hidden, "n=3 hidden middle bracket")
zero(n3_middle_outer - a * e * (Hp * K + 5 * H * Kp), "n=3 outer Euler")


# ---------------------------------------------------------------------------
# n=2: the square law makes the negative hidden coefficient an H^2 arm.
# ---------------------------------------------------------------------------

sigma = sp.symbols("sigma", nonzero=True)
n2 = sp.Integer(2)
rho2 = 3 * t / (2 * d)
f02 = a * H**4
g02 = c * H**3
f22 = d * K**2
g22 = e * K
g32 = t * K**3
M2 = sigma * H**2
L2 = rho2 * sigma * H**2 * K
S2 = K**2 * M2
n2_upper = wedge(-2, M2, 3, g32) + wedge(2, f22, -1, L2)
n2_middle = wedge(-4, f02, 1, g22) + wedge(-2, M2, -1, L2)
n2_primitive = a * e * H**4 * K**4 + rho2 * S2**2 / 2
zero(n2_upper, "n=2 upper double")
zero(
    n2_middle * K**3 - derivative(n2_primitive),
    "n=2 square primitive",
)
zero(
    n2_middle.subs(sigma**2, -2 * a * e / rho2),
    "n=2 square-law middle row",
)
n2_scalar = wedge(-4, f02, 3, g32)
n2_scalar += wedge(-2, M2, 1, g22)
n2_scalar += wedge(2, f22, -3, g02)
n2_arm_quotient = sp.cancel(n2_scalar / H)
n2_denominator = sp.denom(n2_arm_quotient)
require(
    not any(n2_denominator.has(variable) for variable in (H, K, Hp, Kp)),
    "n=2 scalar quotient is polynomial in source variables",
)
zero(n2_scalar - H * n2_arm_quotient, "n=2 scalar arm factor")


print("THM-3717 exact W003 scalar-triple audit")
print("placement = P(-n-2,-2,2n-2), Q(1-2n,1-n,1,n+1), n>=2")
print("generic upper double and square primitive identities PASS")
print("generic lower Euler reduction in all five residue branches PASS")
print("n=3 middle singleton contradiction PASS")
print("n=4 degenerate-coefficient contradiction PASS")
print("n=2 square-arm scalar divisibility PASS")
print("scope = complete residual W003 family C, scalar fibre 03+12+20")
print("ALL CHECKS PASSED")
