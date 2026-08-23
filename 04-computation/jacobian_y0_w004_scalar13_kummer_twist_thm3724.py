#!/usr/bin/env python3
"""Exact companion for THM-3724 (W004 scalar-13 Kummer twist)."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition: object, label: str) -> None:
    if condition is not True and condition != sp.S.true:
        raise ArithmeticError(label)


# ---------------------------------------------------------------------------
# W004 support/fibre and absolute-weight audit.
# ---------------------------------------------------------------------------

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
q_weights = (-2 * n - 2, -n - 2, -2, 2 * n - 2)
require(sp.expand(p_weights[1] + q_weights[3]) == -1, "scalar address 13")
require(sp.expand(p_weights[2] + q_weights[2]) == -1, "scalar address 22")
require(sp.expand(p_weights[2] + q_weights[3]) == 2 * n - 1, "top singleton")

delta_residues = tuple(gcd(r - 3, 8) for r in range(8))
require(delta_residues == (1, 2, 1, 8, 1, 2, 1, 4), "delta residue table")
for nv in range(2, 258):
    dv = gcd(nv - 3, 8)
    require(dv == gcd(3 * nv - 1, 2 * nv + 2), f"singleton gcd n={nv}")
    require((3 * nv - 1) % dv == 0, f"alpha integrality n={nv}")
    require((2 * nv + 2) % dv == 0, f"beta integrality n={nv}")
    require((4 * nv - 4) % dv == 0, f"m integrality n={nv}")
    require(gcd(dv, 2 * nv - 1) == 1, f"Kummer coprimality n={nv}")


# ---------------------------------------------------------------------------
# End-row transport identities.
# ---------------------------------------------------------------------------

delta = sp.symbols("delta", integer=True, positive=True)
alpha = (3 * n - 1) / delta
beta = (2 * n + 2) / delta
R = 2 * n - 1
T = 2 * n - 2
C = n + 2
m = 4 * (n - 1) / delta

a, c, d, t, rho = sp.symbols("a c d t rho", nonzero=True)
H, K, Hp, Kp = sp.symbols("H K Hp Kp", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.diff(expr, H) * Hp + sp.diff(expr, K) * Kp


def bracket(r: sp.Expr, s: sp.Expr, f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(s * deriv(f) * g - r * f * deriv(g))


f0 = a * H**alpha
f2 = d * K
g3 = t * K**T
L = rho * H**alpha * K ** (R - 2)
upper_row = bracket(-(3 * n - 1), T, f0, g3) + bracket(1, -C, f2, L)
require(sp.factor(upper_row.subs(rho, a * t * T / d)) == 0, "upper transport")

require(sp.simplify(delta * alpha - C - (R - 2)) == 0, "low Hp coefficient")
require(sp.simplify(delta * m - R - (R - 2)) == 0, "low resonance denominator")
require(sp.simplify(2 * alpha - beta - m) == 0, "low H exponent")
require(sp.simplify(delta * m - R - 1 - (R - 3)) == 0, "low K exponent")

# Directly reconstruct the two middle source rows across two full residue
# periods.  This is a sign/orientation control; the symbolic Cramer algebra
# below is the all-n identity.
Pi0, Pi0p, Zs, Zsp = sp.symbols("Pi0 Pi0p Zs Zsp", nonzero=True)
for nv in range(2, 18):
    dv = gcd(nv - 3, 8)
    av = (3 * nv - 1) // dv
    bv = (2 * nv + 2) // dv
    rv = 2 * nv - 1
    cv = nv + 2
    Uv = H * K**dv
    Uvp = Hp * K**dv + dv * H * K ** (dv - 1) * Kp
    f0v = a * H**av
    f0vp = a * av * H ** (av - 1) * Hp
    Mv = K ** (-rv) * Pi0
    Mvp = K ** (-rv) * (Pi0p - rv * Pi0 * Kp / K)
    Lv = rho * H**av * K ** (rv - 2)
    Lvp = rho * (
        av * H ** (av - 1) * Hp * K ** (rv - 2)
        + (rv - 2) * H**av * K ** (rv - 3) * Kp
    )
    Nv = K**-2 * Zs
    Nvp = K**-2 * (Zsp - 2 * Zs * Kp / K)
    g0v = c * H**bv
    g0vp = c * bv * H ** (bv - 1) * Hp

    row_12 = (
        -2 * Mvp * Nv
        + rv * Mv * Nvp
        - dv * bv * d * Kp * g0v
        - d * K * g0vp
    )
    reduced_12 = K ** (-rv - 2) * (
        rv * Pi0 * Zsp
        - 2 * Pi0p * Zs
        - c * d * bv * Uv ** (bv - 1) * Uvp
    )
    require(sp.factor(row_12 - reduced_12) == 0, f"middle row 12+20 n={nv}")

    row_02 = (
        -2 * f0vp * Nv
        + dv * av * f0v * Nvp
        - cv * Mvp * Lv
        + rv * Mv * Lvp
    )
    reduced_02 = av * Uv ** (av - 1) * K ** (-3 * nv - 1) * (
        a * (dv * Uv * Zsp - 2 * Uvp * Zs)
        + rho * (rv * Pi0 * Uvp - sp.Rational(cv, av) * Uv * Pi0p)
    )
    require(sp.factor(row_02 - reduced_02) == 0, f"middle row 02+11 n={nv}")


# ---------------------------------------------------------------------------
# Zero Kummer sheet: exact Cramer solution and derivative incompatibility.
# ---------------------------------------------------------------------------

U = sp.symbols("U", positive=True)
A0 = a * rho * alpha / (c * beta)
e = beta - m

z_m = -R * rho * A0 * (n - 3) / (2 * a * alpha * delta)
z_e = -delta * c * d * beta / (2 * A0 * (R - 2))
y_m = -rho * A0 * m * (n - 3) / (a * alpha * delta)
y_e = -c * d * beta / (A0 * (R - 2))
Z0 = z_m * U**m + z_e * U**e
Y0 = y_m * U ** (m - 1) + y_e * U ** (e - 1)

zero_eq_12 = sp.factor(
    R * A0 * U**m * Y0
    - 2 * A0 * m * U ** (m - 1) * Z0
    - c * d * beta * U ** (beta - 1)
)
zero_eq_02 = sp.factor(
    a * alpha * (delta * U * Y0 - 2 * Z0)
    + rho * A0 * (R * alpha - C * m) * U**m
)
require(sp.powsimp(zero_eq_12, force=True) == 0, "zero-sheet row 12+20")
require(sp.powsimp(zero_eq_02, force=True) == 0, "zero-sheet row 02+11")

zero_claim = (
    c * d * beta * (n - 2) / (A0 * (R - 2)) * U ** (e - 1)
    - rho * A0 * m * (n - 3) * (R - 2) / (2 * a * alpha * delta)
    * U ** (m - 1)
)
require(
    sp.factor(sp.powsimp(sp.diff(Z0, U) - Y0 - zero_claim, force=True)) == 0,
    "zero-sheet compatibility",
)
require(sp.simplify(e - m + 2 * (3 * n - 5) / delta) == 0, "zero exponents distinct")
zero_coeff_1 = c * d * beta * (n - 2) / (A0 * (R - 2))
zero_coeff_2 = -rho * A0 * m * (n - 3) * (R - 2) / (2 * a * alpha * delta)
require(sp.simplify(zero_coeff_1.subs({n: 2, delta: 1})) == 0, "zero n=2 first exit")
require(sp.simplify(zero_coeff_2.subs({n: 2, delta: 1})) != 0, "zero n=2 survivor")
require(sp.simplify(zero_coeff_1.subs({n: 3, delta: 8})) != 0, "zero n=3 survivor")
require(sp.simplify(zero_coeff_2.subs({n: 3, delta: 8})) == 0, "zero n=3 second exit")


# ---------------------------------------------------------------------------
# Nonzero Kummer sheet: exact determinant and compatibility numerator.
# ---------------------------------------------------------------------------

v, B = sp.symbols("v B", positive=True, nonzero=True)
Mexp = 4 * n - 4
Pi = A0 * v**Mexp + B * v**R
Pi_v = sp.diff(Pi, v)

rhs_12 = c * d * beta * delta * v ** (2 * n + 1)
rhs_02 = -rho * (R * alpha * delta * Pi - C * v * Pi_v)
determinant = 2 * a * alpha * delta * (v * Pi_v - R * Pi)
determinant_claim = 2 * a * alpha * delta * A0 * (2 * n - 3) * v**Mexp
require(sp.factor(determinant - determinant_claim) == 0, "Kummer determinant")

# Cramer's rule for unknowns (Y,Z).
Z1 = sp.cancel((R * Pi * rhs_02 - a * alpha * delta * v * rhs_12) / determinant)
Y1 = sp.cancel(((-2 * a * alpha * delta) * rhs_12 + 2 * Pi_v * rhs_02) / determinant)

compat_numerator = (
    -3 * B * c * rho**2 * (n - 2) * (n + 1) * (2 * n - 3) ** 2
    * (2 * n - 1) * v ** (4 * n + 3)
    - 2 * a * rho**3 * (n - 3) * (n - 1) * (2 * n - 3) ** 2
    * (3 * n - 1) * v ** (6 * n)
    + 8 * c**3 * d * (n - 2) * (n + 1) ** 3 * v**10
)
compat_denominator = (
    2 * a * c * rho * (n + 1) * (2 * n - 3) * (3 * n - 1)
    * v ** (2 * n + 5)
)
require(
    sp.factor(sp.together(sp.diff(Z1, v) - Y1 - compat_numerator / compat_denominator))
    == 0,
    "nonzero-sheet compatibility",
)

require(sp.factor(compat_numerator.subs(n, 2)) != 0, "nonzero sheet n=2")
require(sp.factor(compat_numerator.subs(n, 3)) != 0, "nonzero sheet n=3")
power_differences = (
    sp.factor((4 * n + 3) - 6 * n),
    sp.factor((4 * n + 3) - 10),
    sp.factor(6 * n - 10),
)
expected_power_differences = (3 - 2 * n, 4 * n - 7, 2 * (3 * n - 5))
require(
    all(sp.simplify(x - y) == 0 for x, y in zip(power_differences, expected_power_differences)),
    "power gaps",
)


# n=1: the singleton weights are (1,0), and W_(1,0)(f,g)=-f g'.
require((sp.sympify(p_weights[2]).subs(n, 1), q_weights[3].subs(n, 1)) == (1, 0), "n=1 weights")
fp, gp, f, g = sp.symbols("fp gp f g")
require(sp.expand(0 * fp * g - 1 * f * gp) == -f * gp, "n=1 scalar deletion")


print("THM-3724 exact W004 scalar-13 Kummer-twist audit")
print("fibre word = " + "|".join(fibre_word))
print("weights = P(1-3n,1-2n,1), Q(-2n-2,-n-2,-2,2n-2)")
print("delta residues n mod 8 = " + ",".join(map(str, delta_residues)))
print("end transports = upper primitive + Kummer power law verified")
print("zero sheet = two-power derivative incompatibility verified")
print("nonzero sheet = three-power compatibility numerator verified")
print("boundary = n=2,n=3 coefficient exits and n=1 deletion verified")
print("scope = complete named W004 scalar fibre 13+22, all scales")
