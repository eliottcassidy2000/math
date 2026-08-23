#!/usr/bin/env python3
"""Exact companion for THM-3727 (dual W004 ternary Kummer twist)."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition: object, label: str) -> None:
    if condition is not True and condition != sp.S.true:
        raise ArithmeticError(label)


# ---------------------------------------------------------------------------
# W004 fibre word and absolute weights.
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

p_weights = (-n - 2, -2, 2 * n - 2)
q_weights = (1 - 4 * n, 1 - 3 * n, 1 - 2 * n, 1)
require(sp.expand(p_weights[1] + q_weights[3]) == -1, "scalar address 13")
require(sp.expand(p_weights[2] + q_weights[2]) == -1, "scalar address 22")

delta_residues = tuple(gcd(r + 2, 9) for r in range(9))
require(delta_residues == (1, 3, 1, 1, 3, 1, 1, 9, 1), "delta residue table")
for nv in range(2, 291):
    dv = gcd(nv + 2, 9)
    require(dv == gcd(nv + 2, 4 * nv - 1), f"singleton gcd n={nv}")
    require((nv + 2) % dv == 0, f"alpha integrality n={nv}")
    require((4 * nv - 1) % dv == 0, f"beta integrality n={nv}")
    require((5 - 2 * nv) % dv == 0, f"charge integrality n={nv}")
    require(dv % 2 == 1 and gcd(dv, 2) == 1, f"ternary coprimality n={nv}")
    if nv >= 3:
        require((5 - 2 * nv) // dv < 0, f"negative charge n={nv}")


# ---------------------------------------------------------------------------
# End-row transport and Kummer-law identities.
# ---------------------------------------------------------------------------

delta = sp.symbols("delta", integer=True, positive=True)
alpha = (n + 2) / delta
beta = (4 * n - 1) / delta
T = 2 * n - 2
C = 3 * n - 1
S = 2 * n - 3
m = (5 - 2 * n) / delta

a, c, d, t, rho = sp.symbols("a c d t rho", nonzero=True)
H, K, Hp, Kp = sp.symbols("H K Hp Kp", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.diff(expr, H) * Hp + sp.diff(expr, K) * Kp


def bracket(r: sp.Expr, s: sp.Expr, f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(s * deriv(f) * g - r * f * deriv(g))


f0 = a * H**alpha
f2 = d * K**T
g3 = t * K
L = rho * H**alpha * K ** (-S)
upper_row = bracket(-(n + 2), 1, f0, g3) + bracket(T, -C, f2, L)
require(sp.factor(upper_row.subs(rho, a * t / (d * T))) == 0, "upper transport")

require(sp.simplify(delta * alpha - C + S) == 0, "upper K exponent")
require(sp.simplify(2 * alpha - beta - m) == 0, "low H exponent")
require(sp.simplify(delta * m - 2 + S) == 0, "low resonance denominator")

# Reconstruct the lowest source row across two residue periods.  This is a
# direct sign/exponent control; the displayed symbolic identities prove the
# all-n reduction.
Pi0, Pi0p = sp.symbols("Pi0 Pi0p", nonzero=True)
for nv in range(2, 20):
    dv = gcd(nv + 2, 9)
    av = (nv + 2) // dv
    bv = (4 * nv - 1) // dv
    cv = 3 * nv - 1
    sv = 2 * nv - 3
    mv = (5 - 2 * nv) // dv
    Uv = H * K**dv
    Uvp = Hp * K**dv + dv * H * K ** (dv - 1) * Kp
    f0v = a * H**av
    f0vp = a * av * H ** (av - 1) * Hp
    Lv = rho * H**av * K ** (-sv)
    Lvp = rho * (
        av * H ** (av - 1) * Hp * K ** (-sv)
        - sv * H**av * K ** (-sv - 1) * Kp
    )
    Mv = K**-2 * Pi0
    Mvp = K**-2 * (Pi0p - 2 * Pi0 * Kp / K)
    g0v = c * H**bv
    g0vp = c * bv * H ** (bv - 1) * Hp
    source_low = (
        -cv * f0vp * Lv
        + dv * av * f0v * Lvp
        - dv * bv * Mvp * g0v
        + 2 * Mv * g0vp
    )
    Dv = -a * rho * av * sv / (c * bv)
    reduced_low = -c * bv * H ** (bv - 1) * K ** (-dv - 2) * (
        dv * Uv * Pi0p
        - 2 * Pi0 * Uvp
        - Dv * H**mv * K ** (dv * mv) * Uvp
    )
    require(sp.factor(source_low - reduced_low) == 0, f"lowest row n={nv}")


# ---------------------------------------------------------------------------
# Exceptional n=2: exact middle rows and derivative compatibility.
# ---------------------------------------------------------------------------

U = sp.symbols("U", positive=True)
A0 = sp.symbols("A0", nonzero=True)
B = sp.symbols("B")
Pi = A0 * U + B * U**2
Pi_U = sp.diff(Pi, U)
Y, Z = sp.symbols("Y Z")

eq_12 = 2 * Pi * Y - 3 * Pi_U * Z - 14 * c * d * U**6
eq_02 = 4 * a * (U * Y - 3 * Z) + rho * (-5 * U * Pi_U + 8 * Pi)
matrix = sp.Matrix([[2 * Pi, -3 * Pi_U], [4 * a * U, -12 * a]])
require(sp.factor(matrix.det() + 12 * a * A0 * U) == 0, "n=2 determinant")

Z_claim = U * (
    3 * A0**2 * rho
    + A0 * B * rho * U
    - 2 * B**2 * rho * U**2
    + 28 * a * c * d * U**5
) / (6 * A0 * a)
Y_claim = (
    3 * A0**2 * rho
    + 4 * A0 * B * rho * U
    - 4 * B**2 * rho * U**2
    + 56 * a * c * d * U**5
) / (4 * A0 * a)
require(sp.factor(eq_12.subs({Y: Y_claim, Z: Z_claim})) == 0, "n=2 row 12+20")
require(sp.factor(eq_02.subs({Y: Y_claim, Z: Z_claim})) == 0, "n=2 row 02+11")

compatibility = (
    -3 * A0**2 * rho
    - 8 * A0 * B * rho * U
    + 168 * a * c * d * U**5
) / (12 * A0 * a)
require(sp.factor(sp.diff(Z_claim, U) - Y_claim - compatibility) == 0, "n=2 compatibility")
require(sp.Poly(sp.together(compatibility).as_numer_denom()[0], U).degree() == 5,
        "n=2 fifth-power survivor")
require(sp.Poly(sp.together(compatibility).as_numer_denom()[0], U).coeff_monomial(1) != 0,
        "n=2 constant survivor")

# Direct source-coordinate check of the two n=2 middle reductions.
Pis, Pisp, Zs, Zsp = sp.symbols("Pis Pisp Zs Zsp", nonzero=True)
Uv = H * K
Uvp = Hp * K + H * Kp
Mv = K**-2 * Pis
Mvp = K**-2 * (Pisp - 2 * Pis * Kp / K)
Nv = K**-3 * Zs
Nvp = K**-3 * (Zsp - 3 * Zs * Kp / K)
f0v = a * H**4
f0vp = 4 * a * H**3 * Hp
f2v = d * K**2
f2vp = 2 * d * K * Kp
g0v = c * H**7
g0vp = 7 * c * H**6 * Hp
Lv = rho * H**4 / K
Lvp = rho * (4 * H**3 * Hp / K - H**4 * Kp / K**2)

source_12 = -3 * Mvp * Nv + 2 * Mv * Nvp - 7 * f2vp * g0v - 2 * f2v * g0vp
reduced_12 = K**-5 * (2 * Pis * Zsp - 3 * Pisp * Zs - 14 * c * d * Uv**6 * Uvp)
require(sp.factor(source_12 - reduced_12) == 0, "direct n=2 row 12+20")

source_02 = -3 * f0vp * Nv + 4 * f0v * Nvp - 5 * Mvp * Lv + 2 * Mv * Lvp
reduced_02 = H**3 * K**-4 * (
    4 * a * (Uv * Zsp - 3 * Uvp * Zs)
    + rho * (-5 * Uv * Pisp + 8 * Pis * Uvp)
)
require(sp.factor(source_02 - reduced_02) == 0, "direct n=2 row 02+11")


# n=1 singleton: W_(0,1)(f,g)=f'g, so the weight-zero piece is scalar.
require((p_weights[2].subs(n, 1), sp.sympify(q_weights[3]).subs(n, 1)) == (0, 1),
        "n=1 singleton weights")
fp, gp, f, g = sp.symbols("fp gp f g")
require(sp.expand(1 * fp * g - 0 * f * gp) == fp * g, "n=1 scalar deletion")


print("THM-3727 exact dual W004 ternary-Kummer audit")
print("fibre word = " + "|".join(fibre_word))
print("weights = P(-n-2,-2,2n-2), Q(1-4n,1-3n,1-2n,1)")
print("delta residues n mod 9 = " + ",".join(map(str, delta_residues)))
print("end transports = upper primitive + negative Kummer charge verified")
print("generic scales = n>=3 Laurent nonpolynomiality gate verified")
print("exceptional scale = n=2 middle-row compatibility verified")
print("boundary = n=1 weight-zero deletion to 2x4 verified")
print("pairing = dyadic delta|8 versus ternary delta|9")
print("scope = complete named dual W004 scalar fibre 13+22, all scales")
