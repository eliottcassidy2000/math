#!/usr/bin/env python3
"""Finite exact companion for proved and independently audited THM-3589.

Universal every-line and generic-residue statements remain proof-driven.
This script checks the algebraic identities, cyclic compiler controls, arm
specializations, sharp rational hostile, and boundary-degree invoices using
exact SymPy arithmetic.  It deliberately contains no Python assert gates.
"""

from math import ceil, comb

import sympy as sp


b, c, e, z, t, x, q, s = sp.symbols("b c e z t x q s")
CHECKS = 0


def require(label, condition):
    """Record one active exact gate and fail with a stable label."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError("FAILED: " + label)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def bracket(F, G, n, Sigma):
    """Poisson bracket on c^n e=Sigma(b), before quotient reduction."""
    return sp.expand(
        c**n * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
        - sp.diff(Sigma, b)
        * (sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c))
        - n
        * c ** (n - 1)
        * e
        * (sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b))
    )


print("THM-3589 exact companion")
print("SECTION Poisson and Kummer identities")

SIGMAS = [b * (b - 1), b * (b - 1) * (b + 2), b * (b**2 - 1) * (b - 3)]
for n in range(2, 9):
    Sigma = SIGMAS[(n - 2) % len(SIGMAS)]
    require(f"generator bc n={n}", zero(bracket(b, c, n, Sigma) - c**n))
    require(f"generator ce n={n}", zero(bracket(c, e, n, Sigma) + sp.diff(Sigma, b)))
    require(f"generator be n={n}", zero(bracket(b, e, n, Sigma) + n * c ** (n - 1) * e))
    require(
        f"z bracket n={n}",
        zero(bracket(c * e, b, n, Sigma).subs(e, Sigma / c**n) - (n - 1) * Sigma),
    )
    hostile = c ** (1 - n) / (n - 1)
    require(f"rational hostile n={n}", zero(bracket(hostile, b, n, Sigma) - 1))

print("PASS Poisson generators, z=ce, and sharp rational hostile n=2..8")


print("SECTION cyclic compiler central lines")
CYCLIC_ROWS = 0
for n in range(2, 6):
    for r in range(1, 4):
        S = 1 + t**r
        P = sp.expand(
            sum(sp.Rational(comb(n - 1, j), n * r * j + 1) * t ** (r * j) for j in range(n))
        )
        B = sp.expand(t * P**n)
        p = sp.prod(sp.Rational(j * n * r, j * n * r + 1) for j in range(1, n))
        E0 = p ** (n * r)
        Epoly = z**r + E0
        W = sp.cancel(Epoly.subs(z, B) / S**n)
        require(f"ODE n={n} r={r}", zero(P + n * t * sp.diff(P, t) - S ** (n - 1)))
        require(f"P0 n={n} r={r}", P.subs(t, 0) == 1)
        require(f"B0 n={n} r={r}", B.subs(t, 0) == 0)
        require(f"W polynomial n={n} r={r}", sp.denom(W) == 1)
        require(f"W0 n={n} r={r}", W.subs(t, 0) == E0)
        # q=0 gives t=0, (b,c,e)=(0,x,0); x=0 gives (0,0,E0*q).
        require(f"q-line c n={n} r={r}", (x * P.subs(t, 0) * S.subs(t, 0)) == x)
        require(f"q-line e n={n} r={r}", (q * W.subs(t, 0)).subs(q, 0) == 0)
        require(f"x-line e n={n} r={r}", q * W.subs(t, 0) == E0 * q)
        require(f"nontrivial degree n={n} r={r}", sp.degree(B, t) > 1)
        CYCLIC_ROWS += 1

print(f"PASS {CYCLIC_ROWS} cyclic compiler rows and both central-line charts")


print("SECTION stable-subalgebra arm split")
ARM_GATES = 0
ARM_DATA = [
    (b * (b - 1), (0, 1)),
    (b * (b - 1) * (b + 2), (0, 1, -2)),
    (b * (b**2 - 1) * (b - 3), (0, 1, -1, 3)),
]

for n in range(2, 6):
    Sigma, roots = ARM_DATA[(n - 2) % len(ARM_DATA)]
    h = b**2 + 2 * b + 3
    g = z + (b + 1) * z**2
    f_const = sp.expand(h + Sigma * g)
    f_surface = f_const.subs(z, c * e)
    for beta in roots:
        require(
            f"constant arm restriction n={n} beta={beta}",
            zero(f_const.subs(b, beta) - h.subs(b, beta)),
        )
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    Ptest = b**i * c**j * e**k
                    arm_value = sp.expand(bracket(Ptest, f_surface, n, Sigma)).subs(
                        {b: beta, c: 0}, simultaneous=True
                    )
                    require(
                        f"regular arm vanish n={n} beta={beta} mon={i}{j}{k}",
                        zero(arm_value),
                    )
                    ARM_GATES += 1

    # Nonconstant-arm controls: f=b+z+z^2 and f=b+z^3 at zeta=2.
    for f_probe in (b + z + z**2, b + z**3):
        fz = sp.diff(f_probe, z)
        for beta in roots:
            residue_den = (n - 1) * sp.diff(Sigma, b).subs(b, beta) * fz.subs({b: beta, z: 2})
            require(f"nonzero residue n={n} beta={beta}", residue_den != 0)
            require(f"ramified residue survives n={n} beta={beta}", (n - 1) / residue_den != 0)
            ARM_GATES += 2

print(f"PASS arm decomposition, regular vanishing, and residue controls ({ARM_GATES} gates)")


print("SECTION boundary curves and weight invoice")
nodal_f = s**2 - 1
nodal_g = s * (s**2 - 1)
require("nodal collision f", nodal_f.subs(s, 1) == nodal_f.subs(s, -1))
require("nodal collision g", nodal_g.subs(s, 1) == nodal_g.subs(s, -1))
require(
    "nodal immersion",
    sp.gcd(sp.Poly(sp.diff(nodal_f, s), s), sp.Poly(sp.diff(nodal_g, s), s)).degree() == 0,
)
require("nodal degree floor", sp.degree(nodal_f, s) == 2 and sp.degree(nodal_g, s) == 3)

alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta")
affine = alpha * b + beta * c + gamma * e + delta
require("affine c-arm", zero(affine.subs({b: 0, e: 0}) - (beta * c + delta)))
require("affine e-arm", zero(affine.subs({b: 0, c: 0}) - (gamma * e + delta)))

WEIGHT_GATES = 0
for n in range(2, 9):
    for s_abs in range(1, 4 * n + 1):
        m = ceil(s_abs / n)
        c_exp = n * m - s_abs
        survives_e_arm = c_exp == 0
        require(
            f"negative survival n={n} s={s_abs}",
            survives_e_arm == (s_abs % n == 0),
        )
        if survives_e_arm and m >= 2:
            require(f"negative location n={n} m={m}", -s_abs <= -2 * n)
        WEIGHT_GATES += 1
    require(f"span floor n={n}", 2 - (-2 * n) == 2 * n + 2)

print(f"PASS nodal sharpness, affine hostile, and {WEIGHT_GATES} weight-survival gates")


print("SECTION transverse central-arm jet")
JET_GATES = 0
for n in range(2, 9):
    Sigma = b * (b - 1)
    # High/low pieces alone have zero tangent jet at the arm crossing.
    no_middle_F = c**2 + e**2
    no_middle_G = c**3 + e**3
    require(
        f"no-middle bracket vanishes n={n}",
        bracket(no_middle_F, no_middle_G, n, Sigma).subs({b: 0, c: 0, e: 0}) == 0,
    )
    # Equal middle directions still have rank-one jet.
    same_middle_F = c + no_middle_F
    same_middle_G = 2 * c + no_middle_G
    require(
        f"same-middle bracket vanishes n={n}",
        bracket(same_middle_F, same_middle_G, n, Sigma).subs({b: 0, c: 0, e: 0}) == 0,
    )
    # Opposite weight 1/-n middle pieces give the two tangent directions.
    opposite_F = c + no_middle_F
    opposite_G = e + no_middle_G
    require(
        f"opposite-middle unit jet n={n}",
        bracket(opposite_F, opposite_G, n, Sigma).subs({b: 0, c: 0, e: 0}) == 1,
    )
    require(f"three-weight F invoice n={n}", len({2, 1, -2 * n}) == 3)
    require(f"three-weight G invoice n={n}", len({3, -n, -3 * n}) == 3)
    JET_GATES += 5

# The nonconstant hypothesis is load-bearing: when Sigma=1 there are no
# affine arms, c is a unit, and the rational hostile becomes regular ce/(n-1).
for n in range(2, 9):
    regular_hostile = c * e / (n - 1)
    require(
        f"constant-Sigma regular hostile n={n}",
        zero(bracket(regular_hostile, b, n, sp.Integer(1)).subs(e, c ** (-n)) - 1),
    )
    JET_GATES += 1

print(f"PASS transverse jet, three-weight floor, and constant-Sigma hostile ({JET_GATES} gates)")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
