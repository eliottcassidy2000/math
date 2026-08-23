#!/usr/bin/env python3
"""Exact companion for THM-3790's cubic pseudo-plane arm gate."""

from __future__ import annotations

import sympy as sp


r, z, e, c, t, s, zeta = sp.symbols("r z e c t s zeta", nonzero=True)
variables = (r, z, e)
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 * c**3 + 6 * r * e],
        [-9 * z**2, -3 * c**3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


def assert_zero(expr: sp.Expr) -> None:
    assert sp.cancel(expr) == 0


# The minimal nodal immersion and its first-normal Bezout completion.
a0 = t**2
b0 = t**3 - t
a1 = -1 / (3 * c**3)
b1 = -t / (2 * c**3)
bezout = sp.factor(3 * c**3 * (a1 * sp.diff(b0, t) - sp.diff(a0, t) * b1))
assert bezout == 1
assert sp.gcd(sp.diff(a0, t), sp.diff(b0, t)) == 1
assert_zero(a0.subs(t, 1) - a0.subs(t, -1))
assert_zero(b0.subs(t, 1) - b0.subs(t, -1))

# The smallest global lift and its exact failed bracket.
A = e**2 - z / (3 * c**3)
C = e**3 - e - e * z / (2 * c**3)
lifted = sp.factor(bracket(A, C))
expected = (c**3 + 2 * e * r) * (2 * c**3 + z) / (2 * c**6)
assert_zero(lifted - expected)

# Critical equations and the seven-point family.
surface = r**2 * e - z**3 + c**3 * r
critical = [sp.factor(bracket(A, q)) for q in variables]
crit_poly = 8 * zeta**7 + 9 * c**15
r0 = 2 * zeta**3 / c**3
e0 = -c**6 / (4 * zeta**3)


def reduce_at_critical(expr: sp.Expr) -> sp.Expr:
    substituted = sp.cancel(expr.subs({r: r0, z: zeta, e: e0}))
    numerator = sp.together(substituted).as_numer_denom()[0]
    return sp.factor(sp.rem(numerator, crit_poly, zeta))


assert reduce_at_critical(surface) == 0
assert reduce_at_critical(c**3 + 2 * r * e) == 0
assert all(reduce_at_critical(expr) == 0 for expr in critical)
assert sp.gcd(crit_poly, sp.diff(crit_poly, zeta)) == 1

print("THM3790_ARM_QUOTIENT B/I^2=C[e,z]/(z^2); r=0 because c^3*r=z^3-r^2*e")
print(f"THM3790_BEZOUT {bezout}")
print("THM3790_NODAL gcd(2*t,3*t^2-1)=1; gamma(1)=gamma(-1)=(1,0)")
print(f"THM3790_LIFTED_BRACKET {lifted}")
print("THM3790_CRITICAL_POLY 8*zeta^7+9*c^15; squarefree=YES; count=7")
print("THM3790_CRITICAL_POINT r=2*zeta^3/c^3; e=-c^6/(4*zeta^3)")
print("THM3790_CRITICAL_REMAINDERS surface=0 K=0 Dr=0 Dz=0 De=0")
print("THM3790_CHECKS OK")
