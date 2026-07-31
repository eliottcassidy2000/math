#!/usr/bin/env python3
"""Exact scratch controls for the THM-2993 reciprocal cross wall.

This is intentionally a scratch referee, not canonical evidence.
"""

import sympy as sp


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


T = sp.symbols("T")
a, b, c = sp.symbols("a b c")
H = (
    a**6 * c**2
    - 2 * a**4 * b * c**2
    - 2 * a**3 * b**3 * c
    - 26 * a**3 * c**3
    + 29 * a**2 * b**2 * c**2
    - 2 * a * b**4 * c
    - 18 * a * b * c**3
    + b**6
    - 26 * b**3 * c**2
    + 189 * c**4
)


def quartic_invariants(p, q, r):
    aa = 2 * p
    bb = p**2 - 4 * r
    cc = -q**2
    ff = T**4 + p * T**2 + q * T + r
    DD = sp.factor(sp.discriminant(ff, T))
    KK = sp.factor(bb**3 - aa**3 * cc)
    HH = sp.factor(H.subs({a: aa, b: bb, c: cc}))
    return DD, KK, HH


# Independent H wall and four-real sign controls.
require(quartic_invariants(-7, 7, 0) == (2401, -16807, 0), "wall hostile")
require(quartic_invariants(-4, 2, 1) == (592, -320, -154368), "negative control")
require(
    quartic_invariants(-10, 1, 0) == (3973, 992000, 980121756189),
    "positive control",
)


# Ordered marked wall: xyz = +/- Vandermonde is a smooth plane cubic.
x, y, z = sp.symbols("x y z")
F = sp.expand(x * y * z - (x - y) * (x - z) * (y - z))
partials = [sp.diff(F, v) for v in (x, y, z)]
for chart, variables in (({x: 1}, (y, z)), ({y: 1}, (x, z)), ({z: 1}, (x, y))):
    polys = [sp.expand(g.subs(chart)) for g in (F, *partials)]
    require(sp.groebner(polys, *variables, order="lex").polys == [sp.Poly(1, *variables)], "smooth chart")

# The even-permutation subgroup C3 acts freely, and elementary symmetric
# functions are its generic quotient on the chosen sign sheet.  A projective
# fixed point has shape [1:lambda:lambda^2], lambda^3=1; the gcd check proves
# that none lies on C_+.
cyclic = {x: y, y: z, z: x}
require(sp.expand(F.subs(cyclic, simultaneous=True) - F) == 0, "C3 invariance")
lam = sp.symbols("lam")
fixed_value = sp.expand(F.subs({x: 1, y: lam, z: lam**2}))
require(sp.gcd(fixed_value, lam**3 - 1) == 1, "free C3 action")
for symmetric in (x + y + z, x * y + y * z + z * x, x * y * z):
    require(
        sp.expand(symmetric.subs(cyclic, simultaneous=True) - symmetric) == 0,
        "symmetric quotient coordinate",
    )


# Projection from [1:0:0] gives the exact genus-one quartic.
tau, eta = sp.symbols("tau eta")
projected = sp.factor(F.subs({x: 1, z: tau * y}) / y)
qa = tau**2 - tau
qb = -tau**2 + tau + 1
qc = tau - 1
Q4 = sp.factor(qb**2 - 4 * qa * qc)
require(Q4 == tau**4 - 6 * tau**3 + 7 * tau**2 - 2 * tau + 1, "quartic projection")
I = 12 * 1 * 1 - 3 * (-6) * (-2) + 7**2
J = 72 * 1 * 7 * 1 + 9 * (-6) * 7 * (-2) - 27 * 1 * (-2) ** 2 - 27 * (-6) ** 2 * 1 - 2 * 7**3
require((I, J, sp.discriminant(Q4, tau)) == (25, -506, -7168), "binary quartic invariants")
require(sp.Rational(256 * I**3, sp.discriminant(Q4, tau)) == -sp.Rational(15625, 28), "ordered-wall j")


# Symmetric C3 quotient: Disc(g)=C0^2 and exact X_0(14) model.
t, s = sp.symbols("t s")
wall = t**2 - 4 * t**3 - 4 * s + 18 * t * s - 28 * s**2
v = 28 * s + 2 - 9 * t
require(
    sp.expand(v**2 - (2 - 7 * t) * (16 * t**2 - 11 * t + 2) + 28 * wall) == 0,
    "double-cover quotient",
)
xx = 32 - 112 * t
yy = -112 * v
require(
    sp.expand(yy**2 - xx * (xx**2 + 13 * xx + 128) + 351232 * wall) == 0,
    "E model",
)

X, Y = sp.symbols("X Y")
xx_sub = 4 * X - 4
yy_sub = 8 * Y + 4 * X + 4
E0 = sp.expand(yy_sub**2 - xx_sub * (xx_sub**2 + 13 * xx_sub + 128))
Emin = Y**2 + X * Y + Y - X**3 - 4 * X + 6
require(sp.factor(E0 - 64 * Emin) == 0, "minimal-model isomorphism")

# Exact two-real collision control after depression.
sqrt3 = sp.sqrt(3)
p2 = -2 - sqrt3
q2 = -2 * (1 + sqrt3)
r2 = (3 + 4 * sqrt3) / 4
D2, K2, H2 = quartic_invariants(p2, q2, r2)
require(sp.simplify(D2 + 4096 * (7 + 4 * sqrt3)) == 0, "two-real Delta")
require(sp.simplify(K2 + 512 * (12 + 7 * sqrt3)) == 0, "two-real K")
require(H2 == 0, "two-real H")

print("ordered_plane_cubic_smooth=1")
print("ordered_to_symmetric_quotient_degree=3")
print("ordered_quartic_model=w^2-(t^4-6t^3+7t^2-2t+1)=0")
print("ordered_j=-15625/28")
print("symmetric_quotient=y^2-x(x^2+13x+128)=0")
print("minimal_model=Y^2+XY+Y-X^3-4X+6=0")
print("minimal_change=x:4X-4,y:8Y+4X+4")
print("wall_hostile_D_K_H=2401,-16807,0")
print("four_real_negative_H=-154368")
print("four_real_positive_H=980121756189")
print("two_real_wall=1")
