"""Exact polynomial certificate for THM-2566.

The script distinguishes a raw pullback from its saturation.  For the
sporadic Keller target coordinates (a,b,c), it certifies two maps to the
cuspidal cubic Y^2=X^3.  Each raw inverse image contains the Jelonek
hypersurface L=0 and one coefficient-zero plane with multiplicity two.

Run with both

    python3 04-computation/keller_two_cusp_atlas_thm2566.py
    python3 -O 04-computation/keller_two_cusp_atlas_thm2566.py

The checks do not use Python ``assert``, so optimized replay is substantive.
"""

import sympy as sp


a, b, c, u = sp.symbols("a b c u")
x, y, z = sp.symbols("x y z")

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2

# The escape-coordinate (x) cusp chart.
T_c = 4 - 3 * b * c
S_c = 27 * a * c**2 - 9 * b * c + 8

# The companion u=1+xy cusp chart from THM-1335.
T_a = b**2 - 12 * a
S_a = 54 * a**2 * c - 18 * a * b + b**3


def require(condition, message):
    """Raise in ordinary and optimized Python if an exact check fails."""
    if not condition:
        raise RuntimeError(message)


def zero(expr, message):
    require(sp.cancel(expr) == 0, message)


def saturation_generator(pullback, h):
    """Eliminate u from (pullback, 1-u*h), the standard h-saturation."""
    basis = sp.groebner([pullback, 1 - u * h], u, a, b, c, order="lex")
    eliminated = [
        sp.Poly(poly.as_expr(), a, b, c, domain=sp.QQ).monic().as_expr()
        for poly in basis.polys
        if not poly.as_expr().has(u)
    ]
    require(len(eliminated) == 1, f"unexpected saturation basis for {h}")
    return sp.factor(eliminated[0])


print("=" * 78)
print("[1] The two cusp identities and their exact coefficient planes")
print("=" * 78)

pullback_c = sp.factor(S_c**2 - T_c**3)
pullback_a = sp.factor(S_a**2 - T_a**3)
require(pullback_c == 27 * c**2 * L, "c-chart cusp identity failed")
require(pullback_a == 108 * a**2 * L, "a-chart cusp identity failed")
p_weierstrass = 3 * (12 * a - b**2)
q_weierstrass = 2 * S_a
disc_weierstrass = sp.factor(-4 * p_weierstrass**3 - 27 * q_weierstrass**2)
require(
    disc_weierstrass == -(2**4) * (3**6) * a**2 * L,
    "Weierstrass discriminant identity failed",
)

print("  S_c^2-T_c^3 =", pullback_c)
print("  S_a^2-T_a^3 =", pullback_a)
print("  -4p^3-27q^2 =", disc_weierstrass)
print("  raw divisors: 2V(c)+V(L) and 2V(a)+V(L)  [PASS]")
print("  raw sets:     V(c) union V(L), V(a) union V(L)  [PASS]")
print()

print("=" * 78)
print("[2] Groebner elimination certifies both saturated pullbacks")
print("=" * 78)

sat_c = saturation_generator(S_c**2 - T_c**3, c)
sat_a = saturation_generator(S_a**2 - T_a**3, a)
L_monic = sp.Poly(L, a, b, c, domain=sp.QQ).monic().as_expr()
zero(sat_c - L_monic, "c-saturation is not (L)")
zero(sat_a - L_monic, "a-saturation is not (L)")

print("  ((S_c^2-T_c^3):c^infinity) = (L)  [PASS]")
print("  ((S_a^2-T_a^3):a^infinity) = (L)  [PASS]")
print("  Hence each chart recovers the whole irreducible Jelonek closure after")
print("  saturation, and it is a literal pullback on D(c) or D(a), respectively.")
print()

print("=" * 78)
print("[3] Parasitic planes and hostile points")
print("=" * 78)

zero(T_c.subs(c, 0) - 4, "c-plane T image failed")
zero(S_c.subs(c, 0) - 8, "c-plane S image failed")
smooth_gradient = (2 * 8, -3 * 4**2)
require(smooth_gradient != (0, 0), "(8,4) should be a smooth cusp point")

hostile_c = {a: 1, b: 0, c: 0}
require(L.subs(hostile_c) == 16, "c-hostile should lie off L")
require((S_c.subs(hostile_c), T_c.subs(hostile_c)) == (8, 4), "c-hostile image")

zero(T_a.subs(a, 0) - b**2, "a-plane T image failed")
zero(S_a.subs(a, 0) - b**3, "a-plane S image failed")
hostile_a = {a: 0, b: 1, c: 2}
require(L.subs(hostile_a) == 1, "a-hostile should lie off L")
require((S_a.subs(hostile_a), T_a.subs(hostile_a)) == (1, 1), "a-hostile image")

double_hostile = {a: 0, b: 1, c: 0}
require(L.subs(double_hostile) == -1, "double hostile should lie off L")
require(S_c.subs(double_hostile) ** 2 == T_c.subs(double_hostile) ** 3,
        "double hostile missed c-cusp")
require(S_a.subs(double_hostile) ** 2 == T_a.subs(double_hostile) ** 3,
        "double hostile missed a-cusp")

print("  c=0 maps constantly to (S_c,T_c)=(8,4), a smooth cusp point.")
print("  hostile (1,0,0): L=16 but Phi_c=(8,4).  [PASS]")
print("  a=0 maps to (S_a,T_a)=(b^3,b^2), the normalized cusp trace.")
print("  hostile (0,1,2): L=1 but Phi_a=(1,1).  [PASS]")
print("  double hostile (0,1,0): L=-1 but BOTH raw cusp equations vanish.")
print("  Therefore intersecting raw pullbacks still leaves V(a,c); saturation")
print("  or localization is essential.  [PASS]")
print()

print("=" * 78)
print("[4] The coefficient planes are affine targets, not boundary charts")
print("=" * 78)

radical = 1 + x * y
F1 = sp.expand(radical**3 * z + y**2 * radical * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * radical**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)

# A polynomial section of the entire target plane c=0.
section_c = {x: 0, y: b, z: a - 4 * b**2}
for actual, expected in zip(
    (F1.subs(section_c), F2.subs(section_c), F3.subs(section_c)),
    (a, b, 0),
):
    zero(actual - expected, "c-plane section failed")

print("  F(0,b,a-4b^2)=(a,b,0) identically.  [PASS]")
print("  Thus c=0 is an ordinary affine target plane.  THM-2473 supplies the")
print("  generic 1+2 fibre there; away from L=16a-b^2 it has all 3 points.")
print()

print("=" * 78)
print("[5] Punctured two-chart cover and the exceptional origin")
print("=" * 78)

origin_line_restriction = sp.factor(L.subs({a: 0, c: 0}, simultaneous=True))
require(origin_line_restriction == -b**2, "a=c=0 restriction failed")
origin = {a: 0, b: 0, c: 0}
require(L.subs(origin) == 0, "origin must lie on L")
require((S_c.subs(origin), T_c.subs(origin)) == (8, 4), "origin c-image")
require((S_a.subs(origin), T_a.subs(origin)) == (0, 0), "origin a-image")

print("  L|_(a=c=0)=-b^2.  Therefore a point of V(L) with a=c=0 is the origin.")
print("  V(L) minus {0} is covered by the honest cusp charts D(a) and D(c).")
print("  Phi_c(0)=(8,4) (smooth); Phi_a(0)=(0,0) (the cusp).  [PASS]")
print("  Exact atlas:")
print("    V(L) = {0} union (D(c) intersect Phi_c^-1(C))")
print("                 union (D(a) intersect Phi_a^-1(C)).")
print()

print("=" * 78)
print("[6] The c-chart cusp-point fibre is exactly the omitted curve")
print("=" * 78)

# S_c = (27ac^2-4)+3T_c, so (S_c,T_c)=(0,0) is equivalent to
# bc=4/3 and ac^2=4/27.  In particular c is automatically invertible.
zero(S_c - (27 * a * c**2 - 4) - 3 * T_c, "cusp-point reduction failed")
t = sp.symbols("t", nonzero=True)
omitted = {a: sp.Rational(4, 27) / t**2, b: sp.Rational(4, 3) / t, c: t}
require(T_c.subs(omitted) == 0, "omitted curve T_c failed")
require(S_c.subs(omitted) == 0, "omitted curve S_c failed")
require(L.subs(omitted) == 0, "omitted curve L failed")

print("  (S_c,T_c)=(0,0) iff bc=4/3 and ac^2=4/27.")
print("  Hence Phi_c^-1(cusp)={(4/(27t^2),4/(3t),t): t!=0}.")
print("  By THM-2546 this is exactly the empty-fibre curve.  [PASS]")
print()

print("=" * 78)
print("[7] The a-chart cusp point retains its own parasitic axis")
print("=" * 78)

# Modulo T_a, 8S_a=b^3(3bc-4).  The set is the z-axis union E.
zero(
    8 * S_a - b**3 * (3 * b * c - 4)
    - 3 * (12 * a - b**2) * (12 * a * c + b**2 * c - 4 * b),
    "a-chart cusp-point reduction failed",
)
print("  modulo T_a=0: 8S_a=b^3(3bc-4).  [PASS]")
print("  Thus Phi_a^-1(cusp point) is the z-axis {a=b=0} union the omitted")
print("  curve E.  On D(a), only E remains; at the origin Phi_a hits the cusp.")
print()

print("ALL EXACT THM-2566 CUSP-ATLAS CHECKS PASSED")
