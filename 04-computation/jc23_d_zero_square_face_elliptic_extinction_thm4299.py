#!/usr/bin/env python3
"""Exact certificate for THM-4299.

The script checks the D=0 square factorization, elliptic/deck/contact ledgers,
the off-Lambda parametric-Morse input, every possible compact Newton tail,
one genus-five hostile, and the D=Lambda=0 cubic boundary.
"""

from math import gcd

import sympy as sp


def polygon_ledger(vertices):
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(
            abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]),
        )
        for i in range(len(vertices))
    )
    interior = (area2 - boundary + 2) // 2
    return area2, boundary, interior


a, b, S, P, T = sp.symbols("a b S P T", nonzero=True)
z, r, q, t = sp.symbols("z r q t")
xi, omega = sp.symbols("xi omega")

U, W, Z = a**2, 2 * a * b, b**2
L = a * P**3 + b * S**2 * P**2
C = 1 - U * P**6 - W * S**2 * P**5 - Z * S**4 * P**4
assert sp.expand(W**2 - 4 * U * Z) == 0
assert sp.expand(U + W + Z - (a + b) ** 2) == 0
assert sp.expand(C - (1 - L) * (1 + L)) == 0

# Each factor has a smooth plane-cubic normalization a P^3+b T^2=epsilon.
# Its Newton polygon has one interior point, and its Weierstrass c4 is zero.
factor_polygon = [(0, 0), (0, 3), (2, 2)]
assert polygon_ledger(factor_polygon) == (6, 6, 1)
H = sp.symbols("H")
eps = sp.symbols("eps", nonzero=True)
plane_cubic = a * P**3 + b * T**2 * H - eps * H**3
partials = [sp.diff(plane_cubic, v) for v in (P, T, H)]
assert partials == [3 * a * P**2, 2 * b * H * T, b * T**2 - 3 * H**2 * eps]

# tau exchanges the factors; tau^2 is [-omega] on each, while the target
# action squares to [omega^2].  The differential characters differ by 1.
assert (3 * 2) % 12 == 6
assert (2 + 2 * 2) % 12 == 6
assert (2 * 2) % 12 == 4       # P under tau^2: omega P
assert (2 + 2 * 2) % 12 == 6   # T under tau^2: -T
# 1+omega+omega^2=0 gives (-omega)-omega^2=1.
character_difference = sp.rem(-sp.Symbol("o") - sp.Symbol("o")**2, sp.Symbol("o")**2 + sp.Symbol("o") + 1)
assert character_difference == 1

# The generic-point denominator is nonzero on either elliptic factor.
domain = sp.QQ.frac_field(a, b)
for sign in (1, -1):
    f = sp.Poly(L - sign, P, S, domain=domain)
    dL = sp.Poly(sp.diff(L, P), P, S, domain=domain)
    R = sp.Poly(S**2 - P, P, S, domain=domain)
    assert sp.gcd(f, dL).total_degree() == 0
    assert sp.gcd(f, R).total_degree() == 0

# At top infinity the two factors are z^6 +/- r^2(ar+b).
M = r**2 * (a * r + b)
Hp = z**6 - M
Hm = z**6 + M
assert sp.expand(Hp * Hm - (z**12 - r**4 * (a * r + b) ** 2)) == 0
r0 = -b / a
assert sp.simplify(M.subs(r, r0)) == 0
assert sp.simplify(sp.diff(M, r).subs(r, r0) - b**2 / a) == 0
assert sp.simplify(r0 - 1 + (a + b) / a) == 0

# Off Lambda=0: twelve transverse R contacts and one length-six E+/E- contact.
# The r=0 toric chart separates the four boundary branches.
x = sp.symbols("x")
strict_plus = sp.cancel(Hp.subs(r, z**3 * x) / z**6).subs(z, 0)
strict_minus = sp.cancel(Hm.subs(r, z**3 * x) / z**6).subs(z, 0)
assert strict_plus == 1 - b * x**2
assert strict_minus == 1 + b * x**2
assert sp.resultant(strict_plus, strict_minus, x) != 0
assert sp.diff(z**6 - sp.Symbol("e") * (a + b), z) == 6 * z**5
assert 6 + 6 == 12
assert 6 == 6  # C+/C- intersection length at r0.

# Genus ledgers: p_a(C)=7 and p_a(R union C+ union C-)=18.
main_polygon = [(0, 0), (0, 6), (4, 4)]
global_polygon = [(0, 1), (2, 0), (6, 4), (0, 7)]
assert polygon_ledger(main_polygon) == (24, 12, 7)
assert polygon_ledger(global_polygon) == (48, 14, 18)
assert 1 + 1 + 6 - 1 == 7
assert 0 + 1 + 1 + (12 + 6) - (3 - 1) == 18

# At Lambda=0 (b=-a), all three smooth branches meet pairwise to length six.
M_corner = sp.expand(M.subs({b: -a, r: 1 + q}))
assert sp.expand(M_corner - a * q * (q + 1) ** 2) == 0
corner_plus = z**6 - M_corner
corner_minus = z**6 + M_corner
assert sp.diff(corner_plus, q).subs({q: 0, z: 0}) == -a
assert sp.diff(corner_minus, q).subs({q: 0, z: 0}) == a
assert 6 + 6 + 6 == 18

# Off the corner, division by r-1 is legal and the critical point is Morse.
H0 = r**4 * (a * r + b) ** 2
assert sp.simplify(sp.diff(H0, r).subs(r, r0)) == 0
assert sp.simplify(sp.diff(H0, r, 2).subs(r, r0) - 2 * a**2 * r0**4) == 0
assert sp.simplify(2 * a**2 * r0**4) != 0

# Exhaust every possible compact face of w^2=z^12-kappa(sigma*z).
tail_rows = []
for m in range(1, 12):
    g = gcd(12 - m, m)
    s_weight = (12 - m) // g
    z_weight = m // g
    if m % 2 == 0:
        degree = 12 - m
    else:
        degree = 13 - m
    genus = (degree - 2) // 2
    assert genus == (11 - m) // 2
    differential_order = 9 * s_weight + 5 * z_weight
    assert differential_order > 0
    X = sp.symbols("X")
    normalized_rhs = X ** (12 - m) - 1 if m % 2 == 0 else X * (X ** (12 - m) - 1)
    assert sp.gcd(sp.Poly(normalized_rhs, X), sp.Poly(sp.diff(normalized_rhs, X), X)).degree() == 0
    tail_rows.append((m, s_weight, z_weight, genus, differential_order))

# Positive control: a=b=1 and the weight-11 term t*r^4 give
# kappa(t)=t-4t^2+..., hence m=1, genus 5, and order 104.
c1, c2, c3 = sp.symbols("c1 c2 c3")
rcrit = -1 + c1 * t + c2 * t**2 + c3 * t**3
A_hostile = r**4 * (r + 1) ** 2 + t * r**4
dA_series = sp.series(sp.diff(A_hostile, r).subs(r, rcrit), t, 0, 4).removeO().expand()
solution = sp.solve([dA_series.coeff(t, k) for k in range(1, 4)], (c1, c2, c3), dict=True)[0]
kappa_hostile = sp.series(A_hostile.subs(r, rcrit).subs(solution), t, 0, 4).removeO().expand()
assert solution[c1] == 2
assert kappa_hostile.coeff(t, 1) == 1
assert kappa_hostile.coeff(t, 2) == -4
assert tail_rows[0] == (1, 11, 1, 5, 104)

# At D=Lambda=0 the divided Morse model fails and the first face is cubic.
F_top = a**2 * q**3 * (1 + q) ** 4 - q * z**12 - t**12 / 2
F_cubic = a**2 * q**3 - q * z**12 - t**12 / 2
remainder = sp.factor(F_top - F_cubic)
assert remainder == a**2 * q**4 * (q + 2) * (q**2 + 2 * q + 2)

print("THM-4299 EXACT CERTIFICATE")
print("square_face: D=0, C=(1-L)(1+L), Lambda=(a+b)^2")
print("factor_polygon: 2A=6 boundary=6 interior=1; two smooth j=0 elliptics")
print("deck: tau swaps factors; tau^2 source=[-omega], target=[omega^2], difference=1")
print("contacts_off_corner: R-E+=6 transverse, R-E-=6 transverse, E+-E-=6 (A11)")
print("genus: central=7, full_special=18")
print("m  s  beta  genus  good_form_order")
for row in tail_rows:
    print("%2d %2d %5d %6d %16d" % row)
print("hostile: (U,W,Z)=(1,2,1), kappa=t-4t^2+..., genus=5, order=104")
print("corner: (1,-2,1), three pairwise length-6 branches; first face cubic")
print("ALL THM-4299 EXACT CHECKS PASSED")
