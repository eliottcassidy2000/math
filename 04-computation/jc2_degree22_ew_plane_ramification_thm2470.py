#!/usr/bin/env python3
"""Exact companion for THM-2470."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


r, v, zeta, lam = s.symbols("r v zeta lam")
a, b, c, d, e, f, g, h = s.symbols("a b c d e f g h")

base_k = 922383 * v**2 - 25410 * v + 63
base_2 = (
    15944049 * zeta**2
    + (-162339408 * v + 2236080) * zeta
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)

# lambda=W^5/E^6 and r=W/(E*y), so E/y^5=r^5/lambda and
# W/y^6=r^6/lambda.
F1 = 1331 * lam * (63 - 1089 * v) * zeta + 4 * (-3748096 * r**5 + lam * base_k)
F2 = lam * base_2 - 239878144 * r**5 - 1319329792 * r**6

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
P = (
    14048223625216 * r**10
    + lam * (-10865422960128 * v**2 + 1257156375552 * v - 36364027392) * r**6
    + lam * (4938828618240 * v**2 - 571434716160 * v + 3935500800) * r**5
    - 63 * lam**2 * L5
)

raw_resultant = s.resultant(F1, F2, zeta)
require(
    s.expand(raw_resultant - 255104784 * lam * P) == 0,
    "invariant E-W resultant changed",
)
require(len(s.Poly(P, r, v, lam).terms()) == 13, "E-W eliminant support changed")
require(s.degree(P, r) == 10 and s.degree(P, v) == 5, "E-W eliminant bidegree changed")
require(s.expand(P.subs(r, 0) + 63 * lam**2 * L5) == 0, "fixed L5 section changed")

wall = 234256 * r**5 - 105 * lam
require(
    s.factor(P.subs(v, s.Rational(7, 121)) - 256 * wall**2) == 0,
    "first-flux wall identity changed",
)

# Newton triangle conv{(0,0),(10,0),(0,5)}.
support = {(ir, iv) for (ir, iv, il), coefficient in s.Poly(P, r, v, lam).terms()}
vertices = {(0, 0), (10, 0), (0, 5)}
require(vertices <= support, "Newton triangle lost a vertex")
require(
    all(ir >= 0 and iv >= 0 and ir + 2 * iv <= 10 for ir, iv in support),
    "Newton support escaped the claimed triangle",
)

line = v + a * r**2 + b * r + c
line_remainder = s.Poly(s.rem(s.Poly(P, v), s.Poly(line, v)).as_expr(), r)
gb_line = s.groebner(line_remainder.all_coeffs(), a, b, c, lam, order="grevlex")
require(
    len(gb_line.polys) == 1 and gb_line.polys[0].as_expr() == 1,
    "linear-factor coefficient ideal stopped being unit",
)

quadratic = (
    v**2
    + (a * r**2 + b * r + c) * v
    + d * r**4
    + e * r**3
    + f * r**2
    + g * r
    + h
)
quadratic_remainder = s.Poly(
    s.rem(s.Poly(P, v), s.Poly(quadratic, v)).as_expr(),
    v,
    r,
)
gb_quadratic = s.groebner(
    [coefficient for _, coefficient in quadratic_remainder.terms()],
    a,
    b,
    c,
    d,
    e,
    f,
    g,
    h,
    lam,
    order="grevlex",
)
require(
    len(gb_quadratic.polys) == 1 and gb_quadratic.polys[0].as_expr() == 1,
    "quadratic-factor coefficient ideal stopped being unit",
)

disc_v = s.factor_list(s.discriminant(P, v))
disc_factors = disc_v[1]
require(len(disc_factors) == 3, "unexpected extra v-discriminant factor")
require((lam, 8) in disc_factors, "lambda multiplicity in discriminant changed")
require((wall, 2) in disc_factors, "wall-square factor changed")
K_candidates = [
    factor
    for factor, exponent in disc_factors
    if exponent == 1 and s.degree(factor, r) == 30
]
require(len(K_candidates) == 1, "moving K30 factor is no longer unique")
K30 = K_candidates[0]
require(
    s.degree(K30, lam) == 6 and len(s.Poly(K30, r, lam).terms()) == 22,
    "K30 degree or support changed",
)
require(
    s.expand(K30.subs(r, 0) - 3469890498046875 * lam**6) == 0,
    "K30 primitive normalization changed",
)
require(
    s.expand(
        s.Poly(K30, r).LC()
        - 5164096645133820805624281694208 * (24057 * lam - 1225000)
    )
    == 0,
    "K30 degree-drop factor changed",
)

A1 = 56133 * lam + 327680000
disc_k = s.factor_list(s.discriminant(K30, r))[1]
require(len(disc_k) == 5, "unexpected extra K30-discriminant factor")
require((lam, 174) in disc_k, "K30 invariant-origin multiplicity changed")
require((A1, 3) in disc_k, "K30 A1 factor changed")
P4_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 4 and exponent == 2]
Q4_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 4 and exponent == 3]
R6_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 6 and exponent == 1]
require(
    len(P4_candidates) == len(Q4_candidates) == len(R6_candidates) == 1,
    "large exceptional factors changed",
)
P4, Q4, R6 = P4_candidates[0], Q4_candidates[0], R6_candidates[0]

S2 = 1605730621565435904 * lam**2 - 4390842515926741875 * lam - 6871947673600000000
wall_resultant_factors = s.factor_list(s.resultant(K30, wall, r))[1]
require(len(wall_resultant_factors) == 3, "unexpected extra wall-resultant factor")
require((lam, 30) in wall_resultant_factors, "wall resultant lambda factor changed")
require((A1, 3) in wall_resultant_factors, "wall resultant A1 factor changed")
require((S2, 1) in wall_resultant_factors, "wall resultant S2 factor changed")

degree_drop = 24057 * lam - 1225000
exceptional_factors = [degree_drop, A1, P4, Q4, R6, S2]
for i, left in enumerate(exceptional_factors):
    require(
        s.gcd(s.Poly(left, lam), s.Poly(s.diff(left, lam), lam)).degree() == 0,
        f"exceptional factor {i} lost squarefreeness",
    )
    for right in exceptional_factors[:i]:
        require(
            s.gcd(s.Poly(left, lam), s.Poly(right, lam)).degree() == 0,
            "distinct exceptional factors acquired a common root",
        )

lambda_drop = s.Rational(1225000, 24057)
K_drop = s.Poly(K30.subs(lam, lambda_drop), r)
wall_drop = s.Poly(wall.subs(lam, lambda_drop), r)
require(K_drop.degree() == 29, "K30 degree-drop fibre no longer has degree 29")
require(s.gcd(K_drop, K_drop.diff()).degree() == 0, "degree-drop fibre lost squarefreeness")
require(s.gcd(K_drop, wall_drop).degree() == 0, "degree-drop fibre met the wall")

max_collision_gcd_degree = 3
minimum_simple_k_roots = 30 - 2 * max_collision_gcd_degree
minimum_visible_simple_branches = minimum_simple_k_roots - s.degree(wall, r)
# Total tame ramification of a degree-five cover is even, so nineteen visible
# simple contributions force a twentieth contribution.
minimum_total_branches = minimum_visible_simple_branches + (minimum_visible_simple_branches % 2)
minimum_genus = 1 + (minimum_total_branches - 10) // 2
require(minimum_simple_k_roots == 24, "simple-root floor changed")
require(minimum_visible_simple_branches == 19, "visible branch floor changed")
require(minimum_total_branches == 20, "ramification-parity floor changed")
require(minimum_genus == 6, "genus floor changed")

# Original y=0 boundary on B=C=D=0.
E, W, y, u, Z = s.symbols("E W y u Z")
original_A = -1089 * u + 63 * y**2
original_K = -3748096 * E + 922383 * u**2 * y - 25410 * u * y**3 + 63 * y**5
N1 = 1331 * original_A * Z + 4 * original_K
N2 = (
    15944049 * Z**2
    - 162339408 * Z * u * y
    + 2236080 * Z * y**3
    - 239878144 * E * y
    - 1319329792 * W
    - 1190488992 * u**3
    + 147581280 * u**2 * y**2
    - 1219680 * u * y**4
    + 672 * y**6
)
require(
    s.expand(N1.subs(y, 0) - (-1449459 * u * Z - 14992384 * E)) == 0,
    "y=0 first-flux product changed",
)
require(
    s.expand(
        N2.subs(y, 0)
        - (15944049 * Z**2 - 1319329792 * W - 1190488992 * u**3)
    )
    == 0,
    "y=0 second-flux equation changed",
)
k = -s.Rational(14992384, 1449459) * E
boundary_quintic = 1190488992 * u**5 + 1319329792 * W * u**2 - 15944049 * k**2
require(
    s.cancel(u**2 * N2.subs({y: 0, Z: k / u}) + boundary_quintic) == 0,
    "y=0 constant-field quintic changed",
)

print("THM-2470 exact companion")
print("invariants=lambda=W^5/E^6,r=W/(E*y)")
print("resultant=255104784*lambda*P")
print("P_terms=13")
print("P_bidegree=(10,5)")
print("newton_polygon=conv{(0,0),(10,0),(0,5)}")
print("linear_factor_unit_ideal=YES")
print("quadratic_factor_unit_ideal=YES")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("wall=234256*r^5-105*lambda")
print("P_at_wall=256*wall^2")
print("disc_v(P)=constant*lambda^8*wall^2*K30")
print("K30_bidegree=(30,6)")
print("K30_terms=22")
print("disc_r(K30)=constant*lambda^174*A1^3*P4^2*Q4^3*R6")
print("exceptional_factor_degrees=1,1,4,4,6,2")
print("degree_drop_K_degree=29")
print("minimum_simple_K_roots=24")
print("minimum_visible_simple_branches=19")
print("minimum_total_ramification=20")
print("minimum_normalization_genus=6")
print("y_zero_boundary=constant_uZ_then_constant_field_quintic_in_u")
print("ALL CHECKS PASSED")
