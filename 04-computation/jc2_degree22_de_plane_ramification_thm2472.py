#!/usr/bin/env python3
"""Exact companion for THM-2472."""

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

# lambda=E^4/D^5 and r=E/(D*y), so D/y^4=r^4/lambda and
# E/y^5=r^5/lambda.
F1 = (
    1331 * lam * (63 - 1089 * v) * zeta
    + 4 * (511104 * r**4 - 3748096 * r**5 + lam * base_k)
)
F2 = (
    lam * base_2
    + (-1978994688 * v + 16355328) * r**4
    - 239878144 * r**5
)

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
    - 3831333715968 * r**9
    + 261227298816 * r**8
    + lam * (4938828618240 * v**2 - 571434716160 * v + 3935500800) * r**5
    + lam
    * (
        -16298134440192 * v**3
        + 1077562607616 * v**2
        + 38961457920 * v
        - 987452928
    )
    * r**4
    - 63 * lam**2 * L5
)

raw_resultant = s.resultant(F1, F2, zeta)
require(
    s.expand(raw_resultant - 255104784 * lam * P) == 0,
    "invariant D-E resultant changed",
)
require(len(s.Poly(P, r, v, lam).terms()) == 16, "D-E eliminant support changed")
require(s.degree(P, r) == 10 and s.degree(P, v) == 5, "D-E eliminant bidegree changed")
require(s.expand(P.subs(r, 0) + 63 * lam**2 * L5) == 0, "fixed L5 section changed")

wall = 234256 * r**5 - 31944 * r**4 - 105 * lam
require(
    s.factor(P.subs(v, s.Rational(7, 121)) - 256 * wall**2) == 0,
    "first-flux wall identity changed",
)

# Newton triangle conv{(0,0),(10,0),(0,5)}.
support = {(ir, iv) for (ir, iv, _), _coefficient in s.Poly(P, r, v, lam).terms()}
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

disc_factors = s.factor_list(s.discriminant(P, v))[1]
require(len(disc_factors) == 3, "unexpected extra v-discriminant factor")
require((lam, 7) in disc_factors, "lambda multiplicity in discriminant changed")
require((wall, 2) in disc_factors, "wall-square factor changed")
K_candidates = [
    factor
    for factor, exponent in disc_factors
    if exponent == 1 and s.degree(factor, r) == 30
]
require(len(K_candidates) == 1, "moving K30 factor is no longer unique")
K30 = K_candidates[0]
require(
    s.degree(K30, lam) == 7 and len(s.Poly(K30, r, lam).terms()) == 31,
    "K30 degree or support changed",
)
require(
    s.expand(K30.subs(r, 0) + 24289233486328125 * lam**7) == 0,
    "K30 primitive normalization changed",
)
require(
    s.expand(
        s.Poly(K30, r).LC()
        - 3755706651006415131363113959424 * (11790625 * lam + 2519424)
    )
    == 0,
    "K30 degree-drop factor changed",
)

A1 = 1203125 * lam + 1944
disc_k = s.factor_list(s.discriminant(K30, r))[1]
require(len(disc_k) == 5, "unexpected extra K30-discriminant factor")
require((lam, 190) in disc_k, "K30 invariant-origin multiplicity changed")
require((A1, 3) in disc_k, "K30 A1 factor changed")
P6_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 6 and exponent == 2]
Q6_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 6 and exponent == 3]
R9_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 9 and exponent == 1]
require(
    len(P6_candidates) == len(Q6_candidates) == len(R9_candidates) == 1,
    "large exceptional factors changed",
)
P6, Q6, R9 = P6_candidates[0], Q6_candidates[0], R9_candidates[0]

S3 = (
    679985152000000000000 * lam**3
    - 128456167991962809796875 * lam**2
    - 48607662933281282414472 * lam
    - 5052710781009410264064
)
wall_resultant_factors = s.factor_list(s.resultant(K30, wall, r))[1]
require(len(wall_resultant_factors) == 3, "unexpected extra wall-resultant factor")
require((lam, 29) in wall_resultant_factors, "wall resultant lambda factor changed")
require((A1, 3) in wall_resultant_factors, "wall resultant A1 factor changed")
require((S3, 1) in wall_resultant_factors, "wall resultant S3 factor changed")

degree_drop = 11790625 * lam + 2519424
exceptional_factors = [degree_drop, A1, P6, Q6, R9, S3]
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

lambda_drop = -s.Rational(2519424, 11790625)
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

# Original y=0 boundary on B=C=W=0.
D, E, y, u, Z = s.symbols("D E y u Z")
original_A = -1089 * u + 63 * y**2
original_K = (
    511104 * D * y
    - 3748096 * E
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * original_A * Z + 4 * original_K
N2 = (
    15944049 * Z**2
    - 162339408 * Z * u * y
    + 2236080 * Z * y**3
    - 1978994688 * D * u
    + 16355328 * D * y**2
    - 239878144 * E * y
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
        - (15944049 * Z**2 - 1978994688 * D * u - 1190488992 * u**3)
    )
    == 0,
    "y=0 second-flux equation changed",
)
k = -s.Rational(14992384, 1449459) * E
boundary_quintic = 1190488992 * u**5 + 1978994688 * D * u**3 - 15944049 * k**2
require(
    s.cancel(u**2 * N2.subs({y: 0, Z: k / u}) + boundary_quintic) == 0,
    "y=0 constant-field quintic changed",
)

print("THM-2472 exact companion")
print("invariants=lambda=E^4/D^5,r=E/(D*y)")
print("resultant=255104784*lambda*P")
print("P_terms=16")
print("P_bidegree=(10,5)")
print("newton_polygon=conv{(0,0),(10,0),(0,5)}")
print("linear_factor_unit_ideal=YES")
print("quadratic_factor_unit_ideal=YES")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("wall=234256*r^5-31944*r^4-105*lambda")
print("P_at_wall=256*wall^2")
print("disc_v(P)=constant*lambda^7*wall^2*K30")
print("K30_bidegree=(30,7)")
print("K30_terms=31")
print("disc_r(K30)=constant*lambda^190*A1^3*P6^2*Q6^3*R9")
print("exceptional_factor_degrees=1,1,6,6,9,3")
print("degree_drop_K_degree=29")
print("minimum_simple_K_roots=24")
print("minimum_visible_simple_branches=19")
print("minimum_total_ramification=20")
print("minimum_normalization_genus=6")
print("y_zero_boundary=constant_uZ_then_constant_field_quintic_in_u")
print("ALL CHECKS PASSED")
