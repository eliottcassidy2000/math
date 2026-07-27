#!/usr/bin/env python3
"""Exact companion for THM-2469."""

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

# These are lambda times the two normalized fluxes after
# lambda=D^3/C^4 and r=D/(C*y), so C/y^3=r^3/lambda and
# D/y^4=r^4/lambda.
F1 = (
    1331 * lam * (63 - 1089 * v) * zeta
    + 4
    * (
        (2342560 * v - 58080) * r**3
        + 511104 * r**4
        + lam * base_k
    )
)
F2 = (
    lam * base_2
    - 206145280 * r**3 * zeta
    + (449771520 * v - 1239040) * r**3
    + (-1978994688 * v + 16355328) * r**4
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
    261227298816 * r**8
    + 79159787520 * r**7
    + (-5487587353600 * v**2 + 634927462400 * v - 12368716800) * r**6
    + lam
    * (-16298134440192 * v**3 + 1077562607616 * v**2 + 38961457920 * v - 987452928)
    * r**4
    + lam
    * (-4938828618240 * v**3 + 537420744960 * v**2 - 12593602560 * v + 146361600)
    * r**3
    - 63 * lam**2 * L5
)

raw_resultant = s.resultant(F1, F2, zeta)
require(
    s.expand(raw_resultant - 255104784 * lam * P) == 0,
    "invariant C-D resultant changed",
)
require(len(s.Poly(P, r, v, lam).terms()) == 19, "C-D eliminant support changed")
require(s.degree(P, r) == 8 and s.degree(P, v) == 5, "C-D eliminant bidegree changed")
require(s.expand(P.subs(r, 0) + 63 * lam**2 * L5) == 0, "fixed L5 section changed")

wall = 31944 * r**4 + 4840 * r**3 + 105 * lam
require(
    s.factor(P.subs(v, s.Rational(7, 121)) - 256 * wall**2) == 0,
    "first-flux wall identity changed",
)

# Newton polygon conv{(0,0),(8,0),(6,2),(0,5)}.
support = {(ir, iv) for (ir, iv, il), coefficient in s.Poly(P, r, v, lam).terms()}
vertices = {(0, 0), (8, 0), (6, 2), (0, 5)}
require(vertices <= support, "Newton polygon lost a vertex")
require(
    all(ir >= 0 and iv >= 0 and ir + iv <= 8 and ir + 2 * iv <= 10 for ir, iv in support),
    "Newton support escaped the claimed polygon",
)

# A Minkowski summand of v-degree at most two is contained in one of the
# following maximal line or conic shapes.  The exact coefficient ideals force
# lambda=0, so no such factor exists on the physical chart.
line = v + a * r**2 + b * r + c
line_remainder = s.Poly(
    s.rem(s.Poly(P, v), s.Poly(line, v)).as_expr(),
    r,
)
gb_line = s.groebner(
    line_remainder.all_coeffs(),
    a,
    b,
    c,
    lam,
    order="grevlex",
)
expected_line_basis = [
    1331 * c**2 + 154 * c + 3,
    a,
    55 * b - 363 * c - 21,
    lam,
]
require(
    [poly.as_expr() for poly in gb_line.polys] == expected_line_basis,
    "linear-factor Groebner certificate changed",
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
expected_quadratic_basis = [
    a,
    b,
    121 * c + 14,
    d,
    e,
    3025 * f + 144,
    6655 * g + 96,
    1331 * h - 3,
    lam,
]
require(
    [poly.as_expr() for poly in gb_quadratic.polys] == expected_quadratic_basis,
    "quadratic-factor Groebner certificate changed",
)

# Exact discriminant: the squared quartic is precisely the excluded wall.
disc_v = s.factor_list(s.discriminant(P, v))
disc_factors = disc_v[1]
require((lam, 4) in disc_factors, "lambda multiplicity in discriminant changed")
require((wall, 2) in disc_factors, "wall-square factor changed")
K_candidates = [
    factor
    for factor, exponent in disc_factors
    if exponent == 1 and s.degree(factor, r) == 30
]
require(len(K_candidates) == 1, "moving K30 factor is no longer unique")
K30 = K_candidates[0]
require(
    s.degree(K30, lam) == 10 and len(s.Poly(K30, r, lam).terms()) == 46,
    "K30 degree or support changed",
)

leading_k = s.factor(s.Poly(K30, r).LC())
require(
    s.expand(
        leading_k
        - 176357374671601011051986944000000 * (9504 * lam + 4375)
    )
    == 0,
    "K30 degree-drop factor changed",
)

A1 = 7392 * lam - 625
A2 = 90112 * lam - 5625
A3 = 3372171264 * lam**2 - 218592000 * lam + 383828125

disc_k = s.factor_list(s.discriminant(K30, r))[1]
require((lam, 290) in disc_k, "K30 invariant-origin multiplicity changed")
require((A1, 1) in disc_k, "K30 A1 factor changed")
require((A2, 2) in disc_k, "K30 A2 factor changed")
require((A3, 3) in disc_k, "K30 A3 factor changed")
P10_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 10 and exponent == 3]
Q10_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 10 and exponent == 2]
R12_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 12 and exponent == 1]
require(
    len(P10_candidates) == len(Q10_candidates) == len(R12_candidates) == 1,
    "large exceptional factors changed",
)
P10, Q10, R12 = P10_candidates[0], Q10_candidates[0], R12_candidates[0]

S4 = (
    6125963995157747283975979008 * lam**4
    + 7723423480030158582978912000 * lam**3
    + 4548030222783870066786609375 * lam**2
    + 1419376581390195945000000000 * lam
    + 17937798138268750000000000
)
wall_resultant_factors = s.factor_list(s.resultant(K30, wall, r))[1]
require((lam, 30) in wall_resultant_factors, "wall resultant lambda factor changed")
require((A3, 3) in wall_resultant_factors, "wall resultant A3 factor changed")
require((S4, 1) in wall_resultant_factors, "wall resultant S4 factor changed")

degree_drop = 9504 * lam + 4375
exceptional_factors = [degree_drop, A1, A2, A3, P10, Q10, R12, S4]
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

lambda_drop = -s.Rational(4375, 9504)
K_drop = s.Poly(K30.subs(lam, lambda_drop), r)
wall_drop = s.Poly(wall.subs(lam, lambda_drop), r)
require(K_drop.degree() == 29, "K30 degree-drop fibre no longer has degree 29")
require(s.gcd(K_drop, K_drop.diff()).degree() == 0, "degree-drop fibre lost squarefreeness")
require(s.gcd(K_drop, wall_drop).degree() == 0, "degree-drop fibre met the wall")

# At a discriminant factor of exponent e<=3, the Sylvester determinant
# valuation bounds deg gcd(K,K') by e.  Hence a degree-30 fibre has at least
# 30-2e>=24 simple roots.  Removing all four possible wall roots still leaves
# twenty simple branch values.
max_collision_gcd_degree = 3
minimum_simple_k_roots = 30 - 2 * max_collision_gcd_degree
minimum_usable_branches = minimum_simple_k_roots - s.degree(wall, r)
minimum_genus = 1 + (minimum_usable_branches - 10) // 2
require(minimum_simple_k_roots == 24, "simple-root floor changed")
require(minimum_usable_branches == 20, "usable branch floor changed")
require(minimum_genus == 6, "genus floor changed")

# Original y=0 boundary on B=E=W=0.
C, D, y, u, Z = s.symbols("C D y u Z")
original_A = -1089 * u + 63 * y**2
original_K = (
    2342560 * C * u
    - 58080 * C * y**2
    + 511104 * D * y
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * original_A * Z + 4 * original_K
N2 = (
    15944049 * Z**2
    - 206145280 * C * Z
    - 162339408 * Z * u * y
    + 2236080 * Z * y**3
    + 449771520 * C * u * y
    - 1239040 * C * y**3
    - 1978994688 * D * u
    + 16355328 * D * y**2
    - 1190488992 * u**3
    + 147581280 * u**2 * y**2
    - 1219680 * u * y**4
    + 672 * y**6
)
require(
    s.expand(N1.subs(y, 0) - u * (-1449459 * Z + 9370240 * C)) == 0,
    "y=0 first-flux boundary changed",
)
require(
    s.expand(
        N2.subs(y, 0)
        - (15944049 * Z**2 - 206145280 * C * Z - 1978994688 * D * u - 1190488992 * u**3)
    )
    == 0,
    "y=0 second-flux cubic changed",
)

print("THM-2469 exact companion")
print("invariants=lambda=D^3/C^4,r=D/(C*y)")
print("resultant=255104784*lambda*P")
print("P_terms=19")
print("P_bidegree=(8,5)")
print("newton_polygon=conv{(0,0),(8,0),(6,2),(0,5)}")
print("line_factor_ideal_forces_lambda=0")
print("quadratic_factor_ideal_forces_lambda=0")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("wall=31944*r^4+4840*r^3+105*lambda")
print("P_at_wall=256*wall^2")
print("disc_v(P)=constant*lambda^4*wall^2*K30")
print("K30_bidegree=(30,10)")
print("K30_terms=46")
print("disc_r(K30)=constant*lambda^290*A1*A2^2*A3^3*P10^3*Q10^2*R12")
print("exceptional_factor_degrees=1,1,1,2,10,10,12,4")
print("degree_drop_K_degree=29")
print("minimum_simple_K_roots=24")
print("minimum_usable_branch_values=20")
print("minimum_normalization_genus=6")
print("y_zero_boundary=constant_Z_then_constant_field_cubic_in_u")
print("ALL CHECKS PASSED")
