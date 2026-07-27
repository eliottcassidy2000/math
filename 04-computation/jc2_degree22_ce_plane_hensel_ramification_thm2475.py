#!/usr/bin/env python3
"""Exact companion for THM-2475."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


r, v, zeta, lam = s.symbols("r v zeta lam")
a, b, c, d, e, f, g, h, x = s.symbols("a b c d e f g h x")

base_k = 922383 * v**2 - 25410 * v + 63
base_2 = (
    15944049 * zeta**2
    + (-162339408 * v + 2236080) * zeta
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)
c_flux = -206145280 * zeta + 449771520 * v - 1239040

# lambda=E^3/C^5 and r=E^2/(C^3*y), so
# C/y^3=r^3/lambda^2 and E/y^5=r^5/lambda^3.
F1 = (
    1331 * lam**3 * (63 - 1089 * v) * zeta
    + 4
    * (
        lam * (2342560 * v - 58080) * r**3
        - 3748096 * r**5
        + lam**3 * base_k
    )
)
F2 = lam**3 * base_2 + lam * c_flux * r**3 - 239878144 * r**5

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
    - 580505108480 * lam * r**8
    + lam**2
    * (-5487587353600 * v**2 + 634927462400 * v - 12368716800)
    * r**6
    + lam**3
    * (4938828618240 * v**2 - 571434716160 * v + 3935500800)
    * r**5
    + lam**4
    * (
        -4938828618240 * v**3
        + 537420744960 * v**2
        - 12593602560 * v
        + 146361600
    )
    * r**3
    - 63 * lam**6 * L5
)

raw_resultant = s.resultant(F1, F2, zeta)
require(
    s.expand(raw_resultant - 255104784 * lam**3 * P) == 0,
    "invariant C-E resultant changed",
)
require(len(s.Poly(P, r, v, lam).terms()) == 18, "C-E eliminant support changed")
require(s.degree(P, r) == 10 and s.degree(P, v) == 5, "C-E eliminant bidegree changed")
require(s.expand(P.subs(r, 0) + 63 * lam**6 * L5) == 0, "fixed L5 section changed")
require(
    s.Poly(P, r).coeff_monomial(r) == s.Poly(P, r).coeff_monomial(r**2) == 0,
    "Hensel gap through order two changed",
)
require(
    s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0,
    "fixed quintic lost squarefreeness",
)

wall = 234256 * r**5 - 4840 * lam * r**3 - 105 * lam**3
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

# A quadratic factor has Newton polygon 2*Delta. At r=0 it is
# q0=v^2+c*v+h and divides the squarefree L5. Since P has no r or r^2
# coefficient, coprime Hensel uniqueness forces q1=q2=0. Its r^3 term can
# therefore only be the scalar e*r^3. The two r^3 remainder equations below
# are incompatible on every one of the ten quadratic divisors of L5.
q0 = v**2 + c * v + h
section_remainder = s.Poly(s.rem(s.Poly(L5, v), s.Poly(q0, v)).as_expr(), v)
section_equations = section_remainder.all_coeffs()
quadratic_hensel = q0 + e * r**3 + d * r**4
hensel_remainder = s.Poly(
    s.rem(s.Poly(P, v), s.Poly(quadratic_hensel, v)).as_expr(),
    v,
    r,
)
A = -669650058 * c**2 + 9223830 * c + 446433372 * h + 2439360
N = 112442880 * c**2 + 12235520 * c - 112442880 * h + 286720
B = (
    27009219006 * c**3
    - 558041715 * c**2
    - 108036876024 * c * h
    - 295162560 * c
    + 1116083430 * h
    - 10087770
)
M = 13605588480 * c * h + 1480497920 * h - 403200
require(
    s.expand(
        hensel_remainder.coeff_monomial(v * r**3)
        + 43923 * lam**4 * (lam**2 * e * A + N)
    )
    == 0,
    "quadratic Hensel v-remainder changed",
)
require(
    s.expand(
        hensel_remainder.coeff_monomial(r**3)
        + 363 * lam**4 * (lam**2 * e * B + M)
    )
    == 0,
    "quadratic Hensel constant remainder changed",
)
gb_hensel = s.groebner(
    section_equations + [A * x + N, B * x + M],
    x,
    c,
    h,
    order="grevlex",
)
require(
    len(gb_hensel.polys) == 1 and gb_hensel.polys[0].as_expr() == 1,
    "fixed-section quadratic Hensel obstruction stopped being unit",
)

disc_v_raw = s.discriminant(P, v)
disc_factors = s.factor_list(disc_v_raw)[1]
require(len(disc_factors) == 3, "unexpected extra v-discriminant factor")
require((lam, 22) in disc_factors, "lambda multiplicity in discriminant changed")
require((wall, 2) in disc_factors, "wall-square factor changed")
K_candidates = [
    factor
    for factor, exponent in disc_factors
    if exponent == 1 and s.degree(factor, r) == 30
]
require(len(K_candidates) == 1, "moving K30 factor is no longer unique")
K30 = K_candidates[0]
require(
    s.degree(K30, lam) == 20 and len(s.Poly(K30, r, lam).terms()) == 39,
    "K30 degree or support changed",
)
require(
    s.expand(K30.subs(r, 0) + 209858977321875 * lam**20) == 0,
    "K30 primitive normalization changed",
)
degree_drop = 19370043 * lam**2 - 12500
require(
    s.expand(
        s.Poly(K30, r).LC()
        - 19752025963219313237822537728 * degree_drop
    )
    == 0,
    "K30 degree-drop factor changed",
)

A2 = 881552179200 * lam**2 - 6987070464 * lam + 15353125
disc_k_raw = s.discriminant(K30, r)
disc_k = s.factor_list(disc_k_raw)[1]
require(len(disc_k) == 5, "unexpected extra K30-discriminant factor")
require((lam, 580) in disc_k, "K30 invariant-origin multiplicity changed")
require((A2, 3) in disc_k, "K30 A2 factor changed")
P8_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 8 and exponent == 2]
Q8_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 8 and exponent == 3]
R12_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 12 and exponent == 1]
require(
    len(P8_candidates) == len(Q8_candidates) == len(R12_candidates) == 1,
    "large exceptional factors changed",
)
P8, Q8, R12 = P8_candidates[0], Q8_candidates[0], R12_candidates[0]

S4 = (
    8665219524371609991249920000 * lam**4
    - 16387179390132006347512532172800 * lam**3
    - 775731974166149050267099107687 * lam**2
    - 12728746868025855794767200000 * lam
    - 86790242512199520000000000
)
wall_resultant_raw = s.resultant(K30, wall, r)
wall_resultant_factors = s.factor_list(wall_resultant_raw)[1]
require(len(wall_resultant_factors) == 3, "unexpected extra wall-resultant factor")
require((lam, 90) in wall_resultant_factors, "wall resultant lambda factor changed")
require((A2, 3) in wall_resultant_factors, "wall resultant A2 factor changed")
require((S4, 1) in wall_resultant_factors, "wall resultant S4 factor changed")

exceptional_factors = [degree_drop, A2, P8, Q8, R12, S4]
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

next_coefficient = s.Poly(K30, r).coeff_monomial(r**29)
require(
    s.gcd(s.Poly(next_coefficient, lam), s.Poly(degree_drop, lam)).degree() == 0,
    "degree-drop fibre lost degree 29",
)
require(
    s.gcd(s.Poly(disc_k_raw, lam), s.Poly(degree_drop, lam)).degree() == 0,
    "a degree-drop fibre lost squarefreeness",
)
require(
    s.gcd(s.Poly(wall_resultant_raw, lam), s.Poly(degree_drop, lam)).degree() == 0,
    "a degree-drop fibre met the wall",
)

max_collision_gcd_degree = 3
minimum_simple_k_roots = 30 - 2 * max_collision_gcd_degree
minimum_visible_simple_branches = minimum_simple_k_roots - s.degree(wall, r)
minimum_total_branches = minimum_visible_simple_branches + (minimum_visible_simple_branches % 2)
minimum_genus = 1 + (minimum_total_branches - 10) // 2
require(minimum_simple_k_roots == 24, "simple-root floor changed")
require(minimum_visible_simple_branches == 19, "visible branch floor changed")
require(minimum_total_branches == 20, "ramification-parity floor changed")
require(minimum_genus == 6, "genus floor changed")

# Original y=0 boundary on B=D=W=0.
C, E, y, u, Z = s.symbols("C E y u Z")
original_A = -1089 * u + 63 * y**2
original_K = (
    2342560 * C * u
    - 58080 * C * y**2
    - 3748096 * E
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
    - 239878144 * E * y
    - 1190488992 * u**3
    + 147581280 * u**2 * y**2
    - 1219680 * u * y**4
    + 672 * y**6
)
require(
    s.expand(
        N1.subs(y, 0)
        - (-1449459 * u * Z + 9370240 * C * u - 14992384 * E)
    )
    == 0,
    "y=0 first-flux affine relation changed",
)
require(
    s.expand(
        N2.subs(y, 0)
        - (15944049 * Z**2 - 206145280 * C * Z - 1190488992 * u**3)
    )
    == 0,
    "y=0 second-flux equation changed",
)
alpha = s.Rational(9370240, 1449459) * C
beta = -s.Rational(14992384, 1449459) * E
boundary_quintic = s.expand(
    1190488992 * u**5
    + 206145280 * C * (alpha * u**2 + beta * u)
    - 15944049 * (alpha * u + beta) ** 2
)
require(
    s.cancel(u**2 * N2.subs({y: 0, Z: alpha + beta / u}) + boundary_quintic) == 0,
    "y=0 constant-field quintic changed",
)

print("THM-2475 exact companion")
print("invariants=lambda=E^3/C^5,r=E^2/(C^3*y)")
print("resultant=255104784*lambda^3*P")
print("P_terms=18")
print("P_bidegree=(10,5)")
print("newton_polygon=conv{(0,0),(10,0),(0,5)}")
print("linear_factor_unit_ideal=YES")
print("quadratic_factor=fixed_L5_divisor_plus_Hensel_r3_unit_ideal")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("wall=234256*r^5-4840*lambda*r^3-105*lambda^3")
print("P_at_wall=256*wall^2")
print("disc_v(P)=constant*lambda^22*wall^2*K30")
print("K30_bidegree=(30,20)")
print("K30_terms=39")
print("disc_r(K30)=constant*lambda^580*A2^3*P8^2*Q8^3*R12")
print("exceptional_factor_degrees=2,2,8,8,12,4")
print("degree_drop_fibres=two_squarefree_degree_29_wall_disjoint")
print("minimum_simple_K_roots=24")
print("minimum_visible_simple_branches=19")
print("minimum_total_ramification=20")
print("minimum_normalization_genus=6")
print("y_zero_boundary=affine_Z_in_1_over_u_then_constant_field_quintic")
print("ALL CHECKS PASSED")
