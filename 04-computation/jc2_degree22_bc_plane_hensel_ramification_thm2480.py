#!/usr/bin/env python3
"""Exact companion for THM-2480."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


r, v, zeta, lam = s.symbols("r v zeta lam")
a, b, c, d, e, f, g, h = s.symbols("a b c d e f g h")
A, F, E = s.symbols("A F E")

base_k = 922383 * v**2 - 25410 * v + 63
base_2 = (
    15944049 * zeta**2
    + (-162339408 * v + 2236080) * zeta
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)
b_flux = 65591680 * zeta + 1443016960 * v**2 - 71554560 * v + 98560
c_flux = -206145280 * zeta + 449771520 * v - 1239040

# lambda=C^2/B^3 and r=C/(B*y), so
# B/y^2=r^2/lambda and C/y^3=r^3/lambda.
F1 = (
    1331 * (616 * r**2 + lam * (63 - 1089 * v)) * zeta
    + 4
    * (
        (-745360 * v + 6160) * r**2
        + (2342560 * v - 58080) * r**3
        + lam * base_k
    )
)
F2 = r**2 * b_flux + r**3 * c_flux + lam * base_2

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
P = (
    (55873616691200 * v - 1385296281600) * r**8
    + (-24889156526080 * v + 558316380160) * r**7
    + (
        -49388286182400 * lam * v**2
        + 5714347161600 * lam * v
        - 111318451200 * lam
        + 34222590223360 * v**2
        + 3959638538240 * v
        - 44411530240
    )
    * r**6
    + lam
    * (59714927838720 * v**2 - 3653718942720 * v + 64184440320)
    * r**5
    + lam
    * (
        -149234938081152 * v**3
        - 9500102156160 * v**2
        + 695599766400 * v
        - 6017413248
    )
    * r**4
    + lam**2
    * (
        -44449457564160 * v**3
        + 4836786704640 * v**2
        - 113342423040 * v
        + 1317254400
    )
    * r**3
    + lam**2
    * (
        206782580709936 * v**4
        + 6246495741024 * v**3
        - 1509756494400 * v**2
        + 34466937120 * v
        - 193496688
    )
    * r**2
    - 567 * lam**3 * L5
)

raw_resultant = s.resultant(F1, F2, zeta)
require(s.expand(raw_resultant - 28344976 * P) == 0, "invariant B-C resultant changed")
require(len(s.Poly(P, r, v, lam).terms()) == 32, "B-C eliminant support changed")
require(s.degree(P, r) == 8 and s.degree(P, v) == 5, "B-C eliminant bidegree changed")
require(s.expand(P.subs(r, 0) + 567 * lam**3 * L5) == 0, "fixed L5 section changed")
require(s.Poly(P, r).coeff_monomial(r) == 0, "Hensel order-one gap changed")
require(
    s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0,
    "fixed quintic lost squarefreeness",
)

wall = 745360 * r**5 - 71148 * r**4 + 43560 * lam * r**3 + 5082 * lam * r**2 + 945 * lam**2
first_flux_wall_v = (616 * r**2 + 63 * lam) / (1089 * lam)
require(
    s.cancel(P.subs(v, first_flux_wall_v) - 256 * wall**2 / (9 * lam)) == 0,
    "first-flux wall identity changed",
)

# Newton quadrilateral conv{(0,0),(8,0),(8,1),(0,5)}.
support = {(ir, iv) for (ir, iv, _), _coefficient in s.Poly(P, r, v, lam).terms()}
vertices = {(0, 0), (8, 0), (8, 1), (0, 5)}
require(vertices <= support, "Newton quadrilateral lost a vertex")
require(
    all(
        ir >= 0 and iv >= 0 and ir <= 8 and ir + 2 * iv <= 10
        for ir, iv in support
    ),
    "Newton support escaped the claimed quadrilateral",
)
# Edge lengths (8,1,4,5) imply summands (2c,b,c,b+c), b in {0,1}.
# The possible small positive v-degrees are exhausted by these four pairs.
small_summands = [
    (b_edge, c_edge)
    for b_edge in range(2)
    for c_edge in range(5)
    if 1 <= b_edge + c_edge <= 2
]
require(
    small_summands == [(0, 1), (0, 2), (1, 0), (1, 1)],
    "small Minkowski-summand inventory changed",
)

line = v + a * r**2 + b * r + c
line_remainder = s.Poly(s.rem(s.Poly(P, v), s.Poly(line, v)).as_expr(), r)
gb_line = s.groebner(line_remainder.all_coeffs(), a, b, c, lam, order="grevlex")
require(
    len(gb_line.polys) == 1 and gb_line.polys[0].as_expr() == 1,
    "linear-factor coefficient ideal stopped being unit",
)

# A quadratic factor has q0=v^2+c*v+h dividing squarefree L5. The missing
# order one forces q1=0. Write A=lambda*a, F=lambda*f, E=lambda*e for the
# most general order-two and order-three terms allowed by the larger
# triangular summand. The quadrilateral summand is a specialization. After
# this scaling, the six equations through r^3 are independent of lambda.
P_unit = s.expand(P.subs(lam, 1))
P_low_unit = sum(s.Poly(P_unit, r).coeff_monomial(r**i) * r**i for i in range(4))
quadratic_low_unit = v**2 + c * v + h + (A * v + F) * r**2 + E * r**3
low_remainder_unit = s.Poly(
    s.rem(s.Poly(P_low_unit, v), s.Poly(quadratic_low_unit, v)).as_expr(),
    v,
    r,
)
low_equations = [
    coefficient
    for (iv, ir), coefficient in low_remainder_unit.terms()
    if ir in (0, 2, 3)
]
require(len(low_equations) == 6, "low-order quadratic Hensel system changed size")

# Independently check that the lambda rescaling really recovers the original
# low-order equations, rather than merely a lambda=1 specialization.
P_low = sum(s.Poly(P, r).coeff_monomial(r**i) * r**i for i in range(4))
quadratic_low = (
    v**2 + c * v + h + (A * v + F) * r**2 / lam + E * r**3 / lam
)
low_remainder = s.Poly(
    s.cancel(s.rem(s.Poly(P_low, v), s.Poly(quadratic_low, v)).as_expr()),
    v,
    r,
)
for monomial in ((1, 0), (0, 0), (1, 2), (0, 2), (1, 3), (0, 3)):
    scale = lam**3 if monomial[1] == 0 else lam**2
    require(
        s.cancel(
            low_remainder.coeff_monomial(v ** monomial[0] * r ** monomial[1])
            - scale
            * low_remainder_unit.coeff_monomial(v ** monomial[0] * r ** monomial[1])
        )
        == 0,
        "lambda-scaled low-order Hensel equation changed",
    )

gb_quadratic_low = s.groebner(
    low_equations,
    A,
    F,
    E,
    c,
    h,
    order="grevlex",
)
require(
    len(gb_quadratic_low.polys) == 1
    and gb_quadratic_low.polys[0].as_expr() == 1,
    "low-order quadratic Hensel ideal stopped being unit",
)

disc_v_raw = s.discriminant(P, v)
disc_factors = s.factor_list(disc_v_raw)[1]
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
    s.degree(K30, lam) == 13 and len(s.Poly(K30, r, lam).terms()) == 87,
    "K30 degree or support changed",
)
require(
    s.expand(K30.subs(r, 0) - 38724883599623115234375 * lam**13) == 0,
    "K30 primitive normalization changed",
)
degree_drop = (
    49121386875 * lam**3
    + 35843727150 * lam**2
    + 34155999339 * lam
    + 1686616064
)
require(
    s.expand(
        s.Poly(K30, r).LC()
        - 927499340872432064135168000000 * degree_drop
    )
    == 0,
    "K30 degree-drop factor changed",
)

L1 = 1782 * lam + 245
L2 = 22275 * lam + 2744
A5 = (
    134711704729495154443359375 * lam**5
    + 1800011771181150399288281250 * lam**4
    + 6418659057813101346623578125 * lam**3
    + 1620235334410707585796706250 * lam**2
    + 103287017473928512574894775 * lam
    - 384735691349720065722248
)
disc_k_raw = s.discriminant(K30, r)
disc_k = s.factor_list(disc_k_raw)[1]
require(len(disc_k) == 7, "unexpected extra K30-discriminant factor")
require((lam, 325) in disc_k, "K30 invariant-origin multiplicity changed")
require((L1, 2) in disc_k, "K30 L1 factor changed")
require((L2, 4) in disc_k, "K30 L2 factor changed")
require((A5, 3) in disc_k, "K30 A5 factor changed")
R18_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 18 and exponent == 1]
P20_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 20 and exponent == 2]
Q20_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 20 and exponent == 3]
require(
    len(R18_candidates) == len(P20_candidates) == len(Q20_candidates) == 1,
    "large exceptional factors changed",
)
R18, P20, Q20 = R18_candidates[0], P20_candidates[0], Q20_candidates[0]

S8 = (
    12337829835953864392149281489742232500000000000000 * lam**8
    - 44369823212208566246105819275522866775312500000000 * lam**7
    + 16885864307881167661079424840506375193448798828125 * lam**6
    - 21179519183140752458837472696065162224372050000000 * lam**5
    - 9294048284435692495345979933901524446649735250000 * lam**4
    + 131589285074186063043261783346562971089053066250 * lam**3
    + 205953213712224149574057621432756223317422854200 * lam**2
    + 942027416127503933058469688031656737095805128 * lam
    - 1311283001340114972135581249229837964617938539
)
wall_resultant_raw = s.resultant(K30, wall, r)
wall_resultant_factors = s.factor_list(wall_resultant_raw)[1]
require(len(wall_resultant_factors) == 3, "unexpected extra wall-resultant factor")
require((lam, 52) in wall_resultant_factors, "wall resultant lambda factor changed")
require((A5, 3) in wall_resultant_factors, "wall resultant A5 factor changed")
require((S8, 1) in wall_resultant_factors, "wall resultant S8 factor changed")

exceptional_factors = [degree_drop, L1, L2, A5, R18, P20, Q20, S8]
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

minimum_simple_at_L1 = 30 - 2 * 2
minimum_simple_at_L2 = 30 - 2 * 4
minimum_simple_at_A5 = 30 - 2 * 3
minimum_visible_at_A5 = minimum_simple_at_A5 - s.degree(wall, r)
minimum_visible_at_S8 = 30 - s.degree(wall, r)
minimum_visible_branches = min(
    minimum_simple_at_L1,
    minimum_simple_at_L2,
    minimum_visible_at_A5,
    minimum_visible_at_S8,
)
minimum_total_branches = minimum_visible_branches + (minimum_visible_branches % 2)
minimum_genus = 1 + (minimum_total_branches - 10) // 2
require(minimum_simple_at_L1 == 26, "L1 simple-root floor changed")
require(minimum_simple_at_L2 == 22, "L2 simple-root floor changed")
require(minimum_visible_at_A5 == 19, "A5 visible branch floor changed")
require(minimum_visible_at_S8 == 25, "S8 visible branch floor changed")
require(minimum_visible_branches == 19, "uniform visible branch floor changed")
require(minimum_total_branches == 20, "ramification-parity floor changed")
require(minimum_genus == 6, "genus floor changed")

# Original y=0 boundary on D=E=W=0.
Bpar, Cpar, y, u, Z = s.symbols("Bpar Cpar y u Z")
original_A = 616 * Bpar - 1089 * u + 63 * y**2
original_K = (
    -745360 * Bpar * u * y
    + 6160 * Bpar * y**3
    + 2342560 * Cpar * u
    - 58080 * Cpar * y**2
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * original_A * Z + 4 * original_K
N2 = (
    15944049 * Z**2
    + 65591680 * Bpar * Z * y
    - 206145280 * Cpar * Z
    - 162339408 * Z * u * y
    + 2236080 * Z * y**3
    + 1443016960 * Bpar * u**2
    - 71554560 * Bpar * u * y**2
    + 98560 * Bpar * y**4
    + 449771520 * Cpar * u * y
    - 1239040 * Cpar * y**3
    - 1190488992 * u**3
    + 147581280 * u**2 * y**2
    - 1219680 * u * y**4
    + 672 * y**6
)
denominator = 1331 * (616 * Bpar - 1089 * u)
numerator = -9370240 * Cpar * u
require(
    s.expand(N1.subs(y, 0) - (denominator * Z - numerator)) == 0,
    "y=0 first-flux rational reconstruction changed",
)
require(
    s.expand(
        N2.subs(y, 0)
        - (
            15944049 * Z**2
            - 206145280 * Cpar * Z
            + 1443016960 * Bpar * u**2
            - 1190488992 * u**3
        )
    )
    == 0,
    "y=0 second-flux equation changed",
)
cleared_boundary = s.cancel(
    denominator**2 * N2.subs({y: 0, Z: numerator / denominator})
)
require(cleared_boundary.as_numer_denom()[1] == 1, "boundary denominator did not clear")
boundary_quintic = s.Poly(s.expand(cleared_boundary), u)
require(boundary_quintic.degree() == 5, "y=0 constant-field quintic changed degree")
require(boundary_quintic.LC() != 0, "y=0 constant-field quintic lost its leading term")

print("THM-2480 exact companion")
print("invariants=lambda=C^2/B^3,r=C/(B*y)")
print("resultant=28344976*P")
print("P_terms=32")
print("P_bidegree=(8,5)")
print("newton_polygon=conv{(0,0),(8,0),(8,1),(0,5)}")
print("small_minkowski_summands=two_line_shapes_plus_two_quadratic_shapes")
print("linear_factor_unit_ideal=YES")
print("quadratic_factor=fixed_L5_scaled_Hensel_order_three_unit_ideal")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("wall=745360*r^5-71148*r^4+43560*lambda*r^3+5082*lambda*r^2+945*lambda^2")
print("P_at_wall=256*wall^2/(9*lambda)")
print("disc_v(P)=constant*lambda^7*wall^2*K30")
print("K30_bidegree=(30,13)")
print("K30_terms=87")
print("disc_r(K30)=constant*lambda^325*L1^2*L2^4*A5^3*R18*P20^2*Q20^3")
print("exceptional_factor_degrees=3,1,1,5,18,20,20,8")
print("degree_drop_fibres=three_squarefree_degree_29_wall_disjoint")
print("minimum_visible_simple_branches=19")
print("minimum_total_ramification=20")
print("minimum_normalization_genus=6")
print("y_zero_boundary=rational_Z_then_constant_field_quintic")
print("ALL CHECKS PASSED")
