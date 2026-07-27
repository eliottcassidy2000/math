#!/usr/bin/env python3
"""Exact reduction for the degree-22 B-D-W mixed stratum (THM-2617)."""

from __future__ import annotations

from itertools import product

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_unit_basis(groebner_basis: s.GroebnerBasis) -> bool:
    return (
        len(groebner_basis.polys) == 1
        and groebner_basis.polys[0].as_expr() == 1
    )


p, v, zeta, lam, mu = s.symbols("p v zeta lam mu")
a, b, c, d, e, h = s.symbols("a b c d e h")

# B/y^2=p, D/y^4=lam*p^2, W/y^6=mu*p^3, u/y^2=v, Z/y^3=zeta.
f1 = (
    2044416 * lam * p**2
    - 2981440 * p * v
    + 819896 * p * zeta
    + 24640 * p
    + 3689532 * v**2
    - 1449459 * v * zeta
    - 101640 * v
    + 83853 * zeta
    + 252
)
f2 = (
    -1319329792 * mu * p**3
    - 1978994688 * lam * p**2 * v
    + 16355328 * lam * p**2
    + 1443016960 * p * v**2
    - 71554560 * p * v
    + 65591680 * p * zeta
    + 98560 * p
    - 1190488992 * v**3
    + 147581280 * v**2
    - 162339408 * v * zeta
    - 1219680 * v
    + 15944049 * zeta**2
    + 2236080 * zeta
    + 672
)

raw_resultant = s.resultant(f1, f2, zeta)
R = s.Poly(raw_resultant, p, v, lam, mu).primitive()[1].as_expr()
require(s.cancel(raw_resultant / R) == 2**4 * 11**6, "resultant content changed")
require(len(s.Poly(R, p, v, lam, mu).terms()) == 34, "resultant support changed")
require(s.degree(R, p) == 5 and s.degree(R, v) == 5, "resultant bidegree changed")

support_pv = {
    (monomial[0], monomial[1])
    for monomial, _ in s.Poly(R, p, v, lam, mu).terms()
}
expected_triangle = {(i, j) for i in range(6) for j in range(6 - i)}
require(support_pv == expected_triangle, "Newton polygon stopped being the full 5-triangle")

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
p0 = s.Poly(R, p).coeff_monomial(1)
p1 = s.Poly(R, p).coeff_monomial(p)
require(s.expand(p0 + 567 * L5) == 0, "fixed L5 section changed")
require(
    p1
    == 11088
    * (
        18649222647 * v**4
        + 563356398 * v**3
        - 136161300 * v**2
        + 3108490 * v
        - 17451
    ),
    "first moving coefficient changed",
)
require(s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0, "L5 not squarefree")
require(s.gcd(s.Poly(p0, v), s.Poly(p1, v)).degree() == 0, "fixed section lost transversality")

top = s.Integer(0)
for (ip, iv, ilam, imu), coefficient in s.Poly(R, p, v, lam, mu).terms():
    if ip + iv == 5:
        top += coefficient * p**ip * v**iv * lam**ilam * mu**imu
top_cubic = 256 * mu * p**3 + 384 * lam * p**2 * v - 280 * p * v**2 + 231 * v**3
require(
    s.expand(top + 38974342 * (56 * p - 99 * v) ** 2 * top_cubic) == 0,
    "top homogeneous factorization changed",
)

# Every reducible total-degree-five polynomial has a line or quadratic
# factor.  Its top form must select one or two of the five top lines.

# Linear type L1: the repeated wall line.
wall_line_remainder = s.Poly(R.subs(p, s.Rational(99, 56) * v + b), v)
wall_line_gb = s.groebner(
    [coefficient for _, coefficient in wall_line_remainder.terms()],
    b,
    lam,
    mu,
    order="grevlex",
)
require(is_unit_basis(wall_line_gb), "wall-top line factor returned")

# Linear type L2: one root h of the moving cubic top form.
mu_h = (-384 * lam * h**2 + 280 * h - 231) / (256 * h**3)
moving_line = s.cancel(R.subs({p: h * v + b, mu: mu_h}))
moving_line_equations = []
for _, coefficient in s.Poly(moving_line.as_numer_denom()[0], v).terms():
    primitive = s.Poly(coefficient, b, lam, h).primitive()[1].as_expr()
    if primitive != 0:
        moving_line_equations.append(primitive)
moving_line_gb = s.groebner(moving_line_equations, b, lam, h, order="grevlex")
for forced_monomial in (b**3, b**2 * h, b * h**2, h**3):
    require(
        moving_line_gb.reduce(forced_monomial)[1] == 0,
        "moving-top line monomial obstruction changed",
    )
require(s.expand(top_cubic.subs({p: 0, v: 1})) == 231, "h=0 became a top root")

# General monic quadratic factor and its coefficient equations.
quadratic = p**2 + (a * v + b) * p + c * v**2 + d * v + e
quadratic_remainder = s.Poly(
    s.rem(s.Poly(R, p), s.Poly(quadratic, p)).as_expr(),
    p,
    v,
)
quadratic_equations = [coefficient for _, coefficient in quadratic_remainder.terms()]

# Quadratic type Q1: both copies of the wall line.
type1_substitution = {a: -s.Rational(99, 28), c: s.Rational(9801, 3136)}
type1_equations = [s.factor(eq.subs(type1_substitution)) for eq in quadratic_equations]
type1_gb = s.groebner(type1_equations, b, d, e, lam, mu, order="grevlex")
require(is_unit_basis(type1_gb), "wall-square quadratic factor returned")

# Quadratic type Q2: one wall line and one moving-cubic root h.
type2_substitution = {
    a: -(s.Rational(99, 56) + h),
    c: s.Rational(99, 56) * h,
    mu: mu_h,
}
type2_equations = []
for equation in quadratic_equations:
    numerator = s.cancel(equation.subs(type2_substitution)).as_numer_denom()[0]
    primitive = s.Poly(numerator, b, d, e, lam, h).primitive()[1].as_expr()
    if primitive != 0:
        type2_equations.append(primitive)
type2_gb = s.groebner(type2_equations, b, d, e, lam, h, order="grevlex")
require(type2_gb.reduce(h**3)[1] == 0, "wall-moving quadratic no longer forces h=0")

# Quadratic type Q3: the other two roots after selecting one cubic root h.
# This is the sole unresolved type.  Freeze its exact chart and the first
# generic-solving determinant; do not infer its emptiness.
D0 = -384 * lam * h**2 + 280 * h - 231
type3_substitution = {
    a: h * (280 * h - 231) / D0,
    c: -231 * h**2 / D0,
    mu: D0 / (256 * h**3),
}
type3_top = p**2 + type3_substitution[a] * p * v + type3_substitution[c] * v**2
type3_top_remainder = s.cancel(
    s.rem(
        s.Poly(top_cubic.subs(mu, type3_substitution[mu]), p),
        s.Poly(type3_top, p),
    ).as_expr()
)
require(type3_top_remainder == 0, "last quadratic top type stopped dividing the cubic")

type3_by_monomial: dict[tuple[int, int], s.Expr] = {}
for monomial, equation in quadratic_remainder.terms():
    numerator = s.cancel(equation.subs(type3_substitution)).as_numer_denom()[0]
    primitive = s.Poly(numerator, b, d, e, lam, h).primitive()[1].as_expr()
    if primitive != 0:
        type3_by_monomial[monomial] = primitive
eq13 = type3_by_monomial[(1, 3)]
eq04 = type3_by_monomial[(0, 4)]
top_next_matrix = s.Matrix(
    [
        [s.diff(eq13, b), s.diff(eq13, d)],
        [s.diff(eq04, b), s.diff(eq04, d)],
    ]
)
type3_determinant = s.factor(top_next_matrix.det())
expected_type3_determinant = (
    966306
    * h
    * (384 * lam * h**2 - 560 * h + 693)
    * (384 * lam * h**2 - 280 * h + 231) ** 5
    * (114048 * lam * h**2 - 25088 * h**2 - 44352 * h + 68607) ** 2
)
require(type3_determinant == expected_type3_determinant, "last-type determinant changed")

# The two genuine determinant walls are both empty after physical
# saturation.  On the first, the complete ideal contains h^2.
lambda_wall_1 = (560 * h - 693) / (384 * h**2)
wall_1_equations = []
for equation in type3_by_monomial.values():
    numerator = s.cancel(equation.subs(lam, lambda_wall_1)).as_numer_denom()[0]
    primitive = s.Poly(numerator, b, d, e, h).primitive()[1].as_expr()
    if primitive != 0:
        wall_1_equations.append(primitive)
wall_1_gb = s.groebner(wall_1_equations, b, d, e, h, order="grevlex")
require(wall_1_gb.reduce(h**2)[1] == 0, "first determinant wall stopped forcing h=0")

# On the second wall the two top-next linear equations are parallel.  Their
# consistency minors leave only four rational h-values, and every complete
# specialized ideal is the unit ideal.
lambda_wall_2 = (
    25088 * h**2 + 44352 * h - 68607
) / (114048 * h**2)
wall_2_by_monomial: dict[tuple[int, int], s.Expr] = {}
for monomial, equation in type3_by_monomial.items():
    numerator = s.cancel(equation.subs(lam, lambda_wall_2)).as_numer_denom()[0]
    primitive = s.Poly(numerator, b, d, e, h).primitive()[1].as_expr()
    if primitive != 0:
        wall_2_by_monomial[monomial] = primitive
wall_2_top_a = wall_2_by_monomial[(1, 3)]
wall_2_top_b = wall_2_by_monomial[(0, 4)]
wall_2_rows = []
for equation in (wall_2_top_a, wall_2_top_b):
    wall_2_rows.append(
        (
            s.diff(equation, b),
            s.diff(equation, d),
            equation.subs({b: 0, d: 0}),
        )
    )
wall_2_minor_1 = s.expand(
    wall_2_rows[0][0] * wall_2_rows[1][2]
    - wall_2_rows[1][0] * wall_2_rows[0][2]
)
wall_2_minor_2 = s.expand(
    wall_2_rows[0][1] * wall_2_rows[1][2]
    - wall_2_rows[1][1] * wall_2_rows[0][2]
)
wall_2_consistency_gcd = s.Poly(wall_2_minor_1, h).gcd(
    s.Poly(wall_2_minor_2, h)
)
expected_wall_2_gcd = s.Poly(
    h**6
    * (8 * h - 99) ** 3
    * (32 * h - 99)
    * (56 * h - 99) ** 2
    * (64 * h - 99) ** 4,
    h,
)
require(
    wall_2_consistency_gcd.monic() == expected_wall_2_gcd.monic(),
    "second determinant-wall consistency atlas changed",
)
wall_2_candidates = (
    s.Rational(99, 8),
    s.Rational(99, 32),
    s.Rational(99, 56),
    s.Rational(99, 64),
)
for h_value in wall_2_candidates:
    specialized_equations = []
    for equation in wall_2_by_monomial.values():
        specialized = equation.subs(h, h_value)
        if specialized != 0:
            specialized_equations.append(
                s.Poly(specialized, b, d, e).primitive()[1].as_expr()
            )
    specialized_gb = s.groebner(
        specialized_equations,
        b,
        d,
        e,
        order="grevlex",
    )
    require(
        is_unit_basis(specialized_gb),
        f"second determinant wall acquired a factor at h={h_value}",
    )

# Away from the determinant walls, the next two equations are linear in e.
# Their exact compatibility resultant has only the already-inverted physical
# factors, the two closed walls, and one residual bidegree-(7,15) polynomial.
top_a_b = s.diff(eq13, b)
top_a_d = s.diff(eq13, d)
top_a_0 = eq13.subs({b: 0, d: 0})
top_b_b = s.diff(eq04, b)
top_b_d = s.diff(eq04, d)
top_b_0 = eq04.subs({b: 0, d: 0})
det_bd = s.expand(top_a_b * top_b_d - top_b_b * top_a_d)
b_numerator = s.expand(top_a_d * top_b_0 - top_b_d * top_a_0)
d_numerator = s.expand(top_b_b * top_a_0 - top_a_b * top_b_0)


def after_top_solve(equation: s.Expr) -> s.Expr:
    numerator = s.cancel(
        equation.subs({b: b_numerator / det_bd, d: d_numerator / det_bd})
    ).as_numer_denom()[0]
    return s.Poly(numerator, e, lam, h).primitive()[1].as_expr()


next_12 = after_top_solve(type3_by_monomial[(1, 2)])
next_03 = after_top_solve(type3_by_monomial[(0, 3)])
compatibility = s.Poly(
    s.resultant(next_12, next_03, e),
    lam,
    h,
).primitive()[1].as_expr()

known_compatibility_factors = (
    (h, 3),
    (384 * lam * h**2 - 560 * h + 693, 2),
    (384 * lam * h**2 - 280 * h + 231, 3),
    (
        114048 * lam * h**2
        - 25088 * h**2
        - 44352 * h
        + 68607,
        3,
    ),
)
residual_compatibility = s.Poly(compatibility, lam, h)
for factor, expected_exponent in known_compatibility_factors:
    factor_poly = s.Poly(factor, lam, h)
    actual_exponent = 0
    while True:
        quotient, remainder = residual_compatibility.div(factor_poly)
        if not remainder.is_zero:
            break
        residual_compatibility = quotient
        actual_exponent += 1
    require(
        actual_exponent == expected_exponent,
        "generic compatibility multiplicity changed",
    )
residual_compatibility = residual_compatibility.primitive()[1]
require(
    s.degree(residual_compatibility.as_expr(), lam) == 7
    and s.degree(residual_compatibility.as_expr(), h) == 15
    and len(residual_compatibility.terms()) == 69,
    "residual compatibility divisor changed",
)

# The root-scaled chart lowers all nine total degrees from 8--14 to at most
# eight and is the preferred continuation coordinate for the residual curve.
root_u, root_b, root_d, root_e = s.symbols("root_u root_b root_d root_e")
unscaled_degrees = [
    s.Poly(equation, b, d, e, lam, h).total_degree()
    for equation in type3_by_monomial.values()
]
require(
    (min(unscaled_degrees), max(unscaled_degrees)) == (8, 14),
    "unscaled degree range changed",
)
root_scaled_degrees = []
for equation in type3_by_monomial.values():
    numerator = s.cancel(
        equation.subs(
            {
                lam: root_u / h**2,
                b: root_b * h,
                d: root_d * h**2,
                e: root_e * h**2,
            }
        )
    ).as_numer_denom()[0]
    scaled_poly = s.Poly(
        numerator,
        root_b,
        root_d,
        root_e,
        root_u,
        h,
    )
    # We work on h != 0.  Clearing the lambda=u/h^2 denominators can leave a
    # common monomial h-factor, which is a localization artifact rather than
    # part of the equation and must be saturated away before measuring degree.
    common_h_order = min(monomial[-1] for monomial, _ in scaled_poly.terms())
    primitive = s.Poly(
        scaled_poly.as_expr() / h**common_h_order,
        root_b,
        root_d,
        root_e,
        root_u,
        h,
    ).primitive()[1]
    root_scaled_degrees.append(primitive.total_degree())
require(max(root_scaled_degrees) == 8, "root-scaled degree ceiling changed")

# Conditional square-lift closure: if R is irreducible, p has five simple
# zeros, so Y^2=1/p is connected and has at least six branch places.
finite_simple_branches = s.degree(L5, v)
minimum_even_branches = finite_simple_branches + finite_simple_branches % 2
minimum_lift_genus = -1 + minimum_even_branches // 2
require((finite_simple_branches, minimum_even_branches, minimum_lift_genus) == (5, 6, 2), "lift genus floor changed")

# The identically y=0 boundary is already contradictory on the open chart.
B, Dcoef, W, y, u, Z = s.symbols("B Dcoef W y u Z")
first_A = 616 * B - 1089 * u + 63 * y**2
first_K = (
    -745360 * B * u * y
    + 6160 * B * y**3
    + 511104 * Dcoef * y
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * first_A * Z + 4 * first_K
require(
    s.expand(N1.subs(y, 0) - 1331 * (616 * B - 1089 * u) * Z) == 0,
    "y=0 first-flux boundary changed",
)

# Finite hostile sweep only: evidence, never used for the exact reduction.
values = sorted(
    {
        s.Rational(numerator, denominator)
        for denominator in (1, 2, 3)
        for numerator in range(-4, 5)
        if numerator != 0
    }
)
tested_parameters = 0
reducible_parameters = []
for lam_value, mu_value in product(values, repeat=2):
    specialized = s.Poly(R.subs({lam: lam_value, mu: mu_value}), p, v)
    factors = s.factor_list(specialized.as_expr())[1]
    if len(factors) > 1 or any(exponent > 1 for _, exponent in factors):
        reducible_parameters.append((lam_value, mu_value))
    tested_parameters += 1
require(tested_parameters == 324, "finite parameter universe changed")
require(not reducible_parameters, "finite hostile sweep found reducible fibre")

print("THM-2617 exact companion")
print("stratum=B,D,W_nonzero;C=E=0")
print("invariants=lambda=D/B^2,mu=W/B^3")
print("resultant_content=2^4*11^6")
print("resultant_terms=34")
print("resultant_bidegree=(5,5)")
print("newton_polygon=full_triangle_degree_5")
print("fixed_section=-567*L5")
print("fixed_section_simple_points=5")
print("line_top_types=2;excluded=2")
print("quadratic_top_types=3;excluded=2;open=1")
print("open_quadratic_type=two_moving_cubic_roots")
print("open_type_top_next_determinant_factors=4")
print("open_type_determinant_walls_closed=2")
print("generic_compatibility_divisor_bidegree=(7,15)")
print("generic_compatibility_divisor_terms=69")
print("root_scaled_equation_max_total_degree=8")
print("conditional_square_lift_minimum_genus=2")
print("finite_hostile_parameters=324")
print("finite_hostile_reducible=0")
print("uniform_irreducibility=OPEN")
print("ALL CHECKS PASSED")
