#!/usr/bin/env python3
"""Exact companion for THM-2463."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


p, v, zeta, lam = s.symbols("p v zeta lam")
a, b, c, d, e = s.symbols("a b c d e")

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
    -1978994688 * lam * p**2 * v
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
R = s.Poly(raw_resultant, p, v, lam).primitive()[1].as_expr()
require(s.cancel(raw_resultant / R) == 2**4 * 11**6, "resultant content changed")
require(len(s.Poly(R, p, v, lam).terms()) == 28, "resultant support changed")
require(s.degree(R, p) == 4 and s.degree(R, v) == 5, "resultant bidegree changed")

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
constant_in_p = s.Poly(R, p).coeff_monomial(1)
linear_in_p = s.Poly(R, p).coeff_monomial(p)
require(s.expand(constant_in_p + 567 * L5) == 0, "fixed quintic section changed")
require(s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0, "L5 lost squarefreeness")
require(
    s.gcd(s.Poly(constant_in_p, v), s.Poly(linear_in_p, v)).degree() == 0,
    "resultant acquired a vertical factor",
)
require(
    s.expand(s.diff(R, v).subs(p, 0) + 567 * s.diff(L5, v)) == 0,
    "p=0 smooth-section derivative changed",
)

top = s.Integer(0)
for (ip, iv, il), coefficient in s.Poly(R, p, v, lam).terms():
    if ip + iv == 5:
        top += coefficient * p**ip * v**iv * lam**il
expected_top = -38974342 * v * (56 * p - 99 * v) ** 2 * (
    384 * lam * p**2 - 280 * p * v + 231 * v**2
)
require(s.expand(top - expected_top) == 0, "top homogeneous factorization changed")

bottom_right = s.Poly(R, p, v).coeff_monomial(p**4)
require(
    s.factor(bottom_right) == 71243808768 * lam * (33 * lam - 49),
    "Newton exceptional parameter changed",
)
require(
    s.Poly(R.subs(lam, s.Rational(49, 33)), p, v).coeff_monomial(p**3)
    == -627838790656,
    "exceptional Newton polygon lost its (3,0) vertex",
)


def is_unit_basis(groebner_basis: s.GroebnerBasis) -> bool:
    return len(groebner_basis.polys) == 1 and groebner_basis.polys[0].as_expr() == 1


# A linear Minkowski summand has one of the two shapes below.
linear_e0 = p + a * v + b
remainder_e0 = s.Poly(s.rem(s.Poly(R, p), s.Poly(linear_e0, p)).as_expr(), v)
gb_linear_e0 = s.groebner(
    [coefficient for _, coefficient in remainder_e0.terms()],
    a,
    b,
    lam,
    order="grevlex",
)
require(is_unit_basis(gb_linear_e0), "affine linear-factor ideal stopped being the unit ideal")

denominator = v + a
numerator = b * v**2 + c * v + d
cleared_e1 = s.cancel(denominator**4 * R.subs(p, -numerator / denominator))
require(cleared_e1.as_numer_denom()[1] == 1, "rational-root denominator did not clear")
remainder_e1 = s.Poly(s.expand(cleared_e1), v)
gb_linear_e1 = s.groebner(
    [coefficient for _, coefficient in remainder_e1.terms()],
    a,
    b,
    c,
    d,
    lam,
    order="grevlex",
)
require(is_unit_basis(gb_linear_e1), "degree-(2,1) rational-root ideal stopped being the unit ideal")

# In a 2+2 split exactly one factor has triangular Newton polygon Delta(2,0).
quadratic = p**2 + (a * v + b) * p + c * v**2 + d * v + e
quadratic_remainder = s.Poly(
    s.rem(s.Poly(R, p), s.Poly(quadratic, p)).as_expr(),
    p,
    v,
)
quadratic_equations = [coefficient for _, coefficient in quadratic_remainder.terms()]
remainder_by_monomial = dict(quadratic_remainder.terms())

# Type I: (56p-99v)^2.
type_i = [
    s.factor(
        equation.subs(
            {
                a: -s.Rational(99, 28),
                c: s.Rational(9801, 3136),
            }
        )
    )
    for equation in quadratic_equations
]
gb_type_i = s.groebner(type_i, b, d, e, lam, order="grevlex")
require(is_unit_basis(gb_type_i), "quadratic top type I stopped being impossible")

# Type II: the whole quadratic 384*lambda*p^2-280*p*v+231*v^2.
type_ii_substitution = {
    a: -s.Rational(35, 48) / lam,
    c: s.Rational(77, 128) / lam,
}
type_ii = [
    s.factor(equation.subs(type_ii_substitution))
    for equation in quadratic_equations
]
lambda_exception = s.Rational(196, 891)
gb_type_ii_exception = s.groebner(
    [s.factor(equation.subs(lam, lambda_exception)) for equation in type_ii],
    b,
    d,
    e,
    order="grevlex",
)
require(
    is_unit_basis(gb_type_ii_exception),
    "quadratic top type II exceptional fibre stopped being impossible",
)

q13 = s.factor(remainder_by_monomial[(1, 3)].subs(type_ii_substitution))
q04 = s.factor(remainder_by_monomial[(0, 4)].subs(type_ii_substitution))
coefficient_matrix = s.Matrix(
    [
        [s.diff(q13, b), s.diff(q13, d)],
        [s.diff(q04, b), s.diff(q04, d)],
    ]
)
require(
    s.factor(coefficient_matrix.det())
    == 27102253467051586289664 * (891 * lam - 196) ** 2,
    "type-II top-next determinant changed",
)
b_solution = -s.Rational(7) * (
    78682428 * lam**2 - 34767117 * lam + 3841600
) / (23232 * lam * (891 * lam - 196) ** 2)
d_solution = s.Rational(147) * (
    2822688 * lam**2 - 1246014 * lam + 137543
) / (11264 * lam * (891 * lam - 196) ** 2)
require(
    s.factor(q13.subs({b: b_solution, d: d_solution})) == 0
    and s.factor(q04.subs({b: b_solution, d: d_solution})) == 0,
    "type-II top-next solution changed",
)
type_ii_generic_numerators = [
    s.together(equation.subs({b: b_solution, d: d_solution})).as_numer_denom()[0]
    for equation in type_ii
]
gb_type_ii_generic = s.groebner(
    type_ii_generic_numerators,
    e,
    lam,
    order="grevlex",
)
require(
    is_unit_basis(gb_type_ii_generic),
    "quadratic top type II generic branch stopped being impossible",
)

# Type III: one repeated linear root and one root h of the quadratic.
h, x = s.symbols("h x")
lambda_h = (280 * h - 231) / (384 * h**2)
type_iii_substitution = {
    lam: lambda_h,
    a: -(s.Rational(99, 56) + h),
    c: s.Rational(99, 56) * h,
}
type_iii = [
    s.factor(equation.subs(type_iii_substitution))
    for equation in quadratic_equations
]
type_iii_numerators = [
    s.Poly(
        s.together(equation).as_numer_denom()[0],
        b,
        d,
        e,
        h,
    ).primitive()[1].as_expr()
    for equation in type_iii
]
q13_iii = s.factor(remainder_by_monomial[(1, 3)].subs(type_iii_substitution))
q04_iii = s.factor(remainder_by_monomial[(0, 4)].subs(type_iii_substitution))
q13_x = s.factor(q13_iii.subs(d, x - b * h))
q04_x = s.factor(q04_iii.subs(d, x - b * h))
compatibility = s.factor(
    s.Poly(q13_x, x).coeff_monomial(x)
    * s.Poly(q04_x, x).coeff_monomial(1)
    - s.Poly(q04_x, x).coeff_monomial(x)
    * s.Poly(q13_x, x).coeff_monomial(1)
)
expected_compatibility = (
    s.Rational(191848476676752603477, 16)
    * (20 * h - 33)
    * (448 * h**2 - 1320 * h + 1089) ** 2
    / h**5
)
require(
    s.factor(compatibility - expected_compatibility) == 0,
    "type-III compatibility divisor changed",
)
gb_type_iii_rational = s.groebner(
    [equation.subs(h, s.Rational(33, 20)) for equation in type_iii_numerators],
    b,
    d,
    e,
    order="grevlex",
)
require(
    is_unit_basis(gb_type_iii_rational),
    "type-III rational compatibility branch stopped being impossible",
)
type_iii_minpoly = 448 * h**2 - 1320 * h + 1089
gb_type_iii_quadratic = s.groebner(
    type_iii_numerators + [type_iii_minpoly],
    b,
    d,
    e,
    h,
    order="grevlex",
)
require(
    is_unit_basis(gb_type_iii_quadratic),
    "type-III quadratic compatibility branch stopped being impossible",
)

# The y=0 boundary of the original first flux.
B, D, y, u, Z = s.symbols("B D y u Z")
first_flux_A = 616 * B - 1089 * u + 63 * y**2
first_flux_K = (
    -745360 * B * u * y
    + 6160 * B * y**3
    + 511104 * D * y
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * first_flux_A * Z + 4 * first_flux_K
require(
    s.expand(N1.subs(y, 0) - 1331 * (616 * B - 1089 * u) * Z) == 0,
    "y=0 first-flux boundary changed",
)

# Five simple zeros force at least six branch points on the connected
# quadratic lift Y^2=1/p.  Riemann-Hurwitz then gives genus at least two.
finite_simple_branch_points = s.degree(L5, v)
minimum_even_branch_count = finite_simple_branch_points + (
    finite_simple_branch_points % 2
)
minimum_lift_genus = -1 + minimum_even_branch_count // 2
require(finite_simple_branch_points == 5, "square-lift simple-zero count changed")
require(minimum_even_branch_count == 6, "square-lift parity floor changed")
require(minimum_lift_genus == 2, "square-lift genus floor changed")

print("THM-2463 exact companion")
print("resultant_content=2^4*11^6")
print("resultant_terms=28")
print("resultant_bidegree=(4,5)")
print("newton_generic=Delta(4,1)")
print("newton_exception=lambda=49/33")
print("top_factor=v*(56p-99v)^2*(384lambda p^2-280pv+231v^2)")
print("linear_factor_shapes=2")
print("vertical_factors=NONE")
print("linear_factor_unit_ideals=2")
print("quadratic_top_types=3")
print("quadratic_factor_unit_ideals=5")
print("type_III_compatibility=(20h-33)*(448h^2-1320h+1089)^2")
print("p_zero_section=-567*L5")
print("L5_degree=5")
print("L5_squarefree=YES")
print("square_lift_minimum_branch_points=6")
print("square_lift_minimum_genus=2")
print("y_zero_boundary=Z_zero_on_open_first_flux")
print("ALL CHECKS PASSED")
