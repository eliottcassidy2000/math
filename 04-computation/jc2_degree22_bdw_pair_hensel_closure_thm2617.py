#!/usr/bin/env python3
"""Exact fixed-section pair-field Hensel closure for the BDW eliminant."""

from __future__ import annotations

from functools import lru_cache

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


p, v, zeta, lam, mu = s.symbols("p v zeta lam mu")
alpha, beta = s.symbols("alpha beta")

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
R = s.Poly(s.resultant(f1, f2, zeta), p, v, lam, mu).primitive()[1].as_expr()
ris = [s.Poly(R, p).coeff_monomial(p**i) for i in range(6)]
L5 = -ris[0] / 567
require(s.Poly(L5, v).degree() == 5, "fixed section changed")
require(
    s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0,
    "fixed section stopped being squarefree",
)
require(ris[1].has(lam, mu) is False, "first deformation stopped being fixed")

# A monic quadratic divisor qhat=v^2+alpha*v+beta of the squarefree fixed
# quintic is a point of a ten-point pair-of-roots scheme.
qhat_expr = v**2 + alpha * v + beta
pair_remainder = s.Poly(
    s.rem(s.Poly(L5, v), s.Poly(qhat_expr, v)).as_expr(),
    v,
)
pair_equations = [coefficient for _, coefficient in pair_remainder.terms()]
pair_basis = s.groebner(pair_equations, beta, alpha, order="lex")
require(len(pair_basis.polys) == 2, "pair scheme lost its triangular basis")
beta_relation = s.Poly(pair_basis.polys[0].as_expr(), beta, alpha)
alpha_modulus = s.Poly(pair_basis.polys[1].as_expr(), alpha, domain=s.QQ).monic()
beta_as_univariate = s.Poly(beta_relation.as_expr(), beta)
beta_linear_coefficient = beta_as_univariate.coeff_monomial(beta)
require(
    beta_relation.degree(beta) == 1
    and s.degree(beta_linear_coefficient, alpha) == 0
    and beta_linear_coefficient != 0
    and alpha_modulus.degree() == 10
    and alpha_modulus.is_irreducible,
    "fixed-section pair scheme stopped being one degree-ten field",
)
require(
    s.gcd(alpha_modulus, alpha_modulus.diff()).degree() == 0,
    "pair field modulus stopped being reduced",
)


@lru_cache(maxsize=None)
def ared_cached(expression: s.Expr) -> s.Expr:
    expression = s.cancel(expression)
    numerator, denominator = expression.as_numer_denom()
    numerator_poly = s.Poly(numerator, alpha, domain=s.QQ).rem(alpha_modulus)
    denominator_poly = s.Poly(denominator, alpha, domain=s.QQ).rem(alpha_modulus)
    require(not denominator_poly.is_zero, "zero denominator in pair field")
    inverse = s.invert(denominator_poly, alpha_modulus)
    return s.Poly(
        numerator_poly.as_expr() * inverse,
        alpha,
        domain=s.QQ,
    ).rem(alpha_modulus).as_expr()


def ared(expression: s.Expr | int) -> s.Expr:
    return ared_cached(s.sympify(expression))


def azero(expression: s.Expr | int) -> bool:
    return s.Poly(ared(expression), alpha, domain=s.QQ).is_zero


def ainv(expression: s.Expr | int) -> s.Expr:
    expression = ared(expression)
    require(not azero(expression), "attempted to invert zero in pair field")
    return ared(1 / expression)


def field_numerator(expression: s.Expr | int) -> s.Poly:
    return s.Poly(
        s.cancel(ared(expression)).as_numer_denom()[0],
        alpha,
        domain=s.QQ,
    ).primitive()[1]


def coprime_to_pair_modulus(expression: s.Expr | int) -> bool:
    return s.gcd(alpha_modulus, field_numerator(expression)).degree() == 0


beta_constant = beta_relation.as_expr().subs(beta, 0)
beta_value = ared(-beta_constant / beta_linear_coefficient)


# Small polynomial arithmetic in A[v], where
# A=Q[alpha]/(alpha_modulus).  Coefficients are stored low degree first.
def ptrim(poly: list[s.Expr]) -> list[s.Expr]:
    poly = [ared(coefficient) for coefficient in poly]
    while len(poly) > 1 and azero(poly[-1]):
        poly.pop()
    return poly


def padd(left: list[s.Expr], right: list[s.Expr]) -> list[s.Expr]:
    length = max(len(left), len(right))
    return ptrim(
        [
            ared(
                (left[index] if index < len(left) else 0)
                + (right[index] if index < len(right) else 0)
            )
            for index in range(length)
        ]
    )


def pneg(poly: list[s.Expr]) -> list[s.Expr]:
    return [ared(-coefficient) for coefficient in poly]


def pscale(poly: list[s.Expr], scalar: s.Expr | int) -> list[s.Expr]:
    return ptrim([ared(scalar * coefficient) for coefficient in poly])


def pmul(left: list[s.Expr], right: list[s.Expr]) -> list[s.Expr]:
    result = [s.Integer(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] = ared(
                result[i + j] + left_coefficient * right_coefficient
            )
    return ptrim(result)


def pdivmod_monic(
    dividend: list[s.Expr],
    divisor: list[s.Expr],
) -> tuple[list[s.Expr], list[s.Expr]]:
    divisor = ptrim(divisor)
    require(ared(divisor[-1]) == 1, "divisor is not monic")
    remainder = ptrim(dividend[:])
    quotient = [s.Integer(0)] * max(1, len(remainder) - len(divisor) + 1)
    while len(remainder) >= len(divisor) and not (
        len(remainder) == 1 and azero(remainder[0])
    ):
        shift = len(remainder) - len(divisor)
        leading = remainder[-1]
        quotient[shift] = ared(quotient[shift] + leading)
        subtraction = [s.Integer(0)] * shift + pscale(divisor, leading)
        remainder = padd(remainder, pneg(subtraction))
    return ptrim(quotient), ptrim(remainder)


def premainder(dividend: list[s.Expr], divisor: list[s.Expr]) -> list[s.Expr]:
    return pdivmod_monic(dividend, divisor)[1]


def pzero(poly: list[s.Expr]) -> bool:
    return all(azero(coefficient) for coefficient in poly)


def pequal(left: list[s.Expr], right: list[s.Expr]) -> bool:
    return pzero(padd(left, pneg(right)))


def fixed_vpoly(expression: s.Expr) -> list[s.Expr]:
    poly = s.Poly(expression, v, domain=s.QQ)
    return ptrim([ared(poly.nth(index)) for index in range(poly.degree() + 1)])


def evaluated_vpoly(
    expression: s.Expr,
    lambda_value: s.Expr,
    mu_value: s.Expr,
) -> list[s.Expr]:
    poly = s.Poly(expression, v, lam, mu, domain=s.QQ)
    degree_v = max(monomial[0] for monomial, _ in poly.terms())
    result = [s.Integer(0)] * (degree_v + 1)
    for (iv, ilambda, imu), coefficient in poly.terms():
        result[iv] = ared(
            result[iv]
            + coefficient * lambda_value**ilambda * mu_value**imu
        )
    return ptrim(result)


def solve_two_by_two(
    column_1: list[s.Expr],
    column_2: list[s.Expr],
    target: list[s.Expr],
) -> tuple[s.Expr, s.Expr, s.Expr]:
    column_1 = column_1 + [s.Integer(0)] * (2 - len(column_1))
    column_2 = column_2 + [s.Integer(0)] * (2 - len(column_2))
    target = target + [s.Integer(0)] * (2 - len(target))
    determinant = ared(
        column_1[0] * column_2[1] - column_2[0] * column_1[1]
    )
    require(not azero(determinant), "two-by-two Hensel solve became singular")
    first = ared(
        (target[0] * column_2[1] - column_2[0] * target[1])
        * ainv(determinant)
    )
    second = ared(
        (column_1[0] * target[1] - target[0] * column_1[1])
        * ainv(determinant)
    )
    return first, second, determinant


qhat = [beta_value, alpha, s.Integer(1)]
L5_poly = fixed_vpoly(L5)
shat, fixed_remainder = pdivmod_monic(L5_poly, qhat)
require(pzero(fixed_remainder), "pair field does not divide L5")
require(
    len(qhat) - 1 == 2
    and len(shat) - 1 == 3
    and pequal(pmul(qhat, shat), L5_poly),
    "quadratic/cubic factor-cofactor control failed",
)

# First Hensel coefficient.  Put q1=c*g, g=g_a*v+g_b.  The fixed r1
# determines g uniquely because L5 is squarefree.
r1_poly = fixed_vpoly(ris[1])
r1_remainder = premainder(r1_poly, qhat)
ga_column = premainder(pscale(pmul([0, 1], shat), 567), qhat)
gb_column = premainder(pscale(shat, 567), qhat)
ga, gb, first_determinant = solve_two_by_two(
    ga_column,
    gb_column,
    pneg(r1_remainder),
)
g = [gb, ga]
t1_numerator = padd(r1_poly, pscale(pmul(g, shat), 567))
t1, t1_remainder = pdivmod_monic(t1_numerator, qhat)
require(pzero(t1_remainder), "first Hensel quotient has a remainder")

# Second coefficient.  With x=1/c, its two remainder equations are linear
# in x and lambda and determine both over the pair field.
r2_base = fixed_vpoly(ris[2].subs(lam, 0))
r2_lambda = fixed_vpoly(s.diff(ris[2], lam))
n2_base = padd(r2_base, pneg(pmul(g, t1)))
x_column = premainder(pscale(shat, 567), qhat)
lambda_column = premainder(r2_lambda, qhat)
n2_remainder = premainder(n2_base, qhat)
x_value, lambda_value, second_determinant = solve_two_by_two(
    x_column,
    lambda_column,
    pneg(n2_remainder),
)
require(not azero(x_value), "the forced reciprocal scale x=1/c vanished")
n2 = padd(
    padd(n2_base, pscale(shat, 567 * x_value)),
    pscale(r2_lambda, lambda_value),
)
t2, t2_remainder = pdivmod_monic(n2, qhat)
require(pzero(t2_remainder), "second Hensel quotient has a remainder")

# Third coefficient.  The two remainder equations must agree on one mu.
r3_without_mu = evaluated_vpoly(ris[3].subs(mu, 0), lambda_value, 0)
r3_mu = fixed_vpoly(s.diff(ris[3], mu))
n3_base = padd(
    padd(r3_without_mu, pneg(pmul(g, t2))),
    pneg(pscale(t1, x_value)),
)
n3_remainder = premainder(n3_base, qhat)
mu_column = premainder(r3_mu, qhat)
n3_remainder += [s.Integer(0)] * (2 - len(n3_remainder))
mu_column += [s.Integer(0)] * (2 - len(mu_column))
mu_pivot = next(index for index in range(2) if not azero(mu_column[index]))
mu_value = ared(-n3_remainder[mu_pivot] * ainv(mu_column[mu_pivot]))
third_compatibility = [
    ared(n3_remainder[index] + mu_value * mu_column[index])
    for index in range(2)
]

# Hostile scaling control.  Reconstruct q0,q1 and the first three cofactor
# coefficients with c=1/x, then check the product against r0,r1,r2 directly.
c_value = ainv(x_value)
q0 = pscale(qhat, c_value)
q1 = pscale(g, c_value)
s0 = pscale(shat, -567 * x_value)
s1 = pscale(t1, x_value)
s2 = pscale(t2, x_value)
r0_poly = pscale(L5_poly, -567)
r2_poly = padd(r2_base, pscale(r2_lambda, lambda_value))
require(pequal(pmul(q0, s0), r0_poly), "p^0 product scaling check failed")
require(
    pequal(padd(pmul(q0, s1), pmul(q1, s0)), r1_poly),
    "p^1 product scaling check failed",
)
require(
    pequal(
        padd(padd(pmul(q0, s2), pmul(q1, s1)), s0),
        r2_poly,
    ),
    "p^2 product scaling check failed",
)

r3_poly = evaluated_vpoly(ris[3], lambda_value, mu_value)
p3_direct_remainder = premainder(
    padd(r3_poly, pneg(padd(pmul(q1, s2), s1))),
    qhat,
)
require(
    pequal(p3_direct_remainder, third_compatibility),
    "direct p^3 remainder disagrees with the Hensel recursion",
)

nonzero_compatibility = [
    entry for entry in third_compatibility if not azero(entry)
]
require(
    len(nonzero_compatibility) == 1,
    "p^3 compatibility no longer has one exact obstruction",
)
obstruction = ared(nonzero_compatibility[0])
obstruction_numerator = field_numerator(obstruction)
solve_units = (
    first_determinant,
    second_determinant,
    x_value,
    mu_column[mu_pivot],
)
require(
    all(coprime_to_pair_modulus(unit) for unit in solve_units)
    and all(
        field_numerator(unit).degree() == 9
        and len(field_numerator(unit).terms()) == 10
        for unit in solve_units
    ),
    "a Hensel solve unit stopped being invertible on the pair scheme",
)
require(
    obstruction_numerator.degree() == 9
    and len(obstruction_numerator.terms()) == 10
    and coprime_to_pair_modulus(obstruction),
    "p^3 obstruction stopped being nonzero on every root pair",
)

if all(azero(entry) for entry in third_compatibility):
    n3 = padd(n3_base, pscale(r3_mu, mu_value))
    t3, t3_remainder = pdivmod_monic(n3, qhat)
    require(pzero(t3_remainder), "third Hensel quotient has a remainder")
    r4_poly = evaluated_vpoly(ris[4], lambda_value, mu_value)
    r5_poly = evaluated_vpoly(ris[5], lambda_value, mu_value)
    fourth_obstruction = padd(
        padd(r4_poly, pneg(pmul(g, t3))),
        pneg(pscale(t2, x_value)),
    )
    fifth_obstruction = padd(r5_poly, pneg(pscale(t3, x_value)))
else:
    t3 = []
    fourth_obstruction = []
    fifth_obstruction = []

print("BDW pair-field Hensel closure")
print(f"pair_field_degree={alpha_modulus.degree()}")
print(f"pair_field_irreducible={alpha_modulus.is_irreducible}")
print(f"pair_field_reduced={s.gcd(alpha_modulus, alpha_modulus.diff()).degree() == 0}")
print("unordered_root_pairs=10")
print(f"pair_complement_degree={len(shat) - 1}")
print("factor_cofactor_swap=quadratic_side_selected")
print(f"first_solve_determinant_nonzero={not azero(first_determinant)}")
print(f"second_solve_determinant_nonzero={not azero(second_determinant)}")
print(f"forced_reciprocal_scale_nonzero={not azero(x_value)}")
print("hensel_solve_units_coprime_to_pair_modulus=True")
print("product_scaling_checks=p0,p1,p2")
print("third_compatibility_zero=False")
print(f"third_obstruction_degree={obstruction_numerator.degree()}")
print(f"third_obstruction_terms={len(obstruction_numerator.terms())}")
print("third_obstruction_coprime_to_pair_modulus=True")
if t3:
    print(f"fourth_obstruction_zero={pzero(fourth_obstruction)}")
    print(f"fifth_obstruction_zero={pzero(fifth_obstruction)}")
    print(
        "fourth_nonzero_coefficients="
        f"{sum(not azero(entry) for entry in fourth_obstruction)}"
    )
    print(
        "fifth_nonzero_coefficients="
        f"{sum(not azero(entry) for entry in fifth_obstruction)}"
    )
print("quadratic_factor_scheme=EMPTY")
print("ALL CHECKS PASSED")
