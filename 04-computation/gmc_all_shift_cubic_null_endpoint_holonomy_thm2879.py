#!/usr/bin/env python3
"""Exact companion for THM-2879.

For

    L(s^q)=q!,  f_j=s^j/j!,  d_j=f_(j+1)-f_j,

consider the THM-2855 family

    U=d_n+x*d_(n+2),       V=d_(n+1)+y*d_(n+2).

THM-2855 proves that, for every integer n>=1, there is a unique positive
pair (x_n,y_n) at which the binary quadratic moment divides the cubic.
This companion proves that the quartic endpoint-holonomy determinant is
strictly negative at every one of those points.

All truth-bearing checks use explicit ``require`` calls, so optimized and
ordinary Python execute the same audit.
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


n, x, y = sp.symbols("n x y")


def factorial_ratio(base: sp.Expr, shift: int) -> sp.Expr:
    """Return (base+shift)!/base! as a rational product."""

    if shift >= 0:
        return sp.prod(
            (base + index for index in range(1, shift + 1)),
            start=sp.Integer(1),
        )
    return 1 / sp.prod(
        (base - index + 1 for index in range(1, -shift + 1)),
        start=sp.Integer(1),
    )


def direct_adjacent_tensor(indices: tuple[int, ...]) -> sp.Integer:
    """Independent 2^k expansion of L(prod d_i)."""

    value = sp.Integer(0)
    for mask in range(1 << len(indices)):
        degrees: list[int] = []
        sign = 1
        for position, index in enumerate(indices):
            if mask & (1 << position):
                degrees.append(index + 1)
            else:
                degrees.append(index)
                sign = -sign
        term = sp.factorial(sum(degrees))
        for degree in degrees:
            term /= sp.factorial(degree)
        value += sign * term
    require(value.q == 1, "direct adjacent tensor lost integrality")
    return sp.Integer(value)


def normalized_adjacent_tensor(offsets: tuple[int, ...]) -> sp.Expr:
    """Adjacent tensor divided by L(d_(n+2)^k), for offsets in {0,1,2}."""

    order = len(offsets)
    value = sp.Integer(0)
    for mask in range(1 << order):
        increments = tuple((mask >> position) & 1 for position in range(order))
        increment_sum = sum(increments)
        term = sp.Integer((-1) ** (order - increment_sum))
        term *= factorial_ratio(
            order * n + 2 * order,
            sum(offsets) + increment_sum - 2 * order,
        )
        for offset, increment in zip(offsets, increments):
            term *= factorial_ratio(
                n + offset + increment,
                2 - offset - increment,
            )
        value += term
    return sp.cancel(value)


lower = {0: sp.Integer(1), 2: x}
upper = {1: sp.Integer(1), 2: y}


def multilinear(forms: tuple[dict[int, sp.Expr], ...]) -> sp.Expr:
    return sp.factor(
        sum(
            sp.prod((coefficient for _, coefficient in terms), start=sp.Integer(1))
            * normalized_adjacent_tensor(
                tuple(index for index, _ in terms)
            )
            for terms in product(*(tuple(form.items()) for form in forms))
        )
    )


g0 = multilinear((lower, lower))
g1 = multilinear((lower, upper))
g2 = multilinear((upper, upper))

t0 = multilinear((lower, lower, lower))
t1 = multilinear((lower, lower, upper))
t2 = multilinear((lower, upper, upper))
t3 = multilinear((upper, upper, upper))

quartic = tuple(
    multilinear(
        tuple([lower] * (4 - upper_count) + [upper] * upper_count)
    )
    for upper_count in range(5)
)

cubic_one = sp.factor(
    3 * t1 * g0 * g2 - t3 * g0**2 - 2 * t0 * g1 * g2
)
cubic_two = sp.factor(
    3 * t2 * g0 * g2 - 2 * t3 * g1 * g0 - t0 * g2**2
)

endpoint_holonomy = sp.factor(
    (2 * quartic[1] * g0 - quartic[0] * g1) * g2**2
    - (2 * quartic[3] * g2 - quartic[4] * g1) * g0**2
)


def all_coefficients_have_sign(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    sign: int,
) -> bool:
    coefficients = sp.Poly(expression, *variables).coeffs()
    if sign > 0:
        return all(bool(coefficient > 0) for coefficient in coefficients)
    return all(bool(coefficient < 0) for coefficient in coefficients)


def sign_variations(values: tuple[sp.Expr, ...] | list[sp.Expr]) -> int:
    nonzero = [value for value in values if value != 0]
    return sum(
        bool(left > 0) != bool(right > 0)
        for left, right in zip(nonzero, nonzero[1:])
    )


def canonical_digest(polynomial: sp.Poly) -> str:
    records = "\n".join(
        f"{','.join(str(exponent) for exponent in monomial)}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )
    return sha256((records + "\n").encode()).hexdigest()


def direct_moment(
    forms: tuple[dict[int, sp.Expr], ...],
) -> sp.Expr:
    """Literal factorial tensor, independent of the normalized formulas."""

    return sp.expand(
        sum(
            sp.prod((coefficient for _, coefficient in terms), start=sp.Integer(1))
            * direct_adjacent_tensor(tuple(index for index, _ in terms))
            for terms in product(*(tuple(form.items()) for form in forms))
        )
    )


def fixed_direct_exit(depth: int) -> tuple[int, int, int]:
    """Independent fixed-depth elimination and holonomy control."""

    X, Y = sp.symbols("X Y")
    first = {depth: sp.Integer(1), depth + 2: X}
    second = {depth + 1: sp.Integer(1), depth + 2: Y}

    fixed_g = (
        direct_moment((first, first)),
        direct_moment((first, second)),
        direct_moment((second, second)),
    )
    fixed_t = (
        direct_moment((first, first, first)),
        direct_moment((first, first, second)),
        direct_moment((first, second, second)),
        direct_moment((second, second, second)),
    )
    fixed_four = tuple(
        direct_moment(
            tuple([first] * (4 - upper_count) + [second] * upper_count)
        )
        for upper_count in range(5)
    )

    fixed_one = sp.expand(
        3 * fixed_t[1] * fixed_g[0] * fixed_g[2]
        - fixed_t[3] * fixed_g[0] ** 2
        - 2 * fixed_t[0] * fixed_g[1] * fixed_g[2]
    )
    fixed_two = sp.expand(
        3 * fixed_t[2] * fixed_g[0] * fixed_g[2]
        - 2 * fixed_t[3] * fixed_g[1] * fixed_g[0]
        - fixed_t[0] * fixed_g[2] ** 2
    )
    fixed_j = sp.expand(
        (2 * fixed_four[1] * fixed_g[0] - fixed_four[0] * fixed_g[1])
        * fixed_g[2] ** 2
        - (
            2 * fixed_four[3] * fixed_g[2]
            - fixed_four[4] * fixed_g[1]
        )
        * fixed_g[0] ** 2
    )

    fixed_subresultants = sp.subresultants(fixed_one, fixed_two, X)
    require(
        [sp.degree(item, X) for item in fixed_subresultants]
        == [4, 3, 2, 1, 0],
        f"fixed depth {depth}: subresultant profile changed",
    )
    factors = sp.factor_list(sp.Poly(fixed_subresultants[-1], Y))[1]
    require(
        sorted((factor.degree(), exponent) for factor, exponent in factors)
        == [(2, 2), (15, 1)],
        f"fixed depth {depth}: resultant factor profile changed",
    )
    eliminant = next(
        factor.as_expr() for factor, _ in factors if factor.degree() == 15
    )
    if sp.Poly(eliminant, Y).LC() < 0:
        eliminant = -eliminant
    eliminant_coefficients = sp.Poly(eliminant, Y).all_coeffs()
    require(
        sign_variations(eliminant_coefficients) == 1
        and eliminant_coefficients[0] > 0
        and eliminant_coefficients[-1] < 0,
        f"fixed depth {depth}: positive branch is not unique",
    )

    linear = sp.Poly(fixed_subresultants[-2], X)
    linear_coefficient = linear.nth(1)
    constant_coefficient = linear.nth(0)
    content = sp.gcd(linear_coefficient, constant_coefficient)
    selector_a = sp.cancel(linear_coefficient / content)
    selector_n = sp.cancel(-constant_coefficient / content)
    if selector_a.subs(Y, 1) < 0:
        selector_a = -selector_a
        selector_n = -selector_n
    require(
        all_coefficients_have_sign(selector_a, (Y,), 1)
        and all_coefficients_have_sign(selector_n, (Y,), 1),
        f"fixed depth {depth}: selector positivity changed",
    )

    fixed_j_poly = sp.Poly(fixed_j, X)
    degree = fixed_j_poly.degree()
    selector_a_poly = sp.Poly(selector_a, Y, domain=sp.QQ)
    selector_n_poly = sp.Poly(selector_n, Y, domain=sp.QQ)
    powers_a = [sp.Poly(1, Y, domain=sp.QQ)]
    powers_n = [sp.Poly(1, Y, domain=sp.QQ)]
    for _ in range(degree):
        powers_a.append(powers_a[-1] * selector_a_poly)
        powers_n.append(powers_n[-1] * selector_n_poly)
    substituted_numerator = sum(
        (
            sp.Poly(fixed_j_poly.nth(index), Y, domain=sp.QQ)
            * powers_n[index]
            * powers_a[degree - index]
            for index in range(degree + 1)
        ),
        sp.Poly(0, Y, domain=sp.QQ),
    )
    remainder = sp.Poly(
        sp.prem(substituted_numerator.as_expr(), eliminant, Y),
        Y,
        domain=sp.QQ,
    )
    require(
        remainder.degree() == 14
        and all(bool(coefficient < 0) for coefficient in remainder.coeffs()),
        f"fixed depth {depth}: endpoint remainder lost negativity",
    )
    return (
        len(sp.Poly(selector_a, Y).terms()),
        len(sp.Poly(selector_n, Y).terms()),
        len(remainder.terms()),
    )


def main() -> None:
    # Independent literal-tensor controls for every offset profile used.
    tensor_controls = 0
    scale_controls = 0
    for depth in (1, 2, 3):
        for order in (2, 3, 4):
            base = sp.Rational(
                factorial(order * (depth + 2)),
                factorial(depth + 2) ** order,
            )
            for offsets in product(range(3), repeat=order):
                normalized = normalized_adjacent_tensor(offsets).subs(n, depth)
                direct = direct_adjacent_tensor(
                    tuple(depth + offset for offset in offsets)
                )
                require(
                    sp.cancel(normalized * base - direct) == 0,
                    "normalized/direct tensor mismatch",
                )
                tensor_controls += 1

        # The endpoint determinant uses one quartic and three quadratic
        # factors, so its normalization scale is positive B4*B2^3.
        rational_point = {x: sp.Rational(2, 7), y: sp.Rational(1, 9)}
        first = {
            depth: sp.Integer(1),
            depth + 2: rational_point[x],
        }
        second = {
            depth + 1: sp.Integer(1),
            depth + 2: rational_point[y],
        }
        actual_g = (
            direct_moment((first, first)),
            direct_moment((first, second)),
            direct_moment((second, second)),
        )
        actual_four = tuple(
            direct_moment(
                tuple([first] * (4 - upper_count) + [second] * upper_count)
            )
            for upper_count in range(5)
        )
        actual_j = (
            (2 * actual_four[1] * actual_g[0] - actual_four[0] * actual_g[1])
            * actual_g[2] ** 2
            - (
                2 * actual_four[3] * actual_g[2]
                - actual_four[4] * actual_g[1]
            )
            * actual_g[0] ** 2
        )
        base_two = sp.Rational(
            factorial(2 * (depth + 2)),
            factorial(depth + 2) ** 2,
        )
        base_four = sp.Rational(
            factorial(4 * (depth + 2)),
            factorial(depth + 2) ** 4,
        )
        normalized_j = endpoint_holonomy.subs(
            {n: depth, **rational_point}
        )
        require(
            sp.cancel(
                actual_j - base_four * base_two**3 * normalized_j
            )
            == 0,
            "endpoint normalization scale mismatch",
        )
        scale_controls += 1

    invariant_numerators = tuple(
        sp.together(invariant).as_numer_denom()[0]
        for invariant in (cubic_one, cubic_two)
    )
    invariant_denominators = tuple(
        sp.together(invariant).as_numer_denom()[1]
        for invariant in (cubic_one, cubic_two)
    )
    require(
        all(
            all_coefficients_have_sign(denominator, (n,), 1)
            for denominator in invariant_denominators
        ),
        "a cubic-invariant denominator lost positivity",
    )

    # Recover the exact THM-2855 positive branch.
    subresultants = sp.subresultants(
        invariant_numerators[0], invariant_numerators[1], x
    )
    require(
        [sp.degree(item, x) for item in subresultants]
        == [4, 3, 2, 1, 0],
        "universal subresultant profile changed",
    )
    resultant_factors = sp.factor_list(subresultants[-1])[1]
    y_factors = [
        (factor, exponent)
        for factor, exponent in resultant_factors
        if sp.degree(factor, y) > 0
    ]
    require(
        sorted((sp.degree(factor, y), exponent) for factor, exponent in y_factors)
        == [(2, 2), (15, 1)],
        "universal resultant factor profile changed",
    )
    quadratic_factor = next(
        factor for factor, _ in y_factors if sp.degree(factor, y) == 2
    )
    expected_quadratic = (
        n + 2
        + 2 * (2 * n + 3) * y
        + 2 * (2 * n + 3) * y**2
    )
    quadratic_ratio = sp.cancel(
        sp.Poly(quadratic_factor, y).LC()
        / sp.Poly(expected_quadratic, y).LC()
    )
    require(
        sp.expand(quadratic_factor - quadratic_ratio * expected_quadratic)
        == 0
        and sp.discriminant(expected_quadratic, y) == -4 * (2 * n + 3),
        "positive quadratic factor changed",
    )

    eliminant = next(
        factor for factor, _ in y_factors if sp.degree(factor, y) == 15
    )
    if sp.Poly(eliminant, y).LC().subs(n, 1) < 0:
        eliminant = -eliminant
    eliminant_coefficients = sp.Poly(eliminant, y).all_coeffs()
    require(
        len(eliminant_coefficients) == 16
        and all(
            all_coefficients_have_sign(coefficient, (n,), 1)
            for coefficient in eliminant_coefficients[:3]
        ),
        "permanent-positive eliminant coefficients changed",
    )
    first_negative_depths = (
        102,
        43,
        25,
        16,
        11,
        8,
        6,
        5,
        4,
        3,
        2,
        2,
        1,
    )
    require(
        all(
            left >= right
            for left, right in zip(
                first_negative_depths,
                first_negative_depths[1:],
            )
        ),
        "eliminant threshold staircase changed",
    )
    for coefficient, threshold in zip(
        eliminant_coefficients[3:],
        first_negative_depths,
    ):
        polynomial = sp.Poly(coefficient, n)
        require(
            sign_variations(polynomial.all_coeffs()) == 1
            and coefficient.subs(n, threshold - 1) > 0
            and coefficient.subs(n, threshold) < 0,
            "an eliminant coefficient threshold changed",
        )

    # The n=0 boundary has no positive branch at all.
    bottom_eliminant = sp.Poly(eliminant.subs(n, 0), y)
    require(
        all(bool(coefficient > 0) for coefficient in bottom_eliminant.all_coeffs())
        and bottom_eliminant.count_roots(0, sp.oo) == 0,
        "bottom boundary acquired a positive branch",
    )

    linear = sp.Poly(subresultants[-2], x)
    linear_coefficient = linear.nth(1)
    constant_coefficient = linear.nth(0)
    common_content = sp.gcd(linear_coefficient, constant_coefficient)
    for content_factor, _ in sp.factor_list(common_content)[1]:
        oriented = content_factor
        if sp.Poly(oriented, n, y).LC() < 0:
            oriented = -oriented
        require(
            all_coefficients_have_sign(oriented, (n, y), 1),
            "linear-selector content can vanish in the positive quadrant",
        )
    selector_a = sp.cancel(linear_coefficient / common_content)
    selector_n = sp.cancel(-constant_coefficient / common_content)
    if selector_a.subs({n: 1, y: 1}) < 0:
        selector_a = -selector_a
        selector_n = -selector_n
    require(
        sp.degree(selector_a, y) == 10
        and sp.degree(selector_n, y) == 10
        and len(sp.Poly(selector_a, n, y).terms()) == 251
        and len(sp.Poly(selector_n, n, y).terms()) == 253
        and all_coefficients_have_sign(selector_a, (n, y), 1)
        and all_coefficients_have_sign(selector_n, (n, y), 1),
        "positive universal selector changed",
    )

    # Clear the positive denominator and the positive selector power in J.
    endpoint_numerator, endpoint_denominator = sp.together(
        endpoint_holonomy
    ).as_numer_denom()
    expected_endpoint_denominator = (
        1024
        * (n + 3)
        * (2 * n + 1) ** 2
        * (2 * n + 3) ** 4
        * (4 * n + 1)
        * (4 * n + 3)
        * (4 * n + 5)
        * (4 * n + 7)
    )
    require(
        sp.factor(endpoint_denominator - expected_endpoint_denominator) == 0,
        "endpoint denominator changed",
    )
    endpoint_polynomial = sp.Poly(endpoint_numerator, x)
    endpoint_degree = endpoint_polynomial.degree()
    require(endpoint_degree == 5, "endpoint x-degree changed")

    selector_a_poly = sp.Poly(selector_a, n, y, domain=sp.QQ)
    selector_n_poly = sp.Poly(selector_n, n, y, domain=sp.QQ)
    powers_a = [sp.Poly(1, n, y, domain=sp.QQ)]
    powers_n = [sp.Poly(1, n, y, domain=sp.QQ)]
    for _ in range(endpoint_degree):
        powers_a.append(powers_a[-1] * selector_a_poly)
        powers_n.append(powers_n[-1] * selector_n_poly)
    substituted_numerator = sum(
        (
            sp.Poly(endpoint_polynomial.nth(index), n, y, domain=sp.QQ)
            * powers_n[index]
            * powers_a[endpoint_degree - index]
            for index in range(endpoint_degree + 1)
        ),
        sp.Poly(0, n, y, domain=sp.QQ),
    )
    require(
        substituted_numerator.degree(n) == 121
        and substituted_numerator.degree(y) == 55
        and len(substituted_numerator.terms()) == 6820,
        "cleared endpoint numerator changed: "
        f"{substituted_numerator.degree(n)},"
        f"{substituted_numerator.degree(y)},"
        f"{len(substituted_numerator.terms())}",
    )

    eliminant_lead = sp.Poly(eliminant, y).LC()
    require(
        all_coefficients_have_sign(eliminant_lead, (n,), 1),
        "pseudo-division multiplier can lose positivity",
    )
    pseudo_exponent = (
        substituted_numerator.degree(y) - sp.degree(eliminant, y) + 1
    )
    require(pseudo_exponent == 41, "pseudo-division exponent changed")
    remainder = sp.Poly(
        sp.prem(substituted_numerator.as_expr(), eliminant, y),
        n,
        y,
        domain=sp.QQ,
    )
    remainder_digest = canonical_digest(remainder)
    require(
        remainder.degree(n) == 982
        and remainder.degree(y) == 14
        and len(remainder.terms()) == 14745
        and all(bool(coefficient < 0) for coefficient in remainder.coeffs()),
        "universal endpoint remainder changed",
    )
    require(
        remainder_digest
        == "c1ed1aa0ff682a7226f67c752aceb7bb4924708a2973126fe903c62c686d67a2",
        "universal endpoint remainder digest changed",
    )

    fixed_controls = {
        depth: fixed_direct_exit(depth)
        for depth in (1, 2, 8)
    }

    print("THM-2879 ALL-SHIFT ENDPOINT-HOLONOMY EXIT")
    print("status=PROVED-CANDIDATE+VERIFIED-EXACT")
    print(f"normalized_direct_tensor_controls={tensor_controls}")
    print(f"endpoint_scale_controls={scale_controls}")
    print("cubic_resultant=positive_factors*Gram(y)^2*P15(y)")
    print("positive_cubic_branch=unique for every integer n>=1")
    print("bottom_n0_positive_branch_count=0")
    print("selector=A(n,y)*x-N(n,y); terms=251,253; coefficients=positive")
    print("endpoint_denominator=positive")
    print(
        "pseudo_remainder="
        f"degree_n:{remainder.degree(n)},degree_y:{remainder.degree(y)},"
        f"terms:{len(remainder.terms())},all_coefficients:negative"
    )
    print(f"pseudo_multiplier=LC(P15)^{pseudo_exponent}>0")
    print(f"remainder_sha256={remainder_digest}")
    print(f"fixed_direct_controls={fixed_controls}")
    print("consequence=J_n<0 and beta*kappa_U<alpha*kappa_V for every n>=1")
    print("shared_quartic_line_on_shifted_family=EMPTY")
    print(
        "scope=the THM-2855 shifted high-tooth family only; "
        "general mixed shared line remains open"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
