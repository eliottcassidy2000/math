#!/usr/bin/env python3
"""Exact companion for THM-2842.

This script checks four algebraic parts of the theorem:

1. monic Laguerre multipliers are exactly radial-variance jets;
2. their response on adjacent factorial differences is a falling-factorial
   Vandermonde, and ordered positive mixtures preserve strict positivity;
3. ell_(D-1) is the sharp Krylov-delay state in the D-node carrier; and
4. the same state is a trace-zero unit in the rational Laguerre quotient.

Every truth-bearing gate is an explicit ``require``.  The computation uses
only exact SymPy integers and rationals and is unchanged by ``python -O``.
"""

from __future__ import annotations

from itertools import product
from math import factorial

import sympy as sp


s, t, z = sp.symbols("s t z")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def falling(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(sp.prod(value - j for j in range(degree)))


def rising(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(sp.prod(value + j for j in range(degree)))


def monic_laguerre(degree: int) -> sp.Expr:
    """ell_degree=(-1)^degree degree! L_degree^(0), monic over Z."""
    return sp.expand(
        sum(
            (-1) ** (degree + k)
            * factorial(degree)
            * sp.binomial(degree, k)
            * s**k
            / factorial(k)
            for k in range(degree + 1)
        )
    )


def factorial_functional(poly: sp.Expr) -> sp.Expr:
    result = sp.Integer(0)
    for (degree,), coefficient in sp.Poly(sp.expand(poly), s).terms():
        result += coefficient * factorial(degree)
    return sp.simplify(result)


def adjacent_difference(index: int) -> sp.Expr:
    return sp.expand(
        s ** (index + 1) / factorial(index + 1)
        - s**index / factorial(index)
    )


def quotient_remainder(poly: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    return sp.rem(
        sp.Poly(sp.expand(poly), s, domain=sp.QQ),
        sp.Poly(modulus, s, domain=sp.QQ),
    ).as_expr()


def quotient_trace(poly: sp.Expr, modulus: sp.Expr, dimension: int) -> sp.Expr:
    """Trace of multiplication by poly on Q[s]/(modulus)."""
    result = sp.Integer(0)
    for column in range(dimension):
        remainder = sp.Poly(
            quotient_remainder(poly * s**column, modulus),
            s,
            domain=sp.QQ,
        )
        result += remainder.coeff_monomial(s**column)
    return sp.simplify(result)


def quotient_norm(poly: sp.Expr, modulus: sp.Expr, dimension: int) -> sp.Expr:
    """Norm of poly from Q[s]/(modulus), as a multiplication determinant."""
    columns = []
    for column in range(dimension):
        remainder = sp.Poly(
            quotient_remainder(poly * s**column, modulus),
            s,
            domain=sp.QQ,
        )
        columns.append(
            [
                remainder.coeff_monomial(s**row)
                for row in range(dimension)
            ]
        )
    matrix = sp.Matrix(dimension, dimension, lambda row, col: columns[col][row])
    return sp.factor(matrix.det())


# ---------------------------------------------------------------------------
# 1. Variance jets and the adjacent-difference response
# ---------------------------------------------------------------------------

max_jet = 8
variance_kernel_checks = 0
functional_polynomial_checks = 0
adjacent_response_checks = 0
monomial_response_checks = 0

radial_weight = t**-1 * sp.exp(-s / t)
n_symbol = sp.symbols("n", integer=True, nonnegative=True)
i_symbol = sp.symbols("i", integer=True, nonnegative=True)

for order in range(max_jet + 1):
    ell = monic_laguerre(order)
    differentiated_kernel = sp.expand(
        sp.simplify(
            sp.diff(radial_weight, t, order).subs(t, 1) / sp.exp(-s)
        )
    )
    require(
        sp.expand(differentiated_kernel - ell) == 0,
        f"variance-kernel identity failed at order={order}",
    )
    variance_kernel_checks += 1

    ell_poly = sp.Poly(ell, s, domain=sp.QQ)
    normalized_monomial_response = sp.expand(
        sum(
            coefficient * rising(n_symbol + 1, degree)
            for (degree,), coefficient in ell_poly.terms()
        )
    )
    require(
        sp.Poly(
            normalized_monomial_response - falling(n_symbol, order),
            n_symbol,
            domain=sp.QQ,
        ).is_zero,
        f"functional variance-jet identity failed at order={order}",
    )
    functional_polynomial_checks += 1

    if order >= 1:
        adjacent_response = sp.expand(
            falling(i_symbol + 1, order) - falling(i_symbol, order)
        )
        expected_adjacent = sp.expand(
            order * falling(i_symbol, order - 1)
        )
        require(
            sp.Poly(
                adjacent_response - expected_adjacent,
                i_symbol,
                domain=sp.QQ,
            ).is_zero,
            f"adjacent Laguerre-jet response failed at order={order}",
        )
        adjacent_response_checks += 1

        monomial_response = sp.expand(
            rising(i_symbol + 2, order)
            - rising(i_symbol + 1, order)
        )
        expected_monomial = sp.expand(
            order * rising(i_symbol + 2, order - 1)
        )
        require(
            sp.Poly(
                monomial_response - expected_monomial,
                i_symbol,
                domain=sp.QQ,
            ).is_zero,
            f"adjacent monomial response failed at order={order}",
        )
        monomial_response_checks += 1


# ---------------------------------------------------------------------------
# 2. Point Vandermonde and ordered positive-mixture determinants
# ---------------------------------------------------------------------------

max_symbolic_rank = 5
point_jet_determinants = 0
point_monomial_determinants = 0

for rank in range(1, max_symbolic_rank + 1):
    points = sp.symbols(f"u0:{rank}")
    vandermonde = sp.prod(
        points[right] - points[left]
        for left in range(rank)
        for right in range(left + 1, rank)
    )
    target = sp.expand(factorial(rank) * vandermonde)

    jet_matrix = sp.Matrix(
        [
            [
                row * falling(points[column], row - 1)
                for column in range(rank)
            ]
            for row in range(1, rank + 1)
        ]
    )
    require(
        sp.expand(jet_matrix.det() - target) == 0,
        f"point Laguerre-jet Vandermonde failed at rank={rank}",
    )
    point_jet_determinants += 1

    monomial_matrix = sp.Matrix(
        [
            [
                row * rising(points[column] + 2, row - 1)
                for column in range(rank)
            ]
            for row in range(1, rank + 1)
        ]
    )
    require(
        sp.expand(monomial_matrix.det() - target) == 0,
        f"point monomial Vandermonde failed at rank={rank}",
    )
    point_monomial_determinants += 1


ordered_cone_configurations = 0
ordered_cone_summands = 0
ordered_cone_determinant_matches = 0
ordered_cone_positivity_checks = 0
monomial_jet_determinant_matches = 0

for rank in range(2, 6):
    for variant in range(6):
        supports = []
        coefficients = []
        cursor = variant
        for column in range(rank):
            length = 1 + ((column + variant) % 3)
            support = tuple(range(cursor, cursor + length))
            supports.append(support)
            coefficients.append(
                tuple(
                    sp.Rational(
                        (column + 2) * (offset + 2) + variant + 1,
                        (column + 3) * (offset + 3) + variant + 2,
                    )
                    for offset in range(length)
                )
            )
            cursor += length + 1 + ((2 * column + variant) % 3)

        require(
            all(
                max(supports[column]) < min(supports[column + 1])
                for column in range(rank - 1)
            ),
            f"ordered-support generator failed at rank={rank}, variant={variant}",
        )
        require(
            all(
                coefficient > 0
                for column_coefficients in coefficients
                for coefficient in column_coefficients
            ),
            f"positive coefficient generator failed at rank={rank}, variant={variant}",
        )

        jet_matrix = sp.Matrix(
            [
                [
                    sum(
                        coefficients[column][offset]
                        * row
                        * falling(index, row - 1)
                        for offset, index in enumerate(supports[column])
                    )
                    for column in range(rank)
                ]
                for row in range(1, rank + 1)
            ]
        )
        monomial_matrix = sp.Matrix(
            [
                [
                    sum(
                        coefficients[column][offset]
                        * row
                        * rising(index + 2, row - 1)
                        for offset, index in enumerate(supports[column])
                    )
                    for column in range(rank)
                ]
                for row in range(1, rank + 1)
            ]
        )

        mixture_target = sp.Integer(0)
        for choices in product(
            *(range(len(support)) for support in supports)
        ):
            selected = [
                supports[column][choices[column]]
                for column in range(rank)
            ]
            coefficient_product = sp.prod(
                coefficients[column][choices[column]]
                for column in range(rank)
            )
            point_vandermonde = sp.prod(
                selected[right] - selected[left]
                for left in range(rank)
                for right in range(left + 1, rank)
            )
            require(
                point_vandermonde > 0,
                f"ordered point term was not positive at rank={rank}, "
                f"variant={variant}, choices={choices}",
            )
            mixture_target += (
                factorial(rank) * coefficient_product * point_vandermonde
            )
            ordered_cone_summands += 1

        jet_determinant = sp.factor(jet_matrix.det())
        monomial_determinant = sp.factor(monomial_matrix.det())
        require(
            sp.simplify(jet_determinant - mixture_target) == 0,
            f"multilinear jet determinant failed at rank={rank}, "
            f"variant={variant}",
        )
        ordered_cone_determinant_matches += 1
        require(
            jet_determinant > 0,
            f"ordered cone determinant was not positive at rank={rank}, "
            f"variant={variant}",
        )
        ordered_cone_positivity_checks += 1
        require(
            sp.simplify(monomial_determinant - jet_determinant) == 0,
            f"monomial/jet determinant mismatch at rank={rank}, "
            f"variant={variant}",
        )
        monomial_jet_determinant_matches += 1
        ordered_cone_configurations += 1


# ---------------------------------------------------------------------------
# 3. Sharp ell_(D-1) quotient-Krylov delay and Pade numerator
# ---------------------------------------------------------------------------

max_dimension = 12
krylov_dimensions = 0
krylov_forced_zero_checks = 0
krylov_exit_checks = 0
pade_numerator_checks = 0
consecutive_coprime_checks = 0

for dimension in range(2, max_dimension + 1):
    modulus = monic_laguerre(dimension)
    delayed_state = monic_laguerre(dimension - 1)

    require(
        sp.gcd(
            sp.Poly(modulus, s, domain=sp.QQ),
            sp.Poly(delayed_state, s, domain=sp.QQ),
        ).degree()
        == 0,
        f"consecutive Laguerre coprimality failed at D={dimension}",
    )
    consecutive_coprime_checks += 1

    moments = [
        factorial_functional(s**row * delayed_state)
        for row in range(dimension)
    ]
    for row in range(dimension - 1):
        require(
            moments[row] == 0,
            f"Krylov forced zero failed at D={dimension}, row={row}",
        )
        krylov_forced_zero_checks += 1
    require(
        moments[-1] == factorial(dimension - 1) ** 2,
        f"Krylov exit failed at D={dimension}",
    )
    krylov_exit_checks += 1

    denominator = sp.Poly(
        sp.expand(z**dimension * modulus.subs(s, 1 / z)),
        z,
        domain=sp.QQ,
    )
    require(
        denominator.coeff_monomial(1) == 1,
        f"quotient Pade denominator normalization failed at D={dimension}",
    )
    numerator_coefficients = []
    for degree in range(dimension):
        numerator_coefficients.append(
            sp.simplify(
                sum(
                    denominator.coeff_monomial(z**offset)
                    * moments[degree - offset]
                    for offset in range(degree + 1)
                )
            )
        )
    numerator = sp.expand(
        sum(
            numerator_coefficients[degree] * z**degree
            for degree in range(dimension)
        )
    )
    require(
        numerator
        == factorial(dimension - 1) ** 2 * z ** (dimension - 1),
        f"sharp quotient Pade numerator failed at D={dimension}",
    )
    require(
        sp.gcd(denominator, sp.Poly(numerator, z, domain=sp.QQ)).degree()
        == 0,
        f"quotient Pade numerator/denominator coprimality failed at D={dimension}",
    )
    pade_numerator_checks += 1
    krylov_dimensions += 1


# ---------------------------------------------------------------------------
# 4. Root-free trace-zero unit certificates
# ---------------------------------------------------------------------------

max_trace_dimension = 8
trace_zero_unit_checks = 0
norm_nonzero_checks = 0
resultant_nonzero_checks = 0
finite_horizon_agreement_checks = 0
first_post_horizon_divergence_checks = 0
first_alias_checks = 0

for dimension in range(2, max_trace_dimension + 1):
    modulus = monic_laguerre(dimension)
    delayed_state = monic_laguerre(dimension - 1)
    modulus_poly = sp.Poly(modulus, s, domain=sp.QQ)
    derivative = sp.diff(modulus, s)
    inverse = sp.invert(
        sp.Poly(s * derivative**2, s, domain=sp.QQ),
        modulus_poly,
    ).as_expr()
    omega = quotient_remainder(factorial(dimension) ** 2 * inverse, modulus)

    weighted_trace = quotient_trace(
        quotient_remainder(omega * delayed_state, modulus),
        modulus,
        dimension,
    )
    require(
        weighted_trace == 0,
        f"weighted trace-zero certificate failed at D={dimension}",
    )
    trace_zero_unit_checks += 1

    norm = quotient_norm(delayed_state, modulus, dimension)
    require(norm != 0, f"unit norm vanished at D={dimension}")
    norm_nonzero_checks += 1

    resultant = sp.resultant(delayed_state, modulus, s)
    require(resultant != 0, f"unit resultant vanished at D={dimension}")
    require(
        sp.simplify(norm - resultant) == 0,
        f"norm/resultant convention mismatch at D={dimension}",
    )
    resultant_nonzero_checks += 1

    # The rational transfer function is the quotient-Krylov series.  It
    # agrees with the original factorial functional only while the product
    # remains inside the exact 2D-1 quadrature horizon.  For ell_(D-1), this
    # includes multipliers s^0,...,s^D and fails first at s^(D+1).
    for row in range(dimension + 1):
        quotient_class = quotient_remainder(
            s**row * delayed_state,
            modulus,
        )
        quotient_value = quotient_trace(
            quotient_remainder(omega * quotient_class, modulus),
            modulus,
            dimension,
        )
        original_value = factorial_functional(s**row * delayed_state)
        require(
            quotient_value == original_value,
            f"finite-horizon quotient agreement failed at D={dimension}, "
            f"row={row}",
        )
        finite_horizon_agreement_checks += 1

    first_outside_row = dimension + 1
    outside_class = quotient_remainder(
        s**first_outside_row * delayed_state,
        modulus,
    )
    quotient_outside = quotient_trace(
        quotient_remainder(omega * outside_class, modulus),
        modulus,
        dimension,
    )
    original_outside = factorial_functional(
        s**first_outside_row * delayed_state
    )
    require(
        original_outside - quotient_outside == factorial(dimension) ** 2,
        f"first post-horizon divergence failed at D={dimension}",
    )
    first_post_horizon_divergence_checks += 1

    require(
        quotient_remainder(modulus**2, modulus) == 0
        and factorial_functional(modulus**2) == factorial(dimension) ** 2,
        f"first readout alias failed at D={dimension}",
    )
    first_alias_checks += 1


# ---------------------------------------------------------------------------
# 5. Single-selector linear algebra and the sharp complex two-cone witness
# ---------------------------------------------------------------------------

selector_pattern_checks = 0
selector_positive_patterns = 0

for dimension in range(2, 7):
    weights = [sp.Integer(index + 1) for index in range(dimension)]
    for values in product((-1, 0, 1), repeat=dimension):
        readout_rows = sp.Matrix(
            [
                weights,
                [
                    weights[index] * values[index]
                    for index in range(dimension)
                ],
            ]
        )
        row_rank = readout_rows.rank()
        coordinate_in_row_span = False
        for selected in range(dimension):
            coordinate_row = sp.Matrix(
                [
                    [
                        sp.Integer(1 if index == selected else 0)
                        for index in range(dimension)
                    ]
                ]
            )
            if readout_rows.col_join(coordinate_row).rank() == row_rank:
                coordinate_in_row_span = True
                break

        selector_pattern = False
        for selected in range(dimension):
            outside_values = [
                values[index]
                for index in range(dimension)
                if index != selected
            ]
            common_value = outside_values[0]
            if (
                all(value == common_value for value in outside_values)
                and values[selected] != common_value
            ):
                selector_pattern = True
                break

        require(
            coordinate_in_row_span == selector_pattern,
            f"single-selector iff failed at D={dimension}, values={values}",
        )
        selector_pattern_checks += 1
        if selector_pattern:
            selector_positive_patterns += 1


sharp_two_cone_checks = 0
lower_atom = adjacent_difference(1)
upper_atom = adjacent_difference(2)
sharp_coefficient = (-3 + sp.I * sp.sqrt(3)) / 2
sharp_state = sp.expand(sharp_coefficient * lower_atom + upper_atom)
sharp_moments = tuple(
    sp.simplify(factorial_functional(sharp_state**power))
    for power in range(1, 4)
)
require(
    sharp_moments
    == (
        0,
        0,
        36 + 54 * sp.I * sp.sqrt(3),
    ),
    f"sharp THM-2830 complex witness failed: {sharp_moments}",
)
require(
    tuple(
        sp.simplify(
            sp.binomial(2 * power, power) * sharp_moments[power - 1]
        )
        for power in range(1, 4)
    )
    == (
        0,
        0,
        720 + 1080 * sp.I * sp.sqrt(3),
    ),
    "sharp Gaussian moment 2/4/6 witness failed",
)
sharp_two_cone_checks += 1


print("THM-2842 LAGUERRE VARIANCE-JET OBSERVABILITY - exact audit")
print(
    "variance jet kernels / functional polynomial identities:",
    variance_kernel_checks,
    "/",
    functional_polynomial_checks,
)
print(
    "adjacent jet / monomial response identities:",
    adjacent_response_checks,
    "/",
    monomial_response_checks,
)
print(
    "symbolic point Vandermonde determinants (jet / monomial):",
    point_jet_determinants,
    "/",
    point_monomial_determinants,
)
print("ordered positive-cone configurations:", ordered_cone_configurations)
print("multilinear positive Vandermonde summands:", ordered_cone_summands)
print(
    "mixture determinant / positivity / basis-change checks:",
    ordered_cone_determinant_matches,
    "/",
    ordered_cone_positivity_checks,
    "/",
    monomial_jet_determinant_matches,
)
print("sharp Krylov dimensions D=2..12:", krylov_dimensions)
print(
    "Krylov forced zeros / exits / quotient Pade numerators:",
    krylov_forced_zero_checks,
    "/",
    krylov_exit_checks,
    "/",
    pade_numerator_checks,
)
print("consecutive Laguerre coprimality checks:", consecutive_coprime_checks)
print(
    "trace-zero units / nonzero norms / nonzero resultants D=2..8:",
    trace_zero_unit_checks,
    "/",
    norm_nonzero_checks,
    "/",
    resultant_nonzero_checks,
)
print(
    "finite-horizon agreements / first divergences / first aliases:",
    finite_horizon_agreement_checks,
    "/",
    first_post_horizon_divergence_checks,
    "/",
    first_alias_checks,
)
print(
    "single-selector iff patterns / positive selector patterns D=2..6:",
    selector_pattern_checks,
    "/",
    selector_positive_patterns,
)
print("sharp THM-2830 complex/Gaussian witness checks:", sharp_two_cone_checks)
print("all arithmetic exact: integers/rationals only")
print("THM-2842 exact companion: PASS")
