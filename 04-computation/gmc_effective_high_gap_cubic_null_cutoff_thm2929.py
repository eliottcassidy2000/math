#!/usr/bin/env python3
"""Exact companion for THM-2929's effective high-gap cutoff.

The proof has two parts.

* For 13 <= R <= 19, all four oriented faces, the contraction
  derivative, and the normalized quartic endpoint are certified in
  exact Bernstein form over Q(epsilon_R), where
  epsilon_R^3 = 1/T_3(R).
* For R >= 20, every primitive factorial-ratio error atom contracts by
  a factor strictly smaller than 1/2.  Exact rational envelopes at
  R=20 therefore control the whole tail.

R=12 is retained as a hostile boundary: Bernstein slot three on the
lower eta face is strictly below -1/200.  This says that the chosen box
certificate fails, not that the branch is absent.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from hashlib import sha256
from math import factorial
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_eventual_high_gap_cubic_null_branch_thm2914.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(SOURCE_BYTES).hexdigest()
    == "57e796b1d5115080939b9d19f35d95e1b5de584d632335fbcfd6105267270c8c",
    "THM-2914 dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2914_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2914")
thm2914 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2914)


s, z = sp.symbols("s z")
xi, eta, epsilon = sp.symbols("xi eta epsilon")
box_u, box_v, bernstein_t = sp.symbols("box_u box_v bernstein_t")
ratio_r, ratio_n = sp.symbols(
    "ratio_r ratio_n",
    integer=True,
    nonnegative=True,
)

XI_LOW = sp.Rational(217, 100)
XI_HIGH = sp.Rational(221, 100)
ETA_LOW = sp.Rational(245, 100)
ETA_HIGH = sp.Rational(249, 100)
XI_CENTER = sp.Rational(219, 100)
ETA_CENTER = sp.Rational(247, 100)
BOX_RADIUS = sp.Rational(1, 50)


def f(index: int) -> sp.Expr:
    return s**index / factorial(index)


def factorial_readout(polynomial: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), s)
    return sp.expand(
        sum(
            coefficient * factorial(monomial[0])
            for monomial, coefficient in expanded.terms()
        )
    )


def binary_coefficients(polynomial: sp.Expr, order: int) -> tuple[sp.Expr, ...]:
    expanded = sp.Poly(sp.expand(polynomial), z)
    return tuple(
        sp.cancel(expanded.nth(index) / sp.binomial(order, index))
        for index in range(order + 1)
    )


def reduce_epsilon(polynomial: sp.Expr, inverse_t3: sp.Rational) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), epsilon)
    return sp.expand(
        sum(
            coefficient
            * inverse_t3 ** (monomial[0] // 3)
            * epsilon ** (monomial[0] % 3)
            for monomial, coefficient in expanded.terms()
        )
    )


def epsilon_interval(
    expression: sp.Expr,
    lower: sp.Rational,
    upper: sp.Rational,
) -> tuple[sp.Rational, sp.Rational]:
    """Coefficientwise interval for a polynomial in positive epsilon."""
    polynomial = sp.Poly(sp.expand(expression), epsilon)
    lower_value = sp.Rational(0)
    upper_value = sp.Rational(0)
    for (power,), coefficient in polynomial.terms():
        require(coefficient.is_Rational, "nonrational epsilon coefficient")
        if coefficient >= 0:
            lower_value += coefficient * lower**power
            upper_value += coefficient * upper**power
        else:
            lower_value += coefficient * upper**power
            upper_value += coefficient * lower**power
    return sp.factor(lower_value), sp.factor(upper_value)


def bernstein_one(
    expression: sp.Expr,
    variable: sp.Symbol,
    lower: sp.Rational,
    upper: sp.Rational,
) -> tuple[sp.Expr, ...]:
    transformed = sp.Poly(
        sp.expand(
            expression.subs(
                variable,
                lower + (upper - lower) * bernstein_t,
            )
        ),
        bernstein_t,
    )
    degree = transformed.degree()
    powers = tuple(transformed.nth(index) for index in range(degree + 1))
    return tuple(
        sp.expand(
            sum(
                powers[index]
                * sp.binomial(slot, index)
                / sp.binomial(degree, index)
                for index in range(slot + 1)
            )
        )
        for slot in range(degree + 1)
    )


def bernstein_two(expression: sp.Expr) -> tuple[tuple[sp.Expr, ...], ...]:
    transformed = sp.Poly(
        sp.expand(
            expression.subs(
                {
                    xi: XI_LOW + (XI_HIGH - XI_LOW) * box_u,
                    eta: ETA_LOW + (ETA_HIGH - ETA_LOW) * box_v,
                }
            )
        ),
        box_u,
        box_v,
    )
    degree_u = transformed.degree(box_u)
    degree_v = transformed.degree(box_v)
    powers = {
        (index_u, index_v): transformed.coeff_monomial(
            box_u**index_u * box_v**index_v
        )
        for index_u in range(degree_u + 1)
        for index_v in range(degree_v + 1)
    }
    return tuple(
        tuple(
            sp.expand(
                sum(
                    powers[index_u, index_v]
                    * sp.binomial(slot_u, index_u)
                    / sp.binomial(degree_u, index_u)
                    * sp.binomial(slot_v, index_v)
                    / sp.binomial(degree_v, index_v)
                    for index_u in range(slot_u + 1)
                    for index_v in range(slot_v + 1)
                )
            )
            for slot_v in range(degree_v + 1)
        )
        for slot_u in range(degree_u + 1)
    )


def flatten(table: tuple[tuple[sp.Expr, ...], ...]) -> tuple[sp.Expr, ...]:
    return tuple(entry for row in table for entry in row)


def expression_record(expression: sp.Expr) -> str:
    polynomial = sp.Poly(sp.expand(expression), epsilon)
    return ",".join(
        f"{monomial[0]}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )


def record_digest(records: list[str]) -> str:
    return sha256(("\n".join(records) + "\n").encode()).hexdigest()


def t3_value(top: int) -> int:
    return factorial(3 * top) // factorial(top) ** 3


def t4_value(top: int) -> int:
    return factorial(4 * top) // factorial(top) ** 4


LOW_U = f(1) - f(0)
LOW_V = f(2) - f(1)
LOW_BINARY = LOW_U + z * LOW_V
LOW_G = binary_coefficients(factorial_readout(LOW_BINARY**2), 2)
LOW_T = binary_coefficients(factorial_readout(LOW_BINARY**3), 3)
require(LOW_G == (1, 1, 2), "low quadratic tensor changed")
require(LOW_T == (2, 4, 10, 30), "low cubic tensor changed")

LIMIT_F = sp.Matrix(
    (
        -eta**3 + 6 * eta * xi**2 - 4 * xi**3 - 14,
        -2 * (eta**3 - 3 * eta**2 * xi + 2 * xi**3 + 4),
    )
)
INVERSE_JACOBIAN = sp.Matrix(
    (
        (
            sp.Rational(18870800, 285129003),
            -sp.Rational(6982600, 285129003),
        ),
        (
            sp.Rational(13965200, 285129003),
            sp.Rational(1635200, 95043001),
        ),
    )
)
require(
    INVERSE_JACOBIAN
    * LIMIT_F.jacobian((xi, eta)).subs(
        {xi: XI_CENTER, eta: ETA_CENTER}
    )
    == sp.eye(2),
    "fixed inverse Jacobian changed",
)
LIMIT_G = sp.expand(INVERSE_JACOBIAN * LIMIT_F)
LIMIT_ENDPOINT = sp.expand(
    (eta**2 - 2 * xi**2)
    * (eta**2 - 4 * eta * xi + 2 * xi**2)
)


def finite_data(top: int) -> tuple[
    sp.Matrix,
    sp.Matrix,
    sp.Expr,
]:
    """Return G_R, D(id-G_R), and J_R/H_R over Q(epsilon_R)."""
    high = f(top) - f(2)
    polynomial_u = LOW_U + epsilon * xi * high
    polynomial_v = LOW_V + epsilon * eta * high
    binary = polynomial_u + z * polynomial_v
    quadratic = binary_coefficients(factorial_readout(binary**2), 2)
    cubic = binary_coefficients(factorial_readout(binary**3), 3)
    quartic = binary_coefficients(factorial_readout(binary**4), 4)
    g0, g1, g2 = quadratic
    t0, t1, t2, t3 = cubic
    A0, A1, _A2, A3, A4 = quartic
    invariant_one = sp.expand(
        3 * t1 * g0 * g2
        - t3 * g0**2
        - 2 * t0 * g1 * g2
    )
    invariant_two = sp.expand(
        3 * t2 * g0 * g2
        - 2 * t3 * g1 * g0
        - t0 * g2**2
    )
    t3_integer = t3_value(top)
    inverse_t3 = sp.Rational(1, t3_integer)
    finite_f = sp.Matrix(
        (
            reduce_epsilon(invariant_one, inverse_t3),
            reduce_epsilon(invariant_two, inverse_t3),
        )
    )
    finite_g = sp.expand(INVERSE_JACOBIAN * finite_f)
    contraction_derivative = sp.expand(
        sp.eye(2) - finite_g.jacobian((xi, eta))
    )
    endpoint = sp.expand(
        (2 * A1 * g0 - A0 * g1) * g2**2
        - (2 * A3 * g2 - A4 * g1) * g0**2
    )
    # epsilon^-4 = T_3^2 epsilon^2.
    normalized_endpoint = reduce_epsilon(
        endpoint
        * t3_integer**2
        * epsilon**2
        / t4_value(top),
        inverse_t3,
    )
    return finite_g, contraction_derivative, normalized_endpoint


# The brackets are adjacent multiples of 10^-18 and are verified below
# by cubing.  No decimal approximation is used as a truth check.
EPSILON_BRACKETS = {
    12: (6660271211996, 6660271211997),
    13: (2279040051368, 2279040051369),
    14: (778363173909, 778363173910),
    15: (265396784923, 265396784924),
    16: (90361479086, 90361479087),
    17: (30727000369, 30727000370),
    18: (10436823676, 10436823677),
    19: (3541441644, 3541441645),
}


def exact_epsilon_bracket(top: int) -> tuple[sp.Rational, sp.Rational]:
    denominator = 10**18
    lower_numerator, upper_numerator = EPSILON_BRACKETS[top]
    lower = sp.Rational(lower_numerator, denominator)
    upper = sp.Rational(upper_numerator, denominator)
    require(
        lower**3 * t3_value(top) < 1 < upper**3 * t3_value(top),
        f"epsilon bracket failed at R={top}",
    )
    return lower, upper


def finite_box_audit() -> tuple[str, int, int, int]:
    records: list[str] = []
    face_count = 0
    contraction_count = 0
    endpoint_count = 0
    for top in range(12, 20):
        lower_epsilon, upper_epsilon = exact_epsilon_bracket(top)
        finite_g, contraction_derivative, normalized_endpoint = finite_data(top)
        faces = (
            -finite_g[0].subs(xi, XI_LOW),
            finite_g[0].subs(xi, XI_HIGH),
            -finite_g[1].subs(eta, ETA_LOW),
            finite_g[1].subs(eta, ETA_HIGH),
        )
        face_variables = (eta, eta, xi, xi)
        face_intervals = (
            (ETA_LOW, ETA_HIGH),
            (ETA_LOW, ETA_HIGH),
            (XI_LOW, XI_HIGH),
            (XI_LOW, XI_HIGH),
        )
        face_coefficients = tuple(
            coefficient
            for face, variable, interval in zip(
                faces,
                face_variables,
                face_intervals,
            )
            for coefficient in bernstein_one(
                face,
                variable,
                interval[0],
                interval[1],
            )
        )
        require(
            len(face_coefficients) == 20,
            f"face coefficient count changed at R={top}",
        )
        records.extend(
            f"face:{top}:{index}:{expression_record(coefficient)}"
            for index, coefficient in enumerate(face_coefficients)
        )
        face_count += len(face_coefficients)

        if top == 12:
            hostile = bernstein_one(
                -finite_g[1].subs(eta, ETA_LOW),
                xi,
                XI_LOW,
                XI_HIGH,
            )[3]
            hostile_upper = epsilon_interval(
                hostile,
                lower_epsilon,
                upper_epsilon,
            )[1]
            require(
                hostile_upper < -sp.Rational(1, 200),
                "R=12 lower-eta hostile lost its negative sign",
            )
            continue

        require(
            all(
                epsilon_interval(
                    coefficient,
                    lower_epsilon,
                    upper_epsilon,
                )[0]
                > sp.Rational(1, 160)
                for coefficient in face_coefficients
            ),
            f"oriented face floor failed at R={top}",
        )

        row_bounds = []
        for row in range(2):
            row_bound = sp.Rational(0)
            for column in range(2):
                coefficients = flatten(
                    bernstein_two(contraction_derivative[row, column])
                )
                records.extend(
                    f"contraction:{top}:{row}:{column}:{index}:"
                    f"{expression_record(coefficient)}"
                    for index, coefficient in enumerate(coefficients)
                )
                contraction_count += len(coefficients)
                row_bound += max(
                    max(
                        abs(bound)
                        for bound in epsilon_interval(
                            coefficient,
                            lower_epsilon,
                            upper_epsilon,
                        )
                    )
                    for coefficient in coefficients
                )
            row_bounds.append(sp.factor(row_bound))
        require(
            max(row_bounds) < sp.Rational(1, 5),
            f"contraction bound failed at R={top}",
        )

        endpoint_coefficients = flatten(
            bernstein_two(normalized_endpoint)
        )
        records.extend(
            f"endpoint:{top}:{index}:{expression_record(coefficient)}"
            for index, coefficient in enumerate(endpoint_coefficients)
        )
        endpoint_count += len(endpoint_coefficients)
        require(
            all(
                epsilon_interval(
                    coefficient,
                    lower_epsilon,
                    upper_epsilon,
                )[0]
                > 19
                for coefficient in endpoint_coefficients
            ),
            f"endpoint floor failed at R={top}",
        )
    return (
        record_digest(records),
        face_count,
        contraction_count,
        endpoint_count,
    )


LOW_A = s - 1
LOW_B = s**2 / 2 - s


def primitive_types(
    order: int,
    *,
    normalized_quartic: bool,
) -> tuple[tuple[int, int, int], ...]:
    """Return the (k,j,a) atoms from the literal factorial expansion."""
    answer: set[tuple[int, int, int]] = set()
    for v_count in range(order + 1):
        for high_u in range(order - v_count + 1):
            for high_v in range(v_count + 1):
                high_count = high_u + high_v
                if high_count == 0:
                    continue
                low_polynomial = sp.Poly(
                    sp.expand(
                        LOW_A ** (order - v_count - high_u)
                        * LOW_B ** (v_count - high_v)
                    ),
                    s,
                )
                for high_copies in range(high_count + 1):
                    for (degree,), _coefficient in low_polynomial.terms():
                        if (
                            order == 3
                            and high_count == 3
                            and high_copies == 3
                            and degree == 0
                        ):
                            continue
                        if (
                            normalized_quartic
                            and high_count == 4
                            and high_copies == 4
                            and degree == 0
                        ):
                            continue
                        offset = 2 * (high_count - high_copies) + degree
                        answer.add((high_count, high_copies, offset))
    return tuple(sorted(answer))


QUADRATIC_TYPES = primitive_types(2, normalized_quartic=False)
CUBIC_TYPES = primitive_types(3, normalized_quartic=False)
QUARTIC_TYPES = primitive_types(4, normalized_quartic=True)
require(
    (len(QUADRATIC_TYPES), len(CUBIC_TYPES), len(QUARTIC_TYPES))
    == (9, 22, 45),
    "primitive tail type counts changed",
)


def successive_t3() -> sp.Expr:
    return sp.cancel(
        3
        * (3 * ratio_r + 1)
        * (3 * ratio_r + 2)
        / (ratio_r + 1) ** 2
    )


def successive_t4() -> sp.Expr:
    return sp.cancel(
        4
        * (4 * ratio_r + 1)
        * (4 * ratio_r + 2)
        * (4 * ratio_r + 3)
        / (ratio_r + 1) ** 3
    )


def successive_m(high_copies: int, offset: int) -> sp.Expr:
    return sp.cancel(
        sp.prod(
            high_copies * ratio_r + offset + shift
            for shift in range(1, high_copies + 1)
        )
        / (ratio_r + 1) ** high_copies
    )


def ratio_audit() -> str:
    records: list[str] = []
    t3_ratio = successive_t3()
    t4_ratio = successive_t4()
    for label, types, normalized_quartic in (
        ("quadratic", QUADRATIC_TYPES, False),
        ("cubic", CUBIC_TYPES, False),
        ("quartic", QUARTIC_TYPES, True),
    ):
        for high_count, high_copies, offset in types:
            moment_ratio = successive_m(high_copies, offset)
            if normalized_quartic:
                difference = sp.cancel(
                    t4_ratio**3
                    - 8
                    * t3_ratio ** (4 - high_count)
                    * moment_ratio**3
                )
            else:
                difference = sp.cancel(
                    t3_ratio**high_count - 8 * moment_ratio**3
                )
            numerator = sp.Poly(
                sp.expand(
                    difference.as_numer_denom()[0].subs(
                        ratio_r,
                        ratio_n + 20,
                    )
                ),
                ratio_n,
            )
            require(
                all(coefficient > 0 for coefficient in numerator.all_coeffs()),
                "tail ratio lost coefficientwise positivity: "
                f"{label},{high_count},{high_copies},{offset}",
            )
            records.append(
                f"{label}:{high_count}:{high_copies}:{offset}:"
                + ",".join(map(str, numerator.all_coeffs()))
            )
    return record_digest(records)


def coefficient_l1(expression: sp.Expr) -> sp.Rational:
    polynomial = sp.Poly(sp.expand(expression), box_u, box_v)
    return sp.factor(sum(abs(coefficient) for _, coefficient in polynomial.terms()))


def physical_derivative_l1(
    expression: sp.Expr,
    variable: sp.Symbol,
) -> sp.Rational:
    # xi=xi_0+box_u/50 and eta=eta_0+box_v/50.
    return sp.factor(50 * coefficient_l1(sp.diff(expression, variable)))


def bounded_error_atoms(
    order: int,
    index: int,
    *,
    normalized_quartic: bool,
    epsilon_lower: sp.Rational,
    epsilon_upper: sp.Rational,
) -> tuple[sp.Expr, ...]:
    """Absolute-envelope atoms at R=20, retaining their box monomials."""
    top = 20
    high_xi = XI_CENTER + BOX_RADIUS * box_u
    high_eta = ETA_CENTER + BOX_RADIUS * box_v
    answer: list[sp.Expr] = []
    for high_u in range(order - index + 1):
        for high_v in range(index + 1):
            high_count = high_u + high_v
            if high_count == 0:
                continue
            low_polynomial = sp.Poly(
                sp.expand(
                    LOW_A ** (order - index - high_u)
                    * LOW_B ** (index - high_v)
                ),
                s,
            )
            box_monomial = (
                sp.binomial(order - index, high_u)
                * sp.binomial(index, high_v)
                * high_xi**high_u
                * high_eta**high_v
            )
            for high_copies in range(high_count + 1):
                for (degree,), low_coefficient in low_polynomial.terms():
                    if (
                        order == 3
                        and high_count == 3
                        and high_copies == 3
                        and degree == 0
                    ):
                        continue
                    if (
                        normalized_quartic
                        and high_count == 4
                        and high_copies == 4
                        and degree == 0
                    ):
                        continue
                    offset = 2 * (high_count - high_copies) + degree
                    factorial_ratio = sp.Rational(
                        factorial(high_copies * top + offset),
                        factorial(top) ** high_copies
                        * 2 ** (high_count - high_copies),
                    )
                    signed_coefficient = (
                        (-1) ** (high_count - high_copies)
                        * sp.binomial(high_count, high_copies)
                        * low_coefficient
                    )
                    if normalized_quartic:
                        if high_count < 4:
                            epsilon_factor = epsilon_lower ** (high_count - 4)
                        else:
                            epsilon_factor = sp.Integer(1)
                        value_bound = (
                            factorial_ratio
                            * epsilon_factor
                            / t4_value(top)
                        )
                    else:
                        value_bound = (
                            factorial_ratio * epsilon_upper**high_count
                        )
                    answer.append(
                        sp.expand(
                            signed_coefficient
                            * box_monomial
                            * value_bound
                        )
                    )
    return tuple(answer)


def moment_envelopes(
    order: int,
    *,
    normalized_quartic: bool,
    epsilon_lower: sp.Rational,
    epsilon_upper: sp.Rational,
) -> tuple[
    tuple[sp.Rational, ...],
    tuple[sp.Rational, ...],
    tuple[sp.Rational, ...],
]:
    atom_lists = tuple(
        bounded_error_atoms(
            order,
            index,
            normalized_quartic=normalized_quartic,
            epsilon_lower=epsilon_lower,
            epsilon_upper=epsilon_upper,
        )
        for index in range(order + 1)
    )
    values = tuple(
        sp.factor(sum(coefficient_l1(atom) for atom in atoms))
        for atoms in atom_lists
    )
    xi_derivatives = tuple(
        sp.factor(
            sum(
                physical_derivative_l1(atom, box_u)
                for atom in atoms
            )
        )
        for atoms in atom_lists
    )
    eta_derivatives = tuple(
        sp.factor(
            sum(
                physical_derivative_l1(atom, box_v)
                for atom in atoms
            )
        )
        for atoms in atom_lists
    )
    return values, xi_derivatives, eta_derivatives


def product_difference(
    factors: tuple[
        tuple[sp.Rational, sp.Rational, sp.Rational, sp.Rational],
        ...,
    ],
) -> sp.Rational:
    return sp.factor(
        sp.prod(base + error for base, error, _, _ in factors)
        - sp.prod(base for base, _, _, _ in factors)
    )


def derivative_product_difference(
    factors: tuple[
        tuple[sp.Rational, sp.Rational, sp.Rational, sp.Rational],
        ...,
    ],
) -> sp.Rational:
    answer = sp.Rational(0)
    for index, (_base, _error, base_derivative, error_derivative) in enumerate(
        factors
    ):
        answer += (
            (base_derivative + error_derivative)
            * sp.prod(
                factors[other][0] + factors[other][1]
                for other in range(len(factors))
                if other != index
            )
            - base_derivative
            * sp.prod(
                factors[other][0]
                for other in range(len(factors))
                if other != index
            )
        )
    return sp.factor(answer)


def invariant_error_bounds(
    g_error: tuple[sp.Rational, ...],
    t_error: tuple[sp.Rational, ...],
    g_derivative_error: tuple[sp.Rational, ...],
    t_derivative_error: tuple[sp.Rational, ...],
    base_t_derivative: tuple[sp.Rational, ...],
) -> tuple[
    tuple[sp.Rational, sp.Rational],
    tuple[sp.Rational, sp.Rational],
]:
    limit_g_norm = (sp.Integer(1), sp.Integer(1), sp.Integer(2))
    limit_t_polynomials = (
        2 + (XI_CENTER + BOX_RADIUS * box_u) ** 3,
        4
        + (XI_CENTER + BOX_RADIUS * box_u) ** 2
        * (ETA_CENTER + BOX_RADIUS * box_v),
        10
        + (XI_CENTER + BOX_RADIUS * box_u)
        * (ETA_CENTER + BOX_RADIUS * box_v) ** 2,
        30 + (ETA_CENTER + BOX_RADIUS * box_v) ** 3,
    )
    limit_t_norm = tuple(
        coefficient_l1(polynomial)
        for polynomial in limit_t_polynomials
    )

    def factors(t_index: int, *g_indices: int) -> tuple[
        tuple[sp.Rational, sp.Rational, sp.Rational, sp.Rational],
        ...,
    ]:
        return (
            (
                limit_t_norm[t_index],
                t_error[t_index],
                base_t_derivative[t_index],
                t_derivative_error[t_index],
            ),
        ) + tuple(
            (
                limit_g_norm[g_index],
                g_error[g_index],
                sp.Rational(0),
                g_derivative_error[g_index],
            )
            for g_index in g_indices
        )

    first_terms = (
        (3, factors(1, 0, 2)),
        (1, factors(3, 0, 0)),
        (2, factors(0, 1, 2)),
    )
    second_terms = (
        (3, factors(2, 0, 2)),
        (2, factors(3, 1, 0)),
        (1, factors(0, 2, 2)),
    )
    value_bounds = tuple(
        sp.factor(
            sum(
                multiplier * product_difference(term_factors)
                for multiplier, term_factors in terms
            )
        )
        for terms in (first_terms, second_terms)
    )
    derivative_bounds = tuple(
        sp.factor(
            sum(
                multiplier
                * derivative_product_difference(term_factors)
                for multiplier, term_factors in terms
            )
        )
        for terms in (first_terms, second_terms)
    )
    return value_bounds, derivative_bounds


def tail_envelope_audit() -> tuple[str, tuple[sp.Rational, ...]]:
    top = 20
    lower_epsilon = sp.Rational(1200603, 10**15)
    upper_epsilon = sp.Rational(1200604, 10**15)
    require(
        lower_epsilon**3 * t3_value(top)
        < 1
        < upper_epsilon**3 * t3_value(top),
        "R=20 epsilon bracket changed",
    )
    g_error, g_xi_error, g_eta_error = moment_envelopes(
        2,
        normalized_quartic=False,
        epsilon_lower=lower_epsilon,
        epsilon_upper=upper_epsilon,
    )
    t_error, t_xi_error, t_eta_error = moment_envelopes(
        3,
        normalized_quartic=False,
        epsilon_lower=lower_epsilon,
        epsilon_upper=upper_epsilon,
    )
    quartic_error, _quartic_xi_error, _quartic_eta_error = moment_envelopes(
        4,
        normalized_quartic=True,
        epsilon_lower=lower_epsilon,
        epsilon_upper=upper_epsilon,
    )

    limit_t = (
        2 + (XI_CENTER + BOX_RADIUS * box_u) ** 3,
        4
        + (XI_CENTER + BOX_RADIUS * box_u) ** 2
        * (ETA_CENTER + BOX_RADIUS * box_v),
        10
        + (XI_CENTER + BOX_RADIUS * box_u)
        * (ETA_CENTER + BOX_RADIUS * box_v) ** 2,
        30 + (ETA_CENTER + BOX_RADIUS * box_v) ** 3,
    )
    base_t_xi = tuple(
        physical_derivative_l1(polynomial, box_u)
        for polynomial in limit_t
    )
    base_t_eta = tuple(
        physical_derivative_l1(polynomial, box_v)
        for polynomial in limit_t
    )
    invariant_value_error, invariant_xi_error = invariant_error_bounds(
        g_error,
        t_error,
        g_xi_error,
        t_xi_error,
        base_t_xi,
    )
    _unused_value_error, invariant_eta_error = invariant_error_bounds(
        g_error,
        t_error,
        g_eta_error,
        t_eta_error,
        base_t_eta,
    )
    require(
        invariant_value_error == _unused_value_error,
        "invariant value envelope depends on derivative direction",
    )

    transformed_value_error = tuple(
        sp.factor(
            abs(INVERSE_JACOBIAN[row, 0])
            * invariant_value_error[0]
            + abs(INVERSE_JACOBIAN[row, 1])
            * invariant_value_error[1]
        )
        for row in range(2)
    )
    transformed_derivative_row_error = tuple(
        sp.factor(
            abs(INVERSE_JACOBIAN[row, 0])
            * (invariant_xi_error[0] + invariant_eta_error[0])
            + abs(INVERSE_JACOBIAN[row, 1])
            * (invariant_xi_error[1] + invariant_eta_error[1])
        )
        for row in range(2)
    )
    require(
        max(transformed_value_error) < sp.Rational(1, 700),
        "tail map envelope exceeded 1/700",
    )
    require(
        max(transformed_derivative_row_error) < sp.Rational(1, 900),
        "tail derivative envelope exceeded 1/900",
    )

    limit_g_norm = (sp.Integer(1), sp.Integer(1), sp.Integer(2))
    limit_quartic = tuple(
        (XI_CENTER + BOX_RADIUS * box_u) ** (4 - index)
        * (ETA_CENTER + BOX_RADIUS * box_v) ** index
        for index in range(5)
    )
    limit_quartic_norm = tuple(
        coefficient_l1(polynomial)
        for polynomial in limit_quartic
    )

    def plain_difference(
        quartic_index: int,
        *g_indices: int,
    ) -> sp.Rational:
        factors = (
            (
                limit_quartic_norm[quartic_index],
                quartic_error[quartic_index],
                sp.Rational(0),
                sp.Rational(0),
            ),
        ) + tuple(
            (
                limit_g_norm[g_index],
                g_error[g_index],
                sp.Rational(0),
                sp.Rational(0),
            )
            for g_index in g_indices
        )
        return product_difference(factors)

    endpoint_error = sp.factor(
        2 * plain_difference(1, 0, 2, 2)
        + plain_difference(0, 1, 2, 2)
        + 2 * plain_difference(3, 2, 0, 0)
        + plain_difference(4, 1, 0, 0)
    )
    require(
        endpoint_error < sp.Rational(1, 400),
        "tail endpoint envelope exceeded 1/400",
    )
    envelope_values = (
        *transformed_value_error,
        *transformed_derivative_row_error,
        endpoint_error,
    )
    records = [
        f"g:{index}:{value}" for index, value in enumerate(g_error)
    ]
    records += [
        f"t:{index}:{value}" for index, value in enumerate(t_error)
    ]
    records += [
        f"quartic:{index}:{value}"
        for index, value in enumerate(quartic_error)
    ]
    records += [
        f"envelope:{index}:{value}"
        for index, value in enumerate(envelope_values)
    ]
    return record_digest(records), envelope_values


def limit_box_audit() -> tuple[
    sp.Rational,
    sp.Rational,
    sp.Rational,
]:
    faces = (
        -LIMIT_G[0].subs(xi, XI_LOW),
        LIMIT_G[0].subs(xi, XI_HIGH),
        -LIMIT_G[1].subs(eta, ETA_LOW),
        LIMIT_G[1].subs(eta, ETA_HIGH),
    )
    face_variables = (eta, eta, xi, xi)
    face_intervals = (
        (ETA_LOW, ETA_HIGH),
        (ETA_LOW, ETA_HIGH),
        (XI_LOW, XI_HIGH),
        (XI_LOW, XI_HIGH),
    )
    face_minimum = min(
        coefficient
        for face, variable, interval in zip(
            faces,
            face_variables,
            face_intervals,
        )
        for coefficient in bernstein_one(
            face,
            variable,
            interval[0],
            interval[1],
        )
    )
    require(
        face_minimum
        == sp.Rational(13337106413, 712822507500)
        and face_minimum > sp.Rational(1, 54),
        "limit oriented-face floor changed",
    )

    limit_contraction = sp.eye(2) - LIMIT_G.jacobian((xi, eta))
    row_bounds = []
    for row in range(2):
        row_bound = sp.Rational(0)
        for column in range(2):
            row_bound += max(
                abs(coefficient)
                for coefficient in flatten(
                    bernstein_two(limit_contraction[row, column])
                )
            )
        row_bounds.append(sp.factor(row_bound))
    contraction_maximum = max(row_bounds)
    require(
        contraction_maximum
        == sp.Rational(305733792, 2376075025)
        and contraction_maximum < sp.Rational(13, 100),
        "limit contraction floor changed",
    )

    endpoint_minimum = min(
        flatten(bernstein_two(LIMIT_ENDPOINT))
    )
    require(
        endpoint_minimum
        == sp.Rational(1929107681, 100000000)
        and endpoint_minimum > 19,
        "limit endpoint floor changed",
    )
    return face_minimum, contraction_maximum, endpoint_minimum


def selector_audit() -> None:
    selector_u = sp.symbols("selector_u")
    polynomial = (
        2 * selector_u**3
        + 18 * selector_u**2
        - 729 * selector_u
        + 54
    )
    lower = sp.Rational(15069319, 10**6)
    upper = sp.Rational(15069320, 10**6)
    poly = sp.Poly(polynomial, selector_u, domain=sp.QQ)
    require(
        poly.count_roots(lower, upper) == 1
        and poly.count_roots(0, 1) == 1
        and poly.count_roots(-25, -24) == 1
        and poly.count_roots(-sp.oo, sp.oo) == 3
        and polynomial.subs(selector_u, lower) < 0
        < polynomial.subs(selector_u, upper),
        "positive selector root bracket changed",
    )
    eta_lower = sp.Rational(247, 100)
    eta_upper = sp.Rational(248, 100)
    require(
        eta_lower**3 < lower < upper < eta_upper**3,
        "positive selector eta bracket changed",
    )
    selector_ratio = (
        2 * selector_u**2 + 21 * selector_u - 603
    ) / 189
    selector_derivative = sp.diff(selector_ratio, selector_u)
    require(
        selector_derivative == (4 * selector_u + 21) / 189
        and selector_derivative.subs(selector_u, 0) > 0
        and sp.diff(selector_derivative, selector_u) > 0,
        "selector ratio lost monotonicity",
    )
    xi_lower = eta_lower * selector_ratio.subs(selector_u, lower)
    xi_upper = eta_upper * selector_ratio.subs(selector_u, upper)
    require(
        XI_LOW < xi_lower < xi_upper < XI_HIGH,
        "positive selector left the fixed box",
    )
    # The smaller positive eliminant root has negative xi by THM-2914's
    # exact selector, so this is the unique positive-quadrant sheet.
    require(
        selector_ratio.subs(selector_u, 1) < 0,
        "small positive selector became positive-quadrant",
    )


def main() -> None:
    selector_audit()
    finite_digest, face_count, contraction_count, endpoint_count = (
        finite_box_audit()
    )
    ratio_digest = ratio_audit()
    envelope_digest, _envelope_values = tail_envelope_audit()
    face_limit, contraction_limit, endpoint_limit = limit_box_audit()

    tail_face_floor = sp.Rational(1, 54) - sp.Rational(1, 700)
    tail_contraction_ceiling = (
        sp.Rational(13, 100) + sp.Rational(1, 900)
    )
    tail_endpoint_floor = sp.Integer(19) - sp.Rational(1, 400)
    require(
        tail_face_floor == sp.Rational(323, 18900) > 0,
        "combined tail face floor changed",
    )
    require(
        tail_contraction_ceiling == sp.Rational(59, 450) < 1,
        "combined tail contraction bound changed",
    )
    require(
        tail_endpoint_floor == sp.Rational(7599, 400) > 0,
        "combined tail endpoint floor changed",
    )

    print("THM-2929 EFFECTIVE HIGH-GAP CUBIC-NULL CUTOFF")
    print(
        "dependency=THM-2914:"
        "57e796b1d5115080939b9d19f35d95e1b5de584d632335fbcfd6105267270c8c"
    )
    print(
        "box=[217/100,221/100]x[245/100,249/100];"
        "map=G=(D F_inf(219/100,247/100))^-1*F_R"
    )
    print(
        "positive_selector=u_in_(15069319,15069320)/10^6;"
        "unique_positive_quadrant;inside_box"
    )
    print(
        "finite_R=13..19;"
        "oriented_faces>1/160;contraction<1/5;J/H>19"
    )
    print(
        f"finite_counts=faces:{face_count};"
        f"contraction:{contraction_count};endpoint:{endpoint_count}"
    )
    print(f"finite_certificate_digest_sha256={finite_digest}")
    print(
        "R12_boundary=lower_eta_Bernstein_slot3<-1/200;"
        "certificate_failure_only"
    )
    print(
        "tail_types=quadratic:9,cubic:22,quartic:45;"
        "every_successive_ratio<1/2_for_R>=20"
    )
    print(f"tail_ratio_digest_sha256={ratio_digest}")
    print(
        "R20_envelopes=map<1/700;derivative_row<1/900;"
        "endpoint<1/400"
    )
    print(f"tail_envelope_digest_sha256={envelope_digest}")
    print(
        "limit_floors=face:"
        f"{face_limit};contraction:{contraction_limit};"
        f"endpoint:{endpoint_limit}"
    )
    print(
        "combined_tail=face>323/18900;"
        "contraction<59/450;J/H>7599/400"
    )
    print(
        "conclusion=unique_box_zero_and_positive_endpoint_for_every_R>=13;"
        "no_global_branch_uniqueness_or_optimal_cutoff_claim"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
