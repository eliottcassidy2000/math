#!/usr/bin/env python3
"""Exact companion for THM-2921.

The script proves first-window SFC(4) for the three nonconsecutive
four-subsets of a translated five-point exponent window.  It constructs the
ordinary-monomial quadratic, cubic, and quartic moment forms with their
multinomial coefficients, clears their exact depth denominators, and follows
one fixed degree-seven Macaulay maximal minor through all integer depths.

The determinant has degree at most 196.  Its 197 Gregory--Newton
coefficients are strictly positive (at base zero or one, as appropriate).
An independent direct four-variable multinomial constructor checks the
selected determinants modulo a prime at seven hostile depths.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import factorial, isqrt, prod

try:
    from flint import fmpz_mat, fmpz_poly
except ModuleNotFoundError as error:
    raise RuntimeError(
        "THM-2921 exact replay requires python-flint"
    ) from error


OFFSET_FAMILIES = (
    (0, 1, 2, 4),
    (0, 1, 3, 4),
    (0, 2, 3, 4),
)
ORDERS = (2, 3, 4)
TARGET_DEGREE = 7
PRIME = 1_000_003
AUDIT_DEPTHS = (0, 1, 2, 7, 31, 97, 197)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    return all(value % divisor for divisor in range(3, isqrt(value) + 1, 2))


def compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((total,),)
    return tuple(
        (first,) + tail
        for first in range(total + 1)
        for tail in compositions(total - first, parts - 1)
    )


MONOMIALS = {
    degree: compositions(degree, 3)
    for degree in range(TARGET_DEGREE + 1)
}
TARGET_MONOMIALS = MONOMIALS[TARGET_DEGREE]
TARGET_INDEX = {
    monomial: index
    for index, monomial in enumerate(TARGET_MONOMIALS)
}

# THM-2849 row order: 21 quadratic, 15 cubic, then 10 quartic rows.
# This fixed 36-row chart uses 20Q+10C+6F.
SELECTED_ROWS = (
    tuple(range(20))
    + tuple(range(21, 30))
    + (35,)
    + tuple(range(36, 42))
)


def multinomial(exponents: tuple[int, ...]) -> int:
    answer = factorial(sum(exponents))
    for exponent in exponents:
        answer //= factorial(exponent)
    return answer


def rising(start: int, length: int) -> int:
    return prod(range(start, start + length))


def normalized_tensor(
    depth: int,
    offsets: tuple[int, ...],
) -> Fraction:
    """L(prod f_(n+b))/L(f_n^m), with f_a=s^a/a!."""
    numerator = rising(len(offsets) * depth + 1, sum(offsets))
    denominator = prod(rising(depth + 1, offset) for offset in offsets)
    return Fraction(numerator, denominator)


def difference_tensor(
    depth: int,
    directions: tuple[int, ...],
    offsets: tuple[int, int, int, int],
) -> Fraction:
    """L(prod_j(f_(n+b_j)-f_(n+4))) in normalized tensor units."""
    answer = Fraction(0)
    for mask in range(1 << len(directions)):
        selected = tuple(
            offsets[3] if mask & (1 << position) else offsets[direction]
            for position, direction in enumerate(directions)
        )
        sign = -1 if mask.bit_count() % 2 else 1
        answer += sign * normalized_tensor(depth, selected)
    return answer


def moment_form(
    depth: int,
    order: int,
    offsets: tuple[int, int, int, int],
) -> dict[tuple[int, int, int], Fraction]:
    """Ordinary monomial coefficients, including multinomial copies."""
    answer: dict[tuple[int, int, int], Fraction] = {}
    for monomial in MONOMIALS[order]:
        directions = tuple(
            direction
            for direction, count in enumerate(monomial)
            for _ in range(count)
        )
        answer[monomial] = (
            multinomial(monomial)
            * difference_tensor(depth, directions, offsets)
        )
    return answer


def denominator(depth: int, order: int) -> int:
    if order == 2:
        return (depth + 1) * (depth + 2) * (depth + 3)
    if order == 3:
        return (
            (depth + 1) ** 2
            * (depth + 2) ** 2
            * (depth + 3) ** 2
            * (depth + 4)
        )
    if order == 4:
        return (
            (depth + 1) ** 3
            * (depth + 2) ** 3
            * (depth + 3) ** 3
            * (depth + 4) ** 2
        )
    raise ValueError(order)


POLY_X = fmpz_poly([0, 1])


def rising_poly(slope: int, intercept: int, length: int) -> fmpz_poly:
    answer = fmpz_poly(1)
    for offset in range(length):
        answer *= slope * POLY_X + intercept + offset
    return answer


def reduced_fraction(
    numerator: fmpz_poly,
    denominator_value: fmpz_poly,
) -> tuple[fmpz_poly, fmpz_poly]:
    common = numerator.gcd(denominator_value)
    numerator //= common
    denominator_value //= common
    if denominator_value.leading_coefficient() < 0:
        numerator = -numerator
        denominator_value = -denominator_value
    return numerator, denominator_value


def add_fraction(
    left: tuple[fmpz_poly, fmpz_poly],
    right: tuple[fmpz_poly, fmpz_poly],
) -> tuple[fmpz_poly, fmpz_poly]:
    left_num, left_den = left
    right_num, right_den = right
    return reduced_fraction(
        left_num * right_den + right_num * left_den,
        left_den * right_den,
    )


def normalized_tensor_poly(
    selected_offsets: tuple[int, ...],
) -> tuple[fmpz_poly, fmpz_poly]:
    numerator = rising_poly(
        len(selected_offsets),
        1,
        sum(selected_offsets),
    )
    denominator_value = fmpz_poly(1)
    for offset in selected_offsets:
        denominator_value *= rising_poly(1, 1, offset)
    return reduced_fraction(numerator, denominator_value)


def denominator_poly(order: int) -> fmpz_poly:
    if order == 2:
        return (POLY_X + 1) * (POLY_X + 2) * (POLY_X + 3)
    if order == 3:
        return (
            (POLY_X + 1) ** 2
            * (POLY_X + 2) ** 2
            * (POLY_X + 3) ** 2
            * (POLY_X + 4)
        )
    if order == 4:
        return (
            (POLY_X + 1) ** 3
            * (POLY_X + 2) ** 3
            * (POLY_X + 3) ** 3
            * (POLY_X + 4) ** 2
        )
    raise ValueError(order)


def symbolic_scaled_form(
    order: int,
    offsets: tuple[int, int, int, int],
) -> dict[tuple[int, int, int], fmpz_poly]:
    """Construct D_order(n) times the ordinary coefficient form."""
    answer: dict[tuple[int, int, int], fmpz_poly] = {}
    common_denominator = denominator_poly(order)
    for monomial in MONOMIALS[order]:
        directions = tuple(
            direction
            for direction, count in enumerate(monomial)
            for _ in range(count)
        )
        total = (fmpz_poly(0), fmpz_poly(1))
        for mask in range(1 << order):
            selected_offsets = tuple(
                offsets[3] if mask & (1 << position) else offsets[direction]
                for position, direction in enumerate(directions)
            )
            numerator, tensor_denominator = normalized_tensor_poly(
                selected_offsets
            )
            sign = -1 if mask.bit_count() % 2 else 1
            total = add_fraction(
                total,
                (sign * numerator, tensor_denominator),
            )
        numerator, coefficient_denominator = reduced_fraction(
            multinomial(monomial) * total[0],
            total[1],
        )
        require(
            common_denominator % coefficient_denominator == 0,
            f"denominator escape: {offsets}, order={order}, {monomial}",
        )
        answer[monomial] = (
            common_denominator // coefficient_denominator
        ) * numerator
    return answer


def macaulay_rows(
    depth: int,
    offsets: tuple[int, int, int, int],
    *,
    scaled: bool,
) -> list[list[int]]:
    rows: list[list[int]] = []
    for order in ORDERS:
        form = moment_form(depth, order, offsets)
        scale = denominator(depth, order) if scaled else 1
        coefficients = {
            monomial: coefficient * scale
            for monomial, coefficient in form.items()
        }
        require(
            all(
                coefficient.denominator == 1
                for coefficient in coefficients.values()
            ),
            f"nonintegral coefficient: {offsets}, n={depth}, m={order}",
        )
        for multiplier in MONOMIALS[TARGET_DEGREE - order]:
            row = [0] * len(TARGET_MONOMIALS)
            for monomial, coefficient in coefficients.items():
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[TARGET_INDEX[target]] = coefficient.numerator
            rows.append(row)
    return rows


def selected_determinant(
    depth: int,
    offsets: tuple[int, int, int, int],
    *,
    scaled: bool = True,
) -> int:
    rows = macaulay_rows(depth, offsets, scaled=scaled)
    return int(fmpz_mat([rows[index] for index in SELECTED_ROWS]).det())


def direct_four_variable_form_mod(
    depth: int,
    order: int,
    offsets: tuple[int, int, int, int],
) -> dict[tuple[int, int, int], int]:
    """Independent constructor: expand c3=-(x+y+z) after (6) of THM-2849."""
    support = tuple(depth + offset for offset in offsets)
    answer = {
        monomial: 0
        for monomial in MONOMIALS[order]
    }
    factorial_support = tuple(factorial(value) for value in support)
    base_normalization_inverse = (
        pow(factorial(depth) % PRIME, order, PRIME)
        * pow(factorial(order * depth) % PRIME, PRIME - 2, PRIME)
    ) % PRIME
    for counts in compositions(order, 4):
        numerator = factorial(
            sum(count * exponent for count, exponent in zip(counts, support))
        )
        denominator_value = prod(
            factorial_support[index] ** counts[index]
            for index in range(4)
        )
        base = (
            multinomial(counts)
            * (numerator % PRIME)
            * pow(denominator_value % PRIME, PRIME - 2, PRIME)
            * base_normalization_inverse
        ) % PRIME
        for tail in compositions(counts[3], 3):
            target = tuple(
                counts[index] + tail[index]
                for index in range(3)
            )
            value = base * multinomial(tail)
            if counts[3] % 2:
                value = -value
            answer[target] = (answer[target] + value) % PRIME
    return answer


def macaulay_rows_mod(
    forms: tuple[dict[tuple[int, int, int], int], ...],
) -> list[list[int]]:
    rows: list[list[int]] = []
    for order, form in zip(ORDERS, forms):
        for multiplier in MONOMIALS[TARGET_DEGREE - order]:
            row = [0] * len(TARGET_MONOMIALS)
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[TARGET_INDEX[target]] = coefficient % PRIME
            rows.append(row)
    return rows


def determinant_mod(matrix: list[list[int]]) -> int:
    rows = [[value % PRIME for value in row] for row in matrix]
    determinant = 1
    size = len(rows)
    require(all(len(row) == size for row in rows), "nonsquare determinant")
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if rows[row][column] % PRIME
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            rows[pivot], rows[column] = rows[column], rows[pivot]
            determinant = -determinant
        pivot_value = rows[column][column] % PRIME
        determinant = determinant * pivot_value % PRIME
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        for row in range(column + 1, size):
            multiplier_value = rows[row][column] * inverse % PRIME
            if not multiplier_value:
                continue
            for index in range(column, size):
                rows[row][index] = (
                    rows[row][index]
                    - multiplier_value * rows[column][index]
                ) % PRIME
    return determinant % PRIME


def independent_modular_audit(
    offsets: tuple[int, int, int, int],
) -> int:
    checks = 0
    for depth in AUDIT_DEPTHS:
        direct_forms = tuple(
            direct_four_variable_form_mod(depth, order, offsets)
            for order in ORDERS
        )
        for order, direct_form in zip(ORDERS, direct_forms):
            normalized_form = moment_form(depth, order, offsets)
            reduced = {
                monomial: (
                    coefficient.numerator
                    * pow(coefficient.denominator, PRIME - 2, PRIME)
                ) % PRIME
                for monomial, coefficient in normalized_form.items()
            }
            require(
                reduced == direct_form,
                f"independent form mismatch: {offsets}, n={depth}, m={order}",
            )
        direct_rows = macaulay_rows_mod(direct_forms)
        direct_minor = determinant_mod(
            [direct_rows[index] for index in SELECTED_ROWS]
        )
        scale = (
            pow(denominator(depth, 2), 20, PRIME)
            * pow(denominator(depth, 3), 10, PRIME)
            * pow(denominator(depth, 4), 6, PRIME)
        ) % PRIME
        require(
            selected_determinant(depth, offsets) % PRIME
            == direct_minor * scale % PRIME,
            f"independent minor mismatch: {offsets}, n={depth}",
        )
        checks += 1
    return checks


def finite_differences(values: list[int]) -> tuple[int, ...]:
    first_values = []
    row = values
    while row:
        first_values.append(row[0])
        row = [
            row[index + 1] - row[index]
            for index in range(len(row) - 1)
        ]
    return tuple(first_values)


def integer_sequence_digest(values: tuple[int, ...]) -> str:
    return sha256(",".join(map(str, values)).encode()).hexdigest()


def form_digest(
    forms: tuple[dict[tuple[int, int, int], fmpz_poly], ...],
) -> str:
    payload = []
    for order, form in zip(ORDERS, forms):
        for monomial in MONOMIALS[order]:
            payload.append(
                f"{order}:{monomial}:"
                + ",".join(map(str, form[monomial].coeffs()))
            )
    return sha256("|".join(payload).encode()).hexdigest()


def main() -> None:
    require(is_prime(PRIME), "audit modulus is not prime")
    require(
        PRIME > 4 * (max(AUDIT_DEPTHS) + 4),
        "audit factorial denominator met modulus",
    )
    require(len(SELECTED_ROWS) == 36, "selected row count changed")

    degree_bound = 20 * 3 + 10 * 7 + 6 * 11
    require(degree_bound == 196, "determinant degree invoice changed")

    # Freeze the exact bug that invalidated the first scratch certificate:
    # a mixed ordinary coefficient has two ordered tensor copies.
    mixed_tensor = difference_tensor(0, (0, 1), OFFSET_FAMILIES[0])
    mixed_coefficient = moment_form(0, 2, OFFSET_FAMILIES[0])[(1, 1, 0)]
    require(
        mixed_tensor != 0 and mixed_coefficient == 2 * mixed_tensor,
        "multinomial hostile no longer distinguishes the constructors",
    )

    print("THM-2921 DIAMETER-FOUR NONCONSECUTIVE MACAULAY--NEWTON CLOSURE")
    print("constructor=ordinary monomials with explicit multinomial copies")
    print(
        "selected_rows="
        + ",".join(map(str, SELECTED_ROWS))
        + ";allocation=(20Q,10C,6F)"
    )
    print("scaled_row_degrees=(3,7,11);determinant_degree_bound=196")
    print("multinomial_hostile=mixed_quadratic_ratio_2:PASS")

    total_modular_checks = 0
    for offsets in OFFSET_FAMILIES:
        symbolic_forms = tuple(
            symbolic_scaled_form(order, offsets)
            for order in ORDERS
        )
        symbolic_degrees = tuple(
            max(polynomial.degree() for polynomial in form.values())
            for form in symbolic_forms
        )
        require(
            symbolic_degrees == (3, 7, 11),
            f"scaled row degree changed for {offsets}",
        )

        # Exact polynomial-versus-rational checks at the singularly most
        # tempting small depths.
        for depth in range(9):
            for order, symbolic_form in zip(ORDERS, symbolic_forms):
                numeric_form = moment_form(depth, order, offsets)
                scale = denominator(depth, order)
                for monomial in MONOMIALS[order]:
                    expected = numeric_form[monomial] * scale
                    require(
                        expected.denominator == 1,
                        "cross-check lost integrality",
                    )
                    require(
                        int(symbolic_form[monomial](depth))
                        == expected.numerator,
                        f"symbolic/numeric mismatch: {offsets}, n={depth}",
                    )

        total_modular_checks += independent_modular_audit(offsets)

        values = tuple(
            selected_determinant(depth, offsets)
            for depth in range(198)
        )
        raw_zero = selected_determinant(0, offsets, scaled=False)
        scale_zero = (
            denominator(0, 2) ** 20
            * denominator(0, 3) ** 10
            * denominator(0, 4) ** 6
        )
        require(
            values[0] == raw_zero * scale_zero,
            f"n=0 original/scaled audit failed for {offsets}",
        )

        newton_base = 0 if offsets == (0, 2, 3, 4) else 1
        newton = finite_differences(
            list(values[newton_base:newton_base + 197])
        )
        require(len(newton) == 197, "Newton row length changed")
        require(
            all(value > 0 for value in newton),
            f"Newton positivity failed for {offsets}",
        )
        require(values[0] != 0, "depth-zero determinant vanished")

        print(
            f"family={offsets};"
            f"newton_base={newton_base};"
            f"n0_sign={'+' if values[0] > 0 else '-'};"
            f"n0_raw_digits={len(str(abs(raw_zero)))};"
            f"form_digest={form_digest(symbolic_forms)};"
            f"newton_digest={integer_sequence_digest(newton)}"
        )
        print(
            f"family={offsets};"
            f"newton_count={len(newton)};"
            f"first_digits={len(str(newton[0]))};"
            f"last_digits={len(str(newton[-1]))};"
            "all_positive=YES"
        )

    require(total_modular_checks == 21, "modular audit count changed")
    print("symbolic_denominator_division=PASS;families=3")
    print(
        f"independent_direct_constructor_mod_p={PRIME};"
        f"minor_checks={total_modular_checks};PASS"
    )
    print("consequence=fixed maximal minor nonzero for every integer n>=0")
    print(
        "scope=first-window SFC(4) for all three translated "
        "nonconsecutive four-subsets of five consecutive exponents"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
