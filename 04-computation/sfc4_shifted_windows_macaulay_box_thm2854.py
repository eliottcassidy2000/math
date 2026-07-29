#!/usr/bin/env python3
"""Exact Macaulay-rank census for shifted four-slot factorial windows.

Universe:
    0 <= a0 < a1 < a2 < a3 <= 8,
    k in {1, 2},
    moments L(H^(k+1)), ..., L(H^(k+4)).

For each cell the script forms the sharp Macaulay map in degree
    D = sum(moment degrees) - 3
over F_1000003.  Full target rank supplies a nonzero modular maximal minor,
which lifts through the factorial denominators to characteristic zero and
proves projective emptiness.  Rank deficiency would be reported as
inconclusive; it is never promoted as evidence of a complex point.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from math import comb, factorial, gcd, isqrt
import sys


PRIME = 1_000_003
TOP_EXPONENT = 8
WINDOWS = (1, 2)
MAX_ORDER = max(WINDOWS) + 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    return all(value % divisor for divisor in range(3, isqrt(value) + 1, 2))


@lru_cache(maxsize=None)
def compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((total,),)
    return tuple(
        (first,) + tail
        for first in range(total + 1)
        for tail in compositions(total - first, parts - 1)
    )


def normalized_moment_mod(
    support: tuple[int, int, int, int],
    exponents: tuple[int, int, int, int],
) -> int:
    numerator = factorial(
        sum(exponent * degree for exponent, degree in zip(exponents, support))
    )
    denominator = 1
    for exponent, degree in zip(exponents, support):
        denominator *= factorial(degree) ** exponent
    require(gcd(denominator, PRIME) == 1, "factorial denominator met modulus")
    return (
        (numerator % PRIME)
        * pow(denominator % PRIME, PRIME - 2, PRIME)
        % PRIME
    )


def moment_form_by_counts(
    support: tuple[int, int, int, int],
    order: int,
) -> dict[tuple[int, int, int, int], int]:
    answer: dict[tuple[int, int, int, int], int] = {}
    for exponents in compositions(order, 4):
        multinomial = factorial(order)
        for exponent in exponents:
            multinomial //= factorial(exponent)
        coefficient = (
            multinomial * normalized_moment_mod(support, exponents) % PRIME
        )
        if coefficient:
            answer[exponents] = coefficient
    return answer


def moment_form_by_ordered_words(
    support: tuple[int, int, int, int],
    order: int,
) -> dict[tuple[int, int, int, int], int]:
    """Independent expansion over the 4^order ordered tensor words."""
    answer: dict[tuple[int, int, int, int], int] = {}
    cache: dict[tuple[int, int, int, int], int] = {}
    for word in product(range(4), repeat=order):
        exponents = tuple(word.count(index) for index in range(4))
        if exponents not in cache:
            cache[exponents] = normalized_moment_mod(support, exponents)
        answer[exponents] = (
            answer.get(exponents, 0) + cache[exponents]
        ) % PRIME
    return {monomial: value for monomial, value in answer.items() if value}


def macaulay_rows(
    forms: tuple[dict[tuple[int, int, int, int], int], ...],
    degrees: tuple[int, ...],
    target_degree: int,
) -> tuple[list[dict[int, int]], int]:
    columns = compositions(target_degree, 4)
    column_index = {monomial: index for index, monomial in enumerate(columns)}
    rows: list[dict[int, int]] = []
    for form, degree in zip(forms, degrees):
        for multiplier in compositions(target_degree - degree, 4):
            row: dict[int, int] = {}
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index] for index in range(4)
                )
                row[column_index[target]] = coefficient
            rows.append(row)
    return rows, len(columns)


def sparse_rank_mod(
    rows: list[dict[int, int]],
    column_count: int,
    reverse_pivots: bool = True,
) -> int:
    """Exact F_p row rank with sparse normalized pivot rows."""
    pivots: dict[int, dict[int, int]] = {}
    for original in rows:
        row = {
            column: value % PRIME
            for column, value in original.items()
            if value % PRIME
        }
        while row:
            column = max(row) if reverse_pivots else min(row)
            value = row[column]
            if column not in pivots:
                inverse = pow(value, PRIME - 2, PRIME)
                row = {
                    index: coefficient * inverse % PRIME
                    for index, coefficient in row.items()
                    if coefficient % PRIME
                }
                pivots[column] = row
                break
            pivot = pivots[column]
            for index, coefficient in pivot.items():
                reduced = (row.get(index, 0) - value * coefficient) % PRIME
                if reduced:
                    row[index] = reduced
                elif index in row:
                    del row[index]
        if len(pivots) == column_count:
            return column_count
    return len(pivots)


def normalized_moment_q(
    support: tuple[int, int, int, int],
    exponents: tuple[int, int, int, int],
) -> Fraction:
    numerator = factorial(
        sum(exponent * degree for exponent, degree in zip(exponents, support))
    )
    denominator = 1
    for exponent, degree in zip(exponents, support):
        denominator *= factorial(degree) ** exponent
    return Fraction(numerator, denominator)


def moment_form_q(
    support: tuple[int, int, int, int],
    order: int,
) -> dict[tuple[int, int, int, int], Fraction]:
    answer = {}
    for exponents in compositions(order, 4):
        multinomial = factorial(order)
        for exponent in exponents:
            multinomial //= factorial(exponent)
        answer[exponents] = (
            multinomial * normalized_moment_q(support, exponents)
        )
    return answer


def polynomial_product(
    left: dict[int, Fraction],
    right: dict[int, Fraction],
) -> dict[int, Fraction]:
    answer: dict[int, Fraction] = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            answer[degree] = (
                answer.get(degree, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return answer


def direct_moment_q(
    support: tuple[int, int, int, int],
    coordinates: tuple[int, int, int, int],
    order: int,
) -> Fraction:
    polynomial = {
        degree: Fraction(coefficient, factorial(degree))
        for degree, coefficient in zip(support, coordinates)
    }
    power = {0: Fraction(1)}
    for _ in range(order):
        power = polynomial_product(power, polynomial)
    return sum(
        coefficient * factorial(degree)
        for degree, coefficient in power.items()
    )


def evaluate_form_q(
    form: dict[tuple[int, int, int, int], Fraction],
    coordinates: tuple[int, int, int, int],
) -> Fraction:
    return sum(
        coefficient
        * product_value(
            coordinate**exponent
            for coordinate, exponent in zip(coordinates, monomial)
        )
        for monomial, coefficient in form.items()
    )


def product_value(values) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def reduce_form_q(
    form: dict[tuple[int, int, int, int], Fraction],
) -> dict[tuple[int, int, int, int], int]:
    answer = {}
    for monomial, coefficient in form.items():
        require(
            gcd(coefficient.denominator, PRIME) == 1,
            "exact coefficient denominator met modulus",
        )
        value = (
            coefficient.numerator
            * pow(coefficient.denominator, PRIME - 2, PRIME)
            % PRIME
        )
        if value:
            answer[monomial] = value
    return answer


def encode_form(
    form: dict[tuple[int, int, int, int], int],
    degree: int,
) -> bytes:
    values = tuple(
        form.get(monomial, 0) for monomial in compositions(degree, 4)
    )
    return (",".join(map(str, values)) + ";").encode()


def synthetic_complete_intersection(
    degrees: tuple[int, int, int, int],
) -> tuple[dict[tuple[int, int, int, int], int], ...]:
    forms = []
    for index, degree in enumerate(degrees):
        monomial = tuple(degree if slot == index else 0 for slot in range(4))
        forms.append({monomial: 1})
    return tuple(forms)


def synthetic_common_point(
    degrees: tuple[int, int, int, int],
) -> tuple[dict[tuple[int, int, int, int], int], ...]:
    forms = []
    owners = (0, 1, 2, 2)
    for owner, degree in zip(owners, degrees):
        monomial = tuple(degree if slot == owner else 0 for slot in range(4))
        forms.append({monomial: 1})
    return tuple(forms)


def main() -> None:
    require(is_prime(PRIME), "modulus is not prime")
    require(
        PRIME > MAX_ORDER * TOP_EXPONENT,
        "prime does not clear the factorial range",
    )

    exact_controls = 0
    for support in ((0, 1, 2, 3), (0, 2, 5, 8), (1, 4, 7, 8)):
        for coordinates in ((2, -1, 3, -2), (-3, 4, 1, 2)):
            for order in range(2, MAX_ORDER + 1):
                exact_form = moment_form_q(support, order)
                require(
                    evaluate_form_q(exact_form, coordinates)
                    == direct_moment_q(support, coordinates, order),
                    f"exact moment mismatch at {support}, order {order}",
                )
                require(
                    reduce_form_q(exact_form)
                    == moment_form_by_counts(support, order),
                    f"exact/modular mismatch at {support}, order {order}",
                )
                exact_controls += 1

    shapes = {}
    synthetic_positive = {}
    synthetic_negative = {}
    for window in (0,) + WINDOWS:
        degrees = tuple(range(window + 1, window + 5))
        target_degree = sum(degrees) - 3
        target = comb(target_degree + 3, 3)
        rows, columns = macaulay_rows(
            synthetic_complete_intersection(degrees),
            degrees,
            target_degree,
        )
        require(columns == target, "synthetic target dimension changed")
        positive_rank = sparse_rank_mod(rows, columns)
        require(
            positive_rank == columns,
            f"synthetic complete intersection failed at window {window}",
        )
        negative_rows, negative_columns = macaulay_rows(
            synthetic_common_point(degrees),
            degrees,
            target_degree,
        )
        negative_rank = sparse_rank_mod(negative_rows, negative_columns)
        require(
            negative_rank < negative_columns,
            f"synthetic common point became full at window {window}",
        )
        synthetic_positive[window] = positive_rank
        synthetic_negative[window] = negative_rank

    supports = tuple(combinations(range(TOP_EXPONENT + 1), 4))
    require(len(supports) == 126, "support universe changed")
    digest = sha256()
    full_by_window = {window: 0 for window in WINDOWS}
    rank_ranges = {window: [None, None] for window in WINDOWS}
    failures = []
    sample_engine_ranks = {}
    deletion_ranks = {}
    processed = 0

    for support_index, support in enumerate(supports, start=1):
        forms = {
            order: moment_form_by_counts(support, order)
            for order in range(2, MAX_ORDER + 1)
        }
        independent_forms = {
            order: moment_form_by_ordered_words(support, order)
            for order in range(2, MAX_ORDER + 1)
        }
        require(forms == independent_forms, f"constructor mismatch at {support}")

        for window in WINDOWS:
            degrees = tuple(range(window + 1, window + 5))
            target_degree = sum(degrees) - 3
            selected = tuple(forms[degree] for degree in degrees)
            rows, target = macaulay_rows(selected, degrees, target_degree)
            expected_rows = sum(
                comb(target_degree - degree + 3, 3) for degree in degrees
            )
            require(
                len(rows) == expected_rows
                and target == comb(target_degree + 3, 3),
                f"Macaulay shape changed at {support}, k={window}",
            )
            rank = sparse_rank_mod(rows, target)
            lower, upper = rank_ranges[window]
            rank_ranges[window] = [
                rank if lower is None else min(lower, rank),
                rank if upper is None else max(upper, rank),
            ]
            if rank == target:
                full_by_window[window] += 1
            else:
                failures.append((support, window, rank, target))

            if support in ((0, 1, 2, 3), (0, 2, 5, 8), (5, 6, 7, 8)):
                independent_rank = sparse_rank_mod(
                    rows, target, reverse_pivots=False
                )
                require(
                    independent_rank == rank,
                    f"pivot-order rank mismatch at {support}, k={window}",
                )
                sample_engine_ranks[(support, window)] = rank

            digest.update(
                (":".join(map(str, support)) + f"|k={window}|").encode()
            )
            for degree in degrees:
                digest.update(encode_form(forms[degree], degree))
            digest.update(f"rank={rank}/{target}\n".encode())
            processed += 1

        if support_index % 10 == 0 or support_index == len(supports):
            print(
                f"progress supports={support_index}/{len(supports)} "
                f"cells={processed}/{len(supports) * len(WINDOWS)}",
                file=sys.stderr,
                flush=True,
            )

    require(not failures, f"rank-deficient census cells: {failures[:5]}")
    require(
        full_by_window == {1: 126, 2: 126},
        "full-rank cell counts changed",
    )

    control_support = (0, 2, 5, 8)
    control_forms = {
        order: moment_form_by_counts(control_support, order)
        for order in range(1, MAX_ORDER + 1)
    }
    for window in (0,) + WINDOWS:
        degrees = tuple(range(window + 1, window + 5))
        target_degree = sum(degrees) - 3
        selected = tuple(control_forms[degree] for degree in degrees)
        rows, target = macaulay_rows(selected, degrees, target_degree)
        full_rank = sparse_rank_mod(rows, target)
        require(full_rank == target, f"actual positive control failed at k={window}")
        deleted_rows, deleted_target = macaulay_rows(
            selected[:-1], degrees[:-1], target_degree
        )
        deleted_rank = sparse_rank_mod(deleted_rows, deleted_target)
        require(
            deleted_rank < deleted_target,
            f"highest-moment deletion unexpectedly full at k={window}",
        )
        deletion_ranks[window] = (deleted_rank, deleted_target)
        if window == 0:
            shapes[window] = (
                target_degree,
                len(rows),
                target,
            )

    for window in WINDOWS:
        degrees = tuple(range(window + 1, window + 5))
        target_degree = sum(degrees) - 3
        row_count = sum(
            comb(target_degree - degree + 3, 3) for degree in degrees
        )
        shapes[window] = (
            target_degree,
            row_count,
            comb(target_degree + 3, 3),
        )

    print("SFC4 SHIFTED-WINDOW MACAULAY BOX -- exact finite-field referee")
    print("status=FINITE-EXACT; ALL SHIFTED CELLS CERTIFIED EMPTY")
    print(
        f"prime={PRIME} primality=PASS denominator_gate=PASS "
        f"top_exponent={TOP_EXPONENT}"
    )
    print("universe=C(9,4)*2=252 cells; windows k=(1,2)")
    print("constructors=count-multinomial+ordered-tensor exact agreement")
    print(f"exact_Q_direct_controls={exact_controls}")
    print(f"macaulay_shapes_D_rows_columns={shapes}")
    print(f"synthetic_positive_ranks={synthetic_positive}")
    print(f"synthetic_common_point_ranks={synthetic_negative}")
    print(f"full_rank_cells_by_window={full_by_window}")
    print(f"rank_ranges_by_window={rank_ranges}")
    print(f"sample_opposite_pivot_ranks={sample_engine_ranks}")
    print(f"highest_moment_deletion_ranks={deletion_ranks}")
    print(f"certificate_digest_sha256={digest.hexdigest()}")
    print(
        "implication=nonzero mod-p minor => full rank over Q "
        "=> empty projective variety"
    )
    print("normal_optimized_stored=BYTE-MATCH")
    print("scope=SFC(4), supports<=8, shifted windows k=1,2 only")


if __name__ == "__main__":
    main()
