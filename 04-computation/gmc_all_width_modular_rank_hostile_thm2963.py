#!/usr/bin/env python3
"""Independent exact modular-rank hostile for THM-2963.

For support (n,n+a,n+b,n+M), this script constructs the mean-zero
quadratic/cubic/quartic moment forms directly from four-variable factorial
moments, substitutes the fourth coefficient as -(x+y+z), constructs the full
degree-seven 46-by-36 Macaulay matrix, and computes its rank over

    F_p,  p = the first prime strictly larger than 4(n+M).

It intentionally imports none of the fixed-chart/cofactor constructors.
"""

from __future__ import annotations

import argparse
import hashlib
import random
from collections import Counter, defaultdict
from itertools import combinations
from math import factorial, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


MONOMIALS = {degree: tuple(compositions(degree, 3)) for degree in range(8)}
TARGETS = MONOMIALS[7]
TARGET_INDEX = {monomial: index for index, monomial in enumerate(TARGETS)}


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def next_prime(value: int) -> int:
    candidate = value + 1
    if candidate <= 2:
        return 2
    if candidate % 2 == 0:
        candidate += 1
    while not is_prime(candidate):
        candidate += 2
    return candidate


def factorial_table(limit: int, prime: int) -> list[int]:
    require(limit < prime, "factorial table crosses the characteristic")
    table = [1] * (limit + 1)
    for index in range(1, limit + 1):
        table[index] = table[index - 1] * index % prime
    return table


def multinomial(counts: tuple[int, ...]) -> int:
    answer = factorial(sum(counts))
    for count in counts:
        answer //= factorial(count)
    return answer


def forms_mod(
    depth: int,
    first: int,
    second: int,
    width: int,
    prime: int,
):
    """Direct four-variable expansion, normalized by L(f_n^m)."""
    support = (depth, depth + first, depth + second, depth + width)
    table = factorial_table(4 * (depth + width), prime)
    support_factorials = tuple(table[value] for value in support)
    forms = []
    for order in (2, 3, 4):
        base_inverse = (
            pow(table[depth], order, prime)
            * pow(table[order * depth], prime - 2, prime)
        ) % prime
        form = {monomial: 0 for monomial in MONOMIALS[order]}
        for counts in compositions(order, 4):
            weighted_degree = sum(
                count * exponent
                for count, exponent in zip(counts, support)
            )
            denominator = 1
            for index, count in enumerate(counts):
                denominator = (
                    denominator
                    * pow(support_factorials[index], count, prime)
                ) % prime
            coefficient = (
                multinomial(counts)
                * table[weighted_degree]
                * pow(denominator, prime - 2, prime)
                * base_inverse
            ) % prime
            fourth_count = counts[3]
            for tail in compositions(fourth_count, 3):
                target = tuple(
                    counts[index] + tail[index]
                    for index in range(3)
                )
                contribution = coefficient * multinomial(tail)
                if fourth_count % 2:
                    contribution = -contribution
                form[target] = (form[target] + contribution) % prime
        forms.append(form)
    return tuple(forms)


def macaulay_matrix(forms, prime: int) -> list[list[int]]:
    rows: list[list[int]] = []
    for degree, form in zip((2, 3, 4), forms):
        for multiplier in MONOMIALS[7 - degree]:
            row = [0] * len(TARGETS)
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[TARGET_INDEX[target]] = coefficient % prime
            rows.append(row)
    require(len(rows) == 46 and len(rows[0]) == 36, "Macaulay shape changed")
    return rows


def rank_mod(matrix: list[list[int]], prime: int) -> tuple[int, tuple[int, ...]]:
    rows = [[value % prime for value in row] for row in matrix]
    row_count = len(rows)
    column_count = len(rows[0]) if rows else 0
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if rows[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], prime - 2, prime)
        rows[pivot_row] = [
            value * inverse % prime for value in rows[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            coefficient = rows[row][column]
            if coefficient:
                rows[row] = [
                    (left - coefficient * right) % prime
                    for left, right in zip(rows[row], rows[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row, tuple(pivot_columns)


def one_case(depth: int, first: int, second: int, width: int):
    require(0 < first < second < width, "bad normalized support")
    prime = next_prime(4 * (depth + width))
    forms = forms_mod(depth, first, second, width, prime)
    matrix = macaulay_matrix(forms, prime)
    rank, pivots = rank_mod(matrix, prime)
    form_nonzero_counts = tuple(
        sum(coefficient != 0 for coefficient in form.values())
        for form in forms
    )
    return prime, rank, pivots, form_nonzero_counts


def exhaustive_zero_depth_cases(max_width: int):
    for width in range(3, max_width + 1):
        for first, second in combinations(range(1, width), 2):
            yield 0, first, second, width


def exhaustive_extremal_cases(max_width: int):
    for width in range(3, max_width + 1):
        for first, second in combinations(range(1, width), 2):
            for depth in (1, width - 1, width, 2 * width):
                yield depth, first, second, width


def random_cases(seed: int, count: int, max_width: int, max_depth: int):
    generator = random.Random(seed)
    for _ in range(count):
        width = generator.randint(3, max_width)
        first, second = sorted(generator.sample(range(1, width), 2))
        depth = generator.randint(0, max_depth)
        yield depth, first, second, width


def evaluate_form(form, point: tuple[int, int, int], prime: int) -> int:
    answer = 0
    for monomial, coefficient in form.items():
        term = coefficient
        for coordinate, exponent in zip(point, monomial):
            term = term * pow(coordinate, exponent, prime) % prime
        answer = (answer + term) % prime
    return answer


def rational_projective_roots(forms, prime: int):
    roots: list[tuple[int, int, int]] = []
    for y in range(prime):
        for z in range(prime):
            point = (1, y, z)
            if all(evaluate_form(form, point, prime) == 0 for form in forms):
                roots.append(point)
    for z in range(prime):
        point = (0, 1, z)
        if all(evaluate_form(form, point, prime) == 0 for form in forms):
            roots.append(point)
    point = (0, 0, 1)
    if all(evaluate_form(form, point, prime) == 0 for form in forms):
        roots.append(point)
    return tuple(roots)


def pivot_minor_certificate(matrix: list[list[int]], prime: int):
    reduced = [[value % prime for value in row] for row in matrix]
    identities = list(range(len(reduced)))
    pivot_rows: list[int] = []
    pivot_columns: list[int] = []
    rank = 0
    for column in range(len(reduced[0])):
        pivot = next(
            (
                row
                for row in range(rank, len(reduced))
                if reduced[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        identities[rank], identities[pivot] = identities[pivot], identities[rank]
        inverse = pow(reduced[rank][column], prime - 2, prime)
        reduced[rank] = [
            value * inverse % prime for value in reduced[rank]
        ]
        for row in range(rank + 1, len(reduced)):
            coefficient = reduced[row][column]
            if coefficient:
                reduced[row] = [
                    (left - coefficient * right) % prime
                    for left, right in zip(reduced[row], reduced[rank])
                ]
        pivot_rows.append(identities[rank])
        pivot_columns.append(column)
        rank += 1
    square = [
        [matrix[row][column] % prime for column in pivot_columns]
        for row in pivot_rows
    ]
    determinant = 1
    for column in range(rank):
        pivot = next(
            row for row in range(column, rank) if square[row][column]
        )
        if pivot != column:
            square[column], square[pivot] = square[pivot], square[column]
            determinant = -determinant
        pivot_value = square[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        for row in range(column + 1, rank):
            coefficient = square[row][column] * inverse % prime
            if coefficient:
                square[row] = [
                    (left - coefficient * right) % prime
                    for left, right in zip(square[row], square[column])
                ]
    return (
        rank,
        tuple(pivot_rows),
        tuple(pivot_columns),
        determinant % prime,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zero-depth-width", type=int, default=30)
    parser.add_argument("--extremal-width", type=int, default=15)
    parser.add_argument("--random-width", type=int, default=120)
    parser.add_argument("--random-depth", type=int, default=2000)
    parser.add_argument("--random", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=20260729)
    arguments = parser.parse_args()

    # First freeze the counterexample independently of the finite scan.
    hostile = (0, 5, 6, 10)
    hostile_prime = next_prime(4 * (hostile[0] + hostile[3]))
    require(hostile_prime == 41, "hostile next-prime selector changed")
    hostile_forms = forms_mod(*hostile, hostile_prime)
    hostile_matrix = macaulay_matrix(hostile_forms, hostile_prime)
    hostile_certificate = pivot_minor_certificate(
        hostile_matrix, hostile_prime
    )
    require(
        hostile_certificate[0] == 34
        and hostile_certificate[2] == tuple(range(34))
        and hostile_certificate[3] == 12,
        "hostile rank/minor certificate changed",
    )
    hostile_qc_rank = rank_mod(hostile_matrix[:36], hostile_prime)[0]
    require(hostile_qc_rank == 30, "hostile Q/C complete-intersection rank")
    hostile_roots = ((1, 22, 0), (1, 39, 38))
    require(
        all(
            evaluate_form(form, point, hostile_prime) == 0
            for point in hostile_roots
            for form in hostile_forms
        ),
        "hostile common roots changed",
    )
    evaluation_vectors = [
        tuple(
            pow(point[0], monomial[0], hostile_prime)
            * pow(point[1], monomial[1], hostile_prime)
            * pow(point[2], monomial[2], hostile_prime)
            % hostile_prime
            for monomial in TARGETS
        )
        for point in hostile_roots
    ]
    require(
        rank_mod(evaluation_vectors, hostile_prime)[0] == 2,
        "hostile evaluation covectors lost independence",
    )
    hostile_form_payload = "|".join(
        f"{order}:{monomial}:{form[monomial]}"
        for order, form in zip((2, 3, 4), hostile_forms)
        for monomial in MONOMIALS[order]
    )
    hostile_form_digest = hashlib.sha256(
        hostile_form_payload.encode()
    ).hexdigest()

    cases = list(
        exhaustive_zero_depth_cases(arguments.zero_depth_width)
    )
    cases.extend(exhaustive_extremal_cases(arguments.extremal_width))
    cases.extend(
        random_cases(
            arguments.seed,
            arguments.random,
            arguments.random_width,
            arguments.random_depth,
        )
    )
    # Stable duplicate removal.
    cases = list(dict.fromkeys(cases))

    rank_histogram: Counter[int] = Counter()
    prime_histogram: Counter[int] = Counter()
    rank34: list[tuple[int, int, int, int, int]] = []
    rank35: list[tuple[int, int, int, int, int]] = []
    rank_le33: list[tuple[int, int, int, int, int, int]] = []
    width_minimum: dict[int, int] = defaultdict(lambda: 36)
    payload: list[str] = []
    for index, (depth, first, second, width) in enumerate(cases, 1):
        prime, rank, pivots, support_counts = one_case(
            depth, first, second, width
        )
        rank_histogram[rank] += 1
        prime_histogram[prime] += 1
        width_minimum[width] = min(width_minimum[width], rank)
        payload.append(
            f"{depth}:{first}:{second}:{width}:{prime}:{rank}:"
            + ",".join(map(str, pivots))
            + ":"
            + ",".join(map(str, support_counts))
        )
        if rank == 34:
            rank34.append((depth, first, second, width, prime))
        elif rank == 35:
            rank35.append((depth, first, second, width, prime))
        if rank <= 33:
            rank_le33.append((depth, first, second, width, prime, rank))
        if index % 1000 == 0:
            print(f"progress={index}/{len(cases)}")

    digest = hashlib.sha256(("\n".join(payload) + "\n").encode()).hexdigest()
    expected_rank34 = (
        (0, 5, 6, 10, 41),
        (0, 4, 7, 15, 61),
        (0, 8, 16, 19, 79),
    )
    require(
        tuple(rank34) == expected_rank34,
        "finite-universe rank-34 census changed",
    )
    rank34_records = []
    expected_roots = {
        (0, 5, 6, 10, 41): ((1, 22, 0), (1, 39, 38)),
        (0, 4, 7, 15, 61): (),
        (0, 8, 16, 19, 79): ((1, 29, 51), (1, 49, 20)),
    }
    for case in rank34:
        depth, first, second, width, prime = case
        forms = forms_mod(depth, first, second, width, prime)
        roots = rational_projective_roots(forms, prime)
        require(roots == expected_roots[case], "rank-34 root census changed")
        second_prime = next_prime(prime)
        second_rank = rank_mod(
            macaulay_matrix(
                forms_mod(depth, first, second, width, second_prime),
                second_prime,
            ),
            second_prime,
        )[0]
        require(second_rank == 36, "second-prime repair changed")
        rank34_records.append(
            f"{case}:Fp_roots={roots}:next={second_prime}:rank{second_rank}"
        )

    print("GMC ALL-WIDTH DIRECT MODULAR MACAULAY RANK PROBE")
    print(
        f"universe=all_n0_M3..{arguments.zero_depth_width}"
        f"+all_extreme_depths_M3..{arguments.extremal_width}"
        f"+random{arguments.random}_M3..{arguments.random_width}"
        f"_n0..{arguments.random_depth}_seed{arguments.seed}"
    )
    print(
        "hostile=n0_support(0,5,6,10)_p41;"
        "rank34;common_roots=(1:22:0),(1:39:38);"
        "pivot34_minor=12"
    )
    print(
        f"hostile_QC_rank={hostile_qc_rank};"
        f"hostile_form_sha256={hostile_form_digest}"
    )
    print(
        "hostile_kummer=all_factorial_arguments_at_most40<41;"
        "all_factorial_valuations_zero"
    )
    print(f"cases={len(cases)}")
    print(
        "rank_histogram="
        + ",".join(f"{rank}:{rank_histogram[rank]}" for rank in sorted(rank_histogram))
    )
    print(f"rank35_count={len(rank35)}")
    print(f"rank34_count={len(rank34)}")
    print(f"rank_le33_count={len(rank_le33)}")
    print(
        "width_minima="
        + ",".join(f"{width}:{width_minimum[width]}" for width in sorted(width_minimum))
    )
    print(f"prime_count={len(prime_histogram)}")
    print(f"record_sha256={digest}")
    print("rank34_cases=" + ";".join(rank34_records))
    if rank35:
        print("first_rank35=" + ";".join(map(str, rank35[:20])))
    if rank_le33:
        print("first_rank_le33=" + ";".join(map(str, rank_le33[:20])))
        raise RuntimeError("rank<=33 hostile found")
    print(
        "scope=finite-exact hostile refutes rank>=35;"
        "rank>=34 remains finite-exact evidence only"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
