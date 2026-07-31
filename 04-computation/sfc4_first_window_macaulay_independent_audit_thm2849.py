#!/usr/bin/env python3
"""Independent exact audit of the THM-2849 first-window SFC(4) census.

This script does not import the primary companion and does not use NumPy.
It reconstructs the normalized moment forms from symmetric difference
tensors over Q, checks them against direct polynomial powers on hostile
samples, reduces them modulo the declared prime, and reruns all 1,820
Macaulay ranks with a sparse pure-Python eliminator using the opposite
pivot order from the primary implementation.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import factorial, gcd, isqrt


PRIME = 1_000_003
TOP_EXPONENT = 15
TARGET_DEGREE = 7
PRIMARY_DIGEST = "becbc51a7eeeb1ded8cea45d11a497a07427e38b789ee9767f9f17f963e4c10f"


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


MONOMIALS = {
    degree: compositions(degree, 3) for degree in range(TARGET_DEGREE + 1)
}


def normalized_moment_mod(
    support: tuple[int, int, int, int],
    directions: tuple[int, ...],
) -> int:
    numerator = factorial(sum(support[index] for index in directions))
    denominator = 1
    for index in directions:
        denominator *= factorial(support[index])
    require(gcd(denominator, PRIME) == 1, "factorial denominator met modulus")
    return (
        (numerator % PRIME)
        * pow(denominator % PRIME, PRIME - 2, PRIME)
        % PRIME
    )


def difference_tensor_mod(
    support: tuple[int, int, int, int],
    directions: tuple[int, ...],
) -> int:
    """L(prod_j(f_a[d_j]-f_a3)) by an unordered-count reconstruction."""
    answer = 0
    for mask in range(1 << len(directions)):
        selected = tuple(
            3 if mask & (1 << position) else direction
            for position, direction in enumerate(directions)
        )
        term = normalized_moment_mod(support, selected)
        if mask.bit_count() % 2:
            answer -= term
        else:
            answer += term
    return answer % PRIME


def moment_form_mod(
    support: tuple[int, int, int, int],
    order: int,
) -> dict[tuple[int, int, int], int]:
    answer: dict[tuple[int, int, int], int] = {}
    for exponents in compositions(order, 3):
        directions = tuple(
            direction
            for direction, count in enumerate(exponents)
            for _ in range(count)
        )
        multinomial = factorial(order)
        for count in exponents:
            multinomial //= factorial(count)
        coefficient = multinomial * difference_tensor_mod(support, directions)
        coefficient %= PRIME
        if coefficient:
            answer[exponents] = coefficient
    return answer


def macaulay_rows(
    forms: tuple[dict[tuple[int, int, int], int], ...],
    degrees: tuple[int, ...],
) -> list[dict[int, int]]:
    columns = {monomial: index for index, monomial in enumerate(MONOMIALS[7])}
    rows: list[dict[int, int]] = []
    for form, degree in zip(forms, degrees):
        for multiplier in MONOMIALS[7 - degree]:
            row: dict[int, int] = {}
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index] for index in range(3)
                )
                row[columns[target]] = coefficient
            rows.append(row)
    return rows


def sparse_rank_mod(rows: list[dict[int, int]], column_count: int = 36) -> int:
    """Exact row rank over F_p, pivoting from the opposite column end."""
    pivots: dict[int, dict[int, int]] = {}
    for original in rows:
        row = {
            column: value % PRIME
            for column, value in original.items()
            if value % PRIME
        }
        while row:
            column = max(row)
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
    directions: tuple[int, ...],
) -> Fraction:
    denominator = 1
    for index in directions:
        denominator *= factorial(support[index])
    return Fraction(
        factorial(sum(support[index] for index in directions)),
        denominator,
    )


def difference_tensor_q(
    support: tuple[int, int, int, int],
    directions: tuple[int, ...],
) -> Fraction:
    answer = Fraction(0)
    for mask in range(1 << len(directions)):
        selected = tuple(
            3 if mask & (1 << position) else direction
            for position, direction in enumerate(directions)
        )
        sign = -1 if mask.bit_count() % 2 else 1
        answer += sign * normalized_moment_q(support, selected)
    return answer


def moment_form_q(
    support: tuple[int, int, int, int],
    order: int,
) -> dict[tuple[int, int, int], Fraction]:
    answer: dict[tuple[int, int, int], Fraction] = {}
    for exponents in compositions(order, 3):
        directions = tuple(
            direction
            for direction, count in enumerate(exponents)
            for _ in range(count)
        )
        multinomial = factorial(order)
        for count in exponents:
            multinomial //= factorial(count)
        answer[exponents] = multinomial * difference_tensor_q(
            support, directions
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


def direct_moment(
    support: tuple[int, int, int, int],
    coordinates: tuple[int, int, int],
    order: int,
) -> Fraction:
    polynomial = {
        support[index]: Fraction(coordinates[index], factorial(support[index]))
        for index in range(3)
    }
    polynomial[support[3]] = (
        polynomial.get(support[3], Fraction(0))
        - Fraction(sum(coordinates), factorial(support[3]))
    )
    power = {0: Fraction(1)}
    for _ in range(order):
        power = polynomial_product(power, polynomial)
    return sum(
        coefficient * factorial(degree)
        for degree, coefficient in power.items()
    )


def evaluate_form(
    form: dict[tuple[int, int, int], Fraction],
    coordinates: tuple[int, int, int],
) -> Fraction:
    return sum(
        coefficient
        * coordinates[0] ** monomial[0]
        * coordinates[1] ** monomial[1]
        * coordinates[2] ** monomial[2]
        for monomial, coefficient in form.items()
    )


def reduce_form(
    form: dict[tuple[int, int, int], Fraction],
) -> dict[tuple[int, int, int], int]:
    answer = {}
    for monomial, coefficient in form.items():
        require(
            gcd(coefficient.denominator, PRIME) == 1,
            "reduced exact coefficient met modulus",
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
    form: dict[tuple[int, int, int], int],
    degree: int,
) -> bytes:
    values = tuple(form.get(monomial, 0) for monomial in MONOMIALS[degree])
    return (",".join(map(str, values)) + ";").encode()


def main() -> None:
    require(is_prime(PRIME), "modulus is not prime")
    require(PRIME > 4 * TOP_EXPONENT, "factorial denominator gate failed")

    exact_controls = 0
    for support in ((0, 1, 2, 3), (0, 5, 10, 15), (3, 7, 11, 15)):
        for coordinates in ((2, -1, 3), (-3, 4, 1)):
            require(
                direct_moment(support, coordinates, 1) == 0,
                f"mean elimination failed at {support}, {coordinates}",
            )
            for order in (2, 3, 4):
                exact_form = moment_form_q(support, order)
                require(
                    evaluate_form(exact_form, coordinates)
                    == direct_moment(support, coordinates, order),
                    f"exact-Q constructor failed at {support}, order {order}",
                )
                require(
                    reduce_form(exact_form) == moment_form_mod(support, order),
                    f"modular reduction failed at {support}, order {order}",
                )
                exact_controls += 1

    supports = tuple(combinations(range(TOP_EXPONENT + 1), 4))
    require(len(supports) == 1820, "support universe changed")
    digest = sha256()
    minimum_rank = 36
    maximum_rank = 0
    failures: list[tuple[tuple[int, int, int, int], int]] = []

    for support in supports:
        forms = tuple(moment_form_mod(support, order) for order in (2, 3, 4))
        rank = sparse_rank_mod(macaulay_rows(forms, (2, 3, 4)))
        minimum_rank = min(minimum_rank, rank)
        maximum_rank = max(maximum_rank, rank)
        if rank != 36:
            failures.append((support, rank))
        digest.update((":".join(map(str, support)) + "|").encode())
        for form, degree in zip(forms, (2, 3, 4)):
            digest.update(encode_form(form, degree))
        digest.update(f"rank={rank}\n".encode())

    require(not failures, f"independent rank failures: {failures[:5]}")
    require(minimum_rank == maximum_rank == 36, "rank range changed")
    require(digest.hexdigest() == PRIMARY_DIGEST, "primary digest mismatch")

    control_support = (1, 2, 3, 4)
    control_forms = tuple(
        moment_form_mod(control_support, order) for order in (2, 3, 4)
    )
    deletion_rank = sparse_rank_mod(
        macaulay_rows(control_forms[:2], (2, 3))
    )
    require(deletion_rank == 30, "quartic-deletion rank changed")

    fake_forms = (
        {(2, 0, 0): 1},
        {(3, 0, 0): 1},
        {(4, 0, 0): 1},
    )
    fake_rank = sparse_rank_mod(macaulay_rows(fake_forms, (2, 3, 4)))
    require(fake_rank == 21, "common-line hostile rank changed")

    print("THM-2849 INDEPENDENT MACAULAY AUDIT -- exact referee")
    print("status=INDEPENDENTLY-VERIFIED-EXACT")
    print("engine=standard-library sparse F_p elimination; opposite pivot order")
    print(f"prime={PRIME} primality=PASS denominator_gate=PASS")
    print(f"exact_Q_direct_controls={exact_controls} mean_controls=6")
    print(
        f"full_census_cells={len(supports)} "
        f"rank_range=({minimum_rank},{maximum_rank})"
    )
    print(f"independent_digest_sha256={digest.hexdigest()}")
    print(
        f"primary_digest_match=PASS deletion_rank={deletion_rank} "
        f"common_line_rank={fake_rank}"
    )
    print("mod_p_lifting_and_projective_emptiness_audit=PASS")


if __name__ == "__main__":
    main()
