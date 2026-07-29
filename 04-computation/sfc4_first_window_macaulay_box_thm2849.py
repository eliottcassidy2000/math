#!/usr/bin/env python3
"""Exact finite-field referee for the THM-2849 first-window SFC(4) box."""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations, product
from math import factorial, isqrt

import numpy as np


PRIME = 1_000_003
TOP_EXPONENT = 15
TARGET_DEGREE = 7


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


FACT = tuple(factorial(index) % PRIME for index in range(4 * TOP_EXPONENT + 1))
INV_FACT = tuple(pow(value, PRIME - 2, PRIME) for value in FACT)


def factorial_moment(degrees: tuple[int, ...]) -> int:
    value = FACT[sum(degrees)]
    for degree in degrees:
        value = value * INV_FACT[degree] % PRIME
    return value


def moment_form_by_substitution(
    support: tuple[int, int, int, int], order: int
) -> dict[tuple[int, int, int], int]:
    """Expand H=sum_(i<3)x_i f_ai-(sum x_i)f_a3 directly."""
    answer: dict[tuple[int, int, int], int] = {}
    for exponents in compositions(order, 4):
        multinomial = FACT[order]
        degrees = []
        for exponent, slot in zip(exponents, support):
            multinomial = multinomial * INV_FACT[exponent] % PRIME
            degrees.extend([slot] * exponent)
        coefficient = multinomial * factorial_moment(tuple(degrees)) % PRIME
        fourth = exponents[3]
        if fourth % 2:
            coefficient = -coefficient % PRIME
        for tail in compositions(fourth, 3):
            tail_multinomial = FACT[fourth]
            for exponent in tail:
                tail_multinomial = tail_multinomial * INV_FACT[exponent] % PRIME
            monomial = tuple(exponents[index] + tail[index] for index in range(3))
            answer[monomial] = (
                answer.get(monomial, 0) + coefficient * tail_multinomial
            ) % PRIME
    return {monomial: value for monomial, value in answer.items() if value}


def difference_tensor(
    support: tuple[int, int, int, int], directions: tuple[int, ...]
) -> int:
    """Compute L(prod_j(f_a[d_j]-f_a3)) by its independent 2^m expansion."""
    answer = 0
    for mask in range(1 << len(directions)):
        degrees = []
        sign = 1
        for position, direction in enumerate(directions):
            if mask & (1 << position):
                degrees.append(support[3])
                sign = -sign
            else:
                degrees.append(support[direction])
        answer = (answer + sign * factorial_moment(tuple(degrees))) % PRIME
    return answer


def moment_form_by_ordered_tensors(
    support: tuple[int, int, int, int], order: int
) -> dict[tuple[int, int, int], int]:
    """Independent ordered-tensor construction of the same homogeneous form."""
    answer: dict[tuple[int, int, int], int] = {}
    for directions in product(range(3), repeat=order):
        monomial = tuple(directions.count(index) for index in range(3))
        answer[monomial] = (
            answer.get(monomial, 0) + difference_tensor(support, directions)
        ) % PRIME
    return {monomial: value for monomial, value in answer.items() if value}


MONOMIALS = {
    degree: compositions(degree, 3) for degree in range(TARGET_DEGREE + 1)
}


def macaulay_rows(
    forms: tuple[dict[tuple[int, int, int], int], ...],
    degrees: tuple[int, ...],
) -> list[list[int]]:
    columns = MONOMIALS[TARGET_DEGREE]
    column_index = {monomial: index for index, monomial in enumerate(columns)}
    rows: list[list[int]] = []
    for form, degree in zip(forms, degrees):
        for multiplier in MONOMIALS[TARGET_DEGREE - degree]:
            row = [0] * len(columns)
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index] for index in range(3)
                )
                row[column_index[target]] = coefficient
            rows.append(row)
    return rows


def rank_mod_numpy(rows: list[list[int]]) -> int:
    matrix = np.asarray(rows, dtype=np.int64)
    row_count, column_count = matrix.shape
    rank = 0
    for column in range(column_count):
        candidates = np.flatnonzero(matrix[rank:, column])
        if len(candidates) == 0:
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            matrix[[rank, pivot]] = matrix[[pivot, rank]]
        inverse = pow(int(matrix[rank, column]), PRIME - 2, PRIME)
        matrix[rank, column:] = matrix[rank, column:] * inverse % PRIME
        if rank + 1 < row_count:
            factors = matrix[rank + 1 :, column].copy()
            active = np.flatnonzero(factors)
            if len(active):
                target_rows = rank + 1 + active
                matrix[target_rows, column:] = (
                    matrix[target_rows, column:]
                    - factors[active, None] * matrix[rank, column:]
                ) % PRIME
        rank += 1
        if rank == row_count or rank == column_count:
            break
    return rank


def rank_mod_python(rows: list[list[int]]) -> int:
    matrix = [row[:] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if matrix[row][column]), None
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], PRIME - 2, PRIME)
        matrix[rank][column:] = [
            value * inverse % PRIME for value in matrix[rank][column:]
        ]
        for row in range(rank + 1, row_count):
            factor = matrix[row][column]
            if factor:
                matrix[row][column:] = [
                    (left - factor * right) % PRIME
                    for left, right in zip(
                        matrix[row][column:], matrix[rank][column:]
                    )
                ]
        rank += 1
        if rank == row_count or rank == column_count:
            break
    return rank


def encode_form(form: dict[tuple[int, int, int], int], degree: int) -> bytes:
    values = tuple(form.get(monomial, 0) for monomial in MONOMIALS[degree])
    return (",".join(map(str, values)) + ";").encode()


def main() -> None:
    require(is_prime(PRIME), "modulus is no longer prime")
    require(PRIME > 4 * TOP_EXPONENT, "prime no longer clears denominators")
    supports = tuple(combinations(range(TOP_EXPONENT + 1), 4))
    require(len(supports) == 1820, "support universe changed")

    digest = sha256()
    minimum_rank = len(MONOMIALS[TARGET_DEGREE])
    maximum_rank = 0
    rank_failures: list[tuple[int, int, int, int]] = []
    sample_ranks: dict[tuple[int, int, int, int], int] = {}

    for support in supports:
        forms = tuple(moment_form_by_substitution(support, order) for order in (2, 3, 4))
        independent_forms = tuple(
            moment_form_by_ordered_tensors(support, order) for order in (2, 3, 4)
        )
        require(forms == independent_forms, f"moment-form mismatch at {support}")
        require(
            all(
                all(sum(monomial) == degree for monomial in form)
                for form, degree in zip(forms, (2, 3, 4))
            ),
            f"homogeneity failure at {support}",
        )

        rows = macaulay_rows(forms, (2, 3, 4))
        require(len(rows) == 46 and len(rows[0]) == 36, "Macaulay shape changed")
        rank = rank_mod_numpy(rows)
        minimum_rank = min(minimum_rank, rank)
        maximum_rank = max(maximum_rank, rank)
        if rank != 36:
            rank_failures.append(support)
        if support in ((0, 1, 2, 3), (1, 2, 3, 4), (0, 5, 10, 15)):
            require(rank_mod_python(rows) == rank, f"rank-engine mismatch at {support}")
            sample_ranks[support] = rank

        digest.update((":".join(map(str, support)) + "|").encode())
        for form, degree in zip(forms, (2, 3, 4)):
            digest.update(encode_form(form, degree))
        digest.update(f"rank={rank}\n".encode())

    require(not rank_failures, f"rank-deficient cells: {rank_failures[:5]}")
    require(minimum_rank == maximum_rank == 36, "full-rank census changed")

    control_support = (1, 2, 3, 4)
    control_forms = tuple(
        moment_form_by_substitution(control_support, order) for order in (2, 3, 4)
    )
    qc_rows = macaulay_rows(control_forms[:2], (2, 3))
    qc_rank = rank_mod_numpy(qc_rows)
    require(qc_rank < 36, "quartic-deletion hostile unexpectedly became full rank")
    require(
        rank_mod_python(qc_rows) == qc_rank,
        "quartic-deletion rank-engine mismatch",
    )

    print("SFC4 FIRST-WINDOW MACAULAY BOX -- exact finite-field referee")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
    print(f"prime={PRIME} top_exponent={TOP_EXPONENT} target_degree={TARGET_DEGREE}")
    print("universe=C(16,4)=1820 supports; moments=(1,2,3,4)")
    print("mean_elimination=H=sum_(i<3)x_i(f_ai-f_a3)")
    print("constructors=substitution+ordered-difference-tensor exact agreement")
    print("macaulay_rows=Q*R5(21)+C*R4(15)+F*R3(10)=46 columns=R7=36")
    print(f"full_rank_cells={len(supports)} rank_range=({minimum_rank},{maximum_rank})")
    print(f"sample_ranks={sample_ranks}")
    print(f"quartic_deleted_control=(1,2,3,4) rank={qc_rank}<36")
    print(f"certificate_digest_sha256={digest.hexdigest()}")
    print("implication=nonzero mod-p minor => rank36 over Q => empty projective variety")
    print("scope=SFC(4) k=0 and max support exponent15 only")


if __name__ == "__main__":
    main()
