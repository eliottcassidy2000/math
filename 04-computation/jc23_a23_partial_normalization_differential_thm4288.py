#!/usr/bin/env python3
"""Exact finite-jet audit for THM-4288.

The theorem is proved formally in its canon file.  This script supplies an
independent finite-linear-algebra control: it constructs the truncated pair
rings A_m and A_ell inside (F_p[b]/b^N)^2, checks the relative-differential
quotient and annihilator, checks the dual-number closed fibre, and recovers
the unique first-order circuit at the ordinary plane triple point.
"""

from __future__ import annotations

from collections.abc import Iterable


P = 1_000_000_007


def rank_mod(rows: Iterable[list[int]]) -> int:
    matrix = [[entry % P for entry in row] for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (r for r in range(pivot_row, len(matrix)) if matrix[r][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], P - 2, P)
        matrix[pivot_row] = [value * inverse % P for value in matrix[pivot_row]]
        for r in range(len(matrix)):
            if r == pivot_row or matrix[r][column] == 0:
                continue
            scale = matrix[r][column]
            matrix[r] = [
                (left - scale * right) % P
                for left, right in zip(matrix[r], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def add_sparse(n: int, left: dict[int, int], right: dict[int, int]) -> list[int]:
    vector = [0] * (2 * n)
    for degree, coefficient in left.items():
        if degree < n:
            vector[degree] = coefficient % P
    for degree, coefficient in right.items():
        if degree < n:
            vector[n + degree] = coefficient % P
    return vector


def diagonal(n: int, degree: int, coefficient: int = 1) -> list[int]:
    return add_sparse(n, {degree: coefficient}, {degree: coefficient})


def branch(n: int, degree: int, coefficient: int = 1) -> list[int]:
    return add_sparse(n, {}, {degree: coefficient})


def pair_basis(n: int, contact: int) -> list[list[int]]:
    return [diagonal(n, i) for i in range(n)] + [
        branch(n, j) for j in range(contact, n)
    ]


def multiply(n: int, x: list[int], y: list[int]) -> list[int]:
    out = [0] * (2 * n)
    for offset in (0, n):
        for i in range(n):
            if x[offset + i] == 0:
                continue
            for j in range(n - i):
                if y[offset + j]:
                    out[offset + i + j] = (
                        out[offset + i + j]
                        + x[offset + i] * y[offset + j]
                    ) % P
    return out


def ideal_span(n: int, basis: list[list[int]], generators: list[list[int]]) -> list[list[int]]:
    return [multiply(n, element, generator) for element in basis for generator in generators]


def in_span(vector: list[int], rows: list[list[int]]) -> bool:
    return rank_mod(rows + [vector]) == rank_mod(rows)


def audit_partial_normalization(m: int, ell: int) -> tuple[int, int, int]:
    if not 1 <= ell < m:
        raise ValueError((m, ell))
    n = 2 * m + 5
    base = pair_basis(n, m)
    upper = pair_basis(n, ell)
    if rank_mod(base) != 2 * n - m or rank_mod(upper) != 2 * n - ell:
        raise AssertionError("pair-ring basis rank changed")

    t = branch(n, ell)
    generated = base + [multiply(n, a, t) for a in base]
    if rank_mod(generated) != rank_mod(upper):
        raise AssertionError("A_m[t] does not span A_ell")

    exponent = m - ell
    differential_generators = [
        diagonal(n, exponent),
        add_sparse(n, {ell: -1}, {ell: 1}),
    ]
    differential_ideal = ideal_span(n, upper, differential_generators)
    omega_length = rank_mod(upper) - rank_mod(differential_ideal)
    expected = min(m - ell, 2 * ell)
    if omega_length != expected:
        raise AssertionError((m, ell, omega_length, expected))

    # Ann_{A_m}(B/I)=A_m intersect I.  Compare that intersection with the
    # independently generated candidate ideal (b^d,q).
    intersection_dimension = (
        rank_mod(base)
        + rank_mod(differential_ideal)
        - rank_mod(base + differential_ideal)
    )
    q = branch(n, m)
    candidate_annihilator = ideal_span(
        n, base, [diagonal(n, expected), q]
    )
    if not all(in_span(row, differential_ideal) for row in candidate_annihilator):
        raise AssertionError("candidate annihilator is not contained in I")
    if rank_mod(candidate_annihilator) != intersection_dimension:
        raise AssertionError("annihilator dimension mismatch")
    if expected and in_span(diagonal(n, expected - 1), differential_ideal):
        raise AssertionError("annihilator exponent is not sharp")

    closed_fibre_ideal = ideal_span(n, upper, [diagonal(n, 1), q])
    closed_fibre_dimension = rank_mod(upper) - rank_mod(closed_fibre_ideal)
    if closed_fibre_dimension != 2:
        raise AssertionError("closed fibre is not two-dimensional")
    if in_span(t, closed_fibre_ideal):
        raise AssertionError("dual-number generator vanished")
    if not in_span(multiply(n, t, t), closed_fibre_ideal):
        raise AssertionError("dual-number generator did not square to zero")

    conductor_loss = rank_mod(upper) - rank_mod(base)
    if conductor_loss != m - ell:
        raise AssertionError("conductor loss changed")
    return conductor_loss, omega_length, closed_fibre_dimension


def triple_vector(jets: int, i: int, j: int) -> list[int]:
    """Restrictions of b^i z^j to E:(b=0), R:(z=0), C:(z=b)."""
    out = [0] * (3 * jets)
    if i == 0 and j < jets:
        out[j] = 1
    if j == 0 and i < jets:
        out[jets + i] = 1
    if i + j < jets:
        out[2 * jets + i + j] = 1
    return out


def audit_triple_point(jets: int = 5) -> tuple[int, int, int]:
    image = [
        triple_vector(jets, i, j)
        for i in range(jets)
        for j in range(jets - i)
    ]
    image_rank = rank_mod(image)
    equal_value_dimension = 3 * jets - 2
    expected_image_rank = equal_value_dimension - 1
    if image_rank != expected_image_rank:
        raise AssertionError((image_rank, expected_image_rank))

    # lambda(e,r,c)=c_1-r_1-e_1.
    functional = [0] * (3 * jets)
    functional[1] = -1
    functional[jets + 1] = -1
    functional[2 * jets + 1] = 1
    for row in image:
        if sum(a * b for a, b in zip(functional, row)) % P:
            raise AssertionError("triple circuit does not annihilate image")

    hostile = [0] * (3 * jets)
    hostile[2 * jets + 1] = 1  # (e,r,c)=(0,0,b)
    if in_span(hostile, image):
        raise AssertionError("nonzero-lambda hostile descended")
    if rank_mod(image + [hostile]) != equal_value_dimension:
        raise AssertionError("lambda is not the sole equal-value obstruction")
    return image_rank, equal_value_dimension, 1


def main() -> None:
    print("THM4288_A23_PARTIAL_NORMALIZATION_DIFFERENTIAL_AUDIT_V1")
    print("FIELD F_1000000007 TRUNCATION N=2M+5")
    for m, ell in ((2, 1), (12, 1), (12, 2), (12, 4), (12, 11)):
        loss, omega, fibre = audit_partial_normalization(m, ell)
        print(
            f"PAIR M {m} ELL {ell} CONDUCTOR_LOSS {loss} "
            f"OMEGA_LENGTH {omega} ANN_A (b^{omega},q) "
            f"CLOSED_FIBRE_DIM {fibre} CLOSED_FIBRE k[t]/(t^2)"
        )
    image_rank, equal_dimension, defect = audit_triple_point()
    print(
        f"TRIPLE JETS 5 IMAGE_RANK {image_rank} "
        f"EQUAL_VALUE_DIM {equal_dimension} DEFECT {defect}"
    )
    print("TRIPLE_FUNCTIONAL lambda=c1-r1-e1 HOSTILE (0,0,b)")
    print("VERDICT PASS FINITE_JET_LINEAR_ALGEBRA_CONTROLS")


if __name__ == "__main__":
    main()
