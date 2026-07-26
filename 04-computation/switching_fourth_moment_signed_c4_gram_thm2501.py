#!/usr/bin/env python3
"""Exact referee for THM-2501.

Checks the Rademacher second/fourth moments of a signed complete-graph
quadratic form, the signed-C4 expansion, the symmetric Gram-energy
identity, switching invariance, and the parity floor for off-diagonal
row correlations.  Every truth-bearing check raises explicitly under
``python -O``.
"""

from __future__ import annotations

from itertools import combinations, product
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


def sign_matrix(n: int, signs: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    matrix = [[0] * n for _ in range(n)]
    for (i, j), sign in zip(edges(n), signs, strict=True):
        matrix[i][j] = sign
        matrix[j][i] = sign
    return tuple(tuple(row) for row in matrix)


def quadratic(
    edge_bank: tuple[tuple[int, int], ...],
    signs: tuple[int, ...],
    switch: tuple[int, ...],
) -> int:
    return sum(
        sign * switch[i] * switch[j]
        for (i, j), sign in zip(edge_bank, signs, strict=True)
    )


def signed_c4_sum(
    n: int,
    matrix: tuple[tuple[int, ...], ...],
) -> int:
    total = 0
    for a, b, c, d in combinations(range(n), 4):
        total += matrix[a][b] * matrix[b][c] * matrix[c][d] * matrix[d][a]
        total += matrix[a][b] * matrix[b][d] * matrix[d][c] * matrix[c][a]
        total += matrix[a][c] * matrix[c][b] * matrix[b][d] * matrix[d][a]
    return total


def gram_energy(matrix: tuple[tuple[int, ...], ...]) -> int:
    n = len(matrix)
    return sum(
        sum(matrix[i][k] * matrix[j][k] for k in range(n)) ** 2
        for i in range(n)
        for j in range(n)
    )


def switched_matrix(
    matrix: tuple[tuple[int, ...], ...],
    switch: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    n = len(matrix)
    return tuple(
        tuple(switch[i] * matrix[i][j] * switch[j] for j in range(n))
        for i in range(n)
    )


def audit_one(n: int, signs: tuple[int, ...]) -> int:
    edge_bank = edges(n)
    matrix = sign_matrix(n, signs)
    switches = tuple(product((-1, 1), repeat=n))
    values = tuple(quadratic(edge_bank, signs, switch) for switch in switches)
    denominator = 1 << n
    moment_2_numerator = sum(value * value for value in values)
    moment_4_numerator = sum(value**4 for value in values)

    m = comb(n, 2)
    c4 = signed_c4_sum(n, matrix)
    gram = gram_energy(matrix)
    gram_base = n * (n - 1) * (2 * n - 3)
    fourth_formula = 3 * m * m - 2 * m + 24 * c4

    require(moment_2_numerator == denominator * m, "second moment failed")
    require(
        moment_4_numerator == denominator * fourth_formula,
        "fourth-moment C4 formula failed",
    )
    require(gram == gram_base + 8 * c4, "Gram/C4 identity failed")
    require(
        moment_4_numerator
        == denominator
        * (3 * m * m - 2 * m + 3 * (gram - gram_base)),
        "fourth-moment Gram formula failed",
    )
    require(
        max(abs(value) for value in values) ** 4 >= fourth_formula,
        "fourth-moment maximum bound failed",
    )

    parity_floor = n * (n - 1) ** 2
    if n % 2:
        parity_floor += n * (n - 1)
    require(gram >= parity_floor, "Gram parity floor failed")

    for switch in switches[: min(16, len(switches))]:
        switched = switched_matrix(matrix, switch)
        require(
            gram_energy(switched) == gram
            and signed_c4_sum(n, switched) == c4,
            "switching invariant failed",
        )
        upper_sum = sum(switched[i][j] for i, j in edge_bank)
        require(
            upper_sum == quadratic(edge_bank, signs, switch),
            "switching imbalance identity failed",
        )

    # A skew orientation matrix has identically zero x^T K x; it is
    # not the symmetric signed graph used by the quadratic problem.
    skew = tuple(
        tuple(
            0
            if i == j
            else (matrix[i][j] if i < j else -matrix[i][j])
            for j in range(n)
        )
        for i in range(n)
    )
    for switch in switches[: min(16, len(switches))]:
        skew_form = sum(
            switch[i] * skew[i][j] * switch[j]
            for i in range(n)
            for j in range(n)
        )
        require(skew_form == 0, "skew tournament hostile failed")

    return len(switches)


exhaustive_matrices = 0
switch_checks = 0
for n in range(2, 6):
    m = comb(n, 2)
    for signs in product((-1, 1), repeat=m):
        switch_checks += audit_one(n, signs)
        exhaustive_matrices += 1

# Deterministic larger controls, including all-positive, alternating,
# and a reproducible nonlinear sign stream.
large_controls = 0
state = 0x5EED
for n in range(6, 10):
    m = comb(n, 2)
    controls = [
        (1,) * m,
        tuple(-1 if index % 2 else 1 for index in range(m)),
    ]
    for _ in range(6):
        row = []
        for _ in range(m):
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            row.append(1 if state & 1 else -1)
        controls.append(tuple(row))
    for signs in controls:
        switch_checks += audit_one(n, signs)
        large_controls += 1

print("THM-2501 exact switching-moment referee")
print(f"exhaustive_signed_matrices={exhaustive_matrices}")
print(f"larger_deterministic_controls={large_controls}")
print(f"rademacher_switch_evaluations={switch_checks}")
print("second_moment=m_PASS")
print("fourth_moment=3m^2-2m+24_signed_C4_PASS")
print("gram_energy=base+8_signed_C4_PASS")
print("switching_invariance_and_parity_floor=PASS")
print("skew_tournament_quadratic_hostile=PASS")
print("ALL CHECKS PASSED")
