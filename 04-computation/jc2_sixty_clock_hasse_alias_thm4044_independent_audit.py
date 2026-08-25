#!/usr/bin/env python3
"""Independent matrix audit for THM-4044 over F_61."""

from __future__ import annotations

import hashlib
import math
import sys


sys.stdout.reconfigure(newline="\n")
P = 61
PHI = 44
CHECKS = 0
SEMANTIC: list[str] = []


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)
    print(f"PASS  {label}")


def rank_mod(matrix: list[list[int]]) -> int:
    if not matrix:
        return 0
    a = [[entry % P for entry in row] for row in matrix]
    rows = len(a)
    columns = len(a[0])
    rank = 0
    for column in range(columns):
        pivot = next((r for r in range(rank, rows) if a[r][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][column], -1, P)
        a[rank] = [(inverse * value) % P for value in a[rank]]
        for r in range(rows):
            if r == rank or not a[r][column]:
                continue
            factor = a[r][column]
            a[r] = [(left - factor * right) % P for left, right in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def hasse_entry(degree: int, order: int, point: int) -> int:
    if degree < order:
        return 0
    return math.comb(degree, order) * pow(point, degree - order, P) % P


def observer_matrix(k: int, maximum_degree: int, add_boundary_two: bool = False) -> list[list[int]]:
    degrees = list(range(2, maximum_degree + 1))
    nodes = [pow(PHI, r, P) for r in range(60)]
    matrix = [
        [hasse_entry(degree, order, point) for degree in degrees]
        for point in nodes
        for order in range(k)
    ]
    # Boundary orders 0 and 1 vanish identically on P^2 K[P].  The next row
    # is the first informative sidecar.
    if add_boundary_two:
        matrix.append([hasse_entry(degree, 2, 0) for degree in degrees])
    return matrix


def multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % P
    return result


def polynomial_power(base: list[int], exponent: int) -> list[int]:
    result = [1]
    for _ in range(exponent):
        result = multiply(result, base)
    return result


print("STATUS=THM-4044;INDEPENDENT_CONFLUENT_MATRIX_AUDIT")
for k in (1, 2, 3):
    below = 60 * k + 1
    first_failure = 60 * k + 2
    rank_below = rank_mod(observer_matrix(k, below))
    rank_failure = rank_mod(observer_matrix(k, first_failure))
    rank_repaired = rank_mod(observer_matrix(k, first_failure, add_boundary_two=True))
    require(rank_below == 60 * k, f"O_{k} injective on P^2 degrees through 60k+1")
    require(rank_failure == 60 * k, f"O_{k} has a one-dimensional first kernel at degree 60k+2")
    require(rank_repaired == 60 * k + 1, f"boundary Hasse order two repairs first O_{k} kernel")

    delta = [0, 0] + polynomial_power([P - 1] + [0] * 59 + [1], k)
    degrees = list(range(2, first_failure + 1))
    vector = [delta[d] if d < len(delta) else 0 for d in degrees]
    matrix = observer_matrix(k, first_failure)
    require(
        all(sum(entry * coefficient for entry, coefficient in zip(row, vector)) % P == 0 for row in matrix),
        f"explicit Delta_{k} spans the first matrix kernel",
    )
    require(vector[0] == (-1) ** k % P and vector[-1] == 1, f"Delta_{k} endpoint coefficients")
    SEMANTIC.append(f"k={k}:ranks={rank_below},{rank_failure},{rank_repaired}")

print("RESULT independent_rank_boundary=" + ";".join(SEMANTIC))
semantic_hash = hashlib.sha256("\n".join(SEMANTIC).encode()).hexdigest()
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"CHECKS={CHECKS}")
print("ALL THM-4044 INDEPENDENT CHECKS PASSED")

