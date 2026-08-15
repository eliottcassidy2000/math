#!/usr/bin/env python3
"""Exact controls for THM-3439's near-identity Grassmannian orbit law.

This standard-library companion enumerates small finite Grassmannian fibres
as graph charts and computes every orbit.  The theorem itself is analytic;
the finite controls attack its dimension, depth, invariant-plane, and p=2
boundaries independently.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "a446d26482fdcc5ac913ca882efea17fc2cd0603c27223532884d87d6b4cc51a"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def eye(size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def mmul(left, right, modulus: int):
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    require(len(left[0]) == middle, ("matrix shape", rows, middle, columns))
    return tuple(
        tuple(
            sum(left[i][h] * right[h][j] for h in range(middle)) % modulus
            for j in range(columns)
        )
        for i in range(rows)
    )


def mpow(matrix, exponent: int, modulus: int):
    answer = eye(len(matrix))
    base = matrix
    while exponent:
        if exponent & 1:
            answer = mmul(answer, base, modulus)
        base = mmul(base, base, modulus)
        exponent //= 2
    return answer


def minv(matrix, modulus: int):
    size = len(matrix)
    work = [
        [entry % modulus for entry in row]
        + [int(i == j) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size)
             if gcd(work[row][column], modulus) == 1),
            None,
        )
        require(pivot is not None, ("nonunit pivot", matrix, modulus, column))
        work[column], work[pivot] = work[pivot], work[column]
        inverse = pow(work[column][column], -1, modulus)
        work[column] = [(inverse * value) % modulus for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [
                (work[row][j] - factor * work[column][j]) % modulus
                for j in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def mvec(matrix, vector, modulus: int):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector))) % modulus
        for i in range(len(matrix))
    )


def companion(polynomial, prime: int):
    """Column companion for a monic low-to-high coefficient tuple."""
    degree = len(polynomial) - 1
    require(polynomial[-1] % prime == 1, ("nonmonic", polynomial, prime))
    matrix = [[0] * degree for _ in range(degree)]
    for column in range(degree - 1):
        matrix[column + 1][column] = 1
    for row in range(degree):
        matrix[row][-1] = -polynomial[row] % prime
    return tuple(tuple(row) for row in matrix)


def polynomial_remainder(dividend, divisor, prime: int):
    result = [entry % prime for entry in dividend]
    degree = len(divisor) - 1
    for index in range(len(result) - 1, degree - 1, -1):
        factor = result[index]
        if factor:
            for j, entry in enumerate(divisor):
                result[index - degree + j] = (
                    result[index - degree + j] - factor * entry
                ) % prime
    return tuple(result[:degree])


def is_irreducible(polynomial, prime: int) -> bool:
    degree = len(polynomial) - 1
    if polynomial[-1] % prime != 1:
        return False
    for divisor_degree in range(1, degree // 2 + 1):
        for coefficients in product(range(prime), repeat=divisor_degree):
            divisor = tuple(coefficients) + (1,)
            if not any(polynomial_remainder(polynomial, divisor, prime)):
                return False
    return True


def rref(rows, prime: int):
    work = [list(row) for row in rows if any(entry % prime for entry in row)]
    if not work:
        return ()
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(work))
             if work[row][column] % prime),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column] % prime, -1, prime)
        work[pivot_row] = [(inverse * value) % prime for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row:
                continue
            factor = work[row][column] % prime
            work[row] = [
                (work[row][j] - factor * work[pivot_row][j]) % prime
                for j in range(columns)
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    nonzero = [tuple(row) for row in work if any(row)]
    nonzero.sort(key=lambda row: next(i for i, entry in enumerate(row) if entry))
    return tuple(nonzero)


def subspaces(prime: int, dimension: int, rank: int):
    vectors = tuple(
        vector for vector in product(range(prime), repeat=dimension)
        if any(vector)
    )
    spaces = set()
    for generators in product(vectors, repeat=rank):
        basis = rref(generators, prime)
        if len(basis) == rank:
            spaces.add(basis)
    return tuple(sorted(spaces))


def gaussian_binomial(prime: int, dimension: int, rank: int) -> int:
    numerator = 1
    denominator = 1
    for index in range(rank):
        numerator *= prime ** (dimension - index) - 1
        denominator *= prime ** (rank - index) - 1
    require(numerator % denominator == 0, "Gaussian binomial integrality")
    return numerator // denominator


def invariant_subspaces(matrix, prime: int, rank: int):
    spaces = subspaces(prime, len(matrix), rank)
    invariant = []
    for basis in spaces:
        if all(
            len(rref(basis + (mvec(matrix, vector, prime),), prime)) == rank
            for vector in basis
        ):
            invariant.append(basis)
    require(len(spaces) == gaussian_binomial(prime, len(matrix), rank),
            ("subspace count", prime, len(matrix), rank, len(spaces)))
    return spaces, tuple(invariant)


def block(matrix, row_start, row_stop, column_start, column_stop):
    return tuple(
        tuple(matrix[i][j] for j in range(column_start, column_stop))
        for i in range(row_start, row_stop)
    )


def graph_action(matrix, graph, rank: int, modulus: int):
    dimension = len(matrix)
    a = block(matrix, 0, rank, 0, rank)
    b = block(matrix, 0, rank, rank, dimension)
    c = block(matrix, rank, dimension, 0, rank)
    d = block(matrix, rank, dimension, rank, dimension)
    top = tuple(
        tuple((a[i][j] + mmul(b, graph, modulus)[i][j]) % modulus
              for j in range(rank))
        for i in range(rank)
    )
    dx = mmul(d, graph, modulus)
    bottom = tuple(
        tuple((c[i][j] + dx[i][j]) % modulus for j in range(rank))
        for i in range(dimension - rank)
    )
    return mmul(bottom, minv(top, modulus), modulus)


def orbit_partition(mapping):
    unseen = set(mapping)
    cycles = []
    while unseen:
        start = min(unseen)
        cycle = []
        point = start
        while point not in cycle:
            require(point in mapping, ("map escaped fibre", point))
            cycle.append(point)
            point = mapping[point]
        require(point == start, ("tail before cycle", start, point))
        for item in cycle:
            unseen.discard(item)
        cycles.append(tuple(cycle))
    return tuple(sorted(cycles))


def fibre_orbits(prime: int, polynomial, rank: int, carry_depth: int,
                 level: int):
    e = companion(polynomial, prime)
    dimension = len(e)
    modulus = prime ** level
    scale = prime ** carry_depth
    u = tuple(
        tuple((int(i == j) + scale * e[i][j]) % modulus
              for j in range(dimension))
        for i in range(dimension)
    )
    residue_range = range(prime ** (level - carry_depth))
    graphs = tuple(
        tuple(
            tuple(scale * values[i * rank + j] for j in range(rank))
            for i in range(dimension - rank)
        )
        for values in product(residue_range, repeat=rank * (dimension - rank))
    )
    mapping = {
        graph: graph_action(u, graph, rank, modulus)
        for graph in graphs
    }
    require(len(set(mapping.values())) == len(graphs), ("not permutation", prime, polynomial))
    cycles = orbit_partition(mapping)
    expected_length = prime ** (level - carry_depth)
    expected_count = prime ** (
        (rank * (dimension - rank) - 1) * (level - carry_depth)
    )
    require({len(cycle) for cycle in cycles} == {expected_length},
            ("orbit length", prime, polynomial, rank, carry_depth, level,
             tuple(sorted({len(cycle) for cycle in cycles}))))
    require(len(cycles) == expected_count,
            ("orbit count", prime, polynomial, rank, carry_depth, level, len(cycles)))
    return (
        prime, dimension, rank, carry_depth, level,
        len(graphs), expected_length, len(cycles),
    )


def invariant_plane_hostile():
    prime = 3
    carry_depth = 1
    level = 2
    modulus = prime ** level
    e = ((1, 0), (0, 2))
    u = tuple(
        tuple((int(i == j) + prime * e[i][j]) % modulus for j in range(2))
        for i in range(2)
    )
    graphs = tuple(((prime * value,),) for value in range(prime))
    mapping = {graph: graph_action(u, graph, 1, modulus) for graph in graphs}
    require(all(mapping[graph] == graph for graph in graphs), ("hostile moved", mapping))
    return (prime, carry_depth, level, len(graphs), tuple(mapping.values()))


def prime_two_hostile():
    prime = 2
    polynomial = (1, 1, 1)
    e = companion(polynomial, prime)
    modulus = 8
    u = tuple(
        tuple((int(i == j) + 2 * e[i][j]) % modulus for j in range(2))
        for i in range(2)
    )
    square = mpow(u, 2, modulus)
    require(square in (
        ((5, 0), (0, 5)),
        ((1, 0), (0, 1)),
    ), ("unexpected p2 square", square))
    graphs = tuple(((2 * value,),) for value in range(4))
    mapping = {graph: graph_action(u, graph, 1, modulus) for graph in graphs}
    cycles = orbit_partition(mapping)
    require(tuple(sorted(len(cycle) for cycle in cycles)) == (2, 2),
            ("p2 hostile cycles", cycles))
    return (e, u, square, tuple(sorted(len(cycle) for cycle in cycles)))


def first_carry_checks():
    checks = []
    matrices = (
        (3, companion((1, 0, 1), 3)),
        (5, companion((2, 0, 1), 5)),
    )
    for prime, e in matrices:
        for carry_depth in (1, 2, 3):
            for extra_depth in (1, 2, 3):
                modulus = prime ** (carry_depth + extra_depth + 1)
                u = tuple(
                    tuple((int(i == j) + prime ** carry_depth * e[i][j]) % modulus
                          for j in range(2))
                    for i in range(2)
                )
                power = prime ** extra_depth
                powered = mpow(u, power, modulus)
                scale = prime ** (carry_depth + extra_depth)
                carry = tuple(
                    tuple(((powered[i][j] - int(i == j)) // scale) % prime
                          for j in range(2))
                    for i in range(2)
                )
                require(carry == e, ("first carry", prime, carry_depth, extra_depth, carry, e))
                checks.append((prime, carry_depth, extra_depth, carry))
    return tuple(checks)


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert node found")

    polynomials = (
        (3, (1, 0, 1)),
        (3, (1, 2, 0, 1)),
        (3, (1, 0, 1, 1, 1)),
        (5, (2, 0, 1)),
    )
    irreducible_rows = []
    for prime, polynomial in polynomials:
        require(is_irreducible(polynomial, prime), ("reducible control", prime, polynomial))
        e = companion(polynomial, prime)
        for rank in range(1, len(e)):
            spaces, invariant = invariant_subspaces(e, prime, rank)
            require(not invariant, ("irreducible invariant subspace", prime, polynomial, rank))
            irreducible_rows.append((prime, polynomial, rank, len(spaces), len(invariant)))

    orbit_cases = (
        (3, (1, 0, 1), 1, 1, 2),
        (3, (1, 0, 1), 1, 1, 4),
        (5, (2, 0, 1), 1, 2, 4),
        (3, (1, 2, 0, 1), 1, 1, 3),
        (3, (1, 2, 0, 1), 2, 1, 3),
        (3, (1, 0, 1, 1, 1), 2, 1, 2),
        (3, (1, 0, 1, 1, 1), 2, 1, 3),
    )
    orbit_rows = tuple(fibre_orbits(*case) for case in orbit_cases)
    carry_rows = first_carry_checks()
    fixed_hostile = invariant_plane_hostile()
    two_hostile = prime_two_hostile()

    semantic_payload = (
        tuple(irreducible_rows),
        orbit_rows,
        carry_rows,
        fixed_hostile,
        two_hostile,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("near-identity Grassmannian Hensel orbit controls")
    print(f"irreducible_subspace_rows={len(irreducible_rows)}")
    print(f"irreducible_rows={tuple(irreducible_rows)}")
    print(f"orbit_cases={len(orbit_rows)}")
    for row in orbit_rows:
        print(f"orbit_row={row}")
    print(f"first_carry_checks={len(carry_rows)}")
    print(f"invariant_plane_hostile={fixed_hostile}")
    print(f"prime_two_hostile={two_hostile}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
