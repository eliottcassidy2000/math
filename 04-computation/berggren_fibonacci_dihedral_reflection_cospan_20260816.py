#!/usr/bin/env python3
"""Exact dihedral cospan for Berggren, Farey/Fibonacci, and C4/C6 hostiles.

The two matrix factorizations are exact.  The Keller comparison belongs in
the companion reflection because its old-L section permutations are currently
verified numerically rather than proved at all levels.
"""

from __future__ import annotations

from hashlib import sha256
import json
from math import gcd


Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]

EXPECTED_SEMANTIC_SHA256 = "2ed27bf257b908ab33b5ad5909f46d2b171b0b1b14e2f31cf2205ea6e2825114"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def identity(size: int) -> Matrix:
    return tuple(
        tuple(1 if row == column else 0 for column in range(size))
        for row in range(size)
    )


def matmul(left: Matrix, right: Matrix) -> Matrix:
    require(len(left[0]) == len(right), ("matrix product", left, right))
    return tuple(
        tuple(
            sum(left[row][pivot] * right[pivot][column]
                for pivot in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    require(len(matrix[0]) == len(vector), ("matrix vector", matrix, vector))
    return tuple(
        sum(matrix[row][column] * vector[column]
            for column in range(len(vector)))
        for row in range(len(matrix))
    )


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def add(left: Matrix, right: Matrix, right_scale: int = 1) -> Matrix:
    require((len(left), len(left[0])) == (len(right), len(right[0])),
            ("matrix sum", left, right))
    return tuple(
        tuple(left[row][column] + right_scale * right[row][column]
              for column in range(len(left[0])))
        for row in range(len(left))
    )


def scale(scalar: int, matrix: Matrix) -> Matrix:
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def power(matrix: Matrix, exponent: int) -> Matrix:
    require(exponent >= 0 and len(matrix) == len(matrix[0]),
            ("matrix power", matrix, exponent))
    answer = identity(len(matrix))
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = matmul(answer, base)
        base = matmul(base, base)
        remaining //= 2
    return answer


def determinant_two(matrix: Matrix) -> int:
    require((len(matrix), len(matrix[0])) == (2, 2), matrix)
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def fibonacci(index: int) -> int:
    a, b = 0, 1
    for _ in range(index):
        a, b = b, a + b
    return a


def preserves_form(matrix: Matrix, form: Matrix) -> bool:
    return matmul(matmul(transpose(matrix), form), matrix) == form


def tournament_completion_census(size: int) -> dict[str, int]:
    cycle_arcs = {(vertex, (vertex + 1) % size) for vertex in range(size)}
    cycle_pairs = {tuple(sorted(arc)) for arc in cycle_arcs}
    all_pairs = tuple(
        (left, right)
        for left in range(size)
        for right in range(left + 1, size)
    )
    missing_pairs = tuple(pair for pair in all_pairs if pair not in cycle_pairs)
    rotation_invariant = 0
    for mask in range(1 << len(missing_pairs)):
        arcs = set(cycle_arcs)
        for index, (left, right) in enumerate(missing_pairs):
            arcs.add((right, left) if (mask >> index) & 1 else (left, right))
        rotated = {((left + 1) % size, (right + 1) % size) for left, right in arcs}
        if rotated == arcs:
            rotation_invariant += 1
    return {
        "vertices": size,
        "cycle_pairs": len(cycle_pairs),
        "missing_pairs": len(missing_pairs),
        "tournament_completions": 1 << len(missing_pairs),
        "rotation_invariant_completions": rotation_invariant,
    }


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def main() -> None:
    lorentz = ((1, 0, 0), (0, 1, 0), (0, 0, -1))
    berggren_u = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
    berggren_a = ((-1, 0, 0), (0, 1, 0), (0, 0, 1))
    berggren_b = matmul(berggren_u, berggren_a)
    eye_three = identity(3)
    require(matmul(berggren_a, berggren_a) == eye_three, "Berggren A involution")
    require(matmul(berggren_b, berggren_b) == eye_three, "Berggren B involution")
    require(matmul(berggren_b, berggren_a) == berggren_u, "Berggren BA=U")
    require(preserves_form(berggren_u, lorentz), "U Lorentz")
    require(preserves_form(berggren_a, lorentz), "A Lorentz")
    require(preserves_form(berggren_b, lorentz), "B Lorentz")
    require(matmul(matmul(matmul(berggren_a, berggren_u), berggren_a), berggren_u)
            == eye_three, "A reverses U")
    nilpotent = add(berggren_u, eye_three, -1)
    require(power(nilpotent, 2) != ((0, 0, 0),) * 3, "N squared nonzero")
    require(power(nilpotent, 3) == ((0, 0, 0),) * 3, "N cubed zero")

    parameter_t = ((0, 1), (-1, 2))
    parameter_a = ((0, 1), (1, 0))
    parameter_b = matmul(parameter_t, parameter_a)
    eye_two = identity(2)
    require(matmul(parameter_a, parameter_a) == eye_two, "Farey A involution")
    require(matmul(parameter_b, parameter_b) == eye_two, "Farey B involution")
    require(matmul(parameter_b, parameter_a) == parameter_t, "Farey BA=T")
    require(determinant_two(parameter_t) == 1, "Farey determinant")
    parameter_nilpotent = add(parameter_t, eye_two, -1)
    require(power(parameter_nilpotent, 2) == ((0, 0), (0, 0)),
            "Farey parabolic nilpotent")

    berggren_rows = []
    for depth in range(9):
        triple = matvec(power(berggren_u, depth), (3, 4, 5))
        expected = (2 * depth + 3,
                    2 * depth * depth + 6 * depth + 4,
                    2 * depth * depth + 6 * depth + 5)
        require(triple == expected, ("Berggren U spine", depth, triple, expected))
        parameters = matvec(power(parameter_t, depth), (1, 2))
        require(parameters == (depth + 1, depth + 2),
                ("Farey parabolic ray", depth, parameters))
        require(gcd(*parameters) == 1, parameters)
        q_value = 2 * triple[2] + 1
        require(q_value == triple[0] ** 2 + 2, (depth, triple, q_value))
        berggren_rows.append({
            "depth": depth,
            "fraction": f"{parameters[0]}/{parameters[1]}",
            "q": q_value,
            "triple": triple,
        })

    for left, right in zip(berggren_rows, berggren_rows[1:]):
        p, q = (int(value) for value in left["fraction"].split("/"))
        r, s = (int(value) for value in right["fraction"].split("/"))
        require(abs(p * s - r * q) == 1, ("Berggren Farey edge", left, right))

    fibonacci_q = ((1, 1), (1, 0))
    fibonacci_m = matmul(fibonacci_q, fibonacci_q)
    fibonacci_a = ((0, -1), (1, 0))
    fibonacci_b = scale(-1, matmul(fibonacci_m, fibonacci_a))
    minus_eye_two = scale(-1, eye_two)
    require(matmul(fibonacci_a, fibonacci_a) == minus_eye_two,
            "Fibonacci A projective involution")
    require(matmul(fibonacci_b, fibonacci_b) == minus_eye_two,
            "Fibonacci B projective involution")
    require(matmul(fibonacci_b, fibonacci_a) == fibonacci_m,
            "Fibonacci BA=M")
    require(determinant_two(fibonacci_m) == 1, "Fibonacci determinant")
    require(fibonacci_m[0][0] + fibonacci_m[1][1] == 3,
            "Fibonacci hyperbolic trace")
    fibonacci_m_inverse = ((1, -1), (-1, 2))
    require(matmul(fibonacci_m, fibonacci_m_inverse) == eye_two,
            "Fibonacci inverse")
    require(matmul(matmul(fibonacci_a, fibonacci_m), scale(-1, fibonacci_a))
            == fibonacci_m_inverse, "Fibonacci projective reversal")

    fibonacci_rows = []
    for depth in range(1, 10):
        vector = matvec(power(fibonacci_m, depth), (1, 0))
        expected = (fibonacci(2 * depth + 1), fibonacci(2 * depth))
        require(vector == expected, ("Fibonacci ray", depth, vector, expected))
        require(gcd(*vector) == 1, vector)
        fibonacci_rows.append({
            "depth": depth,
            "fraction": f"{vector[1]}/{vector[0]}",
            "vector": vector,
        })
    for left, right in zip(fibonacci_rows, fibonacci_rows[1:]):
        p, q = (int(value) for value in left["fraction"].split("/"))
        r, s = (int(value) for value in right["fraction"].split("/"))
        require(abs(p * s - r * q) == 1, ("Fibonacci Farey edge", left, right))

    tournament_rows = tuple(tournament_completion_census(size) for size in (4, 6))
    require(tournament_rows == (
        {"vertices": 4, "cycle_pairs": 4, "missing_pairs": 2,
         "tournament_completions": 4, "rotation_invariant_completions": 0},
        {"vertices": 6, "cycle_pairs": 6, "missing_pairs": 9,
         "tournament_completions": 512, "rotation_invariant_completions": 0},
    ), tournament_rows)

    ledger = {
        "berggren": {
            "A": berggren_a,
            "B": berggren_b,
            "U": berggren_u,
            "nilpotency_index": 3,
            "rows": berggren_rows,
        },
        "farey_parameter": {
            "A": parameter_a,
            "B": parameter_b,
            "T": parameter_t,
        },
        "fibonacci": {
            "A_projective": fibonacci_a,
            "B_projective": fibonacci_b,
            "M": fibonacci_m,
            "rows": fibonacci_rows,
            "trace": 3,
        },
        "tournaments": tournament_rows,
    }
    semantic_sha256 = digest(ledger)
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    print("== Berggren/Fibonacci dihedral reflection cospan ==")
    print(f"berggren_factorization=A^2=B^2=I;B*A=U;A*U*A=U^-1;"
          f"B={berggren_b}")
    print("berggren_dynamics=parabolic;(U-I)^3=0!=(U-I)^2;"
          f"q_values={tuple(row['q'] for row in berggren_rows)}")
    print("berggren_farey_ray="
          + str(tuple(row["fraction"] for row in berggren_rows)))
    print("fibonacci_factorization=A^2=B^2=-I_in_GL2;B*A=Q^2_in_GL2;"
          "A,B_are_involutions_in_PGL2")
    print("fibonacci_dynamics=hyperbolic;trace(Q^2)=3;"
          + "farey_ray=" + str(tuple(row["fraction"] for row in fibonacci_rows)))
    for row in tournament_rows:
        print("cycle_tournament_hostile="
              + ";".join(f"{key}={value}" for key, value in row.items()))
    print(f"semantic_sha256={semantic_sha256}")
    print("status=VERIFIED-EXACT MATRIX/FAREY/TOURNAMENT COSPAN;"
          "no Keller identification or LRC current")


if __name__ == "__main__":
    main()
