#!/usr/bin/env python3
"""Independent Gcontent replay; imports no code from the information lane."""

from functools import reduce
from itertools import product
import json
from math import gcd


MATRICES = {
    "s": [[0, 1, 1], [0, 0, 1], [1, 0, 1]],
    "c": [[1, 0, 4], [1, 0, 3], [3, 1, 6]],
    "q0": [[1, 0, 8], [1, 0, 7], [7, 1, 18]],
    "q1": [[3, 0, 6], [3, 0, 5], [9, 3, 14]],
    "q2": [[4, 1, 7], [3, 1, 6], [6, 1, 14]],
    "q4": [[3, 0, 12], [3, 0, 9], [5, 1, 14]],
    "q5": [[4, 1, 9], [4, 0, 7], [7, 2, 11]],
}
WORDS = {
    "W5": "q2 c c q2 s s s s q2",
    "W9": "q1 s s q1 s s q1",
    "W21": "q5 s q0 c q0 q5 q0 c s s s q5",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def multiply(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]


def determinant(a):
    return sum(a[0][j] * (a[1][(j + 1) % 3] * a[2][(j + 2) % 3]
                         - a[1][(j + 2) % 3] * a[2][(j + 1) % 3])
               for j in range(3))


def content(row):
    return reduce(gcd, map(abs, row), 0)


def main():
    determinants = {letter: determinant(matrix) for letter, matrix in MATRICES.items()}
    require(determinants == {"s": 1, "c": 1, "q0": 1, "q1": 9,
                             "q2": 5, "q4": 9, "q5": 21}, "letter determinants")
    report = {"status": "FINITE-EXACT independent products passed",
              "letter_determinants": determinants, "words": {}}
    for name, word in WORDS.items():
        matrix = [[int(i == j) for j in range(3)] for i in range(3)]
        for letter in word.split():
            matrix = multiply(matrix, MATRICES[letter])
        scalar = int(name[1:])
        require(content(value for row in matrix for value in row) == scalar,
                "exact scalar matrix content")
        normalized = [[value // scalar for value in row] for row in matrix]
        require(determinant(normalized) == 1, "unimodular normalization")
        require(determinant(matrix) == scalar**3, "raw determinant")
        for row in product(range(-3, 4), repeat=3):
            image = [sum(row[k] * matrix[k][j] for k in range(3)) for j in range(3)]
            require(content(image) == scalar * content(row), "universal-law finite controls")
        report["words"][name] = {
            "word": word,
            "product": matrix,
            "content": scalar,
            "determinant": determinant(matrix),
            "normalized": normalized,
            "normalized_determinant": determinant(normalized),
            "vertex_cost": sum(1 if letter == "s" else 3 if letter == "c" else 4
                               for letter in word.split()),
            "finite_rows_checked": 343,
        }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
