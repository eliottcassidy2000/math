#!/usr/bin/env python3
"""Independent integer audit of the strengthened THM-3874 presentation."""

from __future__ import annotations

import hashlib
import math
import sys


sys.stdout.reconfigure(encoding="utf-8", newline="\n")
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def matvec(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def cartan_a(size: int) -> list[list[int]]:
    return [
        [2 if i == j else -1 if abs(i - j) == 1 else 0 for j in range(size)]
        for i in range(size)
    ]


def cartan_e7() -> list[list[int]]:
    # Central vertex 0 with arms of lengths 1, 2, 3, matching the canonical
    # generator order used in the strengthened checker.
    result = [[0] * 7 for _ in range(7)]
    for i in range(7):
        result[i][i] = 2
    for i, j in ((0, 1), (0, 2), (2, 3), (0, 4), (4, 5), (5, 6)):
        result[i][j] = result[j][i] = -1
    return result


def determinant_bareiss(matrix: list[list[int]]) -> int:
    size = len(matrix)
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (i for i in range(pivot_index + 1, size) if work[i][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for i in range(pivot_index + 1, size):
            for j in range(pivot_index + 1, size):
                numerator = (
                    work[i][j] * pivot
                    - work[i][pivot_index] * work[pivot_index][j]
                )
                gate(numerator % previous == 0, "Bareiss exact division")
                work[i][j] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


generators = (
    "O",
    "F",
    "a1",
    "a2",
    "a3",
    "a4",
    "a5",
    "c1",
    "c2",
    "b1",
    "b2",
    "b3",
    "m1",
    "e1",
    "e2",
    "e3",
    "e4",
    "e5",
    "e6",
    "e7",
    "T",
)
index = {name: i for i, name in enumerate(generators)}
gate(len(generators) == 21 and len(index) == 21, "generator universe")

twice_a5 = [1, 2, 3, 2, 1]
twice_a3 = [1, 2, 1]
twice_e7 = [6, 3, 4, 2, 5, 4, 3]
gate(matvec(cartan_a(5), twice_a5) == [0, 0, 2, 0, 0], "A5 weight")
gate(matvec(cartan_a(3), twice_a3) == [0, 2, 0], "A3 weight")
gate(matvec(cartan_e7(), twice_e7) == [0, 0, 0, 0, 0, 0, 2], "E7 weight")

# Twice T-O-2F+omega_3(A5)+omega_2(A3)+omega(E7)=0.
gluing = [0] * 21
gluing[index["O"]] = -2
gluing[index["F"]] = -4
gluing[index["T"]] = 2
for name, coefficient in zip(("a1", "a2", "a3", "a4", "a5"), twice_a5):
    gluing[index[name]] = coefficient
for name, coefficient in zip(("b1", "b2", "b3"), twice_a3):
    gluing[index[name]] = coefficient
for name, coefficient in zip(
    ("e1", "e2", "e3", "e4", "e5", "e6", "e7"), twice_e7
):
    gluing[index[name]] = coefficient

killed = {
    "O",
    "F",
    "T",
    "a1",
    "a2",
    "a4",
    "a5",
    "c1",
    "c2",
    "b1",
    "b3",
    "m1",
    "e1",
    "e2",
    "e3",
    "e4",
    "e5",
    "e6",
    "e7",
}
survivors = [name for name in generators if name not in killed]
gate(len(killed) == 19 and survivors == ["a3", "b2"], "boundary survivors")
gate([gluing[index[name]] for name in survivors] == [3, 2], "reduced relation")
gate(math.gcd(3, 2) == 1, "primitive reduced relation")

rows = [gluing]
for name in generators:
    if name in killed:
        row = [0] * 21
        row[index[name]] = 1
        rows.append(row)
gate(len(rows) == 20 and all(len(row) == 21 for row in rows), "20 by 21 presentation")

maximal_minors = []
for omitted_column in range(21):
    square = [row[:omitted_column] + row[omitted_column + 1 :] for row in rows]
    maximal_minors.append(determinant_bareiss(square))
minor_gcd = 0
for value in maximal_minors:
    minor_gcd = math.gcd(minor_gcd, abs(value))
gate(any(maximal_minors), "presentation rank 20")
gate(minor_gcd == 1, "torsion-free rank-one cokernel")

kernel = [0] * 21
kernel[index["a3"]] = 2
kernel[index["b2"]] = -3
gate(
    all(sum(a * b for a, b in zip(row, kernel)) == 0 for row in rows),
    "primitive kernel",
)
gate(math.gcd(*(abs(value) for value in kernel)) == 1, "kernel saturation")

semantic = (
    "twice torsion gluing has O,F,T,A5,A3,E7 packets",
    "boundary and ADE rows kill 19 named generators",
    "only a3=d and b2=e survive",
    "reduced relation is 3d+2e",
    "all 20-by-20 minors have gcd one and cokernel is Z",
)
semantic_hash = hashlib.sha256(repr(semantic).encode()).hexdigest()
print("THM3874_FULL_PRESENTATION", "rows=20;columns=21;rank=20;torsion=0;free_rank=1")
print("THM3874_SURVIVORS", "d=a3;e=b2;relation=3d+2e")
print("THM3874_MAXIMAL_MINOR_GCD", minor_gcd)
print("SEMANTIC_SHA256", semantic_hash)
print("GATES", GATES)
print("RESULT", "PASS")
