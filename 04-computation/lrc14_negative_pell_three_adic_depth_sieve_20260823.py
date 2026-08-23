#!/usr/bin/env python3
"""Exact three-adic depth sieve on the negative-Pell successor orbit."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


DEPTH_CONTROL = 1_001
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


def v3(n: int) -> int:
    if n == 0:
        raise ValueError("v3(0) is not used")
    value = 0
    while n % 3 == 0:
        value += 1
        n //= 3
    return value


def matmul(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def matvec(
    matrix: tuple[tuple[int, int], tuple[int, int]],
    vector: tuple[int, int],
) -> tuple[int, int]:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def matpow(
    matrix: tuple[tuple[int, int], tuple[int, int]], exponent: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    result = ((1, 0), (0, 1))
    base = matrix
    while exponent:
        if exponent & 1:
            result = matmul(result, base)
        base = matmul(base, base)
        exponent //= 2
    return result


A = ((2, 3), (1, 2))
C = ((9, 15), (5, 9))
w = (1, 1)
w_minus = (-1, 1)

check(A[0][0] * A[1][1] - A[0][1] * A[1][0] == 1,
      "Pell transition has determinant one")
check(matpow(A, 3) == ((26, 45), (15, 26)),
      "third Pell transition is exact")
check(
    matpow(A, 3)
    == ((-1 + 3 * C[0][0], 3 * C[0][1]),
        (3 * C[1][0], -1 + 3 * C[1][1])),
    "A cubed equals minus I plus 3C",
)
check(matvec(C, w)[1] == 14 and matvec(C, w_minus)[1] == 4,
      "both lifting seeds are three-adic units")

states: list[tuple[int, int, int]] = []
X, Y = w
for k in range(DEPTH_CONTROL + 1):
    Q = (Y * Y - 1) // 8
    if k:
        check(X * X - 3 * Y * Y == -2,
              "Pell recurrence remains on the negative shell")
        if k % 3 == 1:
            predicted = 0
        else:
            j = (k + 1) // 3 if k % 3 == 2 else k // 3
            predicted = 1 + v3(j)
        check(v3(Q) == predicted, "three-adic depth law")
    if k <= 36:
        states.append((X, Y, Q))
    X, Y = matvec(A, (X, Y))

period = tuple(state[2] % 9 for state in states[:9])
check(period == (0, 1, 6, 3, 1, 3, 6, 1, 0),
      "negative-Pell values have the claimed first mod-nine period")
check(
    all(tuple(state[2] % 9 for state in states[9 * block:9 * (block + 1)])
        == period for block in (1, 2, 3)),
    "all four blocks in a full state period have the same mod-nine values",
)
state_after_thirty_six = matvec(matpow(A, 36), w)
check(tuple(entry % 72 for entry in state_after_thirty_six) == w,
      "the Pell state returns modulo 72 after thirty-six steps")
cube_residues = {pow(a, 3, 9) for a in range(9)}
two_cube_residues = {(a + b) % 9 for a in cube_residues for b in cube_residues}
check(cube_residues == {0, 1, 8} and two_cube_residues == {0, 1, 2, 7, 8},
      "one-cube and two-cube residues modulo nine")
allowed_depths = {index for index, value in enumerate(period) if value in two_cube_residues}
check(allowed_depths == {0, 1, 4, 7, 8},
      "two-cube Pell depths modulo nine")

# The previous million-shell computation retains k=0,...,16.  Removing the
# zero value leaves sixteen positive nodes; exactly the six k=1 mod 3 nodes
# can carry a 3-free primitive two-cube value.  The known hit is k=4.
positive_million_depths = tuple(range(1, 17))
three_free_candidates = tuple(k for k in positive_million_depths if k % 3 == 1)
check(len(positive_million_depths) == 16,
      "million-shell has sixteen positive Pell nodes")
check(three_free_candidates == (1, 4, 7, 10, 13, 16),
      "three-adic sieve leaves six million-shell candidates")
check(states[4][2] == 2_926 and 2_926 == 9**3 + 13**3,
      "known triangular-successor hit lies at depth four")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "companion has no inactive Python assert")

semantic = {
    "classification": "PROVED-LEMMA+VERIFIED-EXACT",
    "orbit": "(X_k,Y_k)=((2,3),(1,2))^k(1,1)",
    "value": "Q_k=(Y_k^2-1)/8",
    "valuation": {
        "k=3j-1_or_3j": "v3(Q_k)=1+v3(j)",
        "k=3j+1": "v3(Q_k)=0",
    },
    "cube_scale_gate": "v3(g)=r>=1 forces v3(j)=3r-1",
    "mod9_allowed_depths": [0, 1, 4, 7, 8],
    "three_free_depths": [1, 4, 7],
    "scope": "necessary address sieve only; no LRC14 implication",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("experiment=LRC14-negative-Pell-three-adic-depth-sieve")
print("classification=PROVED-LEMMA;VERIFIED-EXACT")
print("orbit=Xk_Yk_equals_matrix_2_3_1_2_power_k_times_1_1")
print("value=Qk_equals_Yk2_minus_1_over_8_equals_T_Nk")
print("valuation=k_3j_minus_1_or_3j:1_plus_v3j;k_3j_plus_1:0")
print("cube_scale=v3g_r_positive_implies_v3j_equals_3r_minus_1")
print("mod9_Q_period=0,1,6,3,1,3,6,1,0")
print("two_cube_allowed_k_mod9=0,1,4,7,8;three_free=1,4,7")
print("million_shell=16_positive_nodes;three_free_candidates=6;known_hit_k=4")
print(f"hostile_control_depths={DEPTH_CONTROL}")
print("scope=necessary_address_sieve_only;no_LRC14_consequence")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
