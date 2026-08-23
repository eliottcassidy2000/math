#!/usr/bin/env python3
"""Hostile audit of the normalized-unit mod-9 negative-Pell lemma.

This checker deliberately avoids the proposed checker's SymPy calculation and
its order-by-order binomial implementation.  It uses a block-lifting identity
for B^(3^s), direct modular Pell powers for the two index offsets, a separate
full-integer recurrence for the valuation law and first excluded depths, and
explicit finite-level cube witnesses over Z/3^n Z.

There are no Python ``assert`` statements, so optimized replay executes every
gate.
"""

from __future__ import annotations

import ast
import hashlib
import json
from math import gcd
from pathlib import Path
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]


def matmul(left: Matrix, right: Matrix, modulus: int | None = None) -> Matrix:
    result = (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )
    if modulus is None:
        return result
    return tuple(
        tuple(entry % modulus for entry in row) for row in result
    )  # type: ignore[return-value]


def matpow(matrix: Matrix, exponent: int, modulus: int | None = None) -> Matrix:
    result: Matrix = ((1, 0), (0, 1))
    base = matrix
    while exponent:
        if exponent & 1:
            result = matmul(result, base, modulus)
        base = matmul(base, base, modulus)
        exponent //= 2
    return result


def matvec(matrix: Matrix, vector: Vector, modulus: int | None = None) -> Vector:
    result = (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )
    if modulus is None:
        return result
    return result[0] % modulus, result[1] % modulus


def matrix_mod(matrix: Matrix, modulus: int) -> Matrix:
    return tuple(
        tuple(entry % modulus for entry in row) for row in matrix
    )  # type: ignore[return-value]


def v3(number: int) -> int:
    require(number != 0, "v3 is never applied to zero")
    valuation = 0
    while number % 3 == 0:
        number //= 3
        valuation += 1
    return valuation


I: Matrix = ((1, 0), (0, 1))
A: Matrix = ((2, 3), (1, 2))
A_INV: Matrix = ((2, -3), (-1, 2))
C: Matrix = ((9, 15), (5, 9))
B: Matrix = ((-26, -45), (-15, -26))
W_PLUS: Vector = (1, 1)
W_MINUS: Vector = (-1, 1)


# Index and sign bookkeeping.  A^3=-B, and B commutes with A^{-1}; hence
# A^(3j)W_PLUS=(-1)^j B^j W_PLUS and
# A^(3j-1)W_PLUS=(-1)^j B^j W_MINUS.  The sign disappears in Y^2.
require(matmul(A, A_INV) == I and matmul(A_INV, A) == I,
        "displayed inverse of A is exact")
require(matpow(A, 3) == ((26, 45), (15, 26)),
        "third Pell transition")
require(
    matpow(A, 3)
    == tuple(
        tuple(-I[row][column] + 3 * C[row][column] for column in range(2))
        for row in range(2)
    ),
    "A cubed equals -I+3C",
)
require(B == tuple(tuple(-entry for entry in row) for row in matpow(A, 3)),
        "B equals -A cubed")
require(matvec(A_INV, W_PLUS) == W_MINUS,
        "minus seed is the one-step inverse seed")
require(matmul(A_INV, B) == matmul(B, A_INV),
        "the inverse shift commutes with the three-step block")


# Independent block-lift proof certificate.  Put
# E_s=(B^(3^s)-I)/3^(s+1).  Exact computation at s=2 gives
# E_2 = [[0,3],[4,0]] mod 9.  Cubing I+3^(s+1)E_s shows
# E_(s+1)=E_s mod 9 for every s>=2.  Raising the block to a unit h is
# linear modulo 3^(s+3), so its second-coordinate action is 4h or 5h.
base_power = matpow(B, 9)
base_scale = 3**3
require(
    all((base_power[row][column] - I[row][column]) % base_scale == 0
        for row in range(2) for column in range(2)),
    "B^9-I is divisible by 27 entrywise",
)
e_two = tuple(
    tuple((base_power[row][column] - I[row][column]) // base_scale
          for column in range(2))
    for row in range(2)
)  # type: ignore[assignment]
expected_e_mod9: Matrix = ((0, 3), (4, 0))
require(matrix_mod(e_two, 9) == expected_e_mod9,
        "base block-lift matrix modulo nine")
require(matvec(expected_e_mod9, W_PLUS, 9)[1] == 4,
        "plus seed block coefficient")
require(matvec(expected_e_mod9, W_MINUS, 9)[1] == 5,
        "minus seed block coefficient")

unit_representatives = (1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17)
block_cases = 0
for s in range(2, 13):
    scale = 3 ** (s + 1)
    modulus = 9 * scale
    block = matpow(B, 3**s, modulus)
    lifted = tuple(
        tuple(((block[row][column] - I[row][column]) % modulus) // scale
              for column in range(2))
        for row in range(2)
    )  # type: ignore[assignment]
    require(lifted == expected_e_mod9,
            "all tested block lifts stabilize modulo nine")
    for h in unit_representatives:
        require(gcd(h, 3) == 1, "h representative is a three-adic unit")
        powered_block = matpow(block, h, modulus)
        for seed, coefficient in ((W_PLUS, 4), (W_MINUS, 5)):
            coordinate = matvec(powered_block, seed, modulus)[1]
            normalized = ((coordinate - 1) % modulus) // scale
            require(normalized == coefficient * h % 9,
                    "block lift is linear in the unit part")
            block_cases += 1


def q_mod(depth: int, modulus: int) -> int:
    """Compute Q_depth=(Y_depth^2-1)/8 modulo ``modulus`` exactly."""

    state_modulus = 8 * modulus
    y_value = matvec(matpow(A, depth, state_modulus), W_PLUS,
                     state_modulus)[1]
    require(y_value % 2 == 1, "Pell Y remains odd")
    require((y_value * y_value - 1) % 8 == 0,
            "Pell triangular quotient is integral")
    return ((y_value * y_value - 1) // 8) % modulus


# Directly challenge both index shifts without using B or its seed moments.
offset_cases = 0
for s in range(2, 12):
    scale = 3 ** (s + 1)
    modulus = 9 * scale
    for h in unit_representatives:
        j = 3**s * h
        for depth, expected in ((3 * j, h), (3 * j - 1, -h)):
            q_residue = q_mod(depth, modulus)
            require(q_residue % scale == 0,
                    "direct Pell power has the predicted scale")
            require(q_residue // scale % 9 == expected % 9,
                    "direct Pell power confirms the sign and index shift")
            offset_cases += 1


# A disjoint full-integer recurrence checks the inherited valuation law and
# locates the first exclusions, rather than assuming either from the block
# congruence.
x_value, y_value = W_PLUS
excluded_depths: list[int] = []
valuation_depth = 1_200
for depth in range(valuation_depth + 1):
    require(x_value * x_value - 3 * y_value * y_value == -2,
            "full recurrence stays on X^2-3Y^2=-2")
    require((y_value * y_value - 1) % 8 == 0,
            "full recurrence gives an integral Q")
    q_value = (y_value * y_value - 1) // 8
    if depth:
        if depth % 3 == 1:
            predicted_valuation = 0
            j = None
        elif depth % 3 == 0:
            j = depth // 3
            predicted_valuation = 1 + v3(j)
        else:
            j = (depth + 1) // 3
            predicted_valuation = 1 + v3(j)
        require(v3(q_value) == predicted_valuation,
                "full recurrence confirms the valuation formula")
        if (j is not None and predicted_valuation >= 3
                and predicted_valuation % 3 == 0):
            primitive_unit = q_value // (3**predicted_valuation) % 9
            if primitive_unit in {4, 5}:
                excluded_depths.append(depth)
    x_value, y_value = matvec(A, (x_value, y_value))

require(tuple(excluded_depths[:4]) == (107, 108, 134, 135),
        "the first four scale-compatible exclusions are exact")


# Exact local completion.  At each finite precision, enumerate the unit-cube
# image, then construct a two-cube witness for every permitted unit target.
# Necessity is already forced modulo nine by cube residues {0,+1,-1}.
allowed_classes = {1, 2, 7, 8}
require({pow(value, 3, 9) for value in range(9)} == {0, 1, 8},
        "all cube residues modulo nine")
require(
    {
        (left + right) % 9
        for left in {0, 1, 8}
        for right in {0, 1, 8}
        if (left + right) % 3
    }
    == allowed_classes,
    "the four necessary unit two-cube classes modulo nine",
)

local_targets = 0
for precision in range(2, 10):
    modulus = 3**precision
    cube_to_root: dict[int, int] = {}
    for root in range(modulus):
        cube_to_root.setdefault(pow(root, 3, modulus), root)
    actual_unit_cubes = {
        cube for cube in cube_to_root if gcd(cube, 3) == 1
    }
    expected_unit_cubes = {
        unit for unit in range(modulus)
        if gcd(unit, 3) == 1 and unit % 9 in {1, 8}
    }
    require(actual_unit_cubes == expected_unit_cubes,
            "unit cubes are exactly the +/-1 mod-nine classes")

    for target in range(modulus):
        if gcd(target, 3) != 1:
            continue
        permitted = target % 9 in allowed_classes
        if not permitted:
            # No pair can exist: reduction modulo nine would contradict the
            # already exhausted cube-residue calculation above.
            local_targets += 1
            continue
        if target % 9 in {1, 8}:
            second_root = 3
        elif target % 9 == 2:
            second_root = 1
        else:
            second_root = -1
        first_cube = (target - second_root**3) % modulus
        require(first_cube in cube_to_root,
                "constructive permitted target has a first cube root")
        first_root = cube_to_root[first_cube]
        require((first_root**3 + second_root**3 - target) % modulus == 0,
                "constructive two-cube witness is exact")
        local_targets += 1


# The scale corollary only uses a sign-stable property.  A unit factor in g
# contributes a cube congruent to +/-1 mod 9, and the two branch signs also
# preserve the same four allowed classes.
require({(-value) % 9 for value in allowed_classes} == allowed_classes,
        "allowed classes are stable under the branch sign")
require(
    {(unit_cube * value) % 9 for unit_cube in {1, 8}
     for value in allowed_classes} == allowed_classes,
    "allowed classes are stable under the unit part of the common cube scale",
)
require(
    {unit for unit in range(1, 9) if unit % 3 and unit not in allowed_classes}
    == {4, 5},
    "exactly two of the six unit classes are excluded",
)


source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "checker has no optimization-sensitive assert")

semantic = {
    "classification": "PROVED-LEMMA+VERIFIED-EXACT+HOSTILE-AUDIT",
    "orbit": "A=((2,3),(1,2)); Q_k=(Y_k^2-1)/8",
    "block_lift": "E_s=(B^(3^s)-I)/3^(s+1)=((0,3),(4,0)) mod 9 for s>=2",
    "normalized_congruence": {
        "Q_3j": "+h mod 9",
        "Q_3j_minus_1": "-h mod 9",
        "conditions": "s=v3(j)>=2; h=j/3^s",
    },
    "first_excluded_depths": [107, 108, 134, 135],
    "local_iff": "unit z in Z_3 is a sum of two cubes iff z mod 9 is in {1,2,7,8}",
    "stopping_scope": (
        "final only for unconstrained two-summand Z_3 existence after exact "
        "common-cube-scale removal; prescribed summand data and all nonlocal "
        "coordinates remain eligible to prune"
    ),
}
semantic_blob = json.dumps(semantic, sort_keys=True,
                           separators=(",", ":")).encode()

print("experiment=negative-pell-normalized-unit-mod9-hostile-audit")
print("classification=PROVED-LEMMA;VERIFIED-EXACT;HOSTILE-AUDIT")
print("indexing=A_3j_uses_plus_seed;A_3j_minus_1_uses_inverse_minus_seed;global_sign_squares_out")
print("block_lift=E_s_mod9_equals_0_3__4_0_for_all_s_at_least_2")
print("normalized_congruence=Q_3j_over_3_splus1:+h;Q_3jminus1_over_3_splus1:-h_mod9")
print("scale_gate=s_equals_3r_minus_1;allowed_h_mod9=1,2,7,8;excluded=4,5")
print("first_excluded_depths=107,108,134,135")
print("local_iff=Z3_unit_two_cube_sum_iff_residue_mod9_in_1,2,7,8")
print("stopping_scope=unconstrained_Z3_existence_only_after_exact_scale_removal")
print("nonlocal_scope=prescribed_summand_data,positivity,global_integrality,other_places,and_runner_sidecars_not_closed")
print(f"block_cases={block_cases}")
print(f"offset_cases={offset_cases}")
print(f"valuation_depth={valuation_depth}")
print(f"local_targets={local_targets}")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
