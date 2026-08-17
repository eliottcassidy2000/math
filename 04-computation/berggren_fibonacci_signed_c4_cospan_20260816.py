#!/usr/bin/env python3
"""Exact signed-C4 cospan between Berggren walls and Fibonacci mod 3.

The two carriers are not identified.  This script verifies the precise
common quotient: a pointed directed four-cycle equipped with the nonzero
class in H^1(C4; F2).  It also freezes the tournament-completion hostile.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


D = 145
EXPECTED_SEMANTIC_SHA256 = "93a9f158e462686b4393163e46fcb211f607282eaaf3b8fd20a9f214156a7f9e"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def neg(value):
    return -value[0], -value[1]


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    return (
        left[0] * right[0] + D * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def inv(value):
    denominator = value[0] ** 2 - D * value[1] ** 2
    require(denominator != 0, ("zero denominator", value))
    return value[0] / denominator, -value[1] / denominator


def div(left, right):
    return mul(left, inv(right))


def rational(value):
    return Fraction(value), Fraction(0)


def sign(value):
    a_value, b_value = value
    if b_value == 0:
        return (a_value > 0) - (a_value < 0)
    if a_value == 0:
        return (b_value > 0) - (b_value < 0)
    if (a_value > 0) == (b_value > 0):
        return 1 if a_value > 0 else -1
    norm = a_value * a_value - D * b_value * b_value
    return ((norm > 0) - (norm < 0)) if a_value > 0 else -((norm > 0) - (norm < 0))


def less(left, right):
    return sign(sub(left, right)) < 0


def inverse_branch(value):
    if less(value, rational(Fraction(1, 3))):
        return "C", div(value, sub(rational(1), mul(rational(2), value)))
    if less(value, rational(Fraction(1, 2))):
        return "B", sub(inv(value), rational(2))
    return "A", sub(rational(2), inv(value))


def wall_word(value):
    letters = []
    current = value
    for _ in range(4):
        letter, current = inverse_branch(current)
        letters.append(letter)
    require(current == value, ("wall did not return", value, letters, current))
    return tuple(letters)


def normalize_projective(vector, prime=3):
    vector = tuple(value % prime for value in vector)
    pivot = next(value for value in vector if value)
    inverse = pow(pivot, -1, prime)
    return tuple(value * inverse % prime for value in vector)


def matrix_vector(matrix, vector, prime=3):
    return normalize_projective(
        tuple(sum(row[index] * vector[index] for index in range(2)) % prime for row in matrix),
        prime,
    )


def determinant_sign(left, right):
    determinant = (left[0] * right[1] - left[1] * right[0]) % 3
    require(determinant in (1, 2), ("coincident projective lines", left, right))
    return 1 if determinant == 1 else -1


def cohomology_bits(signs):
    bits = tuple(0 if value == 1 else 1 for value in signs)
    require(sum(bits) % 2 == 1, ("trivial C4 class", signs, bits))
    return bits


def gauge(signs, vertex_signs):
    return tuple(
        signs[index] * vertex_signs[index] * vertex_signs[(index + 1) % 4]
        for index in range(4)
    )


def tournament_completion_audit():
    cycle_arcs = {(0, 1), (1, 2), (2, 3), (3, 0)}
    rows = []
    for diagonal_02 in (0, 1):
        for diagonal_13 in (0, 1):
            arcs = set(cycle_arcs)
            arcs.add((0, 2) if diagonal_02 else (2, 0))
            arcs.add((1, 3) if diagonal_13 else (3, 1))
            invariant = all(
                ((left + 1) % 4, (right + 1) % 4) in arcs
                for left, right in arcs
            )
            rows.append((diagonal_02, diagonal_13, invariant))
    require(sum(row[2] for row in rows) == 0, rows)
    return tuple(rows)


def semantic_hash(payload):
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def main():
    source = Path(__file__)
    syntax = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(syntax)
        ),
        "floating literal",
    )

    sqrt_d = Fraction(0), Fraction(1)
    alpha = div(sub(sqrt_d, rational(9)), rational(8))
    beta = div(sub(sqrt_d, rational(8)), rational(9))
    alpha_word = wall_word(alpha)
    beta_word = wall_word(beta)
    require(alpha_word == tuple("BABB"), alpha_word)
    require(beta_word == tuple("BCBB"), beta_word)
    branch_sign = {"A": 1, "B": -1, "C": 1}
    wall_signs = tuple(branch_sign[letter] for letter in alpha_word)
    require(tuple(branch_sign[letter] for letter in beta_word) == wall_signs, beta_word)
    require(wall_signs == (-1, 1, -1, -1), wall_signs)

    states = ((1, 0), (0, 1), (1, 1), (1, 2))
    fibonacci = ((0, 1), (1, 1))
    images = tuple(matrix_vector(fibonacci, state) for state in states)
    require(images == states[1:] + states[:1], images)
    fibonacci_signs = tuple(
        determinant_sign(states[index], states[(index + 1) % 4])
        for index in range(4)
    )
    require(fibonacci_signs == (1, -1, 1, 1), fibonacci_signs)

    alternating_section = (1, -1, 1, -1)
    switched_signs = gauge(fibonacci_signs, alternating_section)
    require(switched_signs == wall_signs, (switched_signs, wall_signs))
    wall_class = cohomology_bits(wall_signs)
    fibonacci_class = cohomology_bits(fibonacci_signs)

    gauge_orbit = {
        gauge(
            fibonacci_signs,
            tuple(-1 if (mask >> index) & 1 else 1 for index in range(4)),
        )
        for mask in range(16)
    }
    require(len(gauge_orbit) == 8, len(gauge_orbit))
    require(wall_signs in gauge_orbit, "wall class absent from Fibonacci gauge orbit")
    completion_rows = tournament_completion_audit()

    payload = {
        "wall_words": (alpha_word, beta_word),
        "wall_signs": wall_signs,
        "fibonacci_states": states,
        "fibonacci_images": images,
        "fibonacci_signs": fibonacci_signs,
        "alternating_section": alternating_section,
        "switched_signs": switched_signs,
        "classes": (wall_class, fibonacci_class),
        "gauge_orbit_size": len(gauge_orbit),
        "tournament_completions": completion_rows,
    }
    semantic = semantic_hash(payload)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic))

    print("BERGGREN/FIBONACCI SIGNED-C4 COSPAN EXACT AUDIT")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=VERIFIED_EXACT_ELEMENTARY_CANDIDATE;unnumbered")
    print(f"wall_words=alpha:{''.join(alpha_word)};beta:{''.join(beta_word)}")
    print(f"wall_edge_signs={wall_signs};wall_H1_bits={wall_class}")
    print(f"fibonacci_P1F3_cycle={states};determinant_edge_signs={fibonacci_signs}")
    print(f"alternating_vertex_section={alternating_section};switched_signs={switched_signs}")
    print("common_quotient=pointed_directed_C4_plus_nonzero_H1(C4,F2)_class")
    print(f"odd_sign_gauge_orbit_size={len(gauge_orbit)}")
    print("six_pairs=four_successor_pairs_plus_two_unobserved_antipodal_pairs")
    print("C4_rotation_invariant_tournament_completions=0_of_4")
    print("lost=wall_values,branch_letters,affine_CDF_offsets,Fibonacci_height,Farey_frame,q15_owner")
    print("scope=no_density_transfer_no_LRC_current_no_JC_flux_no_canonical_tournament")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
