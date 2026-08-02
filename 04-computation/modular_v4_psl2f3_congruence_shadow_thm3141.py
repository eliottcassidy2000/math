#!/usr/bin/env python3
"""Exact four-point modular congruence-shadow controls for THM-3141."""

from __future__ import annotations

import hashlib
from collections import deque
from itertools import product


POINTS = ("o", "x", "y", "z")
S = (1, 0, 3, 2)
R = (0, 2, 3, 1)
IDENTITY = tuple(range(4))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[index]] for index in range(4))


def power(permutation, exponent):
    answer = IDENTITY
    for _ in range(exponent):
        answer = compose(permutation, answer)
    return answer


def matrix_multiply(left, right, modulus=None):
    answer = (
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
        return answer
    return tuple(tuple(entry % modulus for entry in row) for row in answer)


def matrix_power(matrix, exponent, modulus=None):
    answer = ((1, 0), (0, 1))
    for _ in range(exponent):
        answer = matrix_multiply(matrix, answer, modulus)
    return answer


def projective_normalize(vector, modulus=3):
    a, b = (entry % modulus for entry in vector)
    require((a, b) != (0, 0), vector)
    if a:
        inverse = pow(a, -1, modulus)
    else:
        inverse = pow(b, -1, modulus)
    return (a * inverse % modulus, b * inverse % modulus)


def matrix_action(matrix, vector, modulus=3):
    return projective_normalize(
        (
            matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
            matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
        ),
        modulus,
    )


def generated_group(generators):
    seen = {IDENTITY}
    queue = deque([IDENTITY])
    while queue:
        current = queue.popleft()
        for generator in generators:
            nxt = compose(generator, current)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return tuple(sorted(seen))


def reduced_words(length):
    """Normal forms in C2*C3, as alternating S and R/R2 syllables."""
    if length == 0:
        return ((),)
    words = []
    for first_kind in ("S", "R"):
        choices = []
        for index in range(length):
            kind = first_kind if index % 2 == 0 else ("R" if first_kind == "S" else "S")
            choices.append(("S",) if kind == "S" else ("R", "r"))
        words.extend(product(*choices))
    return tuple(words)


def evaluate_word(word):
    """Evaluate ``word`` in chronological application order (leftmost first)."""
    current = IDENTITY
    for letter in word:
        generator = S if letter == "S" else (R if letter == "R" else power(R, 2))
        current = compose(generator, current)
    return current


def gate(vector):
    a, b = vector
    return max(13 * abs(b), abs(a - 12 * b))


def main():
    require(power(S, 2) == IDENTITY, "S order")
    require(power(R, 3) == IDENTITY and power(R, 1) != IDENTITY, "R order")
    SR = compose(S, R)
    require(power(SR, 3) == IDENTITY and power(SR, 1) != IDENTITY, "extra relation")
    group = generated_group((S, R))
    require(len(group) == 12, len(group))
    require(all(sum(permutation[index] == index for index in range(4)) != 2 for permutation in group), "A4 cycle census")

    first_relations = []
    for length in range(1, 7):
        identities = tuple(word for word in reduced_words(length) if evaluate_word(word) == IDENTITY)
        if identities:
            first_relations.append((length, identities))
            break
    require(first_relations and first_relations[0][0] == 6, first_relations)
    relation_words = first_relations[0][1]
    require(("S", "R", "S", "R", "S", "R") in relation_words, relation_words)

    matrix_s = ((0, -1), (1, 0))
    matrix_r = ((0, -1), (1, 1))
    matrix_sr = matrix_multiply(matrix_s, matrix_r)
    require(matrix_sr == ((-1, -1), (0, -1)), matrix_sr)
    require(matrix_power(matrix_sr, 3) == ((-1, -3), (0, -1)), "integer T3")
    require(matrix_power(matrix_sr, 3, 3) == ((2, 0), (0, 2)), "mod3 identity")

    point_vectors = ((1, 1), (1, 2), (1, 0), (0, 1))
    vector_index = {projective_normalize(vector): index for index, vector in enumerate(point_vectors)}
    induced_s = tuple(vector_index[matrix_action(matrix_s, vector)] for vector in point_vectors)
    induced_r = tuple(vector_index[matrix_action(matrix_r, vector)] for vector in point_vectors)
    require(induced_s == S and induced_r == R, (induced_s, induced_r))

    farey_neighbor_checks = 0
    for a, b, c, d in product(range(-6, 7), repeat=4):
        if a * d - b * c not in (-1, 1):
            continue
        farey_neighbor_checks += 1
        require(
            projective_normalize((a, b)) != projective_normalize((c, d)),
            ((a, b), (c, d)),
        )
    require(farey_neighbor_checks > 0, "farey controls")
    converse = ((1, 0), (1, 4))
    require(abs(converse[0][0] * converse[1][1] - converse[0][1] * converse[1][0]) == 4, converse)
    require(projective_normalize(converse[0]) != projective_normalize(converse[1]), converse)

    vector = (-100, -7)
    s_vector = (7, -100)
    r_vector = (7, -107)
    require(91 * gate(vector) == 8281 <= sum(entry * entry for entry in vector), "base gate")
    require(91 * gate(s_vector) == 118300 > sum(entry * entry for entry in s_vector), "S hostile")
    require(91 * gate(r_vector) == 126581 > sum(entry * entry for entry in r_vector), "R hostile")

    record = (
        POINTS,
        S,
        R,
        group,
        first_relations,
        matrix_s,
        matrix_r,
        matrix_sr,
        point_vectors,
        farey_neighbor_checks,
        converse,
        vector,
        s_vector,
        r_vector,
    )
    digest = hashlib.sha256(repr(record).encode()).hexdigest()
    print("THM3141 quartic V4 modular congruence shadow")
    print(f"orders=S:2;R:3;SR:3;generated_group:{len(group)}")
    print(f"first_extra_relation_length=6;words:{relation_words}")
    print("identification=V4_semidirect_C3:A4:PSL2_F3;full_affine:S4:PGL2_F3")
    print("matrix_relation=SR_is_projective_T;integer_(SR)^3_nontrivial;mod3_(SR)^3_identity")
    print(f"farey=neighbor_checks:{farey_neighbor_checks};converse_det:4")
    print("gate_hostile=d:(-100,-7):8281<=10049;S:118300>10049;R:126581>11498")
    print("ternary_boundary=C3_rotates_the_labelled_input_triple;the_endpoint_jet_readout_is_rotation_invariant_and_forgets_the_root")
    print("missing_sidecar=Gamma3_Farey_lift_plus_rooted_C3_block_plus_common_physical_atom")
    print(f"record_sha256={digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
