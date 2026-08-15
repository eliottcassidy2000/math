#!/usr/bin/env python3
"""Exact four-state bridge: Farey/Fibonacci mod 3 and q=15 sheet modes.

The probe keeps three carriers typed:

* P^1(F_3), the four-state reduction of primitive rational slopes;
* the q=15 unit sign quotient under multiplication by two and inversion;
* the ternary Berggren branch action reduced to a four-state automaton.

It also constructs the 24-state mod-six frame-line quotient and audits the
noncanonical tournament sections.  No physical LRC owner is transported.
Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter, deque
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED = (
    (
        "THM-3333-script",
        ROOT / "04-computation/gaussian_farey_pythagorean_triangular_thm3333.py",
        "ddeb800a881fdf3576d3d03d269ced2b486fabe8bc48594416a7c343531441df",
    ),
    (
        "THM-3333-output",
        ROOT / "05-knowledge/results/gaussian_farey_pythagorean_triangular_thm3333.out",
        "7e8eaabb73fc3f86b3099d3b43aa3c3da925fb24a3081e7611a9faddccebaf79",
    ),
    (
        "THM-3339-script",
        ROOT / "04-computation/fibonacci_berggren_three_ray_owner_thm3339.py",
        "094fc254dcb7965791e59247a98f60a337725b40d169db6f3862c1da5943149b",
    ),
    (
        "THM-3339-output",
        ROOT / "05-knowledge/results/fibonacci_berggren_three_ray_owner_thm3339.out",
        "88c5f44971df12bbd61f84925f408a253450095ee2de23b5a617efde1e4ecdfa",
    ),
    (
        "THM-3379-script",
        ROOT / "04-computation/fibonacci_ray_t4_mod3_projection_thm3379.py",
        "803065a6ce419554b04403b13537889547da53bd1e30aa35edc5e5caca549ae8",
    ),
    (
        "THM-3379-output",
        ROOT / "05-knowledge/results/fibonacci_ray_t4_mod3_projection_thm3379.out",
        "d3cf7be59fe988171b73184b9cf26974b26bdfbcd1e50bbb3d40075a7d666778",
    ),
    (
        "q15-script",
        ROOT / "04-computation/lrc15_first_effective_triphase_mode_probe_20260814.py",
        "751aff937e422a6d09e7e89677ded2f40ce7c9b72c12191078d691cbc4316c33",
    ),
    (
        "q15-output",
        ROOT / "05-knowledge/results/lrc15_first_effective_triphase_mode_probe_20260814.out",
        "258ced868f662355b890ad0cc486f6375d96c04eb0f226a723974d0b79d95b8d",
    ),
)

P3_STATES = ((1, 0), (0, 1), (1, 1), (1, 2))
G = ((0, 1), (1, 1))
A = ((0, 1), (-1, 2))
B = ((0, 1), (1, 2))
C = ((1, 0), (2, 1))
REFLECTION = ((1, 1), (0, 2))
Q15_STATES = (1, 2, 4, 7)
EXPECTED_SEMANTIC_DIGEST = "f111020b7aebaa3b4108800fe981923e77befc586c8521818cd8aa8600c87baf"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def matrix_apply(matrix, vector, modulus):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(2)) % modulus
        for row in range(2)
    )


def matrix_multiply(left, right, modulus):
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(2)) % modulus
            for column in range(2)
        )
        for row in range(2)
    )


def projective_state(vector):
    x, y = (coordinate % 3 for coordinate in vector)
    require((x, y) != (0, 0), ("zero projective vector", vector))
    if x:
        inverse = 1 if x == 1 else 2
        normalized = (1, y * inverse % 3)
    else:
        normalized = (0, 1)
    return P3_STATES.index(normalized)


def projective_permutation(matrix):
    return tuple(projective_state(matrix_apply(matrix, state, 3)) for state in P3_STATES)


def compose_permutations(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def generated_permutation_group(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    queue = deque((identity,))
    while queue:
        current = queue.popleft()
        for generator in generators:
            child = compose_permutations(generator, current)
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return tuple(sorted(seen))


def canonical_q15(residue):
    residue %= 15
    require(gcd(residue, 15) == 1, ("q15 unit", residue))
    return min(residue, (-residue) % 15)


def determinant(left, right):
    return left[0] * right[1] - left[1] * right[0]


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant_n(matrix):
    size = len(matrix)
    total = 0
    for permutation in permutations(range(size)):
        product_value = permutation_sign(permutation)
        for row, column in enumerate(permutation):
            product_value *= matrix[row][column]
        total += product_value
    return total


def matrix_n_multiply(left, right):
    size = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(size))
            for column in range(size)
        )
        for row in range(size)
    )


def matrix_n_add(*terms):
    size = len(terms[0][1])
    return tuple(
        tuple(
            sum(coefficient * matrix[row][column] for coefficient, matrix in terms)
            for column in range(size)
        )
        for row in range(size)
    )


def matrix_n_power(matrix, exponent):
    size = len(matrix)
    result = tuple(tuple(int(row == column) for column in range(size)) for row in range(size))
    base = matrix
    while exponent:
        if exponent % 2:
            result = matrix_n_multiply(result, base)
        base = matrix_n_multiply(base, base)
        exponent //= 2
    return result


def row_times_matrix(row, matrix):
    return tuple(
        sum(row[middle] * matrix[middle][column] for middle in range(len(row)))
        for column in range(len(row))
    )


def fibonacci_projective_states(limit):
    values = [0, 1]
    while len(values) <= limit + 1:
        values.append(values[-1] + values[-2])
    return tuple(projective_state((values[index], values[index + 1])) for index in range(limit + 1))


def tournament_sections():
    sections = []
    for flips_tail in product((0, 1), repeat=3):
        flips = (0,) + flips_tail
        vectors = tuple(
            tuple(((-1 if flips[index] else 1) * coordinate) % 3 for coordinate in state)
            for index, state in enumerate(P3_STATES)
        )
        arcs = []
        for left, right in combinations(range(4), 2):
            residue = determinant(vectors[left], vectors[right]) % 3
            require(residue in (1, 2), ("projective determinant", left, right))
            arcs.append((left, right) if residue == 1 else (right, left))
        scores = tuple(sum(source == vertex for source, _ in arcs) for vertex in range(4))
        cyclic = 0
        for triple in combinations(range(4), 3):
            local_scores = tuple(
                sum(source == vertex for source, target in arcs if source in triple and target in triple)
                for vertex in triple
            )
            cyclic += tuple(sorted(local_scores)) == (1, 1, 1)
        sections.append((flips, tuple(arcs), tuple(sorted(scores)), cyclic))
    require(len({entry[1] for entry in sections}) == 8, "section tournaments")
    return tuple(sections)


def all_tournaments_invariant_under(rotation):
    pairs = tuple(combinations(range(4), 2))
    invariant = []
    for bits in product((0, 1), repeat=len(pairs)):
        arcs = {
            pair if bit == 0 else (pair[1], pair[0])
            for pair, bit in zip(pairs, bits)
        }
        moved = {(rotation[source], rotation[target]) for source, target in arcs}
        if moved == arcs:
            invariant.append(tuple(sorted(arcs)))
    return tuple(invariant)


def sl2_mod_six_frame_line_packet():
    matrices = tuple(
        ((a, b), (c, d))
        for a, b, c, d in product(range(6), repeat=4)
        if (a * d - b * c) % 6 == 1
    )
    require(len(matrices) == 144, ("SL2 mod6 count", len(matrices)))

    upper = ((1, 1), (0, 1))
    lower = ((1, 0), (1, 1))
    identity = ((1, 0), (0, 1))
    generated = {identity}
    queue = deque((identity,))
    while queue:
        current = queue.popleft()
        for generator in (upper, lower):
            child = matrix_multiply(generator, current, 6)
            if child not in generated:
                generated.add(child)
                queue.append(child)
    require(generated == set(matrices), ("integral generator image", len(generated)))

    frame_states = tuple(
        matrix
        for matrix in product((0, 1), repeat=4)
        if (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % 2 == 1
    )
    require(len(frame_states) == 6, frame_states)
    frame_index = {frame: index for index, frame in enumerate(frame_states)}

    fibres = Counter()
    for matrix in matrices:
        frame = (
            matrix[0][0] % 2,
            matrix[0][1] % 2,
            matrix[1][0] % 2,
            matrix[1][1] % 2,
        )
        line = projective_state((matrix[0][0], matrix[1][0]))
        fibres[(frame_index[frame], line)] += 1
    require(len(fibres) == 24, ("frame-line states", len(fibres)))
    require(set(fibres.values()) == {6}, fibres)
    return len(matrices), tuple(sorted(fibres.items()))


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    permutations_by_name = {
        "G": projective_permutation(G),
        "A": projective_permutation(A),
        "B": projective_permutation(B),
        "C": projective_permutation(C),
        "S": projective_permutation(REFLECTION),
    }
    require(permutations_by_name["G"] == (1, 2, 3, 0), permutations_by_name)
    require(permutations_by_name["A"] == (1, 3, 2, 0), permutations_by_name)
    require(permutations_by_name["B"] == (1, 3, 0, 2), permutations_by_name)
    require(permutations_by_name["C"] == (3, 1, 0, 2), permutations_by_name)
    require(permutations_by_name["S"] == (0, 3, 2, 1), permutations_by_name)

    berggren_group = generated_permutation_group(
        (permutations_by_name["A"], permutations_by_name["B"], permutations_by_name["C"])
    )
    dihedral_group = generated_permutation_group(
        (permutations_by_name["G"], permutations_by_name["S"])
    )
    require(len(berggren_group) == 24, len(berggren_group))
    require(len(dihedral_group) == 8, len(dihedral_group))
    require(
        compose_permutations(
            permutations_by_name["S"],
            compose_permutations(permutations_by_name["G"], permutations_by_name["S"]),
        )
        == (3, 0, 1, 2),
        "S G S=G^-1",
    )

    q15_rotation = tuple(canonical_q15(2 * residue) for residue in Q15_STATES)
    q15_inversion = tuple(canonical_q15(pow(residue, -1, 15)) for residue in Q15_STATES)
    require(q15_rotation == (2, 4, 7, 1), q15_rotation)
    require(q15_inversion == (1, 7, 4, 2), q15_inversion)
    q15_index_rotation = tuple(Q15_STATES.index(value) for value in q15_rotation)
    q15_index_inversion = tuple(Q15_STATES.index(value) for value in q15_inversion)
    require(q15_index_rotation == permutations_by_name["G"], "cycle conjugacy")
    require(q15_index_inversion == permutations_by_name["S"], "reflection conjugacy")

    unit_blocks = tuple(
        tuple(sorted({pow(speed, -1, 15), (-pow(speed, -1, 15)) % 15}))
        for speed in Q15_STATES
    )
    require(unit_blocks == ((1, 14), (7, 8), (4, 11), (2, 13)), unit_blocks)
    block_indices = tuple(
        Q15_STATES.index(canonical_q15(block[0])) for block in unit_blocks
    )
    require(block_indices == q15_index_inversion, (block_indices, q15_index_inversion))

    fibonacci_states = fibonacci_projective_states(80)
    require(fibonacci_states[:9] == (1, 2, 3, 0, 1, 2, 3, 0, 1), fibonacci_states[:9])
    require(all(fibonacci_states[index + 4] == fibonacci_states[index] for index in range(77)), "projective period four")

    primitive = tuple(
        (numerator, denominator)
        for denominator in range(1, 31)
        for numerator in range(0, denominator + 1)
        if gcd(numerator, denominator) == 1
    )
    farey_edges_checked = 0
    for left, right in combinations(primitive, 2):
        if abs(determinant(left, right)) == 1:
            require(projective_state(left) != projective_state(right), (left, right))
            farey_edges_checked += 1
    require(projective_state((0, 1)) != projective_state((2, 5)), "hostile states")
    require(abs(determinant((0, 1), (2, 5))) == 2, "nonconverse hostile")

    transition = tuple(
        tuple(
            sum(permutations_by_name[name][source_state] == target_state for name in ("A", "B", "C"))
            for target_state in range(4)
        )
        for source_state in range(4)
    )
    require(transition == ((0, 2, 0, 1), (0, 1, 0, 2), (2, 0, 1, 0), (1, 0, 2, 0)), transition)
    identity4 = tuple(tuple(int(row == column) for column in range(4)) for row in range(4))
    cayley_hamilton = matrix_n_add(
        (1, matrix_n_power(transition, 4)),
        (-2, matrix_n_power(transition, 3)),
        (-6, transition),
        (-9, identity4),
    )
    require(cayley_hamilton == tuple((0, 0, 0, 0) for _ in range(4)), cayley_hamilton)
    require(determinant_n(transition) == -9, determinant_n(transition))

    level_counts = []
    count = (1, 0, 0, 0)
    for depth in range(21):
        require(sum(count) == 3**depth, (depth, count))
        level_counts.append(count)
        count = row_times_matrix(count, transition)
    require(tuple(level_counts[:9]) == (
        (1, 0, 0, 0),
        (0, 2, 0, 1),
        (1, 2, 2, 4),
        (8, 4, 10, 5),
        (25, 20, 20, 16),
        (56, 70, 52, 65),
        (169, 182, 182, 196),
        (560, 520, 574, 533),
        (1681, 1640, 1640, 1600),
    ), level_counts[:9])

    sections = tournament_sections()
    section_profile = Counter((entry[2], entry[3]) for entry in sections)
    require(section_profile == Counter({((0, 2, 2, 2), 1): 4, ((1, 1, 1, 3), 1): 4}), section_profile)
    require(not all_tournaments_invariant_under(permutations_by_name["G"]), "C4 invariant tournament")

    sl2_count, frame_line_fibres = sl2_mod_six_frame_line_packet()
    require(sl2_count == 144, sl2_count)
    require(len(frame_line_fibres) == 24, len(frame_line_fibres))

    canonical_edge_harmonic_mass = sum((Fraction(1, speed) for speed in (1, 2, 3, 4, 5, 7)), Fraction())
    require(canonical_edge_harmonic_mass == Fraction(1019, 420), canonical_edge_harmonic_mass)

    semantic = ExactDigest()
    semantic.add(("permutations", permutations_by_name))
    semantic.add(("berggren_group", berggren_group))
    semantic.add(("dihedral_group", dihedral_group))
    semantic.add(("q15", q15_rotation, q15_inversion, unit_blocks, block_indices))
    semantic.add(("fibonacci", fibonacci_states))
    semantic.add(("farey_edges_checked", farey_edges_checked))
    semantic.add(("transition", transition, tuple(level_counts)))
    semantic.add(("tournaments", sections, tuple(sorted(section_profile.items()))))
    semantic.add(("frame_line", sl2_count, frame_line_fibres))
    semantic.add(("harmonic_mass", canonical_edge_harmonic_mass))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("FIBONACCI-FAREY MOD3 / Q15 FOUR-STATE EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=VERIFIED-EXACT finite quotients plus elementary recurrence/density proof candidate;unnumbered_and_not_canon")
    print(f"P1F3_states={P3_STATES};actions={tuple(permutations_by_name.items())}")
    print(f"q15_unit_sign_cycle={Q15_STATES};times2={q15_rotation};inversion={q15_inversion};blocked_pairs={unit_blocks}")
    print("pointed_conjugacy=times2_on_q15_classes_equals_Fibonacci_G_on_P1F3;owner_to_block_is_inversion_reflection")
    print(f"groups=<A,B,C>_on_P1F3=S4_order_{len(berggren_group)};<G,S>=D4_order_{len(dihedral_group)}")
    print(f"Farey_edges_checked={farey_edges_checked};one_way_only_hostile=(0/1,2/5)_determinant_2")
    print(f"Berggren_transition_matrix={transition};charpoly=(x-3)(x+1)(x^2+3);level_counts_0_to_8={tuple(level_counts[:9])}")
    print("tree_state_density=1/4_each;heap_harmonic_coefficient=1/4_each;state_subsets_have_coefficient=cardinality/4")
    print(f"tournament_sections=8;profile={tuple(sorted(section_profile.items()))};C4_invariant_tournaments=0")
    print(f"mod6_frame_line_states={len(frame_line_fibres)};SL2Z6_elements={sl2_count};fibre_size=6")
    print(f"canonical_q15_edge_harmonic_mass={canonical_edge_harmonic_mass.numerator}/{canonical_edge_harmonic_mass.denominator}")
    print("typing=four_state_cycle_and_24_state_frame_line_quotient;not_a_physical_owner_not_a_canonical_tournament")
    print("scope=no_Farey_converse_no_full_tree_word_recovery_no_LRC_current_or_JC_flux")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
