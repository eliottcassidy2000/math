#!/usr/bin/env python3
"""Exact q=3 phase-triangle cover clutter for THM-3388.

For a three-sheet fibre, every speed coprime to three blocks at most one
sheet.  A full cover is therefore witnessed by three distinct speeds.  The
analytic criterion records the three signed centre gaps as integral affine
cochains and asks for their weighted circulation to vanish.

The criterion is checked against an independent rational event sweep.  The
literal [14] cover clutter, its independence polynomial, the complete q=3
body atlas, and a ternary dilation-lattice harmonic recurrence are frozen.
Runtime gates use RuntimeError, so python -O retains every decision.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED = (
    (
        "THM-3387",
        ROOT / "01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md",
        "c540255a185efb54c67035d69f3fbd94f4c1ad3e30c4e31738a8800e81198613",
    ),
    (
        "THM-3387-script",
        ROOT / "04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py",
        "9b0b46874a569d674b937b37cf74a8985fca2b77e3e480a75fb4924ea602f25a",
    ),
    (
        "THM-3387-output",
        ROOT / "05-knowledge/results/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.out",
        "b4d9ce439bab4501bfd5e2cf13eb0b0e3685b7364f30e43b7d5ca9138d25cb5c",
    ),
)

LITERAL_VERTICES = (1, 2, 4, 5, 7, 8, 10, 11, 13, 14)
EXPECTED_LITERAL_EDGE_COUNT = 48
EXPECTED_PAIRWISE_TRIANGLES = 82
EXPECTED_INDEPENDENCE_PROFILE = (1, 10, 45, 72, 38, 6, 0, 0, 0, 0, 0)
EXPECTED_MAXIMAL_FIVE_SETS = (
    (1, 2, 4, 7, 14),
    (1, 2, 4, 8, 11),
    (1, 2, 4, 8, 13),
    (1, 2, 5, 10, 13),
    (2, 4, 7, 8, 13),
    (2, 4, 7, 8, 14),
)
EXPECTED_Q3_CANDIDATES = 2793
EXPECTED_Q3_GLOBAL = 585
EXPECTED_Q3_EXACT = 588
EXPECTED_CORE_RESCUES = (
    ((1, 2, 3), (8, 11, 13)),
    ((1, 2, 4), (8, 11, 13)),
    ((2, 3, 4), (8, 11, 13)),
)
EXPECTED_SEMANTIC_DIGEST = "082e97aa25d8019ba7de49c0a76333c7a3a221dd19cb4bc3e8d5b43ef9a42216"


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


def strict_danger(numerator, denominator):
    residue = numerator % denominator
    return 14 * min(residue, denominator - residue) < denominator


def transverse_danger(speed, sample, scale, sheet):
    # source time is sample/(2*scale); denominator six*scale combines k/3.
    return strict_danger(speed * (3 * sample + 2 * scale * sheet), 6 * scale)


def core_danger(clock, sample, scale):
    # Descended base time is three times the source time.
    return strict_danger(3 * clock * sample, 2 * scale)


def circle_event_samples(transverse, core=()):
    scale = 42 * lcm(*transverse, *core)
    events = {0}
    for speed in transverse:
        require(speed % 3 != 0, ("transverse type", speed))
        for sheet in range(3):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // 3
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    for clock in core:
        for tooth in range(3 * clock):
            for sign in (-1, 1):
                events.add(
                    (
                        scale * tooth // (3 * clock)
                        + sign * scale // (42 * clock)
                    )
                    % scale
                )

    ordered = tuple(sorted(events))
    midpoints = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += scale
        midpoints.append((left + right) % (2 * scale))
    endpoint_samples = tuple(2 * event for event in ordered)
    return scale, endpoint_samples + tuple(midpoints)


def full_transverse_cover(transverse, sample, scale):
    return all(
        any(transverse_danger(speed, sample, scale, sheet) for speed in transverse)
        for sheet in range(3)
    )


def full_cover_leaks(transverse, core=()):
    leaks = []
    scale, samples = circle_event_samples(transverse, core)
    for sample in samples:
        if not full_transverse_cover(transverse, sample, scale):
            continue
        if any(core_danger(clock, sample, scale) for clock in core):
            continue
        leaks.append(Q(sample, 2 * scale))
    return tuple(leaks)


def phase_gap_values(left, right):
    """Attainable oriented gaps for sheets differing by +1 mod three."""
    common = gcd(left, right)
    modulus = 3 * common
    residue = (left * right) % modulus
    bound = (3 * (left + right) - 1) // 14
    first_index = -((bound + residue) // modulus)
    last_index = (bound - residue) // modulus
    return tuple(
        residue + modulus * index
        for index in range(first_index, last_index + 1)
    )


def pairwise_feasible(left, right):
    return bool(phase_gap_values(left, right))


def phase_triangle_ordered(speeds):
    left, middle, right = speeds
    first = phase_gap_values(left, middle)
    second = phase_gap_values(middle, right)
    third = set(phase_gap_values(right, left))
    for p in first:
        for q in second:
            numerator = -(right * p + left * q)
            if numerator % middle != 0:
                continue
            r = numerator // middle
            if r in third:
                return p, q, r
    return None


def phase_triangle(speeds):
    return phase_triangle_ordered(tuple(sorted(speeds)))


def crt_pair(left_residue, left_modulus, right_residue, right_modulus):
    common = gcd(left_modulus, right_modulus)
    require((right_residue - left_residue) % common == 0, "CRT incompatibility")
    right_reduced = right_modulus // common
    if right_reduced == 1:
        return left_residue % left_modulus
    step = (
        (right_residue - left_residue)
        // common
        * pow(left_modulus // common, -1, right_reduced)
    ) % right_reduced
    return (left_residue + left_modulus * step) % (
        left_modulus * right_reduced
    )


def reconstruct_centres(speeds, gaps):
    """Glue a closed phase cochain to integer centre numerators a,b,c."""
    left, middle, right = speeds
    p, q, r = gaps
    first_difference = (p - left * middle) // 3
    second_difference = (q - middle * right) // 3

    common_left = gcd(left, middle)
    modulus_left = middle // common_left
    residue_left = 0
    if modulus_left != 1:
        residue_left = (
            -first_difference
            // common_left
            * pow(left // common_left, -1, modulus_left)
        ) % modulus_left

    common_right = gcd(right, middle)
    modulus_right = middle // common_right
    residue_right = 0
    if modulus_right != 1:
        residue_right = (
            second_difference
            // common_right
            * pow(right // common_right, -1, modulus_right)
        ) % modulus_right

    b = crt_pair(residue_left, modulus_left, residue_right, modulus_right)
    a = (first_difference + left * b) // middle
    c = (right * b - second_difference) // middle
    third_difference = (r + 2 * right * left) // 3
    require(
        (
            middle * a - left * b,
            right * b - middle * c,
            left * c - right * a,
        )
        == (first_difference, second_difference, third_difference),
        ("centre reconstruction", speeds, gaps, a, b, c),
    )
    return a, b, c


def independent(subset, edge_sets):
    chosen = set(subset)
    return not any(edge <= chosen for edge in edge_sets)


def multinomial3(depth, i, j):
    k = depth - i - j
    return factorial(depth) // (factorial(i) * factorial(j) * factorial(k))


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

    pair_checks = 0
    pair_pool = tuple(number for number in range(1, 200) if number % 3 != 0)
    for left, right in combinations(pair_pool, 2):
        predicted = 3 * (left + right) > 14 * gcd(left, right)
        require(pairwise_feasible(left, right) == predicted, ("pair threshold", left, right))
        pair_checks += 1

    # Independent exact rational event sweep through all triples below 41.
    audit_pool = tuple(number for number in range(1, 41) if number % 3 != 0)
    triple_checks = 0
    positive_checks = 0
    permutation_checks = 0
    for triple in combinations(audit_pool, 3):
        witness = phase_triangle(triple)
        exact = bool(full_cover_leaks(triple))
        require((witness is not None) == exact, ("phase criterion", triple, witness))
        if witness is not None:
            reconstruct_centres(triple, witness)
            positive_checks += 1
        if max(triple) <= 20:
            decisions = {
                phase_triangle_ordered(order) is not None
                for order in permutations(triple)
            }
            require(decisions == {exact}, ("orientation gauge", triple, decisions))
            permutation_checks += 6
        triple_checks += 1

    literal_edges = tuple(
        triple
        for triple in combinations(LITERAL_VERTICES, 3)
        if phase_triangle(triple) is not None
    )
    literal_edge_sets = tuple(frozenset(edge) for edge in literal_edges)
    pairwise_triangles = tuple(
        triple
        for triple in combinations(LITERAL_VERTICES, 3)
        if all(pairwise_feasible(left, right) for left, right in combinations(triple, 2))
    )
    require(len(literal_edges) == EXPECTED_LITERAL_EDGE_COUNT, len(literal_edges))
    require(len(pairwise_triangles) == EXPECTED_PAIRWISE_TRIANGLES, len(pairwise_triangles))
    require((1, 4, 7) in pairwise_triangles and (1, 4, 7) not in literal_edges, "pairwise hostile")
    require(phase_triangle((1, 4, 5)) == (1, -1, -1), "positive phase triangle")
    require(
        (
            phase_gap_values(1, 4),
            phase_gap_values(4, 41),
            phase_gap_values(41, 1),
            phase_triangle((1, 4, 41)),
        )
        == (
            (1,),
            (-7, -4, -1, 2, 5, 8),
            (-7, -4, -1, 2, 5, 8),
            None,
        ),
        "extended pairwise hostile",
    )
    require(phase_triangle((2, 8, 10)) == (-2, 2, 2), "dilation sign hostile")

    independence_profile = tuple(
        sum(
            independent(subset, literal_edge_sets)
            for subset in combinations(LITERAL_VERTICES, size)
        )
        for size in range(len(LITERAL_VERTICES) + 1)
    )
    maximal_five_sets = tuple(
        subset
        for subset in combinations(LITERAL_VERTICES, 5)
        if independent(subset, literal_edge_sets)
    )
    require(independence_profile == EXPECTED_INDEPENDENCE_PROFILE, independence_profile)
    require(maximal_five_sets == EXPECTED_MAXIMAL_FIVE_SETS, maximal_five_sets)

    # Complete literal q=3 body atlas by a second rational event route.
    candidates = 0
    global_rows = 0
    exact_rows = 0
    rescues = []
    core_pool = (1, 2, 3, 4)
    for core_size in range(1, 5):
        transverse_size = 6 - core_size
        for core in combinations(core_pool, core_size):
            for transverse in combinations(LITERAL_VERTICES, transverse_size):
                candidates += 1
                globally_safe = independent(transverse, literal_edge_sets)
                leaks = full_cover_leaks(transverse, core)
                pointwise_exact = not leaks
                if globally_safe:
                    require(pointwise_exact, ("global row", core, transverse))
                    global_rows += 1
                if pointwise_exact:
                    exact_rows += 1
                    if not globally_safe:
                        rescues.append((core, transverse))

    require(candidates == EXPECTED_Q3_CANDIDATES, candidates)
    require(global_rows == EXPECTED_Q3_GLOBAL, global_rows)
    require(exact_rows == EXPECTED_Q3_EXACT, exact_rows)
    require(tuple(rescues) == EXPECTED_CORE_RESCUES, rescues)
    hostile_sample = 389
    hostile_scale = 1232
    require(
        full_transverse_cover((8, 11, 13), hostile_sample, hostile_scale)
        and tuple(core_danger(clock, hostile_sample, hostile_scale) for clock in (1, 2, 3, 4))
        == (False, True, False, False),
        "missing-core nonrescue hostile",
    )

    # A genuine ternary word ancestry whose integer support abelianizes to a
    # three-dimensional exponent lattice.  Every scaled (1,4,5) remains an
    # edge because common dilation preserves the phase criterion.
    roots = (1, 4, 5)
    multipliers = (7, 11, 13)
    root_harmonic = sum((Q(1, root) for root in roots), Q(0))
    require(root_harmonic == Q(29, 20), root_harmonic)
    support_shells = []
    word_shells = []
    dilation_checks = 0
    for depth in range(21):
        support_values = []
        support_mass = Q(0)
        word_mass = Q(0)
        for i in range(depth + 1):
            for j in range(depth - i + 1):
                k = depth - i - j
                scale = 7**i * 11**j * 13**k
                require(phase_triangle(tuple(scale * root for root in roots)) is not None, ("dilation edge", depth, i, j))
                dilation_checks += 1
                multiplicity = multinomial3(depth, i, j)
                for root in roots:
                    value = root * scale
                    support_values.append(value)
                    support_mass += Q(1, value)
                    word_mass += Q(multiplicity, value)
        require(len(set(support_values)) == 3 * comb(depth + 2, 2), ("support collisions", depth))
        require(sum(multinomial3(depth, i, j) for i in range(depth + 1) for j in range(depth - i + 1)) == 3**depth, ("word count", depth))
        require(word_mass == root_harmonic * Q(311, 1001) ** depth, ("word mass", depth))
        support_shells.append(support_mass)
        word_shells.append(word_mass)

    for depth in range(3, len(support_shells)):
        require(
            1001 * support_shells[depth]
            == 311 * support_shells[depth - 1]
            - 31 * support_shells[depth - 2]
            + support_shells[depth - 3],
            ("support recurrence", depth),
        )
    for depth in range(len(word_shells) - 1):
        require(1001 * word_shells[depth + 1] == 311 * word_shells[depth], ("word recurrence", depth))

    infinite_harmonic_mass = root_harmonic * Q(7, 6) * Q(11, 10) * Q(13, 12)
    require(infinite_harmonic_mass == Q(29029, 14400), infinite_harmonic_mass)

    edge_degrees = tuple(sorted(Counter(vertex for edge in literal_edges for vertex in edge).items()))
    semantic = ExactDigest()
    semantic.add(("pair", pair_checks, "3(u+v)>14gcd(u,v)"))
    semantic.add(("triple", triple_checks, positive_checks, permutation_checks))
    semantic.add(("literal_edges", literal_edges, edge_degrees))
    semantic.add(("independence", independence_profile, maximal_five_sets))
    semantic.add(("q3_atlas", candidates, global_rows, exact_rows, tuple(rescues)))
    semantic.add(
        (
            "hostile_audit_controls",
            (1, 4, 41),
            (-7, -4, -1, 2, 5, 8),
            (2, 8, 10),
            (-2, 2, 2),
            Q(hostile_sample, 2 * hostile_scale),
        )
    )
    semantic.add(("ternary_lattice", roots, multipliers, tuple(support_shells), tuple(word_shells), infinite_harmonic_mass))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("THM-3388 Q3 PHASE-TRIANGLE COVER CLUTTER")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=PROVED analytic phase-triangle criterion plus FINITE-EXACT literal q3 clutter and atlas;independently_hostile_audited")
    print("gap_set=P(u,v)={p:p=uv_mod_3gcd(u,v),14|p|<3(u+v)}")
    print("cover_iff=exists_pqr_in_pair_gap_sets_with_w*p+u*q+v*r=0")
    print(f"pair_checks={pair_checks};triple_event_checks={triple_checks};positive_reconstructions={positive_checks};permutation_checks={permutation_checks}")
    print(f"literal_vertices={LITERAL_VERTICES};cover_edges={len(literal_edges)};pairwise_triangles={len(pairwise_triangles)};pairwise_false_positives={len(pairwise_triangles)-len(literal_edges)}")
    print(f"literal_edges={literal_edges}")
    print(f"edge_degrees={edge_degrees}")
    print(f"independence_profile={independence_profile};maximal_five_sets={maximal_five_sets}")
    print("pairwise_hostile=(1,4,7):all_gap_sets_nonempty_but_no_closed_phase_triangle")
    print("extended_pairwise_hostile=(1,4,41);gap_tail=(-7,-4,-1,2,5,8);no_closed_phase_triangle")
    print("positive_control=(1,4,5):phase_triangle=(1,-1,-1)")
    print("dilation_sign_hostile=2*(1,4,5):(p,q,r)=(-2,2,2),not_naive_(2,-2,-2)")
    print(f"q3_body_candidates={candidates};global_transverse_rows={global_rows};exact_rows={exact_rows};core_rescues={tuple(rescues)}")
    print("core_nonrescue_hostile=C:(1,3,4),U:(8,11,13),t:389/2464;only_omitted_clock_2_is_dangerous")
    print(f"ternary_dilation=roots:{roots},multipliers:{multipliers};word_nodes=3^d;distinct_scales=C(d+2,2);checks={dilation_checks}")
    print("harmonic_support_recurrence=1001H_d=311H_(d-1)-31H_(d-2)+H_(d-3)")
    print("harmonic_word_recurrence=1001W_(d+1)=311W_d")
    print(f"orbit_harmonic_mass={infinite_harmonic_mass};root_mass={root_harmonic}")
    print("typing=3_uniform_cover_clutter_with_integral_phase_cochain;not_tournament;word_order_abelianized_by_integer_support")
    print("scope=classifies_q3_slice_of_T3387;no_new_refined_ledger_decrement;no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
