#!/usr/bin/env python3
"""Exact q=4 typed cover clutter for proved THM-3389.

Odd transverse speeds block singleton sheets; speeds congruent to two modulo
four block antipodal pairs.  Minimal full covers therefore have block-size
types 2+2, 2+1+1, or 1+1+1+1.  A complete affine gap cochain decides whether
the assigned blocker intervals share a common phase.

The cochain criterion is compared with an independent exact event sweep.
Runtime gates use RuntimeError, so python -O retains every decision.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
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

LITERAL_VERTICES = (1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14)
EXPECTED_MINIMAL_EDGES = (
    (2, 14),
    (6, 10),
    (6, 14),
    (10, 14),
    (1, 9, 10),
    (1, 10, 11),
    (1, 10, 13),
    (1, 11, 14),
    (1, 13, 14),
    (2, 7, 9),
    (2, 9, 11),
    (2, 11, 13),
    (3, 7, 10),
    (3, 10, 13),
    (3, 11, 14),
    (5, 6, 11),
    (5, 9, 14),
    (6, 7, 13),
    (9, 13, 14),
    (1, 3, 5, 7),
    (1, 3, 5, 9),
    (1, 3, 7, 11),
    (1, 3, 9, 11),
    (1, 5, 7, 11),
    (1, 5, 7, 13),
    (1, 5, 9, 13),
    (1, 7, 9, 13),
    (1, 7, 11, 13),
    (1, 9, 11, 13),
    (3, 5, 7, 9),
    (3, 5, 7, 11),
    (3, 5, 11, 13),
    (3, 7, 11, 13),
    (5, 7, 9, 11),
    (5, 9, 11, 13),
    (7, 9, 11, 13),
)
EXPECTED_RANK_PROFILE = ((2, 4), (3, 15), (4, 17))
EXPECTED_INDEPENDENCE_PROFILE = (1, 11, 51, 118, 123, 44, 3, 0, 0, 0, 0, 0)
EXPECTED_MAXIMAL_SIX_SETS = (
    (1, 2, 3, 5, 6, 13),
    (1, 2, 3, 6, 9, 13),
    (2, 3, 5, 6, 9, 13),
)
EXPECTED_BODY_CANDIDATES = 2541
EXPECTED_EXACT_ROWS = 619
EXPECTED_SEMANTIC_DIGEST = "5e19b6e083e8b1faf6c1a082c6e321f148b2e5ff85bf93aa1d15c92542bf87b2"


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
    # source time sample/(2S) plus sheet/4 has denominator 8S.
    return strict_danger(speed * (4 * sample + 2 * scale * sheet), 8 * scale)


def core_danger(clock, sample, scale):
    return strict_danger(4 * clock * sample, 2 * scale)


def circle_event_samples(transverse, core=()):
    scale = 56 * lcm(*transverse, *core)
    events = {0}
    for speed in transverse:
        require(speed % 4 != 0, ("transverse type", speed))
        for sheet in range(4):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // 4
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    for clock in core:
        for tooth in range(4 * clock):
            for sign in (-1, 1):
                events.add(
                    (
                        scale * tooth // (4 * clock)
                        + sign * scale // (56 * clock)
                    )
                    % scale
                )
    ordered = tuple(sorted(events))
    endpoints = tuple(2 * event for event in ordered)
    midpoints = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += scale
        midpoints.append((left + right) % (2 * scale))
    return scale, endpoints + tuple(midpoints)


def full_transverse_cover(transverse, sample, scale):
    return all(
        any(transverse_danger(speed, sample, scale, sheet) for speed in transverse)
        for sheet in range(4)
    )


def full_cover_leaks(transverse, core=()):
    scale, samples = circle_event_samples(transverse, core)
    leaks = []
    for sample in samples:
        if not full_transverse_cover(transverse, sample, scale):
            continue
        if any(core_danger(clock, sample, scale) for clock in core):
            continue
        leaks.append(Q(sample, 2 * scale))
    return tuple(leaks)


def phase_gap_values(left, right, left_label, right_label):
    """Cleared gaps 4uv(x_left-x_right) for assigned sheet labels."""
    common = gcd(left, right)
    modulus = 4 * common
    residue = ((right_label - left_label) * left * right) % modulus
    bound = (2 * (left + right) - 1) // 7
    first_index = -((bound + residue) // modulus)
    last_index = (bound - residue) // modulus
    return tuple(
        residue + modulus * index
        for index in range(first_index, last_index + 1)
    )


def potential_witness(speeds, labels):
    """Find a complete affine gap cochain with zero triangle circulation."""
    require(len(speeds) == len(labels) and 2 <= len(speeds) <= 4, (speeds, labels))
    star_sets = tuple(
        phase_gap_values(speeds[0], speeds[index], labels[0], labels[index])
        for index in range(1, len(speeds))
    )
    for star in product(*star_sets):
        edges = {(0, index): star[index - 1] for index in range(1, len(speeds))}
        good = True
        for left in range(1, len(speeds)):
            for right in range(left + 1, len(speeds)):
                numerator = (
                    speeds[left] * star[right - 1]
                    - speeds[right] * star[left - 1]
                )
                if numerator % speeds[0] != 0:
                    good = False
                    break
                gap = numerator // speeds[0]
                if gap not in phase_gap_values(
                    speeds[left], speeds[right], labels[left], labels[right]
                ):
                    good = False
                    break
                edges[(left, right)] = gap
            if not good:
                break
        if good:
            return tuple((left, right, edges[(left, right)]) for left, right in combinations(range(len(speeds)), 2))
    return None


def typed_edge_witness(subset):
    paired = tuple(speed for speed in subset if speed % 4 == 2)
    single = tuple(speed for speed in subset if speed % 2 == 1)
    if len(subset) == 2 and len(paired) == 2:
        return potential_witness(paired, (0, 1))
    if len(subset) == 3 and len(paired) == 1 and len(single) == 2:
        for order in (single, single[::-1]):
            witness = potential_witness((paired[0],) + order, (0, 1, 3))
            if witness is not None:
                return witness
        return None
    if len(subset) == 4 and len(single) == 4:
        for order in permutations(single):
            witness = potential_witness(order, (0, 1, 2, 3))
            if witness is not None:
                return witness
        return None
    return None


def independent(subset, edges):
    chosen = set(subset)
    return not any(edge <= chosen for edge in edges)


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
    paired_pool = tuple(number for number in range(2, 200, 4))
    for left, right in combinations(paired_pool, 2):
        predicted = left + right > 7 * gcd(left, right)
        cochain = potential_witness((left, right), (0, 1)) is not None
        exact = bool(full_cover_leaks((left, right)))
        require(predicted == cochain == exact, ("paired edge", left, right))
        pair_checks += 1

    rank3_checks = 0
    audit_paired = tuple(number for number in range(2, 41, 4))
    audit_odd = tuple(number for number in range(1, 41, 2))
    for paired_speed in audit_paired:
        for left, right in combinations(audit_odd, 2):
            subset = tuple(sorted((paired_speed, left, right)))
            cochain = typed_edge_witness(subset) is not None
            exact = bool(full_cover_leaks(subset))
            require(cochain == exact, ("rank3 cochain", subset))
            rank3_checks += 1

    rank4_checks = 0
    audit_odd4 = tuple(number for number in range(1, 26, 2))
    for subset in combinations(audit_odd4, 4):
        cochain = typed_edge_witness(subset) is not None
        exact = bool(full_cover_leaks(subset))
        require(cochain == exact, ("rank4 cochain", subset))
        rank4_checks += 1

    # Independent event route finds every inclusion-minimal literal cover.
    actual_minimal = []
    for size in range(1, 7):
        for subset in combinations(LITERAL_VERTICES, size):
            if not full_cover_leaks(subset):
                continue
            if any(set(edge) <= set(subset) for edge in actual_minimal):
                continue
            actual_minimal.append(subset)
    actual_minimal = tuple(actual_minimal)
    predicted_minimal = tuple(
        subset
        for size in (2, 3, 4)
        for subset in combinations(LITERAL_VERTICES, size)
        if typed_edge_witness(subset) is not None
    )
    require(actual_minimal == EXPECTED_MINIMAL_EDGES, actual_minimal)
    require(predicted_minimal == EXPECTED_MINIMAL_EDGES, predicted_minimal)
    rank_profile = tuple(sorted(Counter(map(len, actual_minimal)).items()))
    require(rank_profile == EXPECTED_RANK_PROFILE, rank_profile)

    edge_sets = tuple(frozenset(edge) for edge in actual_minimal)
    independence_profile = tuple(
        sum(
            independent(subset, edge_sets)
            for subset in combinations(LITERAL_VERTICES, size)
        )
        for size in range(len(LITERAL_VERTICES) + 1)
    )
    maximal_six_sets = tuple(
        subset
        for subset in combinations(LITERAL_VERTICES, 6)
        if independent(subset, edge_sets)
    )
    require(independence_profile == EXPECTED_INDEPENDENCE_PROFILE, independence_profile)
    require(maximal_six_sets == EXPECTED_MAXIMAL_SIX_SETS, maximal_six_sets)

    # Complete literal q=4 body atlas.  There are no core rescues.
    core_pool = (1, 2, 3)
    candidates = 0
    global_rows = 0
    exact_rows = 0
    rescues = []
    for core_size in range(1, 4):
        transverse_size = 6 - core_size
        for core in combinations(core_pool, core_size):
            for transverse in combinations(LITERAL_VERTICES, transverse_size):
                candidates += 1
                globally_safe = independent(transverse, edge_sets)
                if globally_safe:
                    global_rows += 1
                    exact_rows += 1
                    continue
                leaks = full_cover_leaks(transverse, core)
                if not leaks:
                    exact_rows += 1
                    rescues.append((core, transverse))
    require(candidates == EXPECTED_BODY_CANDIDATES, candidates)
    require(global_rows == EXPECTED_EXACT_ROWS, global_rows)
    require(exact_rows == EXPECTED_EXACT_ROWS, exact_rows)
    require(not rescues, rescues)
    require(
        3 * independence_profile[5]
        + 3 * independence_profile[4]
        + independence_profile[3]
        == EXPECTED_EXACT_ROWS,
        "body independence identity",
    )

    rank2_positive = potential_witness((2, 14), (0, 1))
    rank3_positive = potential_witness((10, 1, 9), (0, 1, 3))
    rank4_positive = potential_witness((1, 3, 7, 5), (0, 1, 2, 3))
    require(rank2_positive == ((0, 1, -4),), rank2_positive)
    require(rank3_positive == ((0, 1, -2), (0, 2, 2), (1, 2, 2)), rank3_positive)
    require(
        rank4_positive
        == ((0, 1, -1), (0, 2, -2), (0, 3, -1), (1, 2, 1), (1, 3, 2), (2, 3, 3)),
        rank4_positive,
    )
    require(
        all(
            phase_gap_values(speed_left, speed_right, label_left, label_right)
            for speed_left, speed_right, label_left, label_right in (
                (2, 7, 0, 1),
                (2, 11, 0, 3),
                (7, 11, 1, 3),
            )
        )
        and typed_edge_witness((2, 7, 11)) is None,
        "rank3 pairwise hostile",
    )
    hostile_order = (1, 3, 11, 5)
    require(
        all(
            phase_gap_values(hostile_order[left], hostile_order[right], left, right)
            for left, right in combinations(range(4), 2)
        )
        and potential_witness(hostile_order, (0, 1, 2, 3)) is None,
        "rank4 pairwise hostile",
    )

    # Odd common dilation preserves blocker species and permutes the sheets.
    roots = (1, 9, 10)
    multipliers = (7, 11, 13)
    root_mass = sum((Q(1, root) for root in roots), Q(0))
    require(root_mass == Q(109, 90), root_mass)
    support_shells = []
    word_shells = []
    dilation_checks = 0
    for depth in range(16):
        values = []
        support_mass = Q(0)
        word_mass = Q(0)
        for i in range(depth + 1):
            for j in range(depth - i + 1):
                k = depth - i - j
                scale = 7**i * 11**j * 13**k
                scaled = tuple(scale * root for root in roots)
                require(typed_edge_witness(scaled) is not None, ("dilation edge", depth, i, j))
                dilation_checks += 1
                multiplicity = multinomial3(depth, i, j)
                for root in roots:
                    value = root * scale
                    values.append(value)
                    support_mass += Q(1, value)
                    word_mass += Q(multiplicity, value)
        require(len(set(values)) == 3 * comb(depth + 2, 2), ("support collisions", depth))
        require(word_mass == root_mass * Q(311, 1001) ** depth, ("word mass", depth))
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
    orbit_mass = root_mass * Q(7, 6) * Q(11, 10) * Q(13, 12)
    require(orbit_mass == Q(109109, 64800), orbit_mass)

    edge_degrees = tuple(sorted(Counter(vertex for edge in actual_minimal for vertex in edge).items()))
    semantic = ExactDigest()
    semantic.add(("audit", pair_checks, rank3_checks, rank4_checks))
    semantic.add(("minimal", actual_minimal, rank_profile, edge_degrees))
    semantic.add(("independence", independence_profile, maximal_six_sets))
    semantic.add(("atlas", candidates, global_rows, exact_rows, tuple(rescues)))
    semantic.add(("controls", rank2_positive, rank3_positive, rank4_positive))
    semantic.add(("dilation", roots, multipliers, tuple(support_shells), tuple(word_shells), orbit_mass))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("THM-3389 Q4 TYPED COVER CLUTTER")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=PROVED analytic complete-cochain criterion plus FINITE-EXACT literal q4 clutter and atlas;independently_hostile_audited")
    print("blocker_species=odd:singleton;2mod4:antipodal_pair;minimal_partitions=2+2,2+1+1,1+1+1+1")
    print("gap_cochain=p_ij=4u_i u_j(x_i-x_j);affine_congruence_and_pair_overlap_plus_zero_triangle_circulation")
    print(f"audit_checks=rank2:{pair_checks},rank3:{rank3_checks},rank4:{rank4_checks}")
    print(f"literal_vertices={LITERAL_VERTICES};minimal_edges={len(actual_minimal)};rank_profile={rank_profile}")
    print(f"minimal_edge_list={actual_minimal}")
    print(f"edge_degrees={edge_degrees}")
    print(f"independence_profile={independence_profile};maximal_six_sets={maximal_six_sets}")
    print("rank3_pairwise_hostile=(2,7,11);rank4_pairwise_hostile_order=(1,3,11,5)")
    print(f"positive_cochains=rank2:{rank2_positive};rank3:{rank3_positive};rank4:{rank4_positive}")
    print(f"q4_body_candidates={candidates};global_transverse_rows={global_rows};exact_rows={exact_rows};core_rescues={tuple(rescues)}")
    print("body_identity=3I5+3I4+I3=3*44+3*123+118=619")
    print(f"typed_ternary_dilation=roots:{roots},multipliers:{multipliers};checks={dilation_checks};orbit_mass={orbit_mass}")
    print("harmonic_support_recurrence=1001H_d=311H_(d-1)-31H_(d-2)+H_(d-3),d>=3")
    print("harmonic_word_recurrence=1001W_(d+1)=311W_d,d>=0")
    print("typing=nonuniform_cover_clutter_with_block_sizes_and_exact_affine_Kr_cochain;not_tournament")
    print("scope=classifies_q4_slice_of_T3387;no_core_rescue;no_new_refined_ledger_decrement;no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
