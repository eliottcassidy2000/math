#!/usr/bin/env python3
"""Exact q=6 typed cover-clutter probe.

This is an unnumbered research artifact because the live theorem namespace is
concurrently occupied.  It proves the finite literal census and checks the
general affine-cochain candidate against an independent event sweep.

On six sheets, a transverse speed blocks one coset of size gcd(u,6): a
singleton, an antipodal pair, or a parity triple.  Inclusion-minimal covers
are drawn from the 23 irredundant coset covers of Z/6.  Common phase is tested
by a complete affine gap cochain.  Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Q3_PATH = ROOT / "04-computation/lrc14_q3_phase_triangle_clutter_thm3388.py"
PINNED = (
    (
        "THM-3387",
        ROOT / "01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md",
        "c540255a185efb54c67035d69f3fbd94f4c1ad3e30c4e31738a8800e81198613",
    ),
    (
        "THM-3388",
        ROOT / "01-canon/theorems/THM-3388-three-sheet-phase-triangle-cover-clutter.md",
        "30804bfc27f9d9d0bc54335816e23a94ae90abe497efb2b3ccebea4d870fc670",
    ),
    (
        "THM-3388-script",
        Q3_PATH,
        "5323346310a9a6b188caa0131b177b2ae8e23c7113808cda8955f89828e62154",
    ),
    (
        "THM-3388-output",
        ROOT / "05-knowledge/results/lrc14_q3_phase_triangle_clutter_thm3388.out",
        "5a32319fb8a91b476d292da292ae3cc9933f5f94aad7eb0e834f49e52252c535",
    ),
)

LITERAL_VERTICES = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14)
EXPECTED_MINIMAL_EDGES = (
    (2, 8, 10),
    (2, 10, 14),
    (4, 10, 14),
    (1, 4, 7, 10),
    (1, 5, 7, 9),
    (1, 5, 9, 11),
    (1, 5, 9, 14),
    (1, 7, 8, 9),
    (1, 7, 9, 13),
    (1, 8, 9, 14),
    (1, 8, 10, 11),
    (1, 8, 10, 13),
    (1, 8, 11, 14),
    (1, 8, 13, 14),
    (1, 9, 10, 11),
    (1, 9, 11, 13),
    (1, 9, 13, 14),
    (2, 3, 5, 7),
    (2, 3, 5, 14),
    (2, 3, 7, 11),
    (2, 3, 11, 13),
    (2, 3, 11, 14),
    (2, 5, 13, 14),
    (2, 7, 10, 11),
    (2, 11, 13, 14),
    (3, 4, 7, 11),
    (3, 4, 11, 13),
    (3, 8, 11, 13),
    (4, 5, 11, 14),
    (4, 7, 10, 13),
    (5, 8, 11, 14),
    (10, 11, 13, 14),
    (1, 4, 5, 7, 11),
    (1, 5, 7, 8, 11),
    (1, 5, 7, 8, 13),
    (1, 7, 10, 11, 13),
    (2, 5, 7, 11, 13),
    (4, 5, 7, 11, 13),
    (5, 7, 8, 11, 13),
)
EXPECTED_RANK_PROFILE = ((3, 3), (4, 29), (5, 7))
EXPECTED_TYPE_PROFILE = (
    ((2, 1, 1, 1, 1), 7),
    ((2, 2, 1, 1), 12),
    ((2, 2, 2), 3),
    ((3, 1, 1, 1), 4),
    ((3, 2, 1, 1), 10),
    ((3, 2, 2, 1), 3),
)
EXPECTED_BLOCK_PATTERN_PROFILE = (
    ((1, 1, 1, 1, 1, 1), 1),
    ((2, 1, 1, 1, 1), 3),
    ((2, 2, 1, 1), 3),
    ((2, 2, 2), 1),
    ((3, 1, 1, 1), 2),
    ((3, 2, 1, 1), 6),
    ((3, 2, 2, 1), 6),
    ((3, 3), 1),
)
EXPECTED_INDEPENDENCE_PROFILE = (1, 12, 66, 217, 441, 515, 304, 76, 5, 0, 0, 0, 0)
EXPECTED_CANDIDATES = 2079
EXPECTED_GLOBAL_ROWS = 1471
EXPECTED_EXACT_ROWS = 1478
EXPECTED_CORE_RESCUES = (
    ((1,), (1, 3, 8, 11, 13)),
    ((1,), (3, 5, 8, 11, 13)),
    ((1,), (3, 7, 8, 11, 13)),
    ((1,), (3, 8, 9, 11, 13)),
    ((1,), (3, 8, 10, 11, 13)),
    ((1,), (3, 8, 11, 13, 14)),
    ((1, 2), (3, 8, 11, 13)),
)
EXPECTED_SEMANTIC_DIGEST = "c329b80567995d2ce6110b5e216695b99faa062d04cb70fdfa5bc0db7edc8019"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def load_q3():
    spec = spec_from_file_location("thm3388_q3", Q3_PATH)
    require(spec is not None and spec.loader is not None, "q3 import spec")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
    return strict_danger(speed * (6 * sample + 2 * scale * sheet), 12 * scale)


def core_danger(clock, sample, scale):
    return strict_danger(6 * clock * sample, 2 * scale)


def circle_event_samples(transverse, core=()):
    scale = 84 * lcm(*transverse, *core)
    events = {0}
    for speed in transverse:
        require(speed % 6 != 0, ("transverse type", speed))
        for sheet in range(6):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // 6
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    for clock in core:
        for tooth in range(6 * clock):
            for sign in (-1, 1):
                events.add(
                    (
                        scale * tooth // (6 * clock)
                        + sign * scale // (84 * clock)
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
        for sheet in range(6)
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


def minimal_block_patterns():
    sheets = set(range(6))
    blocks = tuple(
        (
            size,
            label,
            frozenset(sheet for sheet in range(6) if sheet % (6 // size) == label),
        )
        for size in (3, 2, 1)
        for label in range(6 // size)
    )
    patterns = []
    for count in range(1, 7):
        for chosen in combinations(blocks, count):
            if set().union(*(block[2] for block in chosen)) != sheets:
                continue
            if any(
                set().union(*(block[2] for index, block in enumerate(chosen) if index != removed))
                == sheets
                for removed in range(count)
            ):
                continue
            patterns.append(tuple((block[0], block[1]) for block in chosen))
    return tuple(patterns)


def phase_gap_values(left, right, left_label, right_label):
    common = gcd(left, right)
    modulus = 6 * common
    residue = ((right_label - left_label) * left * right) % modulus
    bound = (3 * (left + right) - 1) // 7
    first_index = -((bound + residue) // modulus)
    last_index = (bound - residue) // modulus
    return tuple(
        residue + modulus * index
        for index in range(first_index, last_index + 1)
    )


def potential_witness(speeds, labels):
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
            return tuple(
                (left, right, edges[(left, right)])
                for left, right in combinations(range(len(speeds)), 2)
            )
    return None


def typed_witness(subset, block_patterns):
    speed_type = tuple(sorted((gcd(speed, 6) for speed in subset), reverse=True))
    for pattern in block_patterns:
        pattern_type = tuple(sorted((item[0] for item in pattern), reverse=True))
        if pattern_type != speed_type:
            continue
        owner_types = tuple(item[0] for item in pattern)
        labels = tuple(item[1] for item in pattern)
        for order in permutations(subset):
            if tuple(gcd(speed, 6) for speed in order) != owner_types:
                continue
            witness = potential_witness(order, labels)
            if witness is not None:
                return pattern, order, witness
    return None


def independent(subset, edges):
    chosen = set(subset)
    return not any(edge <= chosen for edge in edges)


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))
    q3 = load_q3()

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

    block_patterns = minimal_block_patterns()
    block_profile = tuple(
        sorted(
            Counter(
                tuple(sorted((item[0] for item in pattern), reverse=True))
                for pattern in block_patterns
            ).items()
        )
    )
    require(len(block_patterns) == 23, len(block_patterns))
    require(block_profile == EXPECTED_BLOCK_PATTERN_PROFILE, block_profile)

    # Quotient tower: parity triples reduce to q=2; antipodal pairs reduce to q=3.
    q2_quotient_checks = 0
    triple_blockers = tuple(number for number in range(3, 200, 6))
    parity_pattern = ((3, 0), (3, 1))
    for left, right in combinations(triple_blockers, 2):
        reduced_left, reduced_right = left // 3, right // 3
        predicted = reduced_left + reduced_right > 7 * gcd(reduced_left, reduced_right)
        cochain = potential_witness((left, right), (0, 1)) is not None
        exact = bool(full_cover_leaks((left, right)))
        require(predicted == cochain == exact, ("q2 quotient", left, right))
        q2_quotient_checks += 1
    require(parity_pattern in block_patterns, "parity pattern")

    q3_quotient_checks = 0
    reduced_pool = tuple(number for number in range(1, 31) if number % 3 != 0)
    pair_pattern = ((2, 0), (2, 1), (2, 2))
    for reduced in combinations(reduced_pool, 3):
        lifted = tuple(2 * speed for speed in reduced)
        q3_edge = q3.phase_triangle(reduced) is not None
        q6_edge = potential_witness(lifted, (0, 1, 2)) is not None
        require(q3_edge == q6_edge, ("q3 quotient", reduced, lifted))
        q3_quotient_checks += 1
    require(pair_pattern in block_patterns, "pair pattern")
    require(
        all(
            phase_gap_values(left, right, left_label, right_label)
            for left, right, left_label, right_label in (
                (2, 8, 0, 1),
                (2, 14, 0, 2),
                (8, 14, 1, 2),
            )
        )
        and potential_witness((2, 8, 14), (0, 1, 2)) is None,
        "lifted q3 hostile",
    )

    # Independent event route and affine cochain route give the same 39 edges.
    event_minimal = []
    cochain_minimal = []
    for size in range(1, 7):
        for subset in combinations(LITERAL_VERTICES, size):
            if full_cover_leaks(subset) and not any(
                set(edge) <= set(subset) for edge in event_minimal
            ):
                event_minimal.append(subset)
            witness = typed_witness(subset, block_patterns)
            if witness is not None and not any(
                set(edge) <= set(subset) for edge in cochain_minimal
            ):
                cochain_minimal.append(subset)
    event_minimal = tuple(event_minimal)
    cochain_minimal = tuple(cochain_minimal)
    require(event_minimal == EXPECTED_MINIMAL_EDGES, event_minimal)
    require(cochain_minimal == EXPECTED_MINIMAL_EDGES, cochain_minimal)
    rank_profile = tuple(sorted(Counter(map(len, event_minimal)).items()))
    type_profile = tuple(
        sorted(
            Counter(
                tuple(sorted((gcd(speed, 6) for speed in edge), reverse=True))
                for edge in event_minimal
            ).items()
        )
    )
    require(rank_profile == EXPECTED_RANK_PROFILE, rank_profile)
    require(type_profile == EXPECTED_TYPE_PROFILE, type_profile)

    edge_sets = tuple(frozenset(edge) for edge in event_minimal)
    independence_profile = tuple(
        sum(
            independent(subset, edge_sets)
            for subset in combinations(LITERAL_VERTICES, size)
        )
        for size in range(len(LITERAL_VERTICES) + 1)
    )
    require(independence_profile == EXPECTED_INDEPENDENCE_PROFILE, independence_profile)

    candidates = 0
    global_rows = 0
    exact_rows = 0
    rescues = []
    for core_size in (1, 2):
        transverse_size = 6 - core_size
        for core in combinations((1, 2), core_size):
            for transverse in combinations(LITERAL_VERTICES, transverse_size):
                candidates += 1
                globally_safe = independent(transverse, edge_sets)
                if globally_safe:
                    global_rows += 1
                    exact_rows += 1
                    continue
                if not full_cover_leaks(transverse, core):
                    exact_rows += 1
                    rescues.append((core, transverse))
    require(candidates == EXPECTED_CANDIDATES, candidates)
    require(global_rows == EXPECTED_GLOBAL_ROWS, global_rows)
    require(exact_rows == EXPECTED_EXACT_ROWS, exact_rows)
    require(tuple(rescues) == EXPECTED_CORE_RESCUES, rescues)
    require(
        2 * independence_profile[5] + independence_profile[4]
        == EXPECTED_GLOBAL_ROWS,
        "body identity",
    )

    edge_degrees = tuple(sorted(Counter(vertex for edge in event_minimal for vertex in edge).items()))
    semantic = ExactDigest()
    semantic.add(("block_patterns", block_patterns, block_profile))
    semantic.add(("quotients", q2_quotient_checks, q3_quotient_checks, (2, 8, 14)))
    semantic.add(("minimal", event_minimal, rank_profile, type_profile, edge_degrees))
    semantic.add(("independence", independence_profile))
    semantic.add(("atlas", candidates, global_rows, exact_rows, tuple(rescues)))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("Q6 TYPED COVER-CLUTTER EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=FINITE-EXACT literal q6 probe plus analytic affine-cochain candidate;unnumbered_and_not_canon")
    print("blocker_species=gcd(u,6):1_singleton,2_antipodal_pair,3_parity_triple")
    print(f"minimal_Z6_block_patterns={len(block_patterns)};block_type_profile={block_profile}")
    print(f"quotient_checks=q2_parity_pairs:{q2_quotient_checks},q3_pair_triples:{q3_quotient_checks}")
    print("q3_lift_hostile=(2,8,14)_reduces_to_(1,4,7):all_pairs_feasible_but_no_phase_closure")
    print(f"literal_vertices={LITERAL_VERTICES};minimal_edges={len(event_minimal)};rank_profile={rank_profile};type_profile={type_profile}")
    print(f"minimal_edge_list={event_minimal}")
    print(f"edge_degrees={edge_degrees}")
    print(f"independence_profile={independence_profile}")
    print(f"q6_body_candidates={candidates};global_transverse_rows={global_rows};exact_rows={exact_rows};core_rescues={tuple(rescues)}")
    print("body_identity=2I5+I4=2*515+441=1471;seven_core_rescues_give_1478")
    print("typing=nonuniform_Z6_coset_cover_clutter_plus_complete_affine_cochain;not_tournament")
    print("scope=exact_q6_slice_probe_for_T3387;no_theorem_promotion,no_refined_decrement,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
