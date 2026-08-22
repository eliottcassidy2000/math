#!/usr/bin/env python3
"""Exact q=5 singleton cover-clutter probe.

Five is the prime control between the typed q=4 and CRT-mixed q=6 fibres.
Every transverse speed blocks at most one sheet, so every minimal full cover
has rank five.  Common phase is decided by a complete affine gap cochain.

The literal clutter and atlas are reconstructed both by exact event geometry
and by the cochain criterion.  This is unnumbered and noncanonical.  Runtime
gates use RuntimeError and therefore survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, lcm
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

SHEETS = 5
LITERAL_VERTICES = (1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14)
LABELS = (0, 1, 2, 3, 4)
EXPECTED_EDGE_COUNT = 231
EXPECTED_INDEPENDENCE_PROFILE = (1, 12, 66, 220, 495, 561, 268, 45, 1, 0, 0, 0, 0)
EXPECTED_BODY_CANDIDATES = 2079
EXPECTED_GLOBAL_ROWS = 1617
EXPECTED_EXACT_ROWS = 1619
EXPECTED_CORE_RESCUES = (
    ((1,), (6, 11, 12, 13, 14)),
    ((2,), (7, 9, 11, 12, 13)),
)
EXPECTED_UNIQUE_INDEPENDENT_EIGHT = (1, 2, 3, 6, 7, 12, 13, 14)
ROBUST_HOSTILE_ORDER = (13, 16, 18, 19, 17)
EXPECTED_SEMANTIC_DIGEST = "8e76f915ddb63c6e126d57b2e82cfe8e7a3d60835db85059cd3d2c94a2e963c7"


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


def circle_norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def fires(speed, source_time, sheet=0):
    return circle_norm(speed * (source_time + Q(sheet, SHEETS))) < Q(1, 14)


def covers_at(transverse, source_time):
    return all(
        any(fires(speed, source_time, sheet) for speed in transverse)
        for sheet in range(SHEETS)
    )


def transverse_danger(speed, sample, scale, sheet):
    return strict_danger(
        speed * (SHEETS * sample + 2 * scale * sheet),
        2 * SHEETS * scale,
    )


def core_danger(clock, sample, scale):
    return strict_danger(SHEETS * clock * sample, 2 * scale)


def circle_event_samples(transverse, core=()):
    scale = 14 * SHEETS * lcm(*transverse, *core)
    events = {0}
    for speed in transverse:
        require(speed % SHEETS != 0, ("transverse type", speed))
        for sheet in range(SHEETS):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // SHEETS
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    for clock in core:
        for tooth in range(SHEETS * clock):
            for sign in (-1, 1):
                events.add(
                    (
                        scale * tooth // (SHEETS * clock)
                        + sign * scale // (14 * SHEETS * clock)
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
        any(
            transverse_danger(speed, sample, scale, sheet)
            for speed in transverse
        )
        for sheet in range(SHEETS)
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
    common = gcd(left, right)
    modulus = SHEETS * common
    residue = ((right_label - left_label) * left * right) % modulus
    bound = (SHEETS * (left + right) - 1) // 14
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


def pairwise_feasible(speeds, labels):
    return all(
        phase_gap_values(speeds[left], speeds[right], labels[left], labels[right])
        for left, right in combinations(range(len(speeds)), 2)
    )


def singleton_witness(subset):
    for order in permutations(subset):
        witness = potential_witness(order, LABELS)
        if witness is not None:
            return order, witness
    return None


def first_pairwise_hostile(subsets):
    for subset in subsets:
        for order in permutations(subset):
            if not pairwise_feasible(order, LABELS):
                continue
            if potential_witness(order, LABELS) is None:
                return order
    return None


def independent(subset, edges):
    chosen = set(subset)
    return not any(edge <= chosen for edge in edges)


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

    require(
        all(gcd(speed, SHEETS) == 1 for speed in LITERAL_VERTICES),
        "singleton blocker typing",
    )

    five_subsets = tuple(combinations(LITERAL_VERTICES, SHEETS))
    event_edges = tuple(subset for subset in five_subsets if full_cover_leaks(subset))
    cochain_witnesses = tuple(
        (subset, singleton_witness(subset)) for subset in five_subsets
    )
    cochain_edges = tuple(
        subset for subset, witness in cochain_witnesses if witness is not None
    )
    require(event_edges == cochain_edges, ("event/cochain mismatch", event_edges, cochain_edges))
    require(len(event_edges) == EXPECTED_EDGE_COUNT, len(event_edges))

    edge_sets = tuple(frozenset(edge) for edge in event_edges)
    independence_profile = tuple(
        sum(
            independent(subset, edge_sets)
            for subset in combinations(LITERAL_VERTICES, size)
        )
        for size in range(len(LITERAL_VERTICES) + 1)
    )
    require(
        independence_profile[:5] == (1, 12, 66, 220, 495),
        independence_profile,
    )
    require(
        independence_profile == EXPECTED_INDEPENDENCE_PROFILE,
        independence_profile,
    )
    independent_eights = tuple(
        subset
        for subset in combinations(LITERAL_VERTICES, 8)
        if independent(subset, edge_sets)
    )
    require(independent_eights == (EXPECTED_UNIQUE_INDEPENDENT_EIGHT,), independent_eights)

    positive_subset, positive = next(
        (subset, witness)
        for subset, witness in cochain_witnesses
        if witness is not None
    )
    pairwise_hostile = first_pairwise_hostile(five_subsets)
    require(pairwise_hostile is not None, "missing pairwise hostile")
    require(pairwise_hostile == (1, 4, 2, 3, 9), pairwise_hostile)
    forced_star = tuple(
        phase_gap_values(pairwise_hostile[0], pairwise_hostile[index], 0, index)
        for index in range(1, 5)
    )
    require(forced_star == ((-1,), (-1,), (-1,), (1,)), forced_star)
    forced_p_14 = (
        pairwise_hostile[1] * forced_star[3][0]
        - pairwise_hostile[4] * forced_star[0][0]
    ) // pairwise_hostile[0]
    legal_p_14 = phase_gap_values(
        pairwise_hostile[1], pairwise_hostile[4], 1, 4
    )
    require(forced_p_14 == 13 and legal_p_14 == (-2, 3), (forced_p_14, legal_p_14))

    robust_cardinalities = tuple(
        len(
            phase_gap_values(
                ROBUST_HOSTILE_ORDER[left],
                ROBUST_HOSTILE_ORDER[right],
                left,
                right,
            )
        )
        for left, right in combinations(range(5), 2)
    )
    require(
        robust_cardinalities == (4, 4, 5, 4, 3, 5, 5, 6, 5, 5)
        and pairwise_feasible(ROBUST_HOSTILE_ORDER, LABELS)
        and singleton_witness(tuple(sorted(ROBUST_HOSTILE_ORDER))) is None,
        robust_cardinalities,
    )

    require(
        all(
            (label_gap * left_unit * right_unit) % 5 != 0
            for label_gap in range(1, 5)
            for left_unit in range(1, 5)
            for right_unit in range(1, 5)
        ),
        "q5 no-tangency residue gate",
    )

    require(covers_at(positive_subset, Q(129, 1960)), "positive source control")
    require(
        covers_at(EXPECTED_CORE_RESCUES[0][1], Q(109, 9240))
        and fires(5, Q(109, 9240))
        and not fires(10, Q(109, 9240)),
        "q5a wrong-core control",
    )
    require(
        covers_at(EXPECTED_CORE_RESCUES[1][1], Q(93, 980))
        and fires(10, Q(93, 980))
        and not fires(5, Q(93, 980)),
        "q5b wrong-core control",
    )

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
    require(candidates == EXPECTED_BODY_CANDIDATES, candidates)
    require(global_rows == EXPECTED_GLOBAL_ROWS, global_rows)
    require(exact_rows == EXPECTED_EXACT_ROWS, exact_rows)
    require(tuple(rescues) == EXPECTED_CORE_RESCUES, rescues)
    require(
        2 * independence_profile[5] + independence_profile[4]
        == EXPECTED_GLOBAL_ROWS,
        "body identity",
    )

    edge_degrees = tuple(
        sorted(Counter(vertex for edge in event_edges for vertex in edge).items())
    )
    semantic = ExactDigest()
    semantic.add(("literal", event_edges, edge_degrees, independence_profile))
    semantic.add(
        (
            "controls",
            positive_subset,
            positive,
            pairwise_hostile,
            forced_star,
            forced_p_14,
            legal_p_14,
            ROBUST_HOSTILE_ORDER,
            robust_cardinalities,
            independent_eights,
        )
    )
    semantic.add(("atlas", candidates, global_rows, exact_rows, tuple(rescues)))
    digest = semantic.hexdigest()
    require(
        digest == EXPECTED_SEMANTIC_DIGEST,
        (
            "freeze",
            digest,
            independence_profile,
            positive_subset,
            positive,
            pairwise_hostile,
            edge_degrees,
        ),
    )

    print("Q5 SINGLETON COVER-CLUTTER EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=FINITE-EXACT plus INDEPENDENTLY HOSTILE-AUDITED q5 convergence artifact;THM3395_is_current_structural_source")
    print("blocker_species=gcd(u,5)=1:singleton;inclusion_minimal_full_cover_rank=5")
    print("cochain=p_ij=5u_i u_j(x_i-x_j);p_ij==(k_j-k_i)u_i u_j_mod_5gcd;14|p_ij|<5(u_i+u_j);zero_triangle_circulation")
    print(f"literal_vertices={LITERAL_VERTICES};minimal_edges={len(event_edges)};rank_profile=((5, {len(event_edges)}),)")
    print(f"edge_degrees={edge_degrees}")
    print(f"independence_profile={independence_profile}")
    print(f"unique_independent_eight={EXPECTED_UNIQUE_INDEPENDENT_EIGHT}")
    print(f"positive_control=subset:{positive_subset},order:{positive[0]},cochain:{positive[1]}")
    print(f"pairwise_feasible_nonedge=subset:(1,2,3,4,9),order:{pairwise_hostile},forced_star:{forced_star},forced_p_14:{forced_p_14},legal_p_14:{legal_p_14}")
    print(f"robust_pairwise_nonedge=order:{ROBUST_HOSTILE_ORDER},pair_fibre_cardinalities:{robust_cardinalities},no_ordering_covers")
    print("strict_closed_full_cover_edges_coincide=PROVED:no_distinct_q5_sheet_labels_can_be_tangent_mod5")
    print("rescue_specific_controls=q5a_wrong_core2_leaks_at_109/9240;q5b_wrong_core1_leaks_at_93/980")
    print(f"q5_body_candidates={candidates};global_transverse_rows={global_rows};exact_rows={exact_rows};core_rescues={tuple(rescues)}")
    print("body_identity=2I5+I4=2*561+495=1617;two_core_rescues_give_1619")
    print("typing=prime_singleton_rank5_clutter_plus_complete_affine_K5_cochain;not_tournament")
    print("scope=independent_q5_convergence_audit_for_T3387_and_T3395;no_new_refined_decrement,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
