#!/usr/bin/env python3
"""Exact q=8 singleton/domino mode-clutter probe.

At q=8 an odd transverse speed is the first one able to block two adjacent
phase classes.  Expanding each owner into exact singleton or domino modes
restores a one-interval carrier.  A complete affine cochain on those modes is
compared with independent exact event geometry through rank five, the full
literal rank range relevant to the body atlas.  Higher-rank physical edges
are deliberately outside this probe and are not inferred from its truncated
edge-generated profile.

This is an unnumbered FINITE-EXACT probe plus analytic mode-cochain candidate.
Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
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
        "THM-3387-output",
        ROOT / "05-knowledge/results/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.out",
        "b4d9ce439bab4501bfd5e2cf13eb0b0e3685b7364f30e43b7d5ca9138d25cb5c",
    ),
    (
        "THM-3389-script",
        ROOT / "04-computation/lrc14_q4_typed_cover_clutter_thm3389.py",
        "cd963c20ff47c9840222c6bfd95088e3d649ac90edd518353c0f64f4f5ec9bfd",
    ),
    (
        "THM-3389-output",
        ROOT / "05-knowledge/results/lrc14_q4_typed_cover_clutter_thm3389.out",
        "45459e60a69bea9d7e99746fb1b1ad8dc8f79506bb27eb7f287f887c38270adc",
    ),
)

SHEETS = 8
LITERAL_VERTICES = (1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14)
SEARCH_MAX_RANK = 5
EXPECTED_MINIMAL_EDGES = (
    (1, 3, 5, 7),
    (1, 3, 5, 9),
    (1, 3, 7, 11),
    (1, 3, 9, 11),
    (1, 5, 7, 13),
    (1, 5, 9, 13),
    (1, 7, 11, 13),
    (1, 9, 11, 13),
    (1, 10, 11, 12),
    (1, 10, 12, 13),
    (1, 11, 12, 14),
    (1, 12, 13, 14),
    (2, 6, 10, 14),
    (4, 5, 9, 14),
    (4, 9, 13, 14),
    (1, 3, 10, 11, 14),
    (1, 5, 6, 7, 11),
    (1, 5, 6, 11, 14),
    (1, 5, 7, 11, 12),
    (1, 6, 7, 10, 13),
    (1, 6, 9, 10, 14),
    (1, 6, 10, 11, 14),
    (2, 3, 7, 9, 10),
    (2, 3, 10, 11, 13),
    (2, 3, 10, 13, 14),
    (2, 5, 6, 11, 14),
    (2, 5, 9, 11, 14),
    (2, 6, 7, 9, 10),
    (2, 9, 10, 11, 14),
    (2, 10, 11, 13, 14),
    (4, 5, 7, 9, 11),
    (4, 7, 9, 11, 13),
)
EXPECTED_RANK_PROFILE = ((4, 15), (5, 17))
EXPECTED_GCD_TYPE_PROFILE = (
    ((1, 1, 1, 1), 8),
    ((2, 2, 2, 2), 1),
    ((4, 2, 1, 1), 6),
    ((2, 1, 1, 1, 1), 1),
    ((2, 2, 1, 1, 1), 6),
    ((2, 2, 2, 1, 1), 7),
    ((4, 1, 1, 1, 1), 3),
)
EXPECTED_INDEPENDENCE_PROFILE = (
    1,
    13,
    78,
    286,
    700,
    1152,
    1223,
    777,
    266,
    42,
    2,
    0,
    0,
    0,
)
EXPECTED_BODY_CANDIDATES = 1287
EXPECTED_EXACT_ROWS = 1152
EXPECTED_SEMANTIC_DIGEST = "5051166544008df0a96c50e8fc2c293e4e35b6192974417abe9a487a547d98f4"


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
    return tuple(
        Q(sample, 2 * scale)
        for sample in samples
        if full_transverse_cover(transverse, sample, scale)
        and not any(core_danger(clock, sample, scale) for clock in core)
    )


def owner_modes(speed):
    common = gcd(speed, SHEETS)
    phase_count = SHEETS // common
    phase_step = (speed // common) % phase_count
    modes = []
    for phase in range(phase_count):
        block = tuple(
            sheet
            for sheet in range(SHEETS)
            if phase_step * sheet % phase_count == phase
        )
        centre_16 = (-(16 // phase_count) * phase) % 16
        modes.append((block, centre_16, 8, "singleton", (phase,)))
    if phase_count == 8:
        for phase in range(8):
            phase_pair = (phase, (phase + 1) % 8)
            block = tuple(
                sheet
                for sheet in range(SHEETS)
                if phase_step * sheet % 8 in phase_pair
            )
            centre_16 = (-(2 * phase + 1)) % 16
            modes.append((block, centre_16, 1, "domino", phase_pair))
    return tuple(modes)


def gap_values(left_speed, right_speed, left_mode, right_mode):
    left_centre = left_mode[1]
    right_centre = right_mode[1]
    left_width = left_mode[2]
    right_width = right_mode[2]
    common = gcd(left_speed, right_speed)
    modulus = 16 * common
    residue = (
        left_centre * right_speed - right_centre * left_speed
    ) % modulus
    bound = (left_width * right_speed + right_width * left_speed - 1) // 7
    first_index = -((bound + residue) // modulus)
    last_index = (bound - residue) // modulus
    return tuple(
        residue + modulus * index
        for index in range(first_index, last_index + 1)
    )


def potential_witness(speeds, modes):
    star_sets = tuple(
        gap_values(speeds[0], speeds[index], modes[0], modes[index])
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
                if gap not in gap_values(
                    speeds[left], speeds[right], modes[left], modes[right]
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


def mode_search(subset):
    ordered = tuple(sorted(subset, key=lambda speed: (len(owner_modes(speed)), speed)))
    mode_banks = tuple(owner_modes(speed) for speed in ordered)
    max_sizes = tuple(max(len(mode[0]) for mode in bank) for bank in mode_banks)
    first_pairwise_hostile = None

    def visit(index, assigned, covered):
        nonlocal first_pairwise_hostile
        if len(covered) + sum(max_sizes[index:]) < SHEETS:
            return None
        if index == len(ordered):
            if len(covered) != SHEETS:
                return None
            witness = potential_witness(ordered, tuple(assigned))
            if witness is not None:
                return ordered, tuple(assigned), witness
            if first_pairwise_hostile is None:
                first_pairwise_hostile = (ordered, tuple(assigned))
            return None

        speed = ordered[index]
        for mode in mode_banks[index]:
            if any(
                not gap_values(ordered[prior], speed, assigned[prior], mode)
                for prior in range(index)
            ):
                continue
            result = visit(
                index + 1,
                assigned + [mode],
                covered | set(mode[0]),
            )
            if result is not None:
                return result
        return None

    witness = visit(0, [], set())
    return witness, first_pairwise_hostile


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
        tuple((gcd(speed, 8), len(owner_modes(speed))) for speed in LITERAL_VERTICES)
        == tuple(
            (
                gcd(speed, 8),
                16 if gcd(speed, 8) == 1 else 4 if gcd(speed, 8) == 2 else 2,
            )
            for speed in LITERAL_VERTICES
        ),
        "mode bank typing",
    )

    event_minimal = []
    mode_minimal = []
    mode_witnesses = []
    first_pairwise_hostile = None
    for size in range(1, SEARCH_MAX_RANK + 1):
        for subset in combinations(LITERAL_VERTICES, size):
            if any(set(edge) <= set(subset) for edge in event_minimal):
                continue
            event_edge = bool(full_cover_leaks(subset))
            mode_witness, pairwise_hostile = mode_search(subset)
            mode_edge = mode_witness is not None
            require(event_edge == mode_edge, ("event/mode mismatch", subset, event_edge, mode_witness))
            if pairwise_hostile is not None and first_pairwise_hostile is None:
                first_pairwise_hostile = pairwise_hostile
            if event_edge:
                event_minimal.append(subset)
                mode_minimal.append(subset)
                mode_witnesses.append((subset, mode_witness))

    event_minimal = tuple(event_minimal)
    mode_minimal = tuple(mode_minimal)
    require(event_minimal == EXPECTED_MINIMAL_EDGES, event_minimal)
    require(mode_minimal == EXPECTED_MINIMAL_EDGES, mode_minimal)
    require(first_pairwise_hostile is not None, "missing pairwise hostile")

    rank_profile = tuple(sorted(Counter(map(len, event_minimal)).items()))
    gcd_type_profile = tuple(
        sorted(
            Counter(
                tuple(sorted((gcd(speed, 8) for speed in edge), reverse=True))
                for edge in event_minimal
            ).items(),
            key=lambda item: (len(item[0]), item[0]),
        )
    )
    require(rank_profile == EXPECTED_RANK_PROFILE, rank_profile)
    require(gcd_type_profile == EXPECTED_GCD_TYPE_PROFILE, gcd_type_profile)

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
    for transverse in combinations(LITERAL_VERTICES, 5):
        candidates += 1
        globally_safe = independent(transverse, edge_sets)
        if globally_safe:
            global_rows += 1
            exact_rows += 1
            continue
        if not full_cover_leaks(transverse, (1,)):
            exact_rows += 1
            rescues.append(((1,), transverse))
    require(candidates == EXPECTED_BODY_CANDIDATES, candidates)
    require(global_rows == EXPECTED_EXACT_ROWS, global_rows)
    require(exact_rows == EXPECTED_EXACT_ROWS, exact_rows)
    require(not rescues, rescues)

    domino_control = next(item for item in mode_witnesses if item[0] == (1, 3, 5, 7))
    quotient_control = next(item for item in mode_witnesses if item[0] == (2, 6, 10, 14))
    require(
        tuple(mode[3] for mode in domino_control[1][1]) == ("domino",) * 4,
        domino_control,
    )
    require(
        tuple(len(mode[0]) for mode in quotient_control[1][1]) == (2, 2, 2, 2),
        quotient_control,
    )

    canonical_witness_profile = tuple(
        (
            subset,
            witness[0],
            tuple((mode[0], mode[1], mode[2], mode[3]) for mode in witness[1]),
            witness[2],
        )
        for subset, witness in mode_witnesses
    )
    semantic = ExactDigest()
    semantic.add(("minimal", event_minimal, rank_profile, gcd_type_profile))
    semantic.add(("witnesses", canonical_witness_profile))
    semantic.add(("pairwise_hostile", first_pairwise_hostile))
    semantic.add(("independence", independence_profile))
    semantic.add(("atlas", candidates, global_rows, exact_rows, tuple(rescues)))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    hostile_order, hostile_modes = first_pairwise_hostile
    print("Q8 DOMINO MODE-CLUTTER EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=FINITE-EXACT literal q8 rank<=5 probe plus analytic singleton/domino mode-cochain candidate;unnumbered_and_not_canon")
    print("mode_species=odd:8_singletons_plus_8_adjacent_dominoes;gcd2:4_antipodal_pairs;gcd4:2_parity_quadruples")
    print("mode_radii=singleton_or_coset:1/(14u);odd_domino:1/(112u)")
    print("cochain=p_ij=16u_i u_j(x_i-x_j);congruence_mod_16gcd;overlap_7|p|<w_i*u_j+w_j*u_i;zero_triangle_circulation")
    print(f"literal_vertices={LITERAL_VERTICES};minimal_edges_through_rank5={len(event_minimal)};rank_le5_profile={rank_profile};gcd_type_profile={gcd_type_profile}")
    print(f"minimal_edge_list_through_rank5={event_minimal}")
    print(f"profile_generated_by_rank_le5_edges={independence_profile}")
    print(f"domino_control={domino_control}")
    print(f"q4_quotient_lift_control={quotient_control}")
    print(f"pairwise_feasible_mode_nonedge=order:{hostile_order},modes:{tuple((m[0],m[1],m[2],m[3]) for m in hostile_modes)}")
    print(f"q8_body_candidates={candidates};global_transverse_rows={global_rows};exact_rows={exact_rows};core_rescues={tuple(rescues)}")
    print("body_identity=I5=1152;all_135_unsafe_rows_leak_outside_core_clock1")
    print("typing=mode_expanded_cover_clutter_plus_complete_affine_cochain;sign_tournament_is_lossy")
    print("scope=exact_q8_rank_le5_slice_for_T3387;higher_rank_physical_clutter_not_enumerated;no_theorem_promotion,no_refined_decrement,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
