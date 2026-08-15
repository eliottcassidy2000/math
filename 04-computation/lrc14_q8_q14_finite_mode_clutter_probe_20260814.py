#!/usr/bin/env python3
"""Exact q=8,...,14 finite-mode cover-clutter compiler.

For phase count m=q/gcd(u,q), a mode is a consecutive block of
s<=ceil(m/7) phase classes.  It has one rational centre lattice and one
rational interval radius, so compatible mode covers are again decided by a
complete affine star cochain.  The compiler compares this analytic carrier
with exact event geometry on every literal subset through rank five, the full
rank range relevant to the body atlas.  Higher-rank physical edges are outside
this probe and are not inferred from the truncated edge-generated profiles.

This is an unnumbered FINITE-EXACT probe plus analytic theorem candidate.
Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Q8_PATH = ROOT / "04-computation/lrc14_q8_domino_mode_clutter_probe_20260814.py"
PINNED = (
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
    (
        "THM-3395-script",
        ROOT / "04-computation/lrc14_small_sheet_typed_cover_star_cochain_thm3395.py",
        "0953401d98d62fd3118bd4a7bbeb50bd43a459a04dc120578ca6af355cada067",
    ),
    (
        "THM-3395-output",
        ROOT / "05-knowledge/results/lrc14_small_sheet_typed_cover_star_cochain_thm3395.out",
        "f7ed05e16fdd3660741aa8a79600cf9920bbebd8087c8d25a252ecca0dbc1ce5",
    ),
    (
        "q8-script",
        Q8_PATH,
        "3e523a2ff8cbd6329782347c56fae2d8519a161c3d127697ca452f3891890b9c",
    ),
    (
        "q8-output",
        ROOT / "05-knowledge/results/lrc14_q8_domino_mode_clutter_probe_20260814.out",
        "0f5a421205bc559c8f12dce8462b4d570fcffba0e602740d1ea66c52cd84d045",
    ),
)

EXPECTED_EDGES = {
    9: (
        (1, 5, 6, 7),
        (1, 5, 6, 11),
        (1, 5, 7, 12),
        (1, 5, 11, 12),
        (1, 6, 7, 13),
        (1, 6, 11, 13),
        (1, 7, 12, 13),
        (1, 11, 12, 13),
        (2, 10, 12, 14),
        (2, 3, 5, 7, 12),
        (2, 3, 5, 8, 12),
        (2, 3, 7, 10, 12),
        (2, 3, 8, 10, 12),
        (2, 5, 8, 10, 14),
        (2, 5, 12, 13, 14),
        (2, 7, 10, 11, 12),
        (3, 5, 7, 12, 13),
        (3, 5, 8, 12, 13),
        (3, 7, 10, 12, 13),
        (3, 8, 10, 12, 13),
        (3, 8, 10, 13, 14),
        (6, 8, 10, 13, 14),
    ),
    10: (
        (1, 3, 4, 7, 9),
        (1, 3, 4, 7, 11),
        (1, 3, 4, 9, 13),
        (1, 3, 4, 11, 13),
        (1, 3, 7, 8, 9),
        (1, 3, 7, 8, 11),
        (1, 3, 7, 9, 12),
        (1, 3, 7, 11, 12),
        (1, 3, 8, 9, 13),
        (1, 3, 8, 11, 13),
        (1, 3, 9, 12, 13),
        (1, 3, 11, 12, 13),
        (1, 4, 9, 13, 14),
        (1, 6, 7, 11, 13),
        (1, 7, 11, 12, 13),
        (1, 8, 9, 13, 14),
        (2, 4, 6, 8, 14),
        (2, 4, 8, 12, 14),
    ),
    11: (),
    12: (
        (4, 6, 10, 14),
        (1, 5, 7, 8, 11),
        (1, 5, 7, 8, 13),
        (1, 8, 10, 11, 14),
        (1, 8, 10, 13, 14),
        (4, 5, 6, 11, 14),
        (4, 6, 7, 10, 13),
        (5, 6, 8, 11, 14),
    ),
    13: (),
    14: (),
}

EXPECTED_PROFILES = {
    9: (1, 13, 78, 286, 706, 1205, 1401, 1060, 483, 116, 11, 0, 0, 0),
    10: (1, 13, 78, 286, 715, 1269, 1602, 1404, 811, 284, 53, 4, 0, 0),
    11: (1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1),
    12: (1, 13, 78, 286, 714, 1271, 1629, 1477, 907, 348, 71, 5, 0, 0),
    13: (1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1),
    14: (1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1),
}

EXPECTED_ATLAS_ROWS = {
    8: 1152,
    9: 1205,
    10: 1269,
    11: 1287,
    12: 1271,
    13: 1287,
    14: 1287,
}
SEARCH_MAX_RANK = 5
EXPECTED_SEMANTIC_DIGEST = "66a69a30c49b72ff8ecbf7de94f495025518e04b73969f2d970debeb6f113023"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def load_q8():
    spec = spec_from_file_location("q8_probe", Q8_PATH)
    require(spec is not None and spec.loader is not None, "q8 module spec")
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


def transverse_danger(q, speed, sample, scale, sheet):
    return strict_danger(
        speed * (q * sample + 2 * scale * sheet),
        2 * q * scale,
    )


def core_danger(q, clock, sample, scale):
    return strict_danger(q * clock * sample, 2 * scale)


def circle_event_samples(q, transverse, core=()):
    scale = 14 * q * lcm(*transverse, *core)
    events = {0}
    for speed in transverse:
        require(speed % q != 0, ("transverse type", q, speed))
        for sheet in range(q):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // q
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    for clock in core:
        for tooth in range(q * clock):
            for sign in (-1, 1):
                events.add(
                    (
                        scale * tooth // (q * clock)
                        + sign * scale // (14 * q * clock)
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


def full_cover(q, transverse, sample, scale):
    return all(
        any(
            transverse_danger(q, speed, sample, scale, sheet)
            for speed in transverse
        )
        for sheet in range(q)
    )


def has_cover(q, transverse):
    scale, samples = circle_event_samples(q, transverse)
    return any(full_cover(q, transverse, sample, scale) for sample in samples)


def leaks_outside_core(q, transverse):
    scale, samples = circle_event_samples(q, transverse, (1,))
    return any(
        full_cover(q, transverse, sample, scale)
        and not core_danger(q, 1, sample, scale)
        for sample in samples
    )


def owner_modes(q, speed):
    common = gcd(speed, q)
    phase_count = q // common
    phase_step = (speed // common) % phase_count
    modes = []
    for size in range(1, (phase_count + 6) // 7 + 1):
        width = (q // phase_count) * (phase_count - 7 * (size - 1))
        require(width > 0, ("mode width", q, speed, phase_count, size))
        for start in range(phase_count):
            phase_block = tuple(
                (start + offset) % phase_count for offset in range(size)
            )
            block = tuple(
                sheet
                for sheet in range(q)
                if phase_step * sheet % phase_count in phase_block
            )
            centre_2q = (
                -(q // phase_count) * (2 * start + size - 1)
            ) % (2 * q)
            modes.append((block, centre_2q, width, size, phase_block))
    return tuple(modes)


def gap_values(q, left_speed, right_speed, left_mode, right_mode):
    modulus = 2 * q * gcd(left_speed, right_speed)
    residue = (
        left_mode[1] * right_speed - right_mode[1] * left_speed
    ) % modulus
    bound = (
        left_mode[2] * right_speed + right_mode[2] * left_speed - 1
    ) // 7
    first_index = -((bound + residue) // modulus)
    last_index = (bound - residue) // modulus
    return tuple(
        residue + modulus * index
        for index in range(first_index, last_index + 1)
    )


def potential_witness(q, speeds, modes):
    star_sets = tuple(
        gap_values(q, speeds[0], speeds[index], modes[0], modes[index])
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
                    q, speeds[left], speeds[right], modes[left], modes[right]
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


def mode_search(q, subset):
    ordered = tuple(
        sorted(subset, key=lambda speed: (len(owner_modes(q, speed)), speed))
    )
    banks = tuple(owner_modes(q, speed) for speed in ordered)
    max_sizes = tuple(max(len(mode[0]) for mode in bank) for bank in banks)

    def visit(index, assigned, covered):
        if len(covered) + sum(max_sizes[index:]) < q:
            return None
        if index == len(ordered):
            if len(covered) != q:
                return None
            witness = potential_witness(q, ordered, tuple(assigned))
            if witness is None:
                return None
            return ordered, tuple(assigned), witness
        speed = ordered[index]
        for mode in banks[index]:
            if any(
                not gap_values(q, ordered[prior], speed, assigned[prior], mode)
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

    return visit(0, [], set())


def independent(subset, edge_sets):
    chosen = set(subset)
    return not any(edge <= chosen for edge in edge_sets)


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))
    q8 = load_q8()

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

    expected_edges = {8: q8.EXPECTED_MINIMAL_EDGES, **EXPECTED_EDGES}
    expected_profiles = {
        8: q8.EXPECTED_INDEPENDENCE_PROFILE,
        **EXPECTED_PROFILES,
    }

    summaries = []
    canonical_witnesses = []
    for q in range(8, 15):
        vertices = tuple(speed for speed in range(1, 15) if speed % q != 0)
        event_minimal = []
        for size in range(1, SEARCH_MAX_RANK + 1):
            for subset in combinations(vertices, size):
                if any(set(edge) <= set(subset) for edge in event_minimal):
                    continue
                if has_cover(q, subset):
                    event_minimal.append(subset)
        event_minimal = tuple(event_minimal)
        require(event_minimal == expected_edges[q], ("edge list", q, event_minimal))

        witnesses = []
        for edge in event_minimal:
            witness = mode_search(q, edge)
            require(witness is not None, ("missing mode witness", q, edge))
            require(
                set().union(*(set(mode[0]) for mode in witness[1])) == set(range(q)),
                ("mode blocks", q, edge, witness),
            )
            witnesses.append((edge, witness))

        edge_sets = tuple(frozenset(edge) for edge in event_minimal)
        profile = tuple(
            sum(
                independent(subset, edge_sets)
                for subset in combinations(vertices, size)
            )
            for size in range(len(vertices) + 1)
        )
        require(profile == expected_profiles[q], ("profile", q, profile))
        require(profile[5] == EXPECTED_ATLAS_ROWS[q], ("atlas profile", q))

        unsafe_five_sets = tuple(
            subset
            for subset in combinations(vertices, 5)
            if not independent(subset, edge_sets)
        )
        require(
            all(leaks_outside_core(q, subset) for subset in unsafe_five_sets),
            ("core rescue", q),
        )
        rank_profile = tuple(sorted(Counter(map(len, event_minimal)).items()))
        mode_size_profile = tuple(
            sorted(
                Counter(
                    tuple(sorted((len(mode[0]) for mode in witness[1]), reverse=True))
                    for _, witness in witnesses
                ).items()
            )
        )
        summaries.append(
            (
                q,
                vertices,
                event_minimal,
                rank_profile,
                mode_size_profile,
                profile,
                len(unsafe_five_sets),
            )
        )
        canonical_witnesses.append((q, tuple(witnesses)))

    require((2, 6, 10, 14) in expected_edges[8], "q4 to q8 lift")
    require((2, 4, 6, 8, 14) in expected_edges[10], "q5 to q10 lift")
    require((4, 6, 10, 14) in expected_edges[12], "q6 to q12 lift")
    require(
        not expected_edges[11] and not expected_edges[13] and not expected_edges[14],
        "capacity/no-edge controls",
    )

    semantic = ExactDigest()
    semantic.add(("summaries", tuple(summaries)))
    semantic.add(("witnesses", tuple(canonical_witnesses)))
    semantic.add(
        (
            "doubling_lifts",
            (2, 6, 10, 14),
            (2, 4, 6, 8, 14),
            (4, 6, 10, 14),
        )
    )
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("LRC14 Q8-Q14 FINITE-MODE CLUTTER EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=FINITE-EXACT literal q8..14 rank<=5 compiler plus analytic consecutive-mode star-cochain candidate;unnumbered_and_not_canon")
    print("mode_formula=m=q/gcd(u,q);1<=s<=ceil(m/7);block=s_consecutive_phase_classes")
    print("mode_centre=h/(2qu),h=-(q/m)(2r+s-1)_mod_2q;mode_radius=w/(14qu),w=(q/m)(m-7(s-1))")
    print("cochain=p_ij=2q*u_i*u_j(x_i-x_j);congruence_mod_2qgcd;overlap_7|p|<w_i*u_j+w_j*u_i;compatible_star")
    for q, vertices, edges, ranks, mode_sizes, profile, unsafe in summaries:
        print(f"q={q};vertices={vertices};minimal_edges_through_rank5={len(edges)};rank_le5_profile={ranks};mode_size_profile={mode_sizes};I5={profile[5]};unsafe_five_sets={unsafe};core_rescues=0")
        print(f"q={q};minimal_edge_list_through_rank5={edges}")
        print(f"q={q};profile_generated_by_rank_le5_edges={profile}")
    print("doubling_lifts=q4_to_q8:(1,3,5,7)->(2,6,10,14);q5_to_q10:(1,2,3,4,7)->(2,4,6,8,14);q6_to_q12:(2,3,5,7)->(4,6,10,14)")
    print("atlas_identity=q8..14_I5=(1152,1205,1269,1287,1271,1287,1287);all_unsafe_rows_leak_outside_core_clock1")
    print("typing=finite_owner_mode_bank_plus_complete_affine_star;not_tournament")
    print("scope=structural_reconstruction_of_T3387_q8..14_rank_le5_slices;higher_rank_physical_clutters_not_enumerated;no_new_row,no_refined_decrement,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
