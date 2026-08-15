#!/usr/bin/env python3
"""Independent exact full-clutter audit for LRC sheet degrees q=8,...,15.

The earlier finite-mode probes intentionally stopped at the body-atlas rank
cutoffs (five for q<=14 and six for q=15).  This audit uses no owner modes or
cochain formulas.  It constructs every rational boundary and intervening
open cell, stores each speed/sheet danger set as a sample bitset, and tests
the complete 2^13 or 2^14 literal subset universe.

This is an unnumbered FINITE-EXACT artifact.  Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import lcm
from pathlib import Path


EXPECTED_BOUNDARY_COUNTS = {
    8: 945,
    9: 1441,
    10: 1301,
    11: 2025,
    12: 1249,
    13: 2341,
    14: 756,
    15: 2101,
}

EXPECTED_RANK_PROFILES = {
    8: ((4, 15), (5, 17), (6, 6)),
    9: ((4, 9), (5, 13), (6, 54), (7, 2)),
    10: ((5, 18), (6, 70), (7, 4)),
    11: ((6, 23), (7, 91), (8, 9)),
    12: ((4, 1), (5, 7), (6, 17), (7, 14)),
    13: ((7, 22), (8, 29), (9, 2)),
    14: ((7, 3), (8, 60), (9, 9)),
    15: ((6, 157), (7, 16), (8, 6)),
}

EXPECTED_INDEPENDENCE_PROFILES = {
    8: (1, 13, 78, 286, 700, 1152, 1217, 762, 254, 39, 2, 0, 0, 0),
    9: (1, 13, 78, 286, 706, 1205, 1347, 884, 292, 42, 2, 0, 0, 0),
    10: (1, 13, 78, 286, 715, 1269, 1532, 1115, 411, 59, 1, 0, 0, 0),
    11: (1, 13, 78, 286, 715, 1287, 1693, 1520, 816, 224, 26, 1, 0, 0),
    12: (1, 13, 78, 286, 714, 1271, 1612, 1382, 720, 202, 29, 2, 0, 0),
    13: (1, 13, 78, 286, 715, 1287, 1716, 1694, 1166, 488, 106, 10, 0, 0),
    14: (1, 13, 78, 286, 715, 1287, 1716, 1713, 1212, 540, 133, 17, 1, 0),
    15: (1, 14, 91, 364, 1001, 2002, 2846, 2768, 1743, 642, 117, 8, 0, 0, 0),
}

EXPECTED_EDGE_DIGESTS = {
    8: "a129c5ef4a5babaa8e4c9e1e986cdcce33d4b9d174881492efcbca031c048002",
    9: "1e0321520c2360be3f0e0ffdd887e37e4bcb7d8fc0b11f23e68cd5fe83ed57c5",
    10: "d14ad83ff1ff83b55ccca9d4ac803190cf496629a568c589f1e1f10b56671916",
    11: "dc7112501a40296e447326129efadfb1ff8021bc78f87b7cc40fcf5520540076",
    12: "8479a2e2d0788372aaeb53886bc172378fc13d8b86bff72cf80ef463d847d8bf",
    13: "0915d5d6270a5a661080efc8952380d5064f30db2ca5c9ae5c9a36d61a4aba08",
    14: "09c5990ad50c0d0b7c133b16d10f1c72cb1bc85eb9f7754a48bd4073fce6acbf",
    15: "3980d53a1bd43516da447e9ff52f91a423c1cca4e6904cd1e2d8b249e8c39c12",
}

EXPECTED_LOW_RANK_PROFILES = {
    8: ((4, 15), (5, 17)),
    9: ((4, 9), (5, 13)),
    10: ((5, 18),),
    11: (),
    12: ((4, 1), (5, 7)),
    13: (),
    14: (),
    15: ((6, 157),),
}

EXPECTED_SEMANTIC_DIGEST = "2baecf15b57ac7adcfaf41342d74cf12eba363a6d0726566e67eabfcf78b0d22"


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


def strict_danger(q, speed, sample, scale, sheet):
    denominator = 2 * q * scale
    residue = speed * (q * sample + 2 * scale * sheet) % denominator
    return 14 * min(residue, denominator - residue) < denominator


def exact_samples(q, vertices):
    scale = 14 * q * lcm(*vertices)
    boundaries = {0}
    for speed in vertices:
        for sheet in range(q):
            for tooth in range(speed):
                base = scale * tooth // speed - scale * sheet // q
                delta = scale // (14 * speed)
                boundaries.add((base - delta) % scale)
                boundaries.add((base + delta) % scale)
    ordered = tuple(sorted(boundaries))
    endpoints = tuple(2 * boundary for boundary in ordered)
    midpoints = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += scale
        midpoints.append((left + right) % (2 * scale))
    samples = endpoints + tuple(midpoints)
    require(len(set(samples)) == 2 * len(ordered), (q, "sample collision"))
    return scale, ordered, samples


def danger_bitsets(q, vertices, scale, samples):
    bank = []
    for speed in vertices:
        sheets = []
        for sheet in range(q):
            bits = 0
            for index, sample in enumerate(samples):
                if strict_danger(q, speed, sample, scale, sheet):
                    bits |= 1 << index
            sheets.append(bits)
        bank.append(tuple(sheets))
    return tuple(bank)


def physically_covers(q, subset_indices, bank, sample_universe):
    simultaneous = sample_universe
    for sheet in range(q):
        blocked = 0
        for index in subset_indices:
            blocked |= bank[index][sheet]
        simultaneous &= blocked
        if not simultaneous:
            return False
    return True


def full_clutter(q):
    vertices = tuple(speed for speed in range(1, 15) if speed % q != 0)
    scale, boundaries, samples = exact_samples(q, vertices)
    bank = danger_bitsets(q, vertices, scale, samples)
    sample_universe = (1 << len(samples)) - 1

    minimal_masks = []
    minimal_edges = []
    for size in range(1, len(vertices) + 1):
        for subset_indices in combinations(range(len(vertices)), size):
            mask = sum(1 << index for index in subset_indices)
            if any(edge_mask & mask == edge_mask for edge_mask in minimal_masks):
                continue
            if physically_covers(q, subset_indices, bank, sample_universe):
                minimal_masks.append(mask)
                minimal_edges.append(tuple(vertices[index] for index in subset_indices))

    minimal_edges = tuple(minimal_edges)
    rank_profile = tuple(sorted(Counter(map(len, minimal_edges)).items()))
    independence_profile = []
    for size in range(len(vertices) + 1):
        count = 0
        for subset_indices in combinations(range(len(vertices)), size):
            mask = sum(1 << index for index in subset_indices)
            if not any(edge_mask & mask == edge_mask for edge_mask in minimal_masks):
                count += 1
        independence_profile.append(count)
    independence_profile = tuple(independence_profile)

    edge_digest = sha256(repr(minimal_edges).encode("ascii")).hexdigest()
    cutoff = 6 if q == 15 else 5
    low_rank_profile = tuple((rank, count) for rank, count in rank_profile if rank <= cutoff)
    omitted = tuple(edge for edge in minimal_edges if len(edge) > cutoff)

    require(len(boundaries) == EXPECTED_BOUNDARY_COUNTS[q], (q, "boundaries", len(boundaries)))
    require(rank_profile == EXPECTED_RANK_PROFILES[q], (q, "ranks", rank_profile))
    require(independence_profile == EXPECTED_INDEPENDENCE_PROFILES[q], (q, "profile", independence_profile))
    require(edge_digest == EXPECTED_EDGE_DIGESTS[q], (q, "edge digest", edge_digest))
    require(low_rank_profile == EXPECTED_LOW_RANK_PROFILES[q], (q, "low ranks", low_rank_profile))
    require(omitted, (q, "missing hostile higher-rank edge"))

    for edge in minimal_edges:
        indices = tuple(vertices.index(speed) for speed in edge)
        require(
            all(
                not physically_covers(
                    q,
                    indices[:removed] + indices[removed + 1 :],
                    bank,
                    sample_universe,
                )
                for removed in range(len(indices))
            ),
            (q, "nonminimal edge", edge),
        )

    return (
        q,
        vertices,
        len(boundaries),
        len(samples),
        rank_profile,
        independence_profile,
        edge_digest,
        omitted[0],
    )


def main():
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

    records = tuple(full_clutter(q) for q in range(8, 16))
    require(records[0][7] == (1, 3, 5, 11, 13, 14), records[0][7])
    require(records[-1][7] == (1, 2, 3, 5, 7, 8, 13), records[-1][7])

    semantic = ExactDigest()
    semantic.add(("full_physical_clutters", records))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("LRC Q8-Q15 FULL PHYSICAL CLUTTER INDEPENDENT EXACT AUDIT")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=FINITE-EXACT full literal subset universes;independent event-cell route;unnumbered_and_not_canon")
    print("route=all_rational_boundaries_plus_all_open_cell_midpoints;speed_sheet_sample_bitsets;no_mode_or_cochain_import")
    for q, vertices, boundary_count, sample_count, ranks, profile, edge_digest, first_omitted in records:
        print(f"q={q};vertices={vertices};boundaries={boundary_count};samples={sample_count};true_minimal_rank_profile={ranks};edge_sha256={edge_digest}")
        print(f"q={q};true_independence_profile={profile};first_edge_above_old_cutoff={first_omitted}")
    print("repair=old_q8_q14_profiles_are_generated_by_rank_le5_edges;old_q15_profile_is_generated_by_rank_le6_edges;their_low_rank_counts_and_I5_values_survive")
    print("hostile=q8_first_omitted_edge_(1,3,5,11,13,14);q15_first_omitted_edge_(1,2,3,5,7,8,13)")
    print("scope=full_physical_clutters_on_literal_speed_pools_only;no_core_rescue_change;no_refined_decrement;no_LRC14_ledger_change")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
