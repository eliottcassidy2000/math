#!/usr/bin/env python3
"""Exact audit of reflected high-pair transfer to projected k=2,3.

It reconstructs the historical 113-row/13-family wall left by THM-3351,
which contains the current post-THM-3361 110-row wall, and tests the literal
source typing of the THM-3360/3376 reflected-pair theorems.
It also tests the cheapest plausible surrogate: call the raw drift pair
(z1,H) high by its reduced ratio and measure it on the upper-median body-safe
cell.  Finally it freezes a same-denominator/same-unit height witness showing
that a denominator passport cannot determine a located pair overlap.
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE_ROOT
PREFIX = ROOT / "04-computation/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.py"
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
CHROMATIC = ROOT / "04-computation/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"

EXPECTED_HASHES = {
    PREFIX: "cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62",
    MEDIAN: "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276",
    BASE: "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31",
    CHROMATIC: "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8",
}

LOW_RATIOS = frozenset((F(4, 3), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6)))
ORIENTED_LOW_RATIOS = tuple(sorted(LOW_RATIOS | {1 / value for value in LOW_RATIOS}))


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def intersection_mass(first, second):
    i = j = 0
    total = F()
    while i < len(first) and j < len(second):
        total += max(F(), min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def low_connected(row):
    unseen = set(row[1:])
    frontier = {row[0]}
    while frontier:
        value = frontier.pop()
        hit = {
            other for other in unseen
            if max(value, other) / min(value, other) in LOW_RATIOS
        }
        unseen -= hit
        frontier |= hit
    return not unseen


def canonical_parent(row):
    for value in reversed(row[1:]):
        candidate = tuple(x for x in row if x != value)
        if low_connected(candidate):
            return candidate
    raise RuntimeError(row)


def primitive(row):
    denominator = lcm(*(value.denominator for value in row))
    integers = tuple(int(value * denominator) for value in row)
    common = gcd(*integers)
    return tuple(value // common for value in integers)


def connected_five_atlas():
    layer = [(F(1),)]
    for _size in range(2, 6):
        out = []
        for row in layer:
            have = set(row)
            children = {
                value * ratio
                for value in row
                for ratio in ORIENTED_LOW_RATIOS
                if value * ratio >= 1 and value * ratio not in have
            }
            for child in sorted(children):
                candidate = tuple(sorted((*row, child)))
                if canonical_parent(candidate) == row:
                    out.append(candidate)
        require(len(out) == len(set(out)), ("atlas duplicate", _size))
        layer = out
    answer = tuple(sorted(primitive(row) for row in layer))
    require(len(answer) == len(set(answer)) == 19389, len(answer))
    return answer


def high_edge_count(row):
    count = 0
    for left, right in combinations(row, 2):
        common = gcd(left, right)
        count += left // common + right // common >= 8
    return count


def main():
    for path, expected in EXPECTED_HASHES.items():
        require(sha(path) == expected, (path, "dependency changed", sha(path)))

    prefix = load("projected_transfer_prefix", PREFIX)
    median = load("projected_transfer_median", MEDIAN)
    base = load("projected_transfer_base", BASE)
    chromatic = load("projected_transfer_chromatic", CHROMATIC)
    rows, _components, _inherited, ranked, _prefix, _indices, _hash = prefix.reconstruct_queue()
    families = tuple(ranked[prefix.PREFIX_FAMILIES :])
    indices = tuple(index for family in families for index in family[4])
    require((len(families), len(indices)) == (13, 113), (len(families), len(indices)))
    require(len({rows[index][0] for index in indices}) == 113, "body repetition")
    thm3361_closed = (191, 228, 332)
    require(all(index in indices for index in thm3361_closed), thm3361_closed)
    current_indices = tuple(index for index in indices if index not in thm3361_closed)
    require(len(current_indices) == 110, len(current_indices))

    reflected_bodies = {body for body, _ruler in median.body_universe()}
    body_intersection = tuple(index for index in indices if rows[index][0] in reflected_bodies)
    robust_hist = Counter(len(median.LP.robust_edges(rows[index][0])[2]) for index in indices)
    require(not body_intersection, body_intersection)
    require(robust_hist == Counter({15: 113}), robust_hist)
    active_bodies = {body for body, _ruler in median.body_universe()}
    same_level_exceptions = {row[0] for row in chromatic.EXPECTED_EXCEPTIONS}
    require(active_bodies.isdisjoint(same_level_exceptions), active_bodies & same_level_exceptions)

    # A reflected drift has z=qL-e with q>=1 and e one of the six body labels.
    # Every residual row begins with z1=216, so test the source type literally.
    reflected_first = []
    for index in indices:
        body, ruler, _high, _wall = rows[index]
        witnesses = tuple(
            (e, (216 + e) // ruler)
            for e in body
            if (216 + e) % ruler == 0 and (216 + e) // ruler >= 1
        )
        if witnesses:
            reflected_first.append((index, witnesses))
    require(not reflected_first, reflected_first)
    ruler_range = (min(rows[index][1] for index in indices), max(rows[index][1] for index in indices))
    require(ruler_range == (129360, 5045040), ruler_range)

    # Hostile surrogate: treat the raw pair (216,H) as high using its reduced
    # integer ratio and evaluate it on the upper-median body-safe cell.
    records = []
    for index in indices:
        body, ruler, high, wall = rows[index]
        require(wall and high == (13 * ruler) // 132 + 1, (index, high))
        actual_ruler, ranges = base.safe_cell_ranges(body)
        require(actual_ruler == ruler, (index, actual_ruler, ruler))
        cells = tuple(cell for left, right in ranges for cell in range(left, right))
        cell = cells[len(cells) // 2]
        require(base.body_cell_is_safe(ruler, body, cell), (index, "unsafe median"))
        first_arcs = base.direct_multiplier_arcs(ruler, 216, cell)
        high_arcs = base.direct_multiplier_arcs(ruler, high, cell)
        first_mass = base.interval_mass(base.merge_intervals(first_arcs))
        high_mass = base.interval_mass(base.merge_intervals(high_arcs))
        overlap = intersection_mass(first_arcs, high_arcs)
        common = gcd(216, high)
        reduced = (216 // common, high // common)
        require(sum(reduced) >= 8, (index, reduced))
        records.append((index, body, ruler, high, cell, reduced, first_mass, high_mass, overlap))
    records = tuple(records)
    overlap_hist = Counter(record[-1] for record in records)
    require(overlap_hist == Counter({F(0): 108, F(1): 5}), overlap_hist)
    mass_type_hist = Counter(
        (
            "empty" if record[6] == 0 else "full" if record[6] == 1 else "partial",
            "empty" if record[7] == 0 else "full" if record[7] == 1 else "partial",
            "zero" if record[8] == 0 else "positive",
        )
        for record in records
    )
    require(
        mass_type_hist
        == Counter({
            ("empty", "empty", "zero"): 69,
            ("empty", "partial", "zero"): 13,
            ("full", "empty", "zero"): 26,
            ("full", "full", "positive"): 5,
        }),
        mass_type_hist,
    )

    witness = next(record for record in records if record[0] == 94)
    require(
        witness
        == (94, (1, 3, 8, 10, 11, 14), 129360, 12741, 65835,
            (72, 4247), F(1), F(0), F(0)),
        witness,
    )

    # Same denominator and the same residue/unit ray, but a different height.
    # The local overlap changes, so even (denominator,unit) is insufficient.
    index, body, ruler, high, cell, _reduced, _m1, _m2, _overlap = witness
    shifted = high + ruler
    require(ruler // gcd(ruler, high) == ruler // gcd(ruler, shifted) == 43120, "denominator")
    require(high % ruler == shifted % ruler, "unit/residue")
    first_arcs = base.direct_multiplier_arcs(ruler, 216, cell)
    high_arcs = base.direct_multiplier_arcs(ruler, high, cell)
    shifted_arcs = base.direct_multiplier_arcs(ruler, shifted, cell)
    height_pair = (
        intersection_mass(first_arcs, high_arcs),
        intersection_mass(first_arcs, shifted_arcs),
        base.interval_mass(base.merge_intervals(shifted_arcs)),
        shifted_arcs,
    )
    require(
        height_pair
        == (F(0), F(6160, 47367), F(6160, 47367),
            ((F(9625, 15789), F(35035, 47367)),)),
        height_pair,
    )

    # Positive scoped transfer after restoring the reflected affine type.
    # On a body-safe cell a genuine reflected clause z=qL-e has mass
    # 1/7+e/[7(qL-e)].  Four distinct residues are maximized independently at
    # q=1.  Five distinct levels are maximized at q=1,...,5, with the larger
    # body labels assigned to the smaller q (the two-term exchange function
    # x/((qL-x)(rL-x)) is strictly increasing in x).
    d4_max = None
    d5_distinct_max = None
    for body in combinations(range(1, 15), 6):
        ruler = 14 * lcm(*body)
        d4_values = sorted((F(e, 7 * (ruler - e)) for e in body), reverse=True)
        d4 = sum(d4_values[:4], F())
        d4_row = (d4, body, ruler, tuple(d4_values[:4]))
        if d4_max is None or d4_row > d4_max:
            d4_max = d4_row
        for labels in combinations(body, 5):
            ordered = tuple(sorted(labels, reverse=True))
            debt = sum((F(e, 7 * (q * ruler - e)) for q, e in enumerate(ordered, 1)), F())
            row = (debt, body, ruler, ordered)
            if d5_distinct_max is None or row > d5_distinct_max:
                d5_distinct_max = row
    require(
        d4_max
        == (F(123896, 5540535), (1, 2, 3, 4, 6, 12), 168,
            (F(1, 91), F(1, 189), F(1, 287), F(1, 385))),
        d4_max,
    )
    require(
        d5_distinct_max
        == (F(183680141, 11691304625), (1, 2, 3, 4, 6, 12), 168,
            (12, 6, 4, 3, 2)),
        d5_distinct_max,
    )
    k3_margin = F(3, 91) - d4_max[0]
    dmax = F(186636088362, 11773143757375)
    pair_target = dmax / 5
    k2_two_edge_margin = F(1, 91) + 2 * pair_target - d5_distinct_max[0]
    require(k3_margin == F(58759, 5540535) > 0, k3_margin)
    require(k2_two_edge_margin == F(95318697414, 58865718786875) > 0,
            k2_two_edge_margin)

    # For five distinct reflected levels, >=2 high edges close by the preceding
    # two-edge gate.  The <=1-high connected atlas consists of only 16 rays.
    atlas5 = connected_five_atlas()
    high_hist = Counter(high_edge_count(row) for row in atlas5)
    require(
        tuple(sorted(high_hist.items()))
        == ((0, 2), (1, 14), (2, 82), (3, 363), (4, 1400), (5, 4969), (6, 12559)),
        high_hist,
    )
    rare = tuple(row for row in atlas5 if high_edge_count(row) <= 1)
    require(len(rare) == 16, len(rare))
    simple_gate = []
    bodies649 = median.body_universe()
    for scale in (1, 2):
        for shape in rare:
            high_count = high_edge_count(shape)
            best = None
            for body, ruler in bodies649:
                labels = tuple(sorted(body, reverse=True)[:5])
                debt = sum(
                    (F(e, 7 * (scale * q * ruler - e)) for e, q in zip(labels, shape)),
                    F(),
                )
                candidate = (debt, body, ruler, labels)
                if best is None or candidate > best:
                    best = candidate
            margin = F(1, 91) + high_count * pair_target - best[0]
            simple_gate.append((scale, shape, high_count, margin, best))
    simple_gate = tuple(simple_gate)
    failed_scale1 = tuple(row[1] for row in simple_gate if row[0] == 1 and row[3] <= 0)
    require(
        failed_scale1
        == ((1, 2, 3, 4, 6), (1, 2, 3, 4, 12), (1, 2, 3, 6, 12),
            (1, 2, 4, 6, 8), (1, 2, 4, 6, 12)),
        failed_scale1,
    )
    require(sum(row[0] == 1 and row[3] > 0 for row in simple_gate) == 11,
            "scale-one gate")
    require(all(row[3] > 0 for row in simple_gate if row[0] == 2), "scale-two gate")
    scale2_min = min((row for row in simple_gate if row[0] == 2), key=lambda row: row[3])

    family_packet = tuple((f[2], f[3], len(f[4]), f[4]) for f in families)
    print("LRC14 REFLECTED-PAIR PROJECTED TRANSFER EXACT AUDIT")
    print("scope=historical_THM3351_113_row_superset;post_THM3361_current_rows", len(current_indices))
    print("universe_families_rows", len(families), len(indices), "family_sha256", digest(family_packet))
    print("ruler_range", ruler_range)
    print("THM3360_3376_body_domain_intersection", len(body_intersection), "robust_edge_hist", tuple(sorted(robust_hist.items())))
    print("reflected_first_drift_representations", len(reflected_first))
    print("raw_high_pair_rows", len(records), "upper_median_overlap_hist", tuple(sorted(overlap_hist.items())))
    print("upper_median_mass_type_hist", tuple(sorted(mass_type_hist.items())))
    print("zero_overlap_witness", witness)
    print("same_denominator_unit_height_witness", (high, shifted, 43120, *height_pair[:3]))
    print("reflected_k3_four_residue_debt_max", d4_max, "projected_margin", k3_margin)
    print("reflected_k2_five_distinct_debt_max", d5_distinct_max,
          "two_high_edge_margin", k2_two_edge_margin)
    print("reflected_k2_connected_five_atlas", len(atlas5), "high_edge_hist", tuple(sorted(high_hist.items())))
    print("reflected_k2_low_high_edge_rays", len(rare), "scale1_simple_pass", 11,
          "scale1_heads", failed_scale1, "scale2_all_pass", True)
    print("reflected_k2_scale2_min_record", scale2_min)
    print("record_sha256", digest(records))
    print("conclusion=literal_transfer_refuted;body_and_affine_source_types_fail;raw_high_pair_has_no_positive_located_floor;denominator_unit_quotient_loses_height;scoped_reflected_k3_closes_by_singletons;scoped_reflected_k2_reduces_to_five_scale1_heads")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
