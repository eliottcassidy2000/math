#!/usr/bin/env python3
"""Exact global configuration-switching cap for THM-2672.

Keep the common future coordinate, its half-digit, and the slope-seven carry
orbit while allowing the rail, base cell, present factor, clock, sector, and
private edge to vary independently from target label to target label.  The
computed envelope retains exact delayed-sector support, literal private
half-tooth geometry/root, and the THM-2640 unit flag, but existentially unions
the forgotten coordinates.  Consequently every physical chart maps into this
envelope, while the converse need not hold.

The envelope covers at most twelve of the thirteen carry values on every open
future-coordinate cell.  Since ``delta -> source_carry + 7*delta (mod 13)`` is
a permutation, even arbitrary labelwise configuration switching cannot cover
all thirteen target labels at one common point.  Endpoint conventions cannot
create an open thirteen-fold component; the elementary endpoints themselves
are inherited chart boundaries rather than positive support components.
"""

from bisect import bisect_right
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_predecessor_carry_private_root_atlas_thm2640 as m
import lrc14_slope7_fixed_configuration_carry_nerve_thm2672 as fixed


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_intervals(intervals):
    out = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def prefix_intervals(prefix):
    starts, lengths, _ = prefix
    return tuple((left, left + length)
                 for left, length in zip(starts, lengths))


def contains_doubled(intervals, doubled_point):
    """Strict midpoint membership without constructing a Fraction."""
    starts = [left for left, _ in intervals]
    index = bisect_right(starts, doubled_point // 2) - 1
    return (index >= 0
            and 2 * intervals[index][0] < doubled_point
            < 2 * intervals[index][1])


def build_flags_and_metadata():
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(m.shard, m.core.SHARDS))
    content = 0
    metadata = []
    rows = []
    for _, shard_content, _, _, _, shard_metadata, shard_rows in results:
        content = gcd(content, shard_content)
        metadata.extend(shard_metadata)
        rows.extend(shard_rows)
    require(content == 26 and len(rows) == 162,
            "THM-2640 carrier reconstruction changed")
    return content, tuple(metadata), rows, fixed.build_flags(rows, content)


def build_sector_support(prefixes):
    support = {}
    for sector in range(2):
        for h in range(m.P):
            for kappa in range(2):
                support[sector, h, kappa] = merge_intervals(
                    interval
                    for ell5 in range(m.Q7)
                    for interval in prefix_intervals(
                        prefixes[sector][ell5][h][kappa]
                    )
                )
    return support


def edge_geometrically_allowed(h, kappa, edge, doubled_y):
    """Exact private-half membership on one future half-digit interior.

    Write d=2h+kappa=13b+a and y=(d+v)/26.  Edge 1 is available
    throughout unless a=12, when it needs v<1/14.  Edge 0 is available
    throughout unless a=0, when it needs v>13/14.  Endpoint choices are
    immaterial for positive-length support.
    """
    digit = 2 * h + kappa
    residue = digit % m.P
    left = digit * m.T // (2 * m.P)
    right = (digit + 1) * m.T // (2 * m.P)
    require(2 * left < doubled_y < 2 * right,
            "edge test left its future half-digit")
    if edge == 0 and residue == 0:
        threshold = left + 13 * (right - left) // 14
        return doubled_y > 2 * threshold
    if edge == 1 and residue == 12:
        threshold = left + (right - left) // 14
        return doubled_y < 2 * threshold
    return True


def audit_edge_formula():
    """Compare the two-threshold formula with the literal circular teeth."""
    checks = 0
    for h in range(m.P):
        for kappa in range(2):
            digit = 2 * h + kappa
            left = digit * m.T // (2 * m.P)
            right = (digit + 1) * m.T // (2 * m.P)
            cuts = [left, right]
            if digit % m.P == 0:
                cuts.insert(1, left + 13 * (right - left) // 14)
            if digit % m.P == 12:
                cuts.insert(1, left + (right - left) // 14)
            for a, b in zip(cuts, cuts[1:]):
                doubled_y = a + b
                y = Fraction(doubled_y, 2 * m.T)
                base_half = digit // m.P
                for carry in range(m.P):
                    deep_q = Fraction(28) * (carry + y) % 182
                    for edge in range(2):
                        root = (2 * carry + base_half
                                + (edge == 0)) % m.P
                        if edge == 0:
                            lo, hi = ((14 * root - 13, 14 * root)
                                      if root else (169, 182))
                        else:
                            lo, hi = ((14 * root, 14 * root + 13)
                                      if root else (0, 13))
                        literal = Fraction(lo) < deep_q < Fraction(hi)
                        formula = edge_geometrically_allowed(
                            h, kappa, edge, doubled_y
                        )
                        require(literal == formula,
                                "private-edge threshold formula failed")
                        checks += 1
    return checks


def fine_coordinate_envelope(flags, prefixes):
    sector_support = build_sector_support(prefixes)
    unit_exists = {}
    unit_multiplicity = {}
    for sector in range(2):
        for edge in range(2):
            for carry in range(m.P):
                for kappa in range(2):
                    for h in range(m.P):
                        count = sum(
                            flags[j][sector][edge][carry][kappa][h]
                            for j in range(len(flags))
                        )
                        unit_exists[sector, edge, carry, kappa, h] = bool(count)
                        unit_multiplicity[sector, edge, carry, kappa, h] = count

    candidates = []
    state_length = Counter()
    state_components = Counter()
    coverage_length_hist = Counter()
    coverage_component_hist = Counter()
    maximum_coverage = 0
    maximal_components = []
    for h in range(m.P):
        for kappa in range(2):
            digit = 2 * h + kappa
            left = digit * m.T // (2 * m.P)
            right = (digit + 1) * m.T // (2 * m.P)
            endpoints = {left, right}
            for sector in range(2):
                for a, b in sector_support[sector, h, kappa]:
                    if left < a < right:
                        endpoints.add(a)
                    if left < b < right:
                        endpoints.add(b)
            residue = digit % m.P
            if residue == 0:
                endpoints.add(left + 13 * (right - left) // 14)
            if residue == 12:
                endpoints.add(left + (right - left) // 14)
            endpoints = sorted(endpoints)
            for a, b in zip(endpoints, endpoints[1:]):
                if a >= b:
                    continue
                doubled_y = a + b
                active_sectors = tuple(
                    sector for sector in range(2)
                    if contains_doubled(
                        sector_support[sector, h, kappa], doubled_y
                    )
                )
                allowed_edges = tuple(
                    edge for edge in range(2)
                    if edge_geometrically_allowed(
                        h, kappa, edge, doubled_y
                    )
                )
                carry_options = tuple(
                    tuple(
                        (sector, edge)
                        for sector in active_sectors
                        for edge in allowed_edges
                        if unit_exists[sector, edge, carry, kappa, h]
                    )
                    for carry in range(m.P)
                )
                covered_carries = tuple(
                    carry for carry, options in enumerate(carry_options)
                    if options
                )
                missing_carries = tuple(
                    carry for carry, options in enumerate(carry_options)
                    if not options
                )
                coverage = len(covered_carries)
                state = (
                    h, kappa, active_sectors, allowed_edges,
                    covered_carries, missing_carries,
                )
                state_length[state] += b - a
                state_components[state] += 1
                coverage_length_hist[coverage] += b - a
                coverage_component_hist[coverage] += 1
                if coverage > maximum_coverage:
                    maximum_coverage = coverage
                    maximal_components = []
                if coverage == maximum_coverage:
                    maximal_components.append((
                        h, kappa, a, b, active_sectors, allowed_edges,
                        carry_options, missing_carries,
                    ))
                if coverage == m.P:
                    # Every source carry merely permutes these thirteen
                    # carry options among the target labels.
                    for source_carry in range(m.P):
                        label_options = tuple(
                            carry_options[(source_carry + 7 * delta) % m.P]
                            for delta in range(m.P)
                        )
                        candidates.append((
                            source_carry, h, kappa, a, b,
                            active_sectors, allowed_edges, label_options,
                        ))

    total_length = sum(b - a for _, _, _, a, b, *_ in candidates)
    maximal_missing_hist = Counter(
        component[-1] for component in maximal_components
    )
    maximal_missing_label_hist = Counter(
        tuple(
            2 * (carry - source_carry) % m.P
            for carry in component[-1]
        )
        for component in maximal_components
        for source_carry in range(m.P)
    )
    nonzero_state_summary = tuple(sorted(
        (
            len(covered), h, kappa, sectors, edges, missing,
            state_components[state], Fraction(length, m.T),
        )
        for state, length in state_length.items()
        for h, kappa, sectors, edges, covered, missing in (state,)
        if covered
    ))
    return {
        "sector_support": sector_support,
        "unit_multiplicity": unit_multiplicity,
        "candidates": tuple(candidates),
        "state_length": state_length,
        "state_components": state_components,
        "coverage_length_hist": coverage_length_hist,
        "coverage_component_hist": coverage_component_hist,
        "maximum_coverage": maximum_coverage,
        "maximal_components": tuple(maximal_components),
        "maximal_missing_hist": maximal_missing_hist,
        "maximal_missing_label_hist": maximal_missing_label_hist,
        "nonzero_state_summary": nonzero_state_summary,
        "total_length": total_length,
    }


def main():
    edge_formula_checks = audit_edge_formula()
    require(edge_formula_checks == 780,
            "literal private-edge audit universe changed")
    content, metadata, rows, flags = build_flags_and_metadata()
    module = m.old.load_present_module()
    prefixes = m.build_pair_prefixes(module)
    fine = fine_coordinate_envelope(flags, prefixes)
    candidates = fine["candidates"]
    require(fine["maximum_coverage"] == 12 and not candidates,
            "fine-coordinate global cap changed")
    require(fine["coverage_component_hist"]
            == Counter({0: 47_512, 11: 39_948, 12: 7_536}),
            "fine-coordinate elementary-cell census changed")
    require(fine["coverage_length_hist"] == Counter({
        0: 282_114_031_894_006,
        11: 13_227_102_414_396,
        12: 2_495_763_530_078,
    }), "fine-coordinate exact measure census changed")
    require(sum(fine["coverage_length_hist"].values()) == m.T,
            "future-coordinate partition lost length")
    require(fine["maximal_missing_hist"]
            == Counter({(0,): 2067, (6,): 3402, (12,): 2067}),
            "maximal missing-carry census changed")
    require(fine["maximal_missing_label_hist"]
            == Counter({(delta,): 7536 for delta in range(m.P)}),
            "source-carry rebase stopped distributing every missing label")

    # Positive hostile: the canonical physical twelve-chart component of
    # THM-2672 lies strictly inside this relaxed envelope with missing carry
    # six.  Thus the cap is thirteen-specific, not an empty-chart artifact.
    known_y_grid = Fraction(138_281_588_745_300)
    known_controls = tuple(
        component for component in fine["maximal_components"]
        if (component[0], component[1], component[-1]) == (6, 0, (6,))
        and component[2] < known_y_grid < component[3]
    )
    require(len(known_controls) == 1,
            "canonical twelve-chart midpoint left the relaxed envelope")
    print("LRC14 slope-seven global configuration-switching exact probe")
    print(f"carrier=THM2640 rails:{len(metadata)} content:{content}")
    print(f"literal_private_edge_formula_checks={edge_formula_checks}")
    print("fine_filter=retain source carry,future half-digit,delayed sector,"
          "private edge/root,and exact unit flag; forget rail/base-cell/present")
    print(f"fine_maximum_simultaneously_covered_carries="
          f"{fine['maximum_coverage']}/13")
    print(f"fine_candidate_components={len(candidates)}")
    print(f"fine_candidate_total_y_measure_with_carry_sum="
          f"{Fraction(fine['total_length'], m.T)}")
    print("fine_coverage_elementary_interval_hist="
          f"{tuple(sorted(fine['coverage_component_hist'].items()))}")
    print("fine_coverage_y_measure_hist="
          f"{tuple((coverage, Fraction(length,m.T)) for coverage,length in sorted(fine['coverage_length_hist'].items()))}")
    print("maximal_missing_carry_component_hist="
          f"{tuple(sorted(fine['maximal_missing_hist'].items()))}")
    print("maximal_missing_label_component_hist_over_source_carries="
          f"{tuple(sorted(fine['maximal_missing_label_hist'].items()))}")
    print(f"nonzero_fine_state_summary={fine['nonzero_state_summary']}")
    print("known_twelve_chart_midpoint_control="
          f"(yT={known_y_grid},state={known_controls[0][0:6] + known_controls[0][-1:]})")
    print("scope=necessary fine-coordinate envelope; forgotten coordinates "
          "only enlarge support, so cap 12 proves the physical global "
          "configuration-switched positive open thirteen-fold intersection "
          "is empty")


if __name__ == "__main__":
    main()
