#!/usr/bin/env python3
"""Exact referee for THM-2672's carry-rebase facet torsor.

The selected THM-2640 configuration has one nonunit predecessor carry.  For
each of the thirteen possible source carries this script intersects the twelve
unit slope-seven pullbacks honestly.  It then distinguishes the coarse nerve
obtained by forgetting the source carry from a carry-refined selected
subcomplex.  Full connected-arc gain refinement in the sense of THM-2658 is
not computed here.
"""

from fractions import Fraction
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_predecessor_carry_private_root_atlas_thm2640 as m


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def shift_union(intervals, shift, modulus):
    """Pull back a circular interval union under x -> x+shift."""
    out = []
    for left, right in intervals:
        start = (left - shift) % modulus
        stop = start + right - left
        if stop <= modulus:
            out.append((start, stop))
        else:
            out.append((start, modulus))
            out.append((0, stop - modulus))
    out.sort()
    merged = []
    for left, right in out:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def shift_weighted(pieces, shift, modulus):
    """Pull back a weighted circular partition under x -> x+shift."""
    out = []
    for left, right, weight in pieces:
        start = (left - shift) % modulus
        stop = start + right - left
        if stop <= modulus:
            out.append((start, stop, weight))
        else:
            out.append((start, modulus, weight))
            out.append((0, stop - modulus, weight))
    return sorted(out)


def intersect_root_half(pieces, speed, edge, root):
    """Intersect the pulled-back private half, including root-zero wrap."""
    if root:
        lo, hi = ((14 * root - 13, 14 * root)
                  if edge == 0 else (14 * root, 14 * root + 13))
    elif edge == 0:
        lo, hi = 169, 182
    else:
        lo, hi = 0, 13
    return m.old.intersect_weighted_comb(pieces, speed, 182, lo, hi)


def first_delayed_component(base, prefix, carry):
    """Return the lexicographically first exact open delayed component."""
    starts, lengths, _ = prefix
    period = 13 * m.T
    best = None
    for left, right, weight in base:
        if not weight:
            continue
        left_scaled = left * m.R
        right_scaled = right * m.R
        for start, length in zip(starts, lengths):
            stop = start + length
            a = carry * m.T + start
            b = carry * m.T + stop
            k = (left_scaled - b) // period + 1
            while k * period + a < right_scaled:
                lo = max(left_scaled, k * period + a)
                hi = min(right_scaled, k * period + b)
                if lo < hi:
                    candidate = (lo, hi)
                    if best is None or candidate < best:
                        best = candidate
                    break
                k += 1
    return best


def main():
    # The first exact twelve-chart configuration from the THM-2672 sweep.
    j = 0
    sector, edge, kappa, h = 0, 0, 0, 6
    result = m.shard((j, j + 1))
    metadata, rows = result[5], result[6]
    require(result[1] == 26, "selected rail primitive content changed")
    require(metadata == ((1, 0, 12),), "selected rail metadata changed")

    roots = tuple(
        (2 * c + (2 * h + kappa) // m.P + (edge == 0)) % m.P
        for c in range(m.P)
    )
    flags = tuple(
        m.is_unit(rows[0][sector][edge][c][kappa][h], roots[c], 26)
        for c in range(m.P)
    )
    unit_carries = tuple(c for c, flag in enumerate(flags) if flag)
    require(unit_carries == (0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
            "selected twelve-carry unit row changed")
    missing_carry = 6

    module, _, _, _, rails, present, _ = m.core.build_carrier_data()
    prefixes = m.build_pair_prefixes(module)
    pieces = rails[j][3]
    rail_support = [(left, right) for left, right, weight in pieces if weight]
    delayed_cache = [[[[{} for _ in range(m.P)] for _ in range(m.Q7)]
                      for _ in range(2)]][0]
    shifted_rail = {}
    shifted_present = {}
    facet_rows = []

    for source_carry in range(m.P):
        missing_delta = (2 * (missing_carry - source_carry)) % m.P
        deltas = tuple(delta for delta in range(m.P)
                       if delta != missing_delta)
        require(tuple(sorted((source_carry + 7 * delta) % m.P
                             for delta in deltas)) == unit_carries,
                "active labels no longer pull back the unit-carry set")
        require(all((roots[source_carry] + delta) % m.P != 0
                    for delta in deltas),
                "an active label acquired root zero")
        require((roots[source_carry] + missing_delta) % m.P == 0,
                "missing label is not the root-zero label")

        anchor = deltas[0]
        anchor_shift = 7 * anchor * m.T // m.R
        common_rail = shift_weighted(pieces, anchor_shift, m.T)
        for delta in deltas[1:]:
            key = (j, delta)
            if key not in shifted_rail:
                shift = 7 * delta * m.T // m.R
                shifted_rail[key] = shift_union(rail_support, shift, m.T)
            common_rail = m.old.intersect_weighted_union(
                common_rail, shifted_rail[key]
            )
            if not common_rail:
                break
        require(common_rail, "twelve translated rail supports became disjoint")

        total = 0
        first_component = None
        for ell5 in range(m.Q7):
            source_present = present[ell5, (-h) % m.P]
            common = m.old.intersect_weighted_union(
                common_rail,
                shift_union(source_present, anchor_shift, m.T),
            )
            for delta in deltas[1:]:
                if not common:
                    break
                key = (ell5, h, delta)
                if key not in shifted_present:
                    shift = 7 * delta * m.T // m.R
                    shifted_present[key] = shift_union(
                        source_present, shift, m.T
                    )
                common = m.old.intersect_weighted_union(
                    common, shifted_present[key]
                )
            if not common:
                continue
            half = intersect_root_half(
                common, module.C3, edge, roots[source_carry]
            )
            component = first_delayed_component(
                half, prefixes[sector][ell5][h][kappa], source_carry
            )
            if component is not None:
                candidate = (ell5, component)
                if first_component is None or component < first_component[1]:
                    first_component = candidate
            values = m.delayed_carry_pair(
                half, prefixes[sector][ell5][h],
                delayed_cache[sector][ell5][h],
            )
            total += values[source_carry][kappa]

        require(total > 0 and first_component is not None,
                "one carry-rebased twelve-facet lost positive overlap")
        ell5, (left_raw, right_raw) = first_component
        left = Fraction(left_raw, m.R * m.T)
        right = Fraction(right_raw, m.R * m.T)
        facet_rows.append((source_carry, missing_delta, total, ell5,
                           left, right))

    source_to_missing = tuple((row[0], row[1]) for row in facet_rows)
    require(tuple(sorted(row[1] for row in facet_rows)) == tuple(range(m.P)),
            "source-carry/missing-label map stopped being bijective")
    lengths = tuple(row[5] - row[4] for row in facet_rows)
    component_length = Fraction(3, 12545122758259)
    require(set(lengths) == {component_length},
            "canonical first-component lengths stopped agreeing")

    # Coarse disconnected-union labels V_delta forget source_carry.  Their
    # nerve is the boundary of a 12-simplex.  The thirteen selected first
    # component witnesses instead span disjoint filled 11-simplices: at each
    # fixed delta their nonzero private roots are distinct, and same-edge
    # private halves are separated open components.
    coarse_fvector = tuple(comb(13, k + 1) for k in range(12))
    refined_fvector = tuple(13 * comb(12, k + 1) for k in range(12))
    require(coarse_fvector ==
            (13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13),
            "coarse boundary f-vector changed")
    require(refined_fvector ==
            (156, 858, 2860, 6435, 10296, 12012,
             10296, 6435, 2860, 858, 156, 13),
            "carry-refined selected f-vector changed")
    require(sum((-1) ** k * coarse_fvector[k] for k in range(12)) == 0,
            "coarse S11 Euler characteristic changed")
    require(sum((-1) ** k * refined_fvector[k] for k in range(12)) == 13,
            "refined disjoint-simplex Euler characteristic changed")
    coarse_ridges = comb(13, 11)
    refined_ridges = 13 * comb(12, 11)
    require((coarse_ridges, refined_ridges) == (78, 156),
            "codimension-one ridge incidence count changed")
    for delta in range(m.P):
        selected_roots = tuple(sorted(
            (roots[c] + delta) % m.P
            for c in range(m.P)
            if (2 * (missing_carry - c)) % m.P != delta
        ))
        require(selected_roots == tuple(range(1, m.P)),
                "selected facet components stopped separating by root")

    # One slope-seven rebase sends c0 -> c0+7 and m -> m-1.  Thirteen
    # canonical representatives accumulate the nonzero extension-kernel
    # translation 91/R=7/S; this arithmetic does not supply a chart map.
    require(all((2 * (missing_carry - (c + 7))
                 - ((2 * (missing_carry - c)) - 1)) % m.P == 0
                for c in range(m.P)),
            "missing-label torsor step changed")
    require(Fraction(13 * 7, m.R) == Fraction(7, m.S)
            and 0 < Fraction(7, m.S) < 1,
            "odometer kernel accumulation changed")

    totals = tuple(row[2] for row in facet_rows)
    facet0 = next(row for row in facet_rows if row[1] == 0)
    facet12 = next(row for row in facet_rows if row[1] == 12)

    # The canonical missing-12 component sits immediately to the left-only
    # side of the edge overlap.  On it the missing translated label has
    # edge-0 root zero, while edge 1 would require q<13.  Since h,kappa and
    # carry are pointwise fixed and root-zero rows are never units, even
    # switching rail/sector/configuration cannot fill this particular face.
    deep_left = (module.C3 * facet12[4] % 1) * 182
    deep_right = (module.C3 * facet12[5] % 1) * 182
    require((deep_left, deep_right) ==
            (Fraction(13), Fraction(13) + Fraction(12, 371293)),
            "canonical maximal face left the one-edge deep strip")
    missing_deep_left = deep_left + 14 * 12
    missing_deep_right = deep_right + 14 * 12
    require((missing_deep_left, missing_deep_right) ==
            (Fraction(181), Fraction(181) + Fraction(12, 371293)),
            "missing-label translated deep strip changed")

    print("THM2672 SLOPE-SEVEN REBASE FACET TORSOR EXACT REFEREE")
    print(f"config={(j, sector, edge, kappa, h)} metadata={metadata[0]}")
    print(f"unit_carries={unit_carries} missing_carry={missing_carry}")
    print(f"source_to_missing={source_to_missing}")
    print(f"positive_facets={len(facet_rows)} distinct_totals={tuple(sorted(set(totals)))}")
    print(f"first_component_length={component_length}")
    print(f"missing0_component={(facet0[0], facet0[3], facet0[4], facet0[5])}")
    print(f"missing12_component={(facet12[0], facet12[3], facet12[4], facet12[5])}")
    print(f"missing12_pulled_deep_strip=({deep_left},{deep_right}) "
          "edge0:root0 edge1:absent")
    print(f"missing12_translated_deep_strip="
          f"({missing_deep_left},{missing_deep_right})")
    print(f"coarse_union_nerve_fvector={coarse_fvector}")
    print("coarse_union_nerve=reduced_H11:1 reduced_H0:0 "
          "codimension_one_ridges:78_each_twice")
    print(f"selected_refined_subcomplex_fvector={refined_fvector}")
    print("selected_refined_subcomplex=components:13 reduced_H0:12 H11:0 "
          "codimension_one_ridges:156_each_once")
    print("facet_projection=(carry,delta)->delta creates a virtual S11; "
          "it is not a carry-refined nerve equivalence")
    print(f"thirteen_step_lift=91/{m.R}=7/{m.S} cocycle_class=7")
    print("cocycle_role=obstructs equivariant cyclic gluing; "
          "not used to derive either f-vector")
    print("scope=one full-union maximal face and failure of this canonical "
          "top-cycle lift; no connected-arc gain refinement, global all13 "
          "exclusion, endpoint current, or LRC decrement")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
