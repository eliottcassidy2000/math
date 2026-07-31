#!/usr/bin/env python3
"""Exact common-overlap probe for the largest THM-2640 carry configurations.

It first rebuilds the exact THM-2640 unit flags.  For every fixed physical
configuration that is unit on
12 of the 13 carries, it then intersects all 12 slope-seven translates on the
common T-grid.  The delayed word, predecessor carry, future half-digit and
private deep half-tooth are aligned exactly; every rail and present factor is
translated and intersected honestly.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_predecessor_carry_private_root_atlas_thm2640 as m


def shift_union(intervals, shift, modulus):
    """Preimage under x -> x+shift of a sorted circular interval union."""
    out = []
    for left, right in intervals:
        length = right - left
        start = (left - shift) % modulus
        stop = start + length
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


def build_flags(rows, content):
    flags = []
    for rail_rows in rows:
        rail_flags = [[[[[False] * m.P for _ in range(2)] for _ in range(m.P)]
                       for _ in range(2)] for _ in range(2)]
        for sector in range(2):
            for edge in range(2):
                for c in range(m.P):
                    for kappa in range(2):
                        for h in range(m.P):
                            root = (2 * c + (2 * h + kappa) // m.P
                                    + (edge == 0)) % m.P
                            rail_flags[sector][edge][c][kappa][h] = m.is_unit(
                                rail_rows[sector][edge][c][kappa][h], root,
                                content,
                            )
        flags.append(m.freeze(rail_flags))
    return flags


def main():
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(m.shard, m.core.SHARDS))
    content = 0
    metadata = []
    rows = []
    for _, c, _, _, _, meta, raw in results:
        content = gcd(content, c)
        metadata.extend(meta)
        rows.extend(raw)
    m.require(content == 26, "global content changed")
    flags = build_flags(rows, content)

    module, _, _, _, rails_data, present, starts = m.core.build_carrier_data()
    prefixes = m.build_pair_prefixes(module)
    m.require(tuple((s, ell, theta) for s, ell, theta, _ in rails_data)
              == tuple(metadata), "rail order changed")
    bycell = {}
    for j, (s, ell, _) in enumerate(metadata):
        bycell.setdefault((s, ell), []).append(j)

    shifted_rail = {}
    shifted_present = {}
    delayed_cache = [[[[{} for _ in range(m.P)] for _ in range(m.Q7)]
                      for _ in range(2)]][0]
    tested = 0
    positive = []
    missing_hist = Counter()
    component_base_hist = Counter()
    alternate_prefilter = []
    maximal_strata = Counter()
    alternate_strata = Counter()

    for cell, rails in sorted(bycell.items()):
        for j in rails:
            pieces = rails_data[j][3]
            rail_support = [(left, right) for left, right, weight in pieces
                            if weight]
            for sector in range(2):
                for edge in range(2):
                    for kappa in range(2):
                        for h in range(m.P):
                            carries = tuple(
                                c for c in range(m.P)
                                if flags[j][sector][edge][c][kappa][h]
                            )
                            if len(carries) != 12:
                                continue
                            tested += 1
                            missing = next(c for c in range(m.P)
                                           if c not in carries)
                            missing_hist[missing] += 1
                            missing_root = (
                                2 * missing + (2 * h + kappa) // m.P
                                + (edge == 0)
                            ) % m.P
                            m.require(missing_root == 0,
                                      "maximal row did not miss root zero")
                            stratum = (cell, h, kappa)
                            maximal_strata[stratum] += 1
                            opposite = 1 - edge
                            alternates = tuple(
                                (j_alt, sector_alt, opposite,
                                 (2 * missing + (2 * h + kappa) // m.P
                                  + (opposite == 0)) % m.P)
                                for j_alt in rails
                                for sector_alt in range(2)
                                if flags[j_alt][sector_alt][opposite]
                                        [missing][kappa][h]
                            )
                            m.require(all(root in (1, m.P - 1)
                                          for *_, root in alternates),
                                      "opposite-edge alternate root changed")
                            if alternates:
                                alternate_prefilter.append(
                                    (cell, j, sector, edge, kappa, h,
                                     missing, alternates)
                                )
                                alternate_strata[stratum] += 1
                            c0 = carries[0]
                            deltas = tuple((2 * (c - c0)) % m.P
                                           for c in carries)
                            m.require(len(set(deltas)) == 12 and 0 in deltas,
                                      "carry/delta section changed")

                            common_rail = [(left, right, weight)
                                           for left, right, weight in pieces
                                           if weight]
                            for delta in deltas:
                                if delta == 0:
                                    continue
                                key = (j, delta)
                                if key not in shifted_rail:
                                    shift = 7 * delta * m.T // m.R
                                    shifted_rail[key] = shift_union(
                                        rail_support, shift, m.T
                                    )
                                common_rail = m.old.intersect_weighted_union(
                                    common_rail, shifted_rail[key]
                                )
                                if not common_rail:
                                    break
                            if not common_rail:
                                continue

                            root = (2 * c0 + (2 * h + kappa) // m.P
                                    + (edge == 0)) % m.P
                            m.require(root != 0, "reference root vanished")
                            total = 0
                            base_components = 0
                            for ell5 in range(m.Q7):
                                common = m.old.intersect_weighted_union(
                                    common_rail,
                                    present[ell5, (-h) % m.P],
                                    starts[ell5, (-h) % m.P],
                                )
                                for delta in deltas:
                                    if delta == 0 or not common:
                                        continue
                                    key = (ell5, h, delta)
                                    if key not in shifted_present:
                                        shift = 7 * delta * m.T // m.R
                                        shifted_present[key] = shift_union(
                                            present[ell5, (-h) % m.P],
                                            shift, m.T,
                                        )
                                    common = m.old.intersect_weighted_union(
                                        common, shifted_present[key]
                                    )
                                if not common:
                                    continue
                                lo, hi = ((14 * root - 13, 14 * root)
                                          if edge == 0 else
                                          (14 * root, 14 * root + 13))
                                half = m.old.intersect_weighted_comb(
                                    common, module.C3, 182, lo, hi
                                )
                                base_components += len(half)
                                values = m.delayed_carry_pair(
                                    half, prefixes[sector][ell5][h],
                                    delayed_cache[sector][ell5][h],
                                )
                                total += values[c0][kappa]
                            component_base_hist[base_components] += 1
                            if total:
                                positive.append((cell, j, sector, edge, kappa,
                                                 h, carries, missing, total,
                                                 base_components))

    print("THM-2640 common same-configuration physical overlap probe")
    root_orbits = tuple(
        tuple((root + delta) % m.P for delta in range(m.P))
        for root in range(1, m.P)
    )
    m.require(all(tuple(sorted(orbit)) == tuple(range(m.P))
                  and orbit.count(0) == 1 for orbit in root_orbits),
              "pointwise root-zero obstruction changed")
    missing_tuple = tuple(sorted(missing_hist.items()))
    component_tuple = tuple(sorted(component_base_hist.items()))
    totals = tuple(item[8] for item in positive)
    missing_delta_hist = Counter(
        (2 * (item[7] - item[6][0])) % m.P for item in positive
    )
    missing_delta_tuple = tuple(sorted(missing_delta_hist.items()))
    alternate_by_h = tuple(sorted(Counter(
        item[5] for item in alternate_prefilter
    ).items()))
    alternate_by_kappa = tuple(sorted(Counter(
        item[4] for item in alternate_prefilter
    ).items()))
    alternate_by_cell = tuple(sorted(Counter(
        item[0] for item in alternate_prefilter
    ).items()))
    alternate_strata_tuple = tuple(sorted(
        (cell, h, kappa, count, maximal_strata[cell, h, kappa])
        for (cell, h, kappa), count in alternate_strata.items()
    ))
    total_gcd = 0
    for value in totals:
        total_gcd = gcd(total_gcd, value)
    m.require(tested == 534 and len(positive) == tested,
              "twelve-carry physical common-overlap census changed")
    m.require(missing_tuple == ((0, 133), (6, 269), (12, 132)),
              "missing-carry census changed")
    m.require(missing_delta_tuple == ((11, 265), (12, 269)),
              "missing target-label census changed")
    m.require(not alternate_prefilter and not alternate_strata,
              "a maximal fixed configuration gained an opposite-edge unit")
    m.require(sum(count for _, count in component_tuple) == tested,
              "base-component census lost configurations")
    print(f"cells={len(bycell)} content={content}")
    print("configuration=(cell,rail,sector,edge,kappa,h); "
          "carry_set={c: normalized seven-clock row is a THM2640 unit}")
    print("physical_intersection=all translated rail supports AND all "
          "translated present F_(ell5,-h) factors AND aligned private deep "
          "half-tooth AND predecessor-carry/future-half-digit delayed word")
    print(f"twelve_carry_configurations_tested={tested}")
    print(f"missing_carry_hist={missing_tuple}")
    print(f"missing_delta_hist_at_c0_min={missing_delta_tuple}")
    print(f"positive_twelve_carry_common_overlaps={len(positive)}")
    print("opposite_edge_prefilter=preserve (cell,h,kappa), use the "
          "unique root-zero carry, switch edge, and scan every rail/sector "
          "of that cell for a nonzero THM2640 unit")
    print(f"opposite_edge_candidate_configurations="
          f"{len(alternate_prefilter)}/{tested}")
    print(f"opposite_edge_candidate_by_h={alternate_by_h}")
    print(f"opposite_edge_candidate_by_kappa={alternate_by_kappa}")
    print(f"opposite_edge_candidate_by_cell={alternate_by_cell}")
    print("opposite_edge_candidate_strata="
          "((cell,h,kappa,candidates,maximal_configs),...)="
          f"{alternate_strata_tuple}")
    print(f"raw_overlap_min_max_distinct_gcd_sum="
          f"{(min(totals), max(totals), len(set(totals)), total_gcd, sum(totals))}")
    print(f"first_positive={tuple(positive[:12])}")
    print(f"base_component_count_hist={component_tuple}")
    print("fixed_configuration_upper=for fixed edge/root branch and every "
          "nonzero base root r0, "
          "(r0+delta)_(delta in F13)=F13 with one zero; "
          "one fixed configuration therefore has multiplicity <=12")
    print("scope_boundary=the adjacent left/right half-tooths overlap, so "
          "configuration switching can change root at the missing label; "
          "this census does not bound the full union-labelled chart nerve")


if __name__ == "__main__":
    main()
