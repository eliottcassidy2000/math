#!/usr/bin/env python3
"""Exact companion for THM-2623's guard cospan and half-tooth wall."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_cross_time_target_future_diagonal_thm2616 as core

old = core.old
P, Q7 = core.P, core.Q7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def build_guard_danger_prefixes(module):
    word = module.make_comb(module.C2, 182, -13, 13)
    word = module.intersect_comb(word, module.W[module.GUARD], 91, -13, 13)
    for i in module.UNIT_IDX:
        word = module.subtract_comb(word, module.W[i], 182, -13, 13)
    word = module.subtract_comb(word, module.C3, 182, -13, 13)
    prefixes = [[None] * P for _ in range(Q7)]
    for ell in range(Q7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for h in range(P):
            digit = old.sat.intersect_interval(
                qell, h * core.T // P, (h + 1) * core.T // P
            )
            prefixes[ell][h] = module.make_prefix(digit)
    return prefixes


def shard(bounds):
    start, stop = bounds
    module, safe_prefix, _, safe_masses, rails, present, starts = core.build_carrier_data()
    danger_prefix = build_guard_danger_prefixes(module)
    danger_masses = tuple(tuple(prefix[2][-1] for prefix in row)
                          for row in danger_prefix)
    prefixes = (safe_prefix, danger_prefix)
    caches = [[dict() for _ in range(Q7)] for _ in range(2)]
    content = 0
    positive = 0
    sector_content = [0, 0]
    sector_positive = [0, 0]
    diagonals = []
    metadata = []
    partition_checks = 0
    for j in range(start, stop):
        s, ell, t, pieces = rails[j]
        metadata.append((s, ell, t))
        # [sector safe/danger][edge left/right][q][ell5][r]
        raw = [[[[[0] * P for _ in range(Q7)] for _ in range(P)]
                for _ in range(2)] for _ in range(2)]
        for q in range(P):
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, (-q) % P], starts[ell5, (-q) % P]
                )
                for r in range(1, P):
                    whole = old.intersect_weighted_comb(
                        overlap, module.C3, 182, 14 * r - 13, 14 * r + 13
                    )
                    halves = (
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * r - 13, 14 * r
                        ),
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * r, 14 * r + 13
                        ),
                    )
                    for sector in range(2):
                        whole_values = core.delayed_all_digits(
                            whole, prefixes[sector][ell5], caches[sector][ell5]
                        )
                        for value in whole_values:
                            if value:
                                sector_positive[sector] += 1
                                sector_content[sector] = gcd(
                                    sector_content[sector], value
                                )
                        half_values = tuple(
                            core.delayed_all_digits(
                                half, prefixes[sector][ell5], caches[sector][ell5]
                            ) for half in halves
                        )
                        if tuple(half_values[0][h] + half_values[1][h]
                                 for h in range(P)) != whole_values:
                            raise RuntimeError("edge halves do not partition the deep danger tooth")
                        partition_checks += 1
                        for edge in range(2):
                            values = half_values[edge]
                            for value in values:
                                if value:
                                    positive += 1
                                    content = gcd(content, value)
                            raw[sector][edge][q][ell5][r] = values[q]
        diagonals.append(tuple(
            tuple(
                tuple(tuple(tuple(row) for row in qrows) for qrows in edge_rows)
                for edge_rows in sector_rows
            ) for sector_rows in raw
        ))
    return (bounds, content, positive, tuple(sector_content),
            tuple(sector_positive), partition_checks,
            tuple(metadata), tuple(diagonals), safe_masses, danger_masses)


def unit(raw, content):
    y = tuple(
        sum((raw[ell][r] // content) * pow(r, -1, P)
            for r in range(1, P)) % P
        for ell in range(Q7)
    )
    reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
    return bool(old.sat.multiplication_determinant_7(reduced))


def fixed_root_unit(raw, q, content):
    """Unit test on the affine graph r=-h-1 with h=q."""
    r = (-q - 1) % P
    if r == 0:
        return False
    scale = pow(r, -1, P)
    y = tuple((raw[ell][r] // content) * scale % P for ell in range(Q7))
    reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
    return bool(old.sat.multiplication_determinant_7(reduced))


def main():
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(shard, core.SHARDS))
    content = 0
    positive = 0
    checks = 0
    metadata = []
    diagonal = []
    sector_content = [0, 0]
    sector_positive = [0, 0]
    safe_masses = results[0][8]
    danger_masses = results[0][9]
    for bounds, c, p, sc, sp, n, meta, raw, sm, dm in results:
        require(sm == safe_masses and dm == danger_masses,
                "workers disagree on endpoint-sector masses")
        content = gcd(content, c)
        positive += p
        for sector in range(2):
            sector_content[sector] = gcd(sector_content[sector], sc[sector])
            sector_positive[sector] += sp[sector]
        checks += n
        metadata.extend(meta)
        diagonal.extend(raw)
    flags = [[[[False] * P for _ in range(2)] for _ in range(2)]
             for _ in diagonal]
    graph_flags = [[[[False] * P for _ in range(2)] for _ in range(2)]
                   for _ in diagonal]
    unit_counts = [[0, 0], [0, 0]]
    graph_unit_counts = [[0, 0], [0, 0]]
    positive_pairs = [[0, 0], [0, 0]]
    support_hist = [[Counter(), Counter()], [Counter(), Counter()]]
    unsplit_flags = [[[False] * P for _ in range(2)] for _ in diagonal]
    unsplit_units = [0, 0]
    unsplit_pairs = [0, 0]
    unsplit_support = [Counter(), Counter()]
    for j, rows in enumerate(diagonal):
        for sector in range(2):
            for q in range(P):
                joined = tuple(tuple(
                    rows[sector][0][q][ell][r]
                    + rows[sector][1][q][ell][r]
                    for r in range(P)) for ell in range(Q7))
                unsplit_pairs[sector] += int(any(
                    joined[ell][r] for ell in range(Q7) for r in range(1, P)
                ))
                if unit(joined, 26):
                    unsplit_flags[j][sector][q] = True
                    unsplit_units[sector] += 1
                for ell in range(Q7):
                    size = sum(bool(joined[ell][r]) for r in range(1, P))
                    if size:
                        unsplit_support[sector][size] += 1
            for edge in range(2):
                for q in range(P):
                    fine = rows[sector][edge][q]
                    positive_pairs[sector][edge] += int(any(
                        fine[ell][r] for ell in range(Q7) for r in range(1, P)
                    ))
                    if unit(fine, content):
                        flags[j][sector][edge][q] = True
                        unit_counts[sector][edge] += 1
                    if fixed_root_unit(fine, q, content):
                        graph_flags[j][sector][edge][q] = True
                        graph_unit_counts[sector][edge] += 1
                    for ell in range(Q7):
                        size = sum(bool(fine[ell][r]) for r in range(1, P))
                        if size:
                            support_hist[sector][edge][size] += 1
    bycell = {}
    for j, (s, ell, t) in enumerate(metadata):
        bycell.setdefault((s, ell), []).append(j)
    qsets = {}
    unsplit_qsets = {}
    graph_qsets = {}
    for cell, rails in sorted(bycell.items()):
        qsets[cell] = tuple(
            q for q in range(P)
            if any(flags[j][sector][edge][q]
                   for j in rails for sector in range(2) for edge in range(2))
        )
        unsplit_qsets[cell] = tuple(
            q for q in range(P)
            if any(unsplit_flags[j][sector][q]
                   for j in rails for sector in range(2))
        )
        graph_qsets[cell] = tuple(
            q for q in range(P)
            if any(graph_flags[j][sector][edge][q]
                   for j in rails for sector in range(2) for edge in range(2))
        )
    deficient = tuple(
        (cell, tuple(q for q in range(P) if q not in values))
        for cell, values in qsets.items() if len(values) < P
    )
    unsplit_deficient = tuple(
        (cell, tuple(q for q in range(P) if q not in values))
        for cell, values in unsplit_qsets.items() if len(values) < P
    )
    expected_unsplit_deficient = (
        ((2,2),(11,)), ((6,5),(2,)), ((6,6),(10,)),
        ((7,2),(11,)), ((8,2),(12,)), ((11,2),(12,)),
        ((11,5),(2,)), ((11,6),(10,)),
    )
    safe_support = tuple(h for h in range(P)
                         if any(safe_masses[ell][h] for ell in range(Q7)))
    danger_support = tuple(h for h in range(P)
                           if any(danger_masses[ell][h] for ell in range(Q7)))
    expected_deficient = (
        ((2, 2), (11,)), ((6, 5), (2,)),
        ((7, 2), (11,)), ((11, 5), (2,)),
    )
    qset_hist = tuple(sorted(Counter(map(len, qsets.values())).items()))
    common = tuple(sorted(set.intersection(*(set(v) for v in qsets.values()))))
    graph_hist = tuple(sorted(Counter(map(len, graph_qsets.values())).items()))
    graph_common = tuple(sorted(set.intersection(
        *(set(v) for v in graph_qsets.values())
    )))
    graph_extra = tuple(
        (cell, tuple(q for q in range(P - 1) if q not in values))
        for cell, values in graph_qsets.items()
        if any(q not in values for q in range(P - 1))
    )
    require(content == 26 and positive == 1_661_176,
            "refined carrier content/census changed")
    require(sector_content == [26, 86_814]
            and sector_positive == [649_968, 244_992],
            "unrefined sector content/census changed")
    require(safe_support == tuple(range(1,12))
            and danger_support == (0,1,11,12)
            and set(safe_support) | set(danger_support) == set(range(P)),
            "guard cospan did not fill the future support")
    require(unsplit_units == [1483, 533]
            and unsplit_pairs == [1703, 613]
            and unsplit_support == [Counter({12:4156}), Counter({12:1526})],
            "unrefined sector diagonal census changed")
    require(tuple(sorted(Counter(map(len, unsplit_qsets.values())).items()))
            == ((12,8),(13,76))
            and unsplit_deficient == expected_unsplit_deficient,
            "unrefined eight-cell wall changed")
    require(checks == 353_808, "half-tooth partition census changed")
    require(unit_counts == [[1442, 1450], [547, 530]],
            "half-tooth unit census changed")
    require(graph_unit_counts == [[1356, 1458], [276, 393]],
            "fixed graph unit census changed")
    require(positive_pairs == [[1703, 1703], [613, 613]],
            "half-tooth pair census changed")
    require(support_hist[0][0] == support_hist[0][1]
            == Counter({11: 3800, 12: 356}),
            "safe refined deep support changed")
    require(support_hist[1][0] == Counter({11: 1150, 12: 376})
            and support_hist[1][1] == Counter({11: 1130, 12: 396}),
            "danger refined deep support changed")
    require(qset_hist == ((12, 4), (13, 80))
            and deficient == expected_deficient,
            "half-tooth four-cell wall changed")
    require(graph_hist == ((11, 4), (12, 80))
            and all(12 not in values for values in graph_qsets.values())
            and graph_extra == expected_deficient,
            "fixed affine-graph cross-check changed")
    print("THM-2623 THM-2587-style left/right edge-bit probe")
    print(f"sector_full_positive={tuple(sector_positive)} sector_content={tuple(sector_content)} labelled_joint_content=26")
    print(f"future_support_safe={safe_support} danger={danger_support} union={tuple(range(P))}")
    print(f"unsplit_units={tuple(unsplit_units)} positive_pairs={tuple(unsplit_pairs)} deep_support={tuple(tuple(sorted(c.items())) for c in unsplit_support)}")
    print(f"unsplit_qset_hist=((12,8),(13,76)) deficient={unsplit_deficient}")
    print(f"refined_full_positive={positive} refined_global_content={content} partition_checks={checks}")
    print(f"unit_counts_safe_left_right={tuple(unit_counts[0])} danger_left_right={tuple(unit_counts[1])}")
    print(f"fixed_graph_unit_counts_safe_left_right={tuple(graph_unit_counts[0])} danger_left_right={tuple(graph_unit_counts[1])}")
    print(f"positive_pairs_safe_left_right={tuple(positive_pairs[0])} danger_left_right={tuple(positive_pairs[1])}")
    print(f"deep_support_hist_safe_left_right={tuple(tuple(sorted(c.items())) for c in support_hist[0])} danger_left_right={tuple(tuple(sorted(c.items())) for c in support_hist[1])}")
    print(f"qset_size_hist={qset_hist}")
    print(f"full_qset_cells={sum(len(v)==P for v in qsets.values())}/84")
    print(f"common={common}")
    print(f"deficient={deficient}")
    print(f"fixed_graph_qset_size_hist={graph_hist} common={graph_common}")
    print(f"fixed_graph_h12_zero_cells=84 extra_deficient={graph_extra}")


if __name__ == "__main__":
    main()
