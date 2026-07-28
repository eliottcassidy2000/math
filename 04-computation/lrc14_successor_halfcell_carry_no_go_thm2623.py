#!/usr/bin/env python3
"""Exact THM-2623 successor-half-cell carry-loss and no-private-row probe."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_cross_time_target_future_diagonal_thm2616 as core

old = core.old
P, Q7, T = core.P, core.Q7, core.T
N = 2 * P


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sector_words(module):
    safe = module.build_word_Ta()
    danger = module.make_comb(module.C2, 182, -13, 13)
    danger = module.intersect_comb(danger, module.W[module.GUARD], 91, -13, 13)
    for i in module.UNIT_IDX:
        danger = module.subtract_comb(danger, module.W[i], 182, -13, 13)
    danger = module.subtract_comb(danger, module.C3, 182, -13, 13)
    return safe, danger


def build_subprefixes(module):
    result = []
    for word in sector_words(module):
        by_clock = []
        for ell in range(Q7):
            qell = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            by_clock.append(tuple(
                module.make_prefix(old.sat.intersect_interval(
                    qell, j * T // N, (j + 1) * T // N
                )) for j in range(N)
            ))
        result.append(tuple(by_clock))
    return tuple(result)


def subdigit_phi_vector(x, prefixes, lengths, cache):
    if x in cache:
        return cache[x]
    j = min(N - 1, N * x // T)
    values = tuple(
        lengths[k] if k < j
        else old.phi_at(x, *prefixes[k]) if k == j
        else 0
        for k in range(N)
    )
    cache[x] = values
    return values


def delayed_all_subdigits(pieces, prefixes, cache):
    lengths = tuple(prefix[2][-1] for prefix in prefixes)
    weighted_len = 0
    acc_r = 0
    acc = [0] * N
    rred = old.R % T
    for left, right, weight in pieces:
        rleft = left * rred % T
        rright = right * rred % T
        weighted_len += weight * (right - left)
        acc_r += weight * (rright - rleft)
        vleft = subdigit_phi_vector(rleft, prefixes, lengths, cache)
        vright = subdigit_phi_vector(rright, prefixes, lengths, cache)
        for j in range(N):
            acc[j] += weight * (vright[j] - vleft[j])
    floor_numerator = old.R * weighted_len - acc_r
    if floor_numerator % T:
        raise RuntimeError("weighted floor count is not integral")
    quotient = floor_numerator // T
    result = tuple(lengths[j] * quotient + acc[j] for j in range(N))
    if min(result) < 0:
        raise RuntimeError("negative delayed overlap")
    return result


def freeze(value):
    return tuple(freeze(item) for item in value) if isinstance(value, list) else value


def shard(bounds):
    start, stop = bounds
    module, _, _, _, rails, present, starts = core.build_carrier_data()
    prefixes = build_subprefixes(module)
    caches = [[dict() for _ in range(Q7)] for _ in range(2)]
    content = 0
    positive = 0
    checks = 0
    metadata = []
    diagonal = []
    for j in range(start, stop):
        s, ell, theta, pieces = rails[j]
        metadata.append((s, ell, theta))
        # [sector][edge left/right][kappa][h=q][ell5][r]
        raw = [[[[[[0] * P for _ in range(Q7)] for _ in range(P)]
                 for _ in range(2)] for _ in range(2)] for _ in range(2)]
        for h in range(P):
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, (-h) % P], starts[ell5, (-h) % P]
                )
                for r in range(1, P):
                    halves = (
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * r - 13, 14 * r
                        ),
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * r, 14 * r + 13
                        ),
                    )
                    for sector in range(2):
                        for edge, half in enumerate(halves):
                            values = delayed_all_subdigits(
                                half, prefixes[sector][ell5], caches[sector][ell5]
                            )
                            for value in values:
                                if value:
                                    positive += 1
                                    content = gcd(content, value)
                            for kappa in range(2):
                                raw[sector][edge][kappa][h][ell5][r] = values[2*h+kappa]
                    checks += 1
        diagonal.append(freeze(raw))
    return bounds, content, positive, checks, tuple(metadata), tuple(diagonal)


def is_unit(fine, content):
    y = tuple(
        sum((fine[ell][r] // content) * pow(r, -1, P)
            for r in range(1, P)) % P
        for ell in range(Q7)
    )
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
    for bounds, c, p, n, meta, raw in results:
        content = gcd(content, c)
        positive += p
        checks += n
        metadata.extend(meta)
        diagonal.extend(raw)

    flags = [[[[[False] * P for _ in range(2)] for _ in range(2)]
              for _ in range(2)] for _ in diagonal]
    units = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
    supports = [[[Counter(), Counter()], [Counter(), Counter()]],
                [[Counter(), Counter()], [Counter(), Counter()]]]
    slice_supports = [[[Counter(), Counter()], [Counter(), Counter()]],
                      [[Counter(), Counter()], [Counter(), Counter()]]]
    private_flags = [[[[[None] * P for _ in range(2)] for _ in range(2)]
                      for _ in range(2)] for _ in diagonal]
    violations = []
    for j, rows in enumerate(diagonal):
        for sector in range(2):
            for edge in range(2):
                for kappa in range(2):
                    for h in range(P):
                        fine = rows[sector][edge][kappa][h]
                        for ell in range(Q7):
                            support = tuple(r for r in range(1, P) if fine[ell][r])
                            if support:
                                supports[sector][edge][kappa][len(support)] += 1
                                predicted = ((2*h+kappa+1) if edge == 0
                                             else (2*h+kappa)) % P
                                if support != ((predicted,) if predicted else ()):
                                    violations.append((j, sector, edge, kappa, h,
                                                       ell, support, predicted))
                        total_support = tuple(
                            r for r in range(1, P)
                            if any(fine[ell][r] for ell in range(Q7))
                        )
                        if total_support:
                            slice_supports[sector][edge][kappa][len(total_support)] += 1
                        if is_unit(fine, content):
                            flags[j][sector][edge][kappa][h] = True
                            units[sector][edge][kappa] += 1
                            private_flags[j][sector][edge][kappa][h] = (
                                total_support[0] if len(total_support) == 1 else None
                            )
    bycell = {}
    for j, (s, ell, theta) in enumerate(metadata):
        bycell.setdefault((s, ell), []).append(j)
    lane_qsets = {}
    for edge in range(2):
        for kappa in range(2):
            lane_qsets[edge, kappa] = {
                cell: tuple(h for h in range(P)
                            if any(flags[j][sector][edge][kappa][h]
                                   for j in rails for sector in range(2)))
                for cell, rails in sorted(bycell.items())
            }
    private_rootsets = {}
    private_witnesses = {}
    for cell, rails in sorted(bycell.items()):
        witnesses = {}
        for r in range(1, P):
            options = []
            for j in rails:
                for sector in range(2):
                    for edge in range(2):
                        for kappa in range(2):
                            for h in range(P):
                                private_root = private_flags[j][sector][edge][kappa][h]
                                if private_root == r:
                                    options.append((j, sector, edge, kappa, h))
            if options:
                witnesses[r] = min(options)
        private_rootsets[cell] = tuple(sorted(witnesses))
        private_witnesses[cell] = tuple(sorted(witnesses.items()))
    unit_tuple = tuple(tuple(tuple(row) for row in sector) for sector in units)
    support_tuple = tuple(tuple(tuple(tuple(sorted(c.items())) for c in edge)
                                for edge in sector) for sector in supports)
    slice_tuple = tuple(tuple(tuple(tuple(sorted(c.items())) for c in edge)
                              for edge in sector) for sector in slice_supports)
    lane_summary = tuple(
        (lane,
         tuple(sorted(Counter(map(len, qsets.values())).items())),
         tuple(sorted(set.intersection(*(set(v) for v in qsets.values())))))
        for lane, qsets in sorted(lane_qsets.items())
    )
    root_hist = tuple(sorted(Counter(map(len, private_rootsets.values())).items()))
    expected_first = (
        (0, 0, 0, 0, 2, 1, (2,3,4,5,6,7,8,9,10,11,12), 5),
        (0, 0, 0, 0, 2, 3, (2,3,4,5,6,7,8,9,10,11,12), 5),
        (0, 0, 0, 0, 3, 1, (2,3,4,5,6,7,8,9,10,11,12), 7),
        (0, 0, 0, 0, 3, 3, (2,3,4,5,6,7,8,9,10,11,12), 7),
    )
    expected_lanes = (
        ((0,0), ((9,2),(11,16),(12,66)), (5,6,7,9,10,12)),
        ((0,1), ((10,2),(11,1),(12,15),(13,66)), (1,6,7,8,9,10)),
        ((1,0), ((10,2),(11,1),(12,13),(13,68)), (1,3,4,5,6,7,8)),
        ((1,1), ((9,2),(11,16),(12,66)), (0,1,3,6,7,8)),
    )
    require(positive == 2_986_852 and content == 26 and checks == 176_904,
            "successor half-cell carrier census changed")
    require(unit_tuple == (((1321,1351),(1334,1354)),
                           ((387,530),(515,385))),
            "successor half-cell unit census changed")
    require(len(violations) == 20_778 and tuple(violations[:4]) == expected_first,
            "false singleton census changed")
    require(all(size in (11,12) for sector in slice_supports
                for edge in sector for counter in edge for size in counter),
            "a successor half-cell slice became private")
    require(lane_summary == expected_lanes,
            "successor half-cell lane atlas changed")
    require(root_hist == ((0,84),),
            "successor half-cell acquired a private unit root")
    print("THM-2623 successor-half-digit/private-deep-row probe")
    print(f"full_positive={positive} global_content={content} fine_checks={checks}")
    print(f"unit_counts={unit_tuple}")
    print(f"support_hist={support_tuple}")
    print(f"slice_support_hist={slice_tuple}")
    print(f"false_singleton_rows={len(violations)} first={tuple(violations[:4])}")
    print(f"lane_summary={lane_summary}")
    print(f"private_unit_root_size_hist={root_hist} failures=all84")
    print("carry_identity=a=floor(13*frac(c3*x))=(2*floor(R*x)+floor(2*frac(R*x))) mod13")


if __name__ == "__main__":
    main()
