#!/usr/bin/env python3
"""Reconstruct the broadcast E3+D6 relative-present private bank exactly."""

from collections import Counter
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_relative_present_semantic_lift_probe_20260728 as relative
import lrc14_two_target_present_semantic_attachment_probe_20260728 as semantic


P = relative.P
H = 6


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def complement(intervals, period):
    out = []
    cursor = 0
    for left, right in intervals:
        if cursor < left:
            out.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < period:
        out.append((cursor, period))
    return tuple(out)


def build_semantic_prefixes(module, fork):
    """[sector][clock][kappa] prefixes on y={13^6 x}."""
    result = []
    for word in relative.private.prior.sector_words(module):
        by_clock = []
        for ell in range(7):
            qell = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            qell = semantic.intersect_sorted(qell, fork)
            pair = []
            for kappa in range(2):
                digit = relative.private.old.sat.intersect_interval(
                    qell,
                    (2 * H + kappa) * module.T // (2 * P),
                    (2 * H + kappa + 1) * module.T // (2 * P),
                )
                pair.append(module.make_prefix(digit))
            by_clock.append(tuple(pair))
        result.append(tuple(by_clock))
    return tuple(result)


def build_fdel(module, ell, shift):
    """Canonical F_(ell,shift) with only unshifted c3-safe deleted."""
    intervals = module.make_comb(
        module.C1, 182, 26 * ell - 13, 26 * ell + 13
    )
    intervals = module.subtract_comb(
        intervals, module.W[module.GUARD], 91, -13, 13
    )
    intervals = module.subtract_comb(
        intervals, module.W[1], 182,
        -14 * shift - 13, -14 * shift + 13,
    )
    for index in (2, 3, 4, 5):
        intervals = module.subtract_comb(
            intervals, module.W[index], 182, -13, 13
        )
    intervals = module.subtract_comb(
        intervals, module.C2, 182,
        14 * shift - 13, 14 * shift + 13,
    )
    return tuple(intervals)


def main():
    module, _prefixes, _, _, rails, present, _starts = (
        relative.lift.m.core.build_carrier_data()
    )
    source = semantic.exclusive_source(module, 3)
    fork = semantic.deepest_fork(module)
    semantic_prefixes = build_semantic_prefixes(module, fork)
    original_prefixes = relative.private.build_pair_prefixes(module)
    fdel = tuple(tuple(build_fdel(module, ell, shift)
                       for shift in range(P))
                 for ell in range(7))

    rows = []
    content = 0
    roots = Counter()
    for rail_index in range(14):
        source_rail = relative.private.old.intersect_weighted_union(
            rails[rail_index][3], source
        )
        for shift in range(P):
            for sector in range(2):
                for edge in range(2):
                    root = 1 if edge == 0 else 12
                    for kappa in range(2):
                        carry = (
                            (root - (2 * H + kappa) // P
                             - int(edge == 0)) * 7
                        ) % P
                        vector = []
                        for ell in range(7):
                            pieces = relative.private.old.intersect_weighted_union(
                                source_rail, fdel[ell][shift]
                            )
                            low = 14 * root - 13 if edge == 0 else 14 * root
                            pieces = relative.private.old.intersect_weighted_comb(
                                pieces, module.C3, 182, low, low + 13
                            )
                            values = relative.private.delayed_carry_pair(
                                pieces,
                                (semantic_prefixes[sector][ell][0],
                                 semantic_prefixes[sector][ell][1]),
                                {},
                            )
                            value = values[carry][kappa]
                            vector.append(value)
                            if value:
                                content = gcd(content, value)
                        vector = tuple(vector)
                        if any(vector):
                            roots[root] += 1
                        rows.append((rail_index, shift, sector, edge,
                                     kappa, carry, root, vector))

    require(len(rows) == 1456, "graph-row universe changed")
    positive = sum(any(row[-1]) for row in rows)
    units = [row for row in rows
             if relative.private.is_unit(row[-1], row[-2], content)]
    unit_roots = Counter(row[-2] for row in units)
    print(f"rows={len(rows)} positive={positive} content={content}")
    print(f"positive_roots={tuple(sorted(roots.items()))}")
    print(f"units={len(units)} unit_roots={tuple(sorted(unit_roots.items()))}")
    extrema = []
    for row in units:
        values = row[-1]
        root = row[-2]
        scale = pow(root, -1, P)
        normalized = tuple((value // content) * scale % P for value in values)
        reduced = tuple((normalized[i] - normalized[-1]) % P for i in range(6))
        if sum(bool(value) for value in reduced) == 1:
            extrema.append((row[:-1], normalized, reduced))
    print(f"one_term_units={tuple(extrema[:30])}")


if __name__ == "__main__":
    main()
