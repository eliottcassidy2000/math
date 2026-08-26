#!/usr/bin/env python3
"""Exact standard-library audit for the all-newcomer multideletion theorem.

This script starts from the 61 exact d=3 failures already proved in THM-4170.
It reconstructs the fixed-pool wall arrangement with Fraction arithmetic,
integrates every newcomer comb on the common integer lattice, constructs the
d=4,5,6 repair hypergraphs needed on the residual filtration, and decides the
existence of a transversal of size at most seven by exhaustive recursion.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, lcm


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(v for v in POOL if v not in ANCHORS)
COMMON = 18_241_159_416_480
THRESHOLD = F(4, 63)

# Exact d=3 failure set from proved THM-4170.
B3 = (
    3, 6, 22, 24, 25, 46, 48, 50, 55, 64, 70, 72, 75, 83, 93,
    96, 100, 103, 105, 110, 122, 127, 128, 140, 147, 153, 158, 166,
    172, 173, 183, 186, 192, 206, 210, 220, 256, 260, 270, 282, 294,
    306, 313, 320, 325, 332, 346, 366, 372, 384, 416, 440, 462, 512,
    516, 520, 550, 567, 744, 768, 924,
)
EXPECTED_B4 = (25, 50, 96, 100, 105, 128, 192, 210, 256)
EXPECTED_B5 = (50,)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def labels(mask):
    return tuple(v for i, v in enumerate(OPTIONAL) if mask >> i & 1)


def safe_prefix_numerator(tick, q):
    """Integral of G_q through tick/COMMON on denominator 14*q*COMMON."""
    whole, remainder = divmod(q * tick, COMMON)
    scaled = 14 * remainder
    if scaled <= COMMON:
        partial = 0
    elif scaled >= 13 * COMMON:
        partial = 12 * COMMON
    else:
        partial = scaled - COMMON
    return 12 * whole * COMMON + partial


def find_cover(edges, budget):
    """Exact bounded transversal search with matching lower-bound pruning."""
    failed = set()

    def search(chosen, remaining):
        key = (chosen, remaining)
        if key in failed:
            return None
        uncovered = 0
        matching_used = 0
        matching_size = 0
        for edge in edges:
            if edge & chosen:
                continue
            if not uncovered:
                uncovered = edge
            if not edge & matching_used:
                matching_used |= edge
                matching_size += 1
                if matching_size > remaining:
                    failed.add(key)
                    return None
        if not uncovered:
            return chosen
        if remaining == 0:
            failed.add(key)
            return None
        branch = uncovered
        while branch:
            bit = branch & -branch
            witness = search(chosen | bit, remaining - 1)
            if witness is not None:
                return witness
            branch ^= bit
        failed.add(key)
        return None

    return search(0, budget), len(failed)


def minimum_cover_through_seven(edges):
    total_states = 0
    for budget in range(8):
        witness, states = find_cover(edges, budget)
        total_states += states
        if witness is not None:
            require(witness.bit_count() <= budget, "cover budget mismatch")
            require(all(edge & witness for edge in edges), "invalid cover")
            return budget, witness, total_states
    return None, None, total_states


def main():
    require(len(OPTIONAL) == 27, "optional size")
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)
    require(common == COMMON, "common lattice")

    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    ticks = tuple(int(point * COMMON) for point in walls)
    require(len(walls) == 7_134, "wall count")
    require(all(F(tick, COMMON) == point for tick, point in zip(ticks, walls)),
            "wall embedding")

    failure_masks = []
    cell_histogram = Counter()
    for left, right in zip(walls, walls[1:]):
        midpoint = (left + right) / 2
        if any(not safe_at(midpoint, anchor) for anchor in ANCHORS):
            failure_masks.append(None)
            cell_histogram["anchor"] += 1
            continue
        failure = 0
        for index, speed in enumerate(OPTIONAL):
            if not safe_at(midpoint, speed):
                failure |= 1 << index
        failure_masks.append(failure)
        cell_histogram[failure.bit_count()] += 1
    require(len(failure_masks) == 7_133, "cell count")

    bucket_cache = {}

    def buckets_for(q):
        if q in bucket_cache:
            return bucket_cache[q]
        buckets = defaultdict(int)
        previous = safe_prefix_numerator(ticks[0], q)
        for index, failure in enumerate(failure_masks):
            current = safe_prefix_numerator(ticks[index + 1], q)
            contribution = current - previous
            require(contribution >= 0, "negative cell contribution")
            if failure is not None:
                buckets[failure] += contribution
            previous = current
        bucket_cache[q] = dict(buckets)
        return bucket_cache[q]

    def hypergraph(q, d):
        buckets = buckets_for(q)
        edges = []
        equalities = 0
        minimum_positive = None
        maximum_negative = None
        for vertices in combinations(range(27), d):
            mask = sum(1 << vertex for vertex in vertices)
            numerator = 0
            subset = mask
            while True:
                numerator += buckets.get(subset, 0)
                if subset == 0:
                    break
                subset = (subset - 1) & mask
            difference = 9 * numerator - 8 * q * COMMON
            if difference >= 0:
                edges.append(mask)
                if difference == 0:
                    equalities += 1
                elif minimum_positive is None or difference < minimum_positive:
                    minimum_positive = difference
            elif maximum_negative is None or difference > maximum_negative:
                maximum_negative = difference
        return tuple(edges), equalities, minimum_positive, maximum_negative

    ledgers = {}

    def classify(q, d):
        edges, equalities, positive, negative = hypergraph(q, d)
        tau, cover, states = minimum_cover_through_seven(edges)
        cover8 = None
        states8 = 0
        if tau is None and (q, d) == (50, 6):
            cover8, states8 = find_cover(edges, 8)
            require(cover8 is not None and cover8.bit_count() == 8,
                    "q50 d6 eight-cover")
            require(all(edge & cover8 for edge in edges),
                    "q50 d6 invalid eight-cover")
        ledgers[(q, d)] = {
            "edges": len(edges),
            "edge_labels": tuple(labels(edge) for edge in edges) if len(edges) <= 10 else None,
            "equalities": equalities,
            "tau_le_7": tau,
            "cover": labels(cover) if cover is not None else None,
            "states": states,
            "cover8": labels(cover8) if cover8 is not None else None,
            "states8": states8,
            "min_positive_numerator": positive,
            "closest_negative_numerator": negative,
        }
        return tau is None

    b4 = []
    for q in B3:
        if not classify(q, 4):
            b4.append(q)
    require(tuple(b4) == EXPECTED_B4, "d=4 failure set")

    b5 = []
    for q in b4:
        if not classify(q, 5):
            b5.append(q)
    require(tuple(b5) == EXPECTED_B5, "d=5 failure set")

    b6 = []
    for q in b5:
        if not classify(q, 6):
            b6.append(q)
    require(not b6, "d=6 residual")

    # Exact minimum cover sizes agree with the independent MILP scout.  The
    # recursion may return a different lexicographic witness of the same size.
    expected_tau4 = {25: 6, 50: 1, 96: 6, 100: 5, 105: 6,
                     128: 6, 192: 7, 210: 5, 256: 6}
    for q, tau in expected_tau4.items():
        require(ledgers[(q, 4)]["tau_le_7"] == tau, f"d4 tau q={q}")
    require(ledgers[(50, 5)]["tau_le_7"] == 5, "d5 tau q=50")

    print("LRC14_MULTIDELETION_ALL_NEWCOMERS_EXACT")
    print("pool_size", len(POOL), "optional_size", len(OPTIONAL))
    print("walls", len(walls), "cells", len(walls) - 1, "common", COMMON)
    print("cell_failure_histogram", tuple(sorted(cell_histogram.items(),
                                                  key=lambda item: str(item[0]))))
    print("threshold", f"{THRESHOLD.numerator}/{THRESHOLD.denominator}",
          "comparison", ">=")
    print("d3_failure_count", len(B3))
    print("d3_failure_set", B3)
    print("d4_failure_count", len(b4))
    print("d4_failure_set", tuple(b4))
    print("d5_failure_count", len(b5))
    print("d5_failure_set", tuple(b5))
    print("d6_failure_count", len(b6))
    print("d6_failure_set", tuple(b6))
    print("d4_residual_rows")
    for q in b4:
        print(q, ledgers[(q, 4)])
    print("d5_residual_rows")
    for q in b5:
        print(q, ledgers[(q, 5)])
    print("d6_positive_rows")
    for q in b5:
        print(q, ledgers[(q, 6)])
    print("d4_rescued_count", len(B3) - len(b4))
    print("d5_rescued_count", len(b4) - len(b5))
    print("d6_rescued_count", len(b5) - len(b6))
    print("threshold_equalities_new", sum(row["equalities"] for row in ledgers.values()))
    print("bodies_per_newcomer", comb(27, 7))
    print("old_thm4156_bodies", comb(27, 8))
    print("PASS")


if __name__ == "__main__":
    main()
