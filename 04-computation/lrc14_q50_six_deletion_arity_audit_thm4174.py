#!/usr/bin/env python3
"""Locate the first deletion arity beyond five that rescues q=50."""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import lcm


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(v for v in POOL if v not in ANCHORS)
COMMON = 18_241_159_416_480
Q = 50


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def safe_prefix_tick(tick):
    # J_q(tick/COMMON), represented on denominator 14*q*COMMON.
    product = Q * tick
    whole, rem = divmod(product, COMMON)
    scaled = 14 * rem
    if scaled <= COMMON:
        partial = 0
    elif scaled >= 13 * COMMON:
        partial = 12 * COMMON
    else:
        partial = scaled - COMMON
    return 12 * whole * COMMON + partial


def labels(mask):
    return tuple(v for i, v in enumerate(OPTIONAL) if mask >> i & 1)


def find_cover_exact(edges, budget):
    failed = set()

    def search(chosen, remaining):
        key = (chosen, remaining)
        if key in failed:
            return None
        uncovered = 0
        matching_used = 0
        matching = 0
        for edge in edges:
            if edge & chosen:
                continue
            if not uncovered:
                uncovered = edge
            if not edge & matching_used:
                matching_used |= edge
                matching += 1
                if matching > remaining:
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
            answer = search(chosen | bit, remaining - 1)
            if answer is not None:
                return answer
            branch ^= bit
        failed.add(key)
        return None

    answer = search(0, budget)
    return answer, len(failed)


def main():
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)
    assert common == COMMON
    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    ticks = tuple(int(w * COMMON) for w in walls)
    assert len(ticks) == 7134

    buckets = defaultdict(int)
    hist = Counter()
    previous = safe_prefix_tick(ticks[0])
    for i, (left, right) in enumerate(zip(walls, walls[1:])):
        current = safe_prefix_tick(ticks[i + 1])
        mid = (left + right) / 2
        if all(safe_at(mid, a) for a in ANCHORS):
            failure = 0
            for j, speed in enumerate(OPTIONAL):
                if not safe_at(mid, speed):
                    failure |= 1 << j
            buckets[failure] += current - previous
            hist[failure.bit_count()] += 1
        else:
            hist["anchor"] += 1
        previous = current

    for d in range(6, 21):
        active = []
        equalities = 0
        for vertices in combinations(range(27), d):
            mask = sum(1 << v for v in vertices)
            numerator = 0
            subset = mask
            while True:
                numerator += buckets.get(subset, 0)
                if subset == 0:
                    break
                subset = (subset - 1) & mask
            difference = 9 * numerator - 8 * Q * COMMON
            if difference >= 0:
                active.append(mask)
            if difference == 0:
                equalities += 1
        active = tuple(active)
        cover7, states = find_cover_exact(active, 7)
        cover8 = None
        states8 = 0
        if cover7 is None:
            cover8, states8 = find_cover_exact(active, 8)
            if cover8 is not None:
                assert all(edge & cover8 for edge in active)
        print(
            "ARITY", d,
            "TOTAL", __import__("math").comb(27, d),
            "EDGES", len(active),
            "EQUALITIES", equalities,
            "COVER7", labels(cover7) if cover7 is not None else None,
            "STATES", states,
            "COVER8", labels(cover8) if cover8 is not None else None,
            "STATES8", states8,
            flush=True,
        )
        if cover7 is None:
            break

    print("CELL_HIST", tuple(sorted(hist.items(), key=lambda x: str(x[0]))))


if __name__ == "__main__":
    main()
