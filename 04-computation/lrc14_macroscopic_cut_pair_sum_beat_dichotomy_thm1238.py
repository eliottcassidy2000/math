#!/usr/bin/env python3
"""Exact referee for THM-1238's pair-sum beat/blocker dichotomy.

The script exhausts a deterministic integer bank using cleared endpoint
inequalities, checks the positive/equality curvature alternatives, and
verifies the infinite parity-locked singleton obstruction through m=1000.
It is dependency-free, exact, and uses no ``assert``.
"""

from itertools import combinations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def beat_block(c: int, k: int, q: int) -> range:
    denominator = 14 * c
    lo = ceil_div((14 * k + 1) * q, denominator)
    hi = ((14 * k + 13) * q) // denominator
    return range(lo, hi + 1)


def least_residue_distance(x: int, p: int, q: int) -> int:
    r = (x * p) % q
    return min(r, q - r)


rows = 0
forced_rows = 0
singleton_rows = 0
positive_rows = 0
zero_rows = 0

for c in range(1, 21):
    for k in range(c):
        for x, y in combinations(range(c + 1, 8 * c + 1), 2):
            q = x + y
            block = tuple(beat_block(c, k, q))
            rows += 1
            require(q > 2 * c, "pair sum exceeds twice carrier")
            require(block, "nonempty pair-sum beat block")

            safe = []
            for p in block:
                D = least_residue_distance(x, p, q)
                require(
                    D == least_residue_distance(y, p, q),
                    "sum-beat defining-mask equality",
                )
                if 14 * D >= q:
                    safe.append((p, D))
                    if 14 * D == q:
                        zero_rows += 1
                        require(q % 14 == 0, "zero curvature requires 14|q")
                    else:
                        positive_rows += 1

            sufficient = 3 * q >= 7 * c or y > 6 * x
            if sufficient:
                forced_rows += 1
                require(safe, "two-speed safe-beat supplier")

            if not safe:
                singleton_rows += 1
                require(2 * c < q and 3 * q < 7 * c, "singleton q horn")
                require(y <= 6 * x, "singleton ratio horn")
                require(len(block) == 1, "singleton beat block")
                D = least_residue_distance(x, block[0], q)
                require(14 * D < q, "negative-curvature singleton")

            require(q // gcd(x, y) >= 3, "nontrivial reduced defining period")

# Infinite parity-locked obstruction to extracting a positive fast-fast edge
# from the macroscopic cut alone.
parity_pair_checks = 0
for m in range(1, 1001):
    c = 420 * m + 1
    k = 210 * m
    speeds = (428 * m + 2, 440 * m + 2, 452 * m + 2,
              464 * m + 2, 476 * m + 2, 500 * m + 2)
    require(tuple(sorted(speeds)) == speeds and all(d % 2 == 0 for d in speeds),
            "parity packet order")
    threshold_num = 7 * c
    require(6 * speeds[4] < threshold_num < 6 * speeds[5], "r=1 threshold")
    invoice = sum(speeds[5] - d for d in speeds)
    require(14 * invoice > speeds[5], "fastest Kakeya invoice")
    gaps = tuple(speeds[i + 1] - speeds[i] for i in range(5))
    require(max(gaps) == 24 * m and 180 * max(gaps) > c, "macroscopic cut")

    for x, y in combinations(speeds, 2):
        q = x + y
        block = tuple(beat_block(c, k, q))
        require(block == (q // 2,), "centered singleton block")
        require(q % 2 == 0 and least_residue_distance(x, q // 2, q) == 0,
                "common-zero parity lock")
        require(2 * c < q and 3 * q < 7 * c and y <= 6 * x,
                "singleton horn inequalities")
        require(q // gcd(x, y) >= 3, "parity reduced period")
        parity_pair_checks += 1

print("THM-1238 MACROSCOPIC-CUT PAIR-SUM BEAT DICHOTOMY EXACT AUDIT")
print(f"two-speed rows checked = {rows}")
print(f"rows in a sufficient safe-beat branch = {forced_rows}")
print(f"negative singleton rows = {singleton_rows}")
print(f"positive safe-beat occurrences = {positive_rows}")
print(f"zero-curvature safe-beat occurrences = {zero_rows}")
print(f"parity-family fast-fast pairs checked = {parity_pair_checks}")
print("parity family range = m=1..1000")
print("RESULT: PASS")
