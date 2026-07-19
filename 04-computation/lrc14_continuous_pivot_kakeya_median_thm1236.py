#!/usr/bin/env python3
"""Exact arithmetic audit for THM-1236's continuous-pivot median law.

The geometric input is the paper lemma obtained by freezing six moving open
arcs on an auxiliary full lap.  This dependency-free referee checks every
rational constant, the constrained tilted-L1 optimizer, and the three
path-weight consumers.  It uses no floating arithmetic and no ``assert`` so
that ordinary and optimized Python execute identical checks.
"""

from fractions import Fraction as F
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def l1(packet: tuple[F, ...], y: F) -> F:
    return sum((abs(d - y) for d in packet), F(0))


def tilted(packet: tuple[F, ...], y: F) -> F:
    return l1(packet, y) - y / 14


require(F(1, 2) - F(1, 14) == F(3, 7), "large-radius cutoff")
require(F(6, 7) + 2 * F(1, 14) == 1, "six-arc critical invoice")
require(F(7, 6) / 14 == F(1, 12), "carrier threshold invoice")
require(F(15, 154) / 7 == F(15, 1078), "cross-check strict-spectrum scale")

# Exhaust every ordered six-packet in a modest exact bank and every breakpoint
# of its tilted L1 functional, together with carrier thresholds and midpoints.
packet_count = 0
threshold_count = 0
sample_count = 0
branch_histogram = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}

for raw in combinations(range(2, 18), 6):
    packet = tuple(F(x) for x in raw)
    d4 = packet[3]
    packet_count += 1
    for c_int in range(1, raw[0]):
        a = F(7 * c_int, 6)
        y0 = max(a, d4)
        threshold_count += 1

        samples = {a, y0, *packet}
        ordered_breaks = sorted({a, *packet})
        for left, right in zip(ordered_breaks, ordered_breaks[1:]):
            samples.add((left + right) / 2)
        samples.add(max(ordered_breaks) + 1)
        samples.add(max(ordered_breaks) + 7)

        target = tilted(packet, y0)
        for y in samples:
            if y >= a:
                sample_count += 1
                require(target <= tilted(packet, y), "constrained optimizer")

        r = sum(d >= a for d in packet)
        if r:
            branch_histogram[r] += 1
            gaps = tuple(packet[i + 1] - packet[i] for i in range(5))
            M = max(gaps)
            if r == 1:
                require(y0 == a, "r=1 optimizer location")
                require(l1(packet, y0) <= 15 * M, "r=1 weight 15")
            elif r == 2:
                require(y0 == a, "r=2 optimizer location")
                require(l1(packet, y0) <= 11 * M, "r=2 weight 11")
            else:
                require(y0 == d4, "r>=3 optimizer location")
                require(l1(packet, y0) <= 9 * M, "r>=3 weight 9")
                require(
                    l1(packet, d4)
                    == packet[3] + packet[4] + packet[5]
                    - packet[0] - packet[1] - packet[2],
                    "upper-median half-sum identity",
                )

# Force every threshold branch independently; the compact combination bank
# above is dense in optimizer chambers but its small ambient range realizes
# only the two fastest-suffix counts.
for r in range(1, 7):
    c = F(60)
    a = F(70)
    below = tuple(F(61 + i) for i in range(6 - r))
    above = tuple(F(71 + i) for i in range(r))
    packet = below + above
    require(len(packet) == 6 and tuple(sorted(packet)) == packet, "branch packet order")
    require(sum(d >= a for d in packet) == r, "branch packet suffix count")
    y0 = max(a, packet[3])
    gaps = tuple(packet[i + 1] - packet[i] for i in range(5))
    M = max(gaps)
    weight = 15 if r == 1 else 11 if r == 2 else 9
    require(l1(packet, y0) <= weight * M, "targeted branch weight")
    branch_histogram[r] += 1

require(F(1, 12) / 15 == F(1, 180), "one-pivot edge constant")
require(F(1, 12) / 11 == F(1, 132), "two-pivot edge constant")
require(F(1, 12) / 9 == F(1, 108), "three-pivot edge constant")

print("THM-1236 CONTINUOUS-PIVOT KAKEYA MEDIAN EXACT AUDIT")
print(f"ordered six-packets checked = {packet_count}")
print(f"carrier thresholds checked = {threshold_count}")
print(f"tilted-L1 breakpoint/midpoint comparisons = {sample_count}")
print("constrained optimizer = max(7c/6,d4)")
print("r=1, r=2, r>=3 path-weight sums = 15, 11, 9")
print("carrier edge floors = c/180, c/132, c/108")
print("branch histogram = " + ", ".join(f"r={r}:{branch_histogram[r]}" for r in range(1, 7)))
print("RESULT: PASS")
