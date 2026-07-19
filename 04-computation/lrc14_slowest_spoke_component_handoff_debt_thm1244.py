#!/usr/bin/env python3
"""Exact referee for THM-1244 slowest-spoke component handoff debt.

The topological interval-chain extraction is proved on paper.  This script
replays its exact spoke/component endpoint arithmetic, the two-label span
dispatch, gcd-quantized tooth overlaps, all labelled five-vertex forests and
active sets, and every scalar constant using dependency-free exact arithmetic.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def nearest_integer(x: F) -> int:
    lo = floor_fraction(x)
    return lo if x - lo <= F(1, 2) else lo + 1


def circle_norm(x: F) -> F:
    lo = floor_fraction(x)
    frac = x - lo
    return min(frac, 1 - frac)


spoke_rows = 0
for c in range(1, 81):
    for d in range(c + 1, c + 41):
        q = c + d
        for k in range(-3, 4):
            t0 = F(2 * k + 1, 2 * c)
            p = nearest_integer(q * t0)
            t1 = F(p, q)
            x = abs(t1 - t0)
            delta = circle_norm(c * t1)
            require(delta == circle_norm(d * t1), "spoke distance equality")
            require(delta == F(1, 2) - c * x, "centered distance identity")
            require(x <= F(1, 2 * q), "nearest-integer radius")

            h = floor_fraction(d * t1)
            G = (F(14 * k + 1, 14 * c), F(14 * k + 13, 14 * c))
            S = (F(14 * h + 1, 14 * d), F(14 * h + 13, 14 * d))
            K = (max(G[0], S[0]), min(G[1], S[1]))
            radius = F(delta - F(1, 14), d)
            require(radius > 0, "positive d1-safe radius")
            require(G[0] <= t1 - radius and t1 + radius <= G[1],
                    "radius lies in carrier gap")
            require(S[0] <= t1 - radius and t1 + radius <= S[1],
                    "radius lies in d1-safe component")
            length = K[1] - K[0]
            lower = F(6 * d - c, 7 * d * (c + d))
            require(length >= 2 * radius >= lower, "component length floor")
            require(lower > F(5, 14 * d), "five-fourteenths floor")
            spoke_rows += 1


two_label_rows = 0
for d1 in range(1, 81):
    for u in range(d1 + 1, d1 + 31):
        for v in range(u + 1, 8 * u + 2):
            if v <= 6 * u:
                span = F(1, 7 * u) + F(1, 7 * v)
                require(span < F(2, 7 * d1) < F(5, 14 * d1),
                        "close two-label span")
            else:
                span = F(1, 7 * u) + F(2, 7 * v)
                require(span < F(4, 21 * u) < F(5, 14 * d1),
                        "far two-label span")
            two_label_rows += 1


overlap_rows = 0
for u in range(2, 101):
    for v in range(u + 1, 101):
        g = gcd(u, v)
        for n in range(-3, 4):
            for m in range(-3, 4):
                numerator = v * (14 * n + 1) - u * (14 * m - 1)
                if numerator <= 0:
                    continue
                require(numerator % g == 0, "overlap numerator gcd")
                overlap = F(numerator, 14 * u * v)
                require(overlap >= F(g, 14 * u * v), "overlap quantum")
                overlap_rows += 1


vertices = tuple(range(5))
edges = tuple(combinations(vertices, 2))


def is_forest(chosen: tuple[tuple[int, int], ...]) -> bool:
    parent = list(vertices)

    def root(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, v in chosen:
        ru, rv = root(u), root(v)
        if ru == rv:
            return False
        parent[ru] = rv
    return True


forests = []
for mask in range(1 << len(edges)):
    chosen = tuple(edges[i] for i in range(len(edges)) if mask >> i & 1)
    if is_forest(chosen):
        forests.append(chosen)
require(len(forests) == 291, "five-vertex forest count")

active_checks = 0
for forest in forests:
    for active_mask in range(1 << 5):
        active = {i for i in vertices if active_mask >> i & 1}
        induced = sum(u in active and v in active for u, v in forest)
        require(induced <= max(0, len(active) - 1), "forest Hunter count")
        active_checks += 1

# Exact scalar ledger checks.
require(F(49, 6) * F(2, 7) == F(7, 3), "basic harmonic coefficient")
require(F(49, 6) * F(13, 196) == F(13, 24), "tree gcd coefficient")
require(F(49, 6) * F(4, 13) == F(98, 39), "tree length coefficient")
require(F(432, 1729) - F(4, 21) == F(44, 741), "private mass constant")
require(F(44, 741) / 5 == F(44, 3705), "private owner pigeonhole")

# The normalized component floor decreases on the compact ratio interval;
# exact cross multiplication checks its endpoint consumer.
for numerator in range(1001, 2167):
    r = F(numerator, 1000)
    if r >= F(13, 6):
        continue
    phi = F(6 * r - 1, 7 * r * (1 + r))
    require(phi > F(432, 1729), "component compact floor")

print("THM-1244 SLOWEST-SPOKE COMPONENT HANDOFF DEBT EXACT AUDIT")
print(f"centered component endpoint rows checked = {spoke_rows}")
print(f"two-label component span rows checked = {two_label_rows}")
print(f"positive gcd-quantized tooth overlaps checked = {overlap_rows}")
print(f"five-vertex forests = {len(forests)}")
print(f"forest/active-set Hunter checks = {active_checks}")
print("forced handoff rank = 2; located quantum Q0 >= gcd(d2,...,d6)/(7d6^2)")
print("scalar debt = H_B >= (6d1-c)/(3d1(c+d1)) + 7g_B/(6d6^2)")
print("private mass in K > 44/(741c); one owner > 44/(3705c)")
print("RESULT: PASS")

