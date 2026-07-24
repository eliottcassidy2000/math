#!/usr/bin/env python3
"""Exact controls for THM-2143 (no floating decisions, no large build)."""

from fractions import Fraction as F
from functools import reduce
from math import gcd
from random import Random


def circle_distance(x):
    x %= 1
    return min(x, 1 - x)


def depth_distribution(speeds, h):
    """Exact p_k from the rational danger-comb endpoint arrangement."""
    endpoints = {F(0), F(1)}
    for v in speeds:
        for m in range(v + 1):
            for sign in (-1, 1):
                x = (F(m) + sign * h) / v
                if 0 < x < 1:
                    endpoints.add(x)
    endpoints = sorted(endpoints)
    p = [F(0) for _ in range(len(speeds) + 1)]
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        depth = sum(circle_distance(v * midpoint) < h for v in speeds)
        p[depth] += right - left
    assert sum(p) == 1
    return p


def partition(p, z):
    return sum(mass * z**k for k, mass in enumerate(p))


def artanh_lower(t):
    return 2 * (t + t**3 / 3 + t**5 / 5)


def artanh_upper(t):
    return 2 * (t + t**3 / 3 + t**5 / (5 * (1 - t**2)))


print("[1] owner-supplied rational certificate")
t_a = F(389, 2181)
t_b = F(5872957, 11821757)
coefficient = F(2457, 6592)
margin = coefficient * artanh_lower(t_b) - artanh_upper(t_a) - F(1, 25)
claimed_margin = F(
    391926968594914200867482400554891567498742649630277,
    82738859282193417287303438726081463937219800938169600,
)
assert margin == claimed_margin > 0
print(f"    exact margin reproduced: {margin}")

print("\n[2] exact n=3 Gibbs hostile/positive pair")
z = F(1, 10)
p_tight = depth_distribution([1, 2], F(1, 3))
p_loose = depth_distribution([1, 3], F(1, 3))
assert p_tight == [F(0), F(2, 3), F(1, 3)]
assert p_loose == [F(1, 9), F(4, 9), F(4, 9)]
z_tight = partition(p_tight, z)
z_loose = partition(p_loose, z)
assert z_tight == F(7, 100) < z
assert z_loose == F(4, 25) > z
cayley_t = (z_loose - z) / (z_loose + z)
assert z_loose / z == F(8, 5) and cayley_t == F(3, 13)
print(f"    tight V=[1,2]: p={p_tight}, Z(1/10)={z_tight} < 1/10")
print(f"    loose V=[1,3]: p={p_loose}, Z(1/10)={z_loose} > 1/10")
print(f"    Cayley parameter=(Z-z)/(Z+z)={cayley_t}; Z/z=8/5")

print("\n[3] deterministic h=3/41 defect-seven temperature battery")
rng = Random(2143)
h = F(3, 41)
z_a = F(896, 1285)
z_b = F(2974400, 8847357)
temperatures = (z_a, z_b, F(1, 10))
rows = []
while len(rows) < 8:
    small = sorted(rng.sample(range(1, 14), 6))
    far = sorted(rng.sample(range(14, 181), 7))
    speeds = small + far
    if reduce(gcd, speeds) != 1:
        continue
    p = depth_distribution(speeds, h)
    fires = tuple(partition(p, q) > q for q in temperatures)
    if fires == (False, False, True):
        rows.append((speeds, p[0]))

for index, (speeds, p0) in enumerate(rows, 1):
    print(f"    row {index}: p0={p0} ({float(p0):.9f})  V={speeds}")
print("    fires z_A, z_B, 1/10: 0/8, 0/8, 8/8")

print("\nAll THM-2143 exact controls passed.")
