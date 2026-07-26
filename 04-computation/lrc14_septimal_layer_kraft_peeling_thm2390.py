#!/usr/bin/env python3
"""Exact companion for THM-2390.

All arithmetic relevant to the quantitative bound is rational/integer.
The finite word exhaustions deliberately use arbitrary guard pairs and
singleton positions, a larger universe than the geometrically realizable
seven-root words.
"""

from fractions import Fraction
from itertools import combinations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mod_one(x):
    return x - (x.numerator // x.denominator)


def circle_norm(x):
    y = mod_one(x)
    return min(y, 1 - y)


def danger_count(a, unit, theta):
    return sum(
        circle_norm(theta + Fraction(unit * k, 7)) < Fraction(a, 14)
        for k in range(7)
    )


def cyclic_midpoints(points):
    pts = sorted(points)
    for left, right in zip(pts, pts[1:] + [pts[:1][0] + 1]):
        yield mod_one((left + right) / 2)


def integer_partitions(total, cap, top=None):
    if total == 0:
        yield ()
        return
    if top is None:
        top = min(total, cap)
    for first in range(min(total, top, cap), 0, -1):
        for tail in integer_partitions(total - first, cap, first):
            yield (first,) + tail


def layer_product(parts):
    ans = Fraction(1)
    for weight in parts:
        ans *= Fraction(7 - weight, 7)
    return ans


print("THM-2390 SEPTIMAL LAYER KRAFT PEELING -- exact companion")

# Exact boundary chambers for every permitted width and every seven-unit
# slope.  An endpoint belongs to both translated endpoint lists, so the
# strict open arc loses one root there.
root_profiles = {}
for a in range(1, 7):
    generic_counts = set()
    endpoint_counts = set()
    for unit in range(1, 7):
        endpoints = {
            mod_one(Fraction(sign * a, 14) - Fraction(unit * k, 7))
            for sign in (-1, 1)
            for k in range(7)
        }
        require(len(endpoints) == 7, f"unexpected endpoint orbit a={a}, u={unit}")
        endpoint_counts.update(danger_count(a, unit, theta) for theta in endpoints)
        generic_counts.update(
            danger_count(a, unit, theta) for theta in cyclic_midpoints(endpoints)
        )
    require(generic_counts == {a}, f"generic root count failed for a={a}")
    require(endpoint_counts == {a - 1}, f"endpoint root count failed for a={a}")
    root_profiles[a] = (a - 1, a)

print(
    "seven-root open-arc counts (aligned,generic): "
    + ", ".join(f"a={a}:{root_profiles[a]}" for a in range(1, 7))
)

# The one-peel capacity inequality is stronger than needed at endpoints.
capacity = []
for weight in range(7):
    min_safe = 7 - weight
    require(min_safe >= 1, "positive peel factor lost below weight seven")
    capacity.append((weight, min_safe, Fraction(min_safe, 7)))
print(
    "one-peel safe capacities W=0..6: "
    + ", ".join(f"{w}->{m} ({f})" for w, m, f in capacity)
)

# Exhaust every unordered integer partition of the lower weight eight.
partitions = sorted(set(integer_partitions(8, 6)), reverse=True)
partition_products = [(parts, layer_product(parts)) for parts in partitions]
minimum = min(value for _, value in partition_products)
minimizers = [parts for parts, value in partition_products if value == minimum]
require(minimum == Fraction(5, 49), "wrong lower partition minimum")
require(minimizers == [(6, 2)], "lower partition minimum is not unique")

print(f"lower-weight-8 partitions checked: {len(partitions)}")
print(f"lower product minimum: {minimum} uniquely at {minimizers[0]}")

full_floor = Fraction(6, 7) ** 2 * minimum
require(full_floor == Fraction(180, 2401), "wrong full common-safe floor")
print(f"two singleton higher layers times lower minimum: {full_floor}")

# The guard has weight two, while the six lower ordinary labels have weight
# one.  A heavy layer therefore contains the guard and five or six ordinary
# labels.
weight7_roles = list(combinations(range(6), 5))
weight8_roles = [tuple(range(6))]
require(len(weight7_roles) == 6, "wrong weight-seven role count")
require(len(weight8_roles) == 1, "wrong weight-eight role count")
print("heavy labelled role choices: W=7 -> 6; W=8 -> 1")


def multiplicities(guard_pair, singleton_positions):
    values = [0] * 7
    for root in guard_pair:
        values[root] += 1
    for root in singleton_positions:
        values[root] += 1
    return tuple(values)


# Exhaust a deliberately enlarged abstract universe.  The guard pair is any
# two-subset, not merely an adjacent pair in one labelled guard coordinate.
w7_total = 0
w7_covers = 0
for guard_pair in combinations(range(7), 2):
    for singleton_positions in product(range(7), repeat=5):
        w7_total += 1
        mult = multiplicities(guard_pair, singleton_positions)
        if min(mult) > 0:
            w7_covers += 1
            require(mult == (1,) * 7, "weight-seven cover is not a partition")

require(w7_total == 21 * 7**5, "wrong W=7 abstract universe")
require(w7_covers == 21 * 120, "wrong W=7 covering count")

w8_total = 0
w8_covers = 0
for guard_pair in combinations(range(7), 2):
    for singleton_positions in product(range(7), repeat=6):
        w8_total += 1
        mult = multiplicities(guard_pair, singleton_positions)
        if min(mult) > 0:
            w8_covers += 1
            require(
                sorted(mult) == [1, 1, 1, 1, 1, 1, 2],
                "weight-eight cover does not have exactly one double",
            )

require(w8_total == 21 * 7**6, "wrong W=8 abstract universe")
require(w8_covers == 21 * 3240, "wrong W=8 covering count")

print(f"W=7 abstract labelled words: {w7_total}; full covers: {w7_covers}")
print(f"W=8 abstract labelled words: {w8_total}; full covers: {w8_covers}")
print("terminal profiles: W=7 -> 1^7; W=8 -> 2,1^6")
print("VERDICT: no-heavy-layer gives floor 180/2401; every survivor has W=7 or W=8")
