#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2559.

The theorem is analytic on rational step packets.  This finite companion
checks its complete root combinatorics: target-informed selection for every
pair of disjoint nonempty 13-bit masks, the fixed-slope Cayley chord sign,
the blocker/word root-constancy laws, the exact ordinary/guard failure
capacities, and the sharp blind-mask control bypassed by the new selector.
"""

from collections import Counter
from fractions import Fraction
from math import comb


P = 13
FULL = (1 << P) - 1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def support(code):
    return {r for r in range(P) if (code >> r) & 1}


def cyclic_word_value(code, tau, start):
    value = 0
    for j in range(P):
        value = 2 * value + ((code >> ((start + j * tau) % P)) & 1)
    return value


def marker(code, tau=1):
    """THM-2531 lexicographic marker for a nonempty proper mask."""
    require(0 < code < FULL, "marker requires a mixed mask")
    return max(range(P), key=lambda a: cyclic_word_value(code, tau, a))


def selected_head(code, tau):
    source = marker(code, tau)
    q = next(
        j for j in range(1, P)
        if not ((code >> ((source + j * tau) % P)) & 1)
    )
    return (source + q * tau) % P


def cayley_profile(code, tau):
    """C_tau e, with (P_tau e)_v=e_(v+tau)."""
    return [
        sum(
            (1 if j % 2 else -1)
            * ((code >> ((v + j * tau) % P)) & 1)
            for j in range(1, P)
        )
        for v in range(P)
    ]


# Precompute the fixed orientation-gauge markers and all Cayley profiles.
markers = [None] * (FULL + 1)
cayley = [None] * (FULL + 1)
for code in range(1, FULL):
    markers[code] = marker(code)
    cayley[code] = {
        tau: cayley_profile(code, tau) for tau in range(1, P)
    }


print("== THM-2559: universal target-informed chord ==")
pair_count = 0
chord_checks = 0
slope_histogram = Counter()
source_weight_histogram = Counter()

for source_code in range(1, FULL):
    source_root = markers[source_code]
    complement = FULL ^ source_code
    target_code = complement
    while target_code:
        target_root = markers[target_code]
        slope = (target_root - source_root) % P
        require(slope != 0, "disjoint markers produced a zero chord")
        require((source_code >> source_root) & 1, "source marker is empty")
        require((target_code >> target_root) & 1, "target marker misses failure")
        require(not ((source_code >> target_root) & 1), "target lies in source")

        transformed = cayley[source_code][slope]
        require(
            transformed[source_root] + transformed[target_root] == -1,
            "Cayley chord has the wrong sign",
        )
        pair_count += 1
        chord_checks += 1
        slope_histogram[slope] += 1
        source_weight_histogram[source_code.bit_count()] += 1
        target_code = (target_code - 1) & complement

expected_pairs = 3**P - 2 * 2**P + 1
require(pair_count == expected_pairs == 1_577_940, "disjoint-pair census changed")
require(sum(slope_histogram.values()) == pair_count, "slope strata lost mass")
require(set(slope_histogram) == set(range(1, P)), "a chord slope disappeared")
for weight in range(1, P):
    expected = comb(P, weight) * (2 ** (P - weight) - 1)
    require(source_weight_histogram[weight] == expected,
            "source-weight pair census changed")

print(f"  disjoint nonempty source/failure mask pairs: {pair_count}")
print(f"  Cayley sign checks: {chord_checks}")
print("  nonzero chord-slope histogram: "
      + ",".join(f"{d}:{slope_histogram[d]}" for d in range(1, P)))
print("  every target-informed head is empty and every chord has p_t-p_s=-1")


print("\n== exact root-constancy and role-capacity interface ==")


def floor_q(value):
    return value.numerator // value.denominator


def fractional_part(value):
    return value - floor_q(value)


def circle_norm(value):
    residue = fractional_part(value)
    return min(residue, 1 - residue)


base_points = (Fraction(1, 17), Fraction(5, 29))
profiles = [
    (1, middle, deepest)
    for deepest in range(5, 20)
    for middle in range(1, deepest)
]
require(len(profiles) == 165, "valuation-profile ledger changed")

blocker_checks = 0
for first, middle, deepest in profiles:
    for depth, unit in zip((first, middle, deepest), (1, 5, 11)):
        coefficient = P**depth * unit
        for z in base_points:
            values = {
                fractional_part(coefficient * (z + root) / P)
                for root in range(P)
            }
            require(len(values) == 1, "13-divisible blocker moved on root fibre")
            blocker_checks += P

word_checks = 0
for clock in range(2, 7):
    for z in base_points:
        values = {
            fractional_part(P**clock * (z + root) / P)
            for root in range(P)
        }
        require(len(values) == 1, "positive-clock word moved on root fibre")
        word_checks += P

# THM-2379's exact translated-tooth identity on every cell of a common
# denominator-182 refinement, away from the finite endpoint walls.
capacity_checks = 0
count_values = {1: set(), 2: set()}
for level in (1, 2):
    for cell in range(182):
        y = Fraction(2 * cell + 1, 364)
        count = sum(
            circle_norm(y - Fraction(root, P)) < Fraction(level, 14)
            for root in range(P)
        )
        right = 2 * level - int(circle_norm(P * y) < Fraction(level, 14))
        require(count == right, "translated-tooth count identity failed")
        count_values[level].add(count)
        capacity_checks += 1
require(count_values[1] == {1, 2}, "ordinary role lost a sharp count")
require(count_values[2] == {3, 4}, "guard role lost a sharp count")

print(f"  blocker root-constancy checks: {blocker_checks}")
print(f"  terminal-word root-constancy checks: {word_checks}")
print(f"  translated-tooth identity checks: {capacity_checks}")
print("  ordinary target failure counts: 1 or 2; guard counts: 3 or 4")


print("\n== blind all-slope control and target-informed bypass ==")
blind_source = sum(1 << r for r in (0, 1, 4))
target_failure = 1 << 3
old_heads = {selected_head(blind_source, tau) for tau in range(1, P)}
require(old_heads == {2, 7, 8, 9, 11, 12}, "THM-2558 blind image changed")
require(not (old_heads & support(target_failure)), "old all-slope selector hit control")
source_root = markers[blind_source]
target_root = markers[target_failure]
target_slope = (target_root - source_root) % P
require((source_root, target_root, target_slope) == (0, 3, 3),
        "target-informed hostile chord changed")
transformed = cayley[blind_source][target_slope]
require(transformed[source_root] + transformed[target_root] == -1,
        "target-informed hostile chord lost its sign")
print("  e={0,1,4}; all-slope heads={2,7,8,9,11,12}; failure root 3 is blind")
print("  target-informed chord: 0 -> 3 at slope 3, Cayley endpoint sum -1")


print("\n== inherited THM-2379 exact floors ==")
denominator = P**2 * (P - 1) ** 2
ordinary_floor = Fraction(P - 2, denominator)
guard_floor = Fraction(P - 4, denominator)
require(ordinary_floor == Fraction(11, 24336), "ordinary repair floor changed")
require(guard_floor == Fraction(1, 2704), "guard repair floor changed")
print(f"  ordinary k_a floor per head mass: {ordinary_floor}")
print(f"  guard k_a floor per head mass: {guard_floor}")
print("  paired blocker remains a static sidecar; future semantic root remains absent")
print("\nall exact checks passed")
