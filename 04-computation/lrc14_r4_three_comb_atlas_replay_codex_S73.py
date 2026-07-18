#!/usr/bin/env python3
"""Independent Fraction replay for THM-1097's frozen certificate facts.

The C++ referee constructs safe sets by successive exact tooth subtraction.
This replay instead classifies every cell of the complete rational breakpoint
arrangement.  It independently verifies the 220-core atlas minimum, derives
the guarded-bank cardinality from the three exact cutoff formulae, and
reconstructs the hardest C++ row using fractions.Fraction throughout.
"""

from fractions import Fraction as F
from itertools import combinations


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def safe_set(core: tuple[int, ...]) -> list[tuple[F, F]]:
    breakpoints = {F(0), F(1)}
    for p in core:
        for j in range(p + 1):
            for sign in (-1, 1):
                point = F(j, p) + sign * F(1, 14 * p)
                if 0 <= point <= 1:
                    breakpoints.add(point)
    ordered = sorted(breakpoints)
    result: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(circle_distance(p * midpoint) >= F(1, 14) for p in core):
            if result and result[-1][1] == left:
                result[-1] = (result[-1][0], right)
            else:
                result.append((left, right))
    return result


def remove_bad(region: list[tuple[F, F]], speed: int) -> list[tuple[F, F]]:
    result: list[tuple[F, F]] = []
    radius = F(1, 14 * speed)
    for left, right in region:
        cursor = left
        j_lo = int(left * speed) - 1
        j_hi = int(right * speed) + 1
        for j in range(j_lo, j_hi + 1):
            bad_left = F(j, speed) - radius
            bad_right = F(j, speed) + radius
            if bad_right <= left or right <= bad_left:
                continue
            bad_left = max(left, bad_left)
            bad_right = min(right, bad_right)
            if cursor < bad_left:
                result.append((cursor, bad_left))
            cursor = max(cursor, bad_right)
        if cursor < right:
            result.append((cursor, right))
    return result


def longest(region: list[tuple[F, F]]) -> F:
    return max(right - left for left, right in region)


cores = list(combinations(range(1, 13), 9))
core_rows = [(longest(safe_set(core)), core) for core in cores]
min_ell = min(ell for ell, _ in core_rows)
min_cores = [core for ell, core in core_rows if ell == min_ell]
min_intervals = [
    interval
    for interval in safe_set(min_cores[0])
    if interval[1] - interval[0] == min_ell
]

# Count exactly the guarded finite domain used by the C++ certificate.
finite_triples = 0
max_k1 = max_k2 = max_k3 = 0
max_counted_k1 = max_counted_k2 = max_counted_k3 = 0
k1_boundaries = k2_boundaries = 0
min_first_omitted_k1_margin: F | None = None
min_first_omitted_k2_margin: F | None = None
min_first_omitted_k3_margin: F | None = None


def tail_margin(ell: F, k1: int, k2: int, k3: int) -> F:
    return (
        21 * ell * k3
        - 7 * ell * (k1 + k2)
        - F(6 * k3, k1)
        - F(6 * k3, k2)
        - 37
    )


for ell, core in core_rows:
    lo = 13 * max(core) + 1
    k1_max = int(F(7, 1) / ell) + 1
    max_k1 = max(max_k1, k1_max)
    first_omitted_k1_margin = ell * (k1_max + 1) - 7
    if (
        min_first_omitted_k1_margin is None
        or first_omitted_k1_margin < min_first_omitted_k1_margin
    ):
        min_first_omitted_k1_margin = first_omitted_k1_margin
    for k1 in range(lo, k1_max + 1):
        k1_boundaries += 1
        k2_denominator = 14 * ell * k1 - 6
        assert k2_denominator > 0
        k2_threshold = F(k1) * (7 * ell * k1 + 43) / k2_denominator
        k2_max = int(k2_threshold) + 1
        max_k2 = max(max_k2, k2_max)
        first_omitted_k2 = k2_max + 1
        first_omitted_k2_margin = tail_margin(
            ell, k1, first_omitted_k2, first_omitted_k2
        )
        assert first_omitted_k2_margin > 0
        if (
            min_first_omitted_k2_margin is None
            or first_omitted_k2_margin < min_first_omitted_k2_margin
        ):
            min_first_omitted_k2_margin = first_omitted_k2_margin
        for k2 in range(k1 + 1, k2_max + 1):
            k2_boundaries += 1
            k3_denominator = 21 * ell - F(6, k1) - F(6, k2)
            assert k3_denominator > 0
            k3_threshold = (7 * ell * (k1 + k2) + 37) / k3_denominator
            k3_max = int(k3_threshold) + 1
            max_k3 = max(max_k3, k3_max)
            first_omitted_k3 = k3_max + 1
            first_omitted_k3_margin = tail_margin(
                ell, k1, k2, first_omitted_k3
            )
            assert first_omitted_k3_margin > 0
            if (
                min_first_omitted_k3_margin is None
                or first_omitted_k3_margin < min_first_omitted_k3_margin
            ):
                min_first_omitted_k3_margin = first_omitted_k3_margin
            count = max(0, k3_max - k2)
            finite_triples += count
            if count:
                max_counted_k1 = max(max_counted_k1, k1)
                max_counted_k2 = max(max_counted_k2, k2)
                max_counted_k3 = max(max_counted_k3, k3_max)

# Independently reconstruct the hardest row reported by the exhaustive C++
# scan.  This does not assume the C++ core-safe endpoint construction.
extremal_core = (1, 2, 3, 6, 7, 8, 10, 11, 12)
k1, k2, k3 = 162, 168, 174
remainder = safe_set(extremal_core + (k1, k2, k3))
extremal_length = longest(remainder)
extremal_intervals = [
    interval for interval in remainder if interval[1] - interval[0] == extremal_length
]
metric = 7 * k3 * extremal_length
ratio = 1 / metric

assert len(cores) == 220
assert min_ell == F(11, 1008)
assert min_cores == [(1, 2, 3, 7, 8, 9, 10, 11, 12)]
assert min_intervals == [
    (F(29, 126), F(27, 112)),
    (F(85, 112), F(97, 126)),
]
assert finite_triples == 39_778_595
assert (max_k1, max_k2, max_k3) == (642, 642, 642)
assert (max_counted_k1, max_counted_k2, max_counted_k3) == (639, 640, 641)
assert (k1_boundaries, k2_boundaries) == (23_589, 1_238_741)
assert min_first_omitted_k1_margin == F(17, 1008)
assert min_first_omitted_k2_margin == F(31, 252)
assert min_first_omitted_k3_margin == F(28_429, 182_160)
assert extremal_length == F(1, 522)
assert metric == F(7, 3)
assert ratio == F(3, 7)
assert extremal_intervals == [
    (F(53, 252), F(517, 2436)),
    (F(1919, 2436), F(199, 252)),
]

print("THM-1097 independent Fraction atlas/cardinality replay")
print(f"cores={len(cores)}")
print(f"finite_triples_from_guards={finite_triples}")
print(f"max_guarded_k1={max_k1}")
print(f"max_guarded_k2={max_k2}")
print(f"max_guarded_k3={max_k3}")
print(f"max_counted_k1={max_counted_k1}")
print(f"max_counted_k2={max_counted_k2}")
print(f"max_counted_k3={max_counted_k3}")
print(f"k1_guard_boundaries={k1_boundaries}")
print(f"k2_guard_boundaries={k2_boundaries}")
print(f"min_first_omitted_k1_margin={min_first_omitted_k1_margin}")
print(f"min_first_omitted_k2_margin={min_first_omitted_k2_margin}")
print(f"min_first_omitted_k3_margin={min_first_omitted_k3_margin}")
print(f"min_core_ell={min_ell}")
for core in min_cores:
    print(f"min_ell_core={list(core)}")
for left, right in min_intervals:
    print(f"min_core_longest_interval=[{left},{right}]")
print(f"finite_bank_extremal_core={list(extremal_core)}")
print(f"finite_bank_extremal_triple=({k1},{k2},{k3})")
print(f"finite_bank_extremal_longest_component={extremal_length}")
print(f"finite_bank_extremal_7k3L={metric}")
print(f"finite_bank_extremal_R={ratio}")
for left, right in extremal_intervals:
    print(f"finite_bank_extremal_interval=[{left},{right}]")
print("replay=PASS")
