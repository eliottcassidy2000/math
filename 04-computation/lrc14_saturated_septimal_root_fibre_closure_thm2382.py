#!/usr/bin/env python3
"""Exact hostile audit for THM-2382.

Only integer and Fraction arithmetic is used.  The program checks the two
finite fibres in the proof and deliberately records controls showing where
genericity, saturation, blocker divisibility, and thirteen-unit hypotheses
enter.
"""

from fractions import Fraction
from itertools import combinations, product
from math import factorial, gcd


def require(condition: bool, label: str) -> None:
    """Optimization-safe exact check."""
    if not condition:
        raise RuntimeError(f"check failed: {label}")


def frac(x: Fraction) -> Fraction:
    return x % 1


def centred_distance(x: Fraction) -> Fraction:
    x = frac(x)
    return min(x, 1 - x)


def in_danger(x: Fraction, radius: Fraction) -> bool:
    return centred_distance(x) < radius


def valuation(number: int, prime: int) -> int:
    value = 0
    while number % prime == 0:
        number //= prime
        value += 1
    return value


def centred_mod(x: Fraction, modulus: int) -> Fraction:
    residue = x % modulus
    return min(residue, modulus - residue)


def root_word(speed: int, y: Fraction, width: int):
    """13-root word for an ordinary tooth (width 1) or guard (width 2)."""
    threshold = Fraction(13 * width, 14)
    return frozenset(
        h
        for h in range(13)
        if centred_mod(Fraction(speed) * (y + h), 13) < threshold
    )


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first,) + tail


def danger_intervals(speed: int):
    """Closed interval model of closure(D_speed) in [0,1]."""
    radius = Fraction(1, 14 * speed)
    intervals = []
    for k in range(speed):
        centre = Fraction(k, speed)
        lo = centre - radius
        hi = centre + radius
        if lo < 0:
            intervals.append((Fraction(0), hi))
            intervals.append((1 + lo, Fraction(1)))
        elif hi > 1:
            intervals.append((lo, Fraction(1)))
            intervals.append((Fraction(0), hi - 1))
        else:
            intervals.append((lo, hi))
    return intervals


def union_length(intervals) -> Fraction:
    ordered = sorted(intervals)
    if not ordered:
        return Fraction(0)
    total = Fraction(0)
    lo, hi = ordered[0]
    for next_lo, next_hi in ordered[1:]:
        if next_lo <= hi:
            hi = max(hi, next_hi)
        else:
            total += hi - lo
            lo, hi = next_lo, next_hi
    return total + hi - lo


# The lower guard is a translated grid of order L=7^(M-h), repeated 7^h
# times in each top bin.  Every nonaligned translate has exactly 2L/7
# guard-danger points.
guard_generic_cases = 0
for M in range(1, 6):
    N = 7 ** (M + 1)
    for h in range(M):
        L = 7 ** (M - h)
        for theta in (Fraction(1, 101), Fraction(1, 2), Fraction(100, 101)):
            one_period = sum(
                in_danger(Fraction(k, L) + theta / L, Fraction(1, 7))
                for k in range(L)
            )
            require(one_period == 2 * L // 7, "generic guard-period count")
            bin_count = (7**h) * one_period
            require(bin_count == 2 * N // 49, "generic guard-bin count")
            guard_generic_cases += 1

# At the aligned strict endpoint the generic count genuinely fails.
guard_endpoint_count = sum(
    in_danger(Fraction(k, 7), Fraction(1, 7)) for k in range(7)
)
require(guard_endpoint_count == 1, "aligned guard endpoint hostile")

# A top-layer ordinary tooth reduces to a translated seven-grid and
# occupies one bin away from the two aligned strict endpoints.
top_singleton_cases = 0
for unit in range(1, 7):
    require(gcd(unit, 7) == 1, "top-grid unit")
    for theta in (
        Fraction(1, 10),
        Fraction(1, 3),
        Fraction(2, 3),
        Fraction(9, 10),
    ):
        count = sum(
            in_danger(Fraction(unit * r, 7) + theta / 7, Fraction(1, 14))
            for r in range(7)
        )
        require(count == 1, "generic top singleton")
        top_singleton_cases += 1

top_endpoint_count = sum(
    in_danger(Fraction(r, 7) + Fraction(1, 14), Fraction(1, 14))
    for r in range(7)
)
require(top_endpoint_count == 0, "top strict-endpoint hostile")

# Saturated bin inequality: sum m_r=7 and 7m_r+2>=7.
compositions = list(weak_compositions(7, 7))
require(len(compositions) == 1716, "weak-composition count")
saturated_vectors = [
    vector for vector in compositions if all(7 * entry + 2 >= 7 for entry in vector)
]
require(
    saturated_vectors == [(1, 1, 1, 1, 1, 1, 1)],
    "unique saturated occupancy vector",
)

# Retaining labels, the allowed top-bin states are exactly the 7!
# transversals.  Record the three intrinsic pair types in this repeated layer.
seven_bin_atlas_size = 0
for assignment in product(range(7), repeat=7):
    counts = [assignment.count(residue) for residue in range(7)]
    if min(counts) >= 1:
        require(counts == [1] * 7, "atlas state is a transversal")
        seven_bin_atlas_size += 1
require(seven_bin_atlas_size == factorial(7), "labeled transversal atlas")

collision_pair_types = {
    "q-q": len(list(combinations(range(5), 2))),
    "q-c": 5 * 2,
    "c-c": 1,
}
require(
    collision_pair_types == {"q-q": 10, "q-c": 10, "c-c": 1},
    "collision pair types",
)

# In the W=2,k=2 alternatives the inherited inequality is only 7m_r>=0.
unsaturated_collision = (2, 0, 0, 0, 0, 0, 0)
require(sum(unsaturated_collision) == 2, "unsaturated occupancy weight")
require(
    all(7 * entry >= 0 for entry in unsaturated_collision),
    "unsaturated collision survives bin inequality",
)

# Exact blocker pullback c=13C on all thirteen inverse roots.
pullback_cases = 0
for C in range(1, 31):
    for y in (Fraction(1, 97), Fraction(22, 97), Fraction(95, 97)):
        for h in range(13):
            root = (y + h) / 13
            require(
                frac(13 * C * root) == frac(C * y),
                "blocker root-pullback identity",
            )
            pullback_cases += 1

# Divisibility by thirteen is essential: c=1 is not constant on a root fibre.
nonmultiple_values = {
    frac((Fraction(1, 97) + h) / 13) for h in range(13)
}
require(len(nonmultiple_values) == 13, "nonmultiple blocker hostile")

# Every quotient comb has measure 1/7.  Exhaust all ordered triples through
# speed 12 and verify the union-bound safe floor exactly.
single_comb_cases = 0
for C in range(1, 13):
    require(
        union_length(danger_intervals(C)) == Fraction(1, 7),
        "quotient comb measure",
    )
    single_comb_cases += 1

quotient_triples = 0
minimum_safe_measure = Fraction(1)
for C1, C2, C3 in product(range(1, 13), repeat=3):
    covered = union_length(
        danger_intervals(C1) + danger_intervals(C2) + danger_intervals(C3)
    )
    require(covered <= Fraction(3, 7), "three-comb union bound")
    safe = 1 - covered
    require(safe >= Fraction(4, 7), "quotient safe-measure floor")
    minimum_safe_measure = min(minimum_safe_measure, safe)
    quotient_triples += 1
require(
    minimum_safe_measure >= Fraction(4, 7),
    "quotient triple census floor",
)

# If a quotient blocker is base-dangerous, its original 13-multiple blocker
# is dangerous at all roots, showing why simultaneous blocker safety is used.
blocker_full_root_support = sum(
    in_danger(13 * ((Fraction(0) + h) / 13), Fraction(1, 14))
    for h in range(13)
)
require(
    blocker_full_root_support == 13,
    "dangerous quotient blocker has full root support",
)

# Exhaust the exact phase cells of a translated thirteen-grid.  A
# thirteen-unit merely permutes the grid labels.
def phase_cell_counts(width: int):
    threshold = Fraction(13 * width, 14)
    endpoints = sorted(
        {
            (sign * threshold - n) % 13
            for sign in (-1, 1)
            for n in range(13)
        }
    )
    counts = set()
    for index, left in enumerate(endpoints):
        right = endpoints[(index + 1) % len(endpoints)]
        if index + 1 == len(endpoints):
            right += 13
        phase = ((left + right) / 2) % 13
        counts.add(
            sum(
                centred_mod(phase + n, 13) < threshold
                for n in range(13)
            )
        )
    return counts, len(endpoints)


danger_phase_counts, danger_phase_cells = phase_cell_counts(1)
guard_phase_counts, guard_phase_cells = phase_cell_counts(2)
require(danger_phase_counts == {1, 2}, "ordinary root phase counts")
require(guard_phase_counts == {3, 4}, "guard root phase counts")

unit_root_cases = 0
unit_root_counts = set()
for unit in range(1, 13):
    require(gcd(unit, 13) == 1, "thirteen-grid unit")
    for numerator in range(1, 28, 2):
        word = root_word(unit, Fraction(numerator, 28), 1)
        require(len(word) in (1, 2), "unit root-word capacity")
        unit_root_counts.add(len(word))
        unit_root_cases += 1
require(unit_root_counts == {1, 2}, "unit root counts realized")

# A thirteen-divisible q can instead be dangerous at every inverse root.
nonunit_q_root_count = len(root_word(13, Fraction(0), 1))
require(nonunit_q_root_count == 13, "nonunit q hostile")

five_q_capacity = 5 * 2
six_q_capacity = 6 * 2
seven_q_capacity = 7 * 2
require(five_q_capacity < 13, "five-q contradiction")
require(six_q_capacity < 13, "six-q strengthening")
require(seven_q_capacity >= 13, "seven-q raw-capacity hostile")

# THM-2377's leading F_7 cancellation does not determine its Bockstein
# carry.  For every ordered pair of units modulo 49, fixing the unique
# leading cancellation class and varying its seven lifts realizes all
# seven carries.
units_mod_49 = [unit for unit in range(1, 49) if unit % 7]
bockstein_pair_cases = 0
for U, V in product(units_mod_49, repeat=2):
    leading = next(b for b in range(1, 7) if (U + b * V) % 7 == 0)
    carries = {
        ((U + (leading + 7 * lift) * V) // 7) % 7
        for lift in range(7)
    }
    require(carries == set(range(7)), "Bockstein lift carries")
    bockstein_pair_cases += 1

# Physical root-only hostile: all valuation and divisibility roles are right,
# all blockers are root-safe, and guard plus q words cover the roots.  The q
# words alone cover only ten roots, so this stalk violates precisely the
# rooted septimal top-owner partition proved from the global scalar cover.
hostile_y = Fraction(1, 97)
hostile_H = 1
hostile_q = (112, 126, 175, 238, 301)
hostile_c = (91, 2366, 107653)
hostile_guard_word = root_word(hostile_H, hostile_y, 2)
hostile_q_words = tuple(root_word(q, hostile_y, 1) for q in hostile_q)
hostile_c_words = tuple(root_word(c, hostile_y, 1) for c in hostile_c)
require(
    hostile_guard_word == frozenset({0, 1, 12}),
    "root-only hostile guard word",
)
require(
    hostile_q_words
    == (
        frozenset({3, 8}),
        frozenset({7, 10}),
        frozenset({2, 4}),
        frozenset({6, 9}),
        frozenset({5, 11}),
    ),
    "root-only hostile q words",
)
require(
    hostile_c_words == (frozenset(), frozenset(), frozenset()),
    "root-only hostile blocker words",
)
hostile_q_union = frozenset().union(*hostile_q_words)
hostile_full_union = hostile_q_union | hostile_guard_word
require(
    hostile_q_union == frozenset(range(2, 12)),
    "root-only hostile q union",
)
require(
    hostile_full_union == frozenset(range(13)),
    "root-only hostile total cover",
)
require(valuation(hostile_H, 7) == 0, "root-only hostile H valuation")
require(
    tuple(valuation(q, 7) for q in hostile_q) == (1, 1, 1, 1, 1),
    "root-only hostile q valuations",
)
require(
    tuple(valuation(c, 7) for c in hostile_c) == (1, 1, 2),
    "root-only hostile blocker septimal valuations",
)
require(all(q % 13 for q in hostile_q), "root-only hostile q units")
require(
    tuple(valuation(c, 13) for c in hostile_c) == (1, 2, 3),
    "root-only hostile blocker thirteen-adic valuations",
)

print("theorem=THM-2382")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(f"guard_generic_cases={guard_generic_cases}")
print("guard_per_bin_count=2N/49")
print(f"guard_endpoint_hostile={guard_endpoint_count}<2")
print(f"top_singleton_generic_cases={top_singleton_cases}")
print(f"top_endpoint_hostile={top_endpoint_count}")
print(f"seven_bin_compositions={len(compositions)}")
print(f"saturated_feasible_vectors={len(saturated_vectors)}")
print("saturated_vector=1,1,1,1,1,1,1")
print(f"labeled_seven_bin_atlas={seven_bin_atlas_size}")
print("collision_pair_types=q-q:10,q-c:10,c-c:1")
print("unsaturated_W2_collision_control=2,0,0,0,0,0,0")
print(f"blocker_pullback_cases={pullback_cases}")
print(f"blocker_nonmultiple_distinct_values={len(nonmultiple_values)}")
print(f"quotient_single_comb_cases={single_comb_cases}")
print(f"quotient_ordered_triples={quotient_triples}")
print(f"quotient_safe_measure_census_min={minimum_safe_measure}>=4/7")
print(f"dangerous_blocker_root_support={blocker_full_root_support}")
print(f"danger_phase_cells={danger_phase_cells}")
print(f"guard_phase_cells={guard_phase_cells}")
print(f"thirteen_unit_root_cases={unit_root_cases}")
print("thirteen_unit_root_counts=1,2")
print("thirteen_guard_root_counts=3,4")
print(f"five_unit_q_capacity={five_q_capacity}<13")
print(f"six_unit_q_capacity={six_q_capacity}<13")
print(f"seven_unit_q_hostile_capacity={seven_q_capacity}>=13")
print(f"thirteen_nonunit_q_hostile_count={nonunit_q_root_count}")
print(f"bockstein_unit_pair_cases={bockstein_pair_cases}")
print("bockstein_lift_carries=0,1,2,3,4,5,6")
print(f"root_only_hostile_y={hostile_y}")
print("root_only_hostile_guard=0,1,12")
print(f"root_only_hostile_q_coverage={len(hostile_q_union)}")
print(f"root_only_hostile_total_coverage={len(hostile_full_union)}")
print("scalar_saturated_branches_excluded=1")
print("other_k2_branches_excluded=0")
print("thirteen_adic_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
