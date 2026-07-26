#!/usr/bin/env python3
"""Exact companion for THM-2430.

The finite universe consists of arithmetic progressions of length 13 in
Z/91Z whose common difference is a unit.  After normalizing the guard
progression to {0,...,25}, the script classifies all exact covers of its
complement by five such progressions.  Two independent enumerations are
used: a point-pivot exact-cover search and a pair/triple meet-in-the-middle
count.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import gcd


N = 91
FULL = (1 << N) - 1
GUARD = sum(1 << x for x in range(26))
COMPLEMENT = FULL ^ GUARD


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def progression(start: int, step: int, length: int = 13) -> int:
    mask = 0
    for j in range(length):
        mask |= 1 << ((start + j * step) % N)
    return mask


def points(mask: int) -> tuple[int, ...]:
    return tuple(x for x in range(N) if (mask >> x) & 1)


def affine_image(mask: int, multiplier: int, translate: int) -> int:
    answer = 0
    for x in points(mask):
        answer |= 1 << ((multiplier * x + translate) % N)
    return answer


def reflect(mask: int) -> int:
    return affine_image(mask, -1, 25)


# Every unoriented length-thirteen unit progression has exactly two
# oriented descriptions.  Retain the unique description with 1<=step<=45.
representations: dict[int, list[tuple[int, int]]] = {}
for step in range(1, N):
    if gcd(step, N) != 1:
        continue
    for start in range(N):
        mask = progression(start, step)
        representations.setdefault(mask, []).append((start, step))

require(len(representations) == 3276, "wrong full progression universe")

canonical_rep: dict[int, tuple[int, int]] = {}
for mask, reps in representations.items():
    forward = sorted((start, step) for start, step in reps if step <= N // 2)
    require(len(reps) == 2 and len(forward) == 1, "unexpected AP stabilizer")
    canonical_rep[mask] = forward[0]

eligible = sorted(
    mask for mask in representations if mask & GUARD == 0
)
require(len(eligible) == 182, "wrong normalized eligible bank")

by_point: dict[int, list[int]] = {x: [] for x in points(COMPLEMENT)}
for mask in eligible:
    for x in points(mask):
        by_point[x].append(mask)


def exact_cover_search() -> set[frozenset[int]]:
    """Point-pivot Algorithm X without any unsafe index-order pruning."""

    solutions: set[frozenset[int]] = set()

    def visit(remainder: int, chosen: tuple[int, ...]) -> None:
        if remainder == 0:
            require(len(chosen) == 5, "wrong cover cardinality")
            solutions.add(frozenset(chosen))
            return
        require(len(chosen) < 5, "positive remainder after five sets")
        remaining_points = points(remainder)
        pivot = min(
            remaining_points,
            key=lambda x: sum(
                (candidate & remainder) == candidate
                for candidate in by_point[x]
            ),
        )
        for candidate in by_point[pivot]:
            if candidate & remainder == candidate:
                visit(remainder ^ candidate, chosen + (candidate,))

    visit(COMPLEMENT, ())
    return solutions


solutions = exact_cover_search()
require(len(solutions) == 62, "wrong exact-cover count")


# Independent meet-in-the-middle count.  A five-set cover is seen once
# for each choice of its two-set side, hence exactly binom(5,2)=10 times.
pair_unions: Counter[int] = Counter()
for i, first in enumerate(eligible):
    for second in eligible[i + 1 :]:
        if first & second == 0:
            pair_unions[first | second] += 1

pair_triple_matches = 0
for i, first in enumerate(eligible):
    for j in range(i + 1, len(eligible)):
        second = eligible[j]
        if first & second:
            continue
        pair = first | second
        for third in eligible[j + 1 :]:
            if pair & third:
                continue
            triple = pair | third
            pair_triple_matches += pair_unions[COMPLEMENT ^ triple]

require(pair_triple_matches == 620, "meet-in-the-middle mismatch")
require(pair_triple_matches // 10 == len(solutions), "cover count mismatch")


# The affine stabilizer of the normalized guard is identity plus reflection.
guard_stabilizer = []
for multiplier in range(N):
    if gcd(multiplier, N) != 1:
        continue
    for translate in range(N):
        if affine_image(GUARD, multiplier, translate) == GUARD:
            guard_stabilizer.append((multiplier, translate))
require(guard_stabilizer == [(1, 0), (90, 25)], "wrong guard stabilizer")


def cover_signature(cover: frozenset[int]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(canonical_rep[mask] for mask in cover))


canonical_orbits: set[tuple[tuple[int, int], ...]] = set()
reflection_fixed = 0
for cover in solutions:
    reflected = frozenset(reflect(mask) for mask in cover)
    require(reflected in solutions, "reflection left solution bank")
    if reflected == cover:
        reflection_fixed += 1
    canonical_orbits.add(
        min(cover_signature(cover), cover_signature(reflected))
    )

require(reflection_fixed == 6, "wrong reflection-fixed count")
require(len(canonical_orbits) == 34, "wrong affine-orbit count")


step_patterns = Counter(
    tuple(sorted(step for _, step in signature))
    for signature in canonical_orbits
)
expected_patterns = Counter(
    {
        (1, 1, 1, 1, 1): 1,
        (1, 1, 1, 1, 30): 2,
        (1, 1, 1, 1, 45): 2,
        (1, 1, 1, 2, 2): 2,
        (1, 1, 2, 2, 30): 3,
        (1, 1, 2, 2, 45): 4,
        (1, 1, 3, 3, 3): 2,
        (1, 1, 44, 44, 44): 2,
        (1, 1, 45, 45, 45): 1,
        (1, 2, 2, 2, 2): 2,
        (1, 2, 2, 4, 4): 2,
        (1, 3, 3, 3, 45): 2,
        (1, 4, 4, 4, 4): 1,
        (2, 2, 2, 2, 30): 2,
        (2, 2, 3, 3, 3): 1,
        (2, 2, 44, 44, 44): 3,
        (2, 2, 45, 45, 45): 1,
        (5, 5, 5, 5, 5): 1,
    }
)
require(step_patterns == expected_patterns, "wrong step-spectrum bank")

step_classes = sorted({step for pattern in step_patterns for step in pattern})
ratio_classes = sorted(
    {
        min(pow(step, -1, N), N - pow(step, -1, N))
        for step in step_classes
    }
)
require(step_classes == [1, 2, 3, 4, 5, 30, 44, 45], "wrong steps")
require(
    ratio_classes == [1, 2, 3, 18, 23, 30, 31, 45],
    "wrong signed speed-ratio classes",
)


# THM-2427's positive-chamber hostile is the first affine orbit.
hostile_guard = sum(1 << x for x in (*range(13), *range(78, 91)))
hostile_words = frozenset(
    sum(1 << x for x in range(13 * j, 13 * (j + 1)))
    for j in range(1, 6)
)
normalized_guard = affine_image(hostile_guard, 1, 13)
normalized_words = frozenset(
    affine_image(mask, 1, 13) for mask in hostile_words
)
require(normalized_guard == GUARD, "hostile guard normalization failed")
require(normalized_words in solutions, "hostile word cover absent")
require(
    cover_signature(normalized_words)
    == ((26, 1), (39, 1), (52, 1), (65, 1), (78, 1)),
    "wrong hostile orbit",
)


def circle_norm(value: Fraction) -> Fraction:
    residue = Fraction(value.numerator % value.denominator, value.denominator)
    return min(residue, 1 - residue)


def root_mask(speed: int, parent: Fraction, threshold: Fraction) -> set[int]:
    return {
        root
        for root in range(13)
        if circle_norm(speed * (parent + root) / 13) < threshold
    }


def check_thirteen_partition(
    guard_speed: int, ordinary_speeds: tuple[int, ...], singleton_count: int
) -> None:
    parent = Fraction(1, 2)
    guard = root_mask(guard_speed, parent, Fraction(1, 7))
    words = [
        root_mask(speed, parent, Fraction(1, 14))
        for speed in ordinary_speeds
    ]
    require(len(guard) in (3, 4), "wrong guard root count")
    require(
        sum(len(word) == 1 for word in words) == singleton_count,
        "wrong singleton count",
    )
    require(
        sum(map(len, words)) + len(guard) == 13,
        "wrong exact-partition incidence",
    )
    require(
        guard.isdisjoint(set().union(*words))
        and sum(map(len, words)) == len(set().union(*words)),
        "thirteen-root masks overlap",
    )
    require(guard | set().union(*words) == set(range(13)), "root gap")


check_thirteen_partition(1, (2, 3, 5, 11, 19), singleton_count=1)
check_thirteen_partition(2, (3, 5, 9, 11, 19), singleton_count=0)


print("THM-2430 exact companion")
print(f"all_unit_ap13={len(representations)}")
print(f"eligible_ap13={len(eligible)}")
print(f"exact_normalized_tilings={len(solutions)}")
print(f"pair_triple_matches={pair_triple_matches}")
print(f"guard_stabilizer={guard_stabilizer}")
print(f"reflection_fixed_tilings={reflection_fixed}")
print(f"affine_tiling_orbits={len(canonical_orbits)}")
print(f"unsigned_step_multisets={len(step_patterns)}")
print("signed_ratio_classes=" + ",".join(map(str, ratio_classes)))
print("thirteen_slice_exact_types=(g,s,E)=(3,0,0),(4,1,0)")
for index, signature in enumerate(sorted(canonical_orbits), 1):
    encoded = " ".join(f"{start}:{step}" for start, step in signature)
    print(f"orbit_{index:02d}={encoded}")
print("thm2427_hostile_orbit=PASS")
print("ALL CHECKS PASSED")
