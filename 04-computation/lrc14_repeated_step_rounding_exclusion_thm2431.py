#!/usr/bin/env python3
"""Exact companion for THM-2431.

The script reconstructs the normalized THM-2430 exact-cover atlas,
extracts every directed centre-difference bank for a repeated unsigned
step, thickens those banks by the nearest-integer addition defect, and
checks the resulting 50/91 versus 4/7 contradiction.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd


L = 91
FULL = (1 << L) - 1
GUARD = sum(1 << x for x in range(26))
COMPLEMENT = FULL ^ GUARD


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def progression(start: int, step: int, length: int = 13) -> int:
    mask = 0
    for j in range(length):
        mask |= 1 << ((start + j * step) % L)
    return mask


def points(mask: int) -> tuple[int, ...]:
    return tuple(x for x in range(L) if (mask >> x) & 1)


def nearest(value: Fraction) -> int:
    shifted = value + Fraction(1, 2)
    return shifted.numerator // shifted.denominator


def circle_norm(value: Fraction) -> Fraction:
    residue = value - (value.numerator // value.denominator)
    return min(residue, 1 - residue)


# Build the unoriented progression universe and retain its unique
# orientation with step in {1,...,45}.
representations: dict[int, list[tuple[int, int]]] = {}
for step in range(1, L):
    if gcd(step, L) != 1:
        continue
    for start in range(L):
        mask = progression(start, step)
        representations.setdefault(mask, []).append((start, step))

require(len(representations) == 3276, "wrong progression universe")

canonical_rep: dict[int, tuple[int, int]] = {}
for mask, reps in representations.items():
    forward = sorted((start, step) for start, step in reps if step <= L // 2)
    require(len(reps) == 2 and len(forward) == 1, "unexpected AP symmetry")
    canonical_rep[mask] = forward[0]

eligible = sorted(mask for mask in representations if mask & GUARD == 0)
require(len(eligible) == 182, "wrong eligible progression bank")

by_point: dict[int, list[int]] = {x: [] for x in points(COMPLEMENT)}
for mask in eligible:
    for x in points(mask):
        by_point[x].append(mask)


def exact_covers() -> set[frozenset[int]]:
    """Point-pivot Algorithm X with no unsafe index-order pruning."""

    answers: set[frozenset[int]] = set()

    def visit(remainder: int, chosen: tuple[int, ...]) -> None:
        if remainder == 0:
            require(len(chosen) == 5, "wrong cover size")
            answers.add(frozenset(chosen))
            return
        require(len(chosen) < 5, "five sets left a positive remainder")
        pivot = min(
            points(remainder),
            key=lambda x: sum(
                (candidate & remainder) == candidate
                for candidate in by_point[x]
            ),
        )
        for candidate in by_point[pivot]:
            if candidate & remainder == candidate:
                visit(remainder ^ candidate, chosen + (candidate,))

    visit(COMPLEMENT, ())
    return answers


solutions = exact_covers()
require(len(solutions) == 62, "wrong normalized cover count")


# An independent pair/triple meet-in-the-middle invoice sees each cover
# once for each of its binom(5,2)=10 choices of a two-set side.
pair_unions: Counter[int] = Counter()
for index, first in enumerate(eligible):
    for second in eligible[index + 1 :]:
        if first & second == 0:
            pair_unions[first | second] += 1

pair_triple_matches = 0
for first_index, first in enumerate(eligible):
    for second_index in range(first_index + 1, len(eligible)):
        second = eligible[second_index]
        if first & second:
            continue
        pair = first | second
        for third in eligible[second_index + 1 :]:
            if pair & third:
                continue
            pair_triple_matches += pair_unions[COMPLEMENT ^ (pair | third)]

require(pair_triple_matches == 620, "independent cover invoice failed")
require(pair_triple_matches == 10 * len(solutions), "cover invoices disagree")


# Direct physical controls for the intrinsic centre formula. On the
# common root fibre (Y+t)/91, the ordinary D_u window is the progression
# with signed step u^{-1} and centre -u^{-1} round(uY).
for speed in (1, 2, 3, 5, 8, 17, 37, 88):
    require(gcd(speed, L) == 1, "physical control speed is not a unit")
    inverse = pow(speed, -1, L)
    for parent in (Fraction(1, 17), Fraction(11, 37), Fraction(123, 257)):
        actual = {
            root
            for root in range(L)
            if circle_norm(Fraction(speed) * (parent + root) / L)
            < Fraction(1, 14)
        }
        predicted = {
            (-inverse * nearest(speed * parent) + offset * inverse) % L
            for offset in range(-6, 7)
        }
        require(len(actual) == 13, "physical control hit an endpoint")
        require(actual == predicted, "intrinsic centre formula failed")


# A length-thirteen progression A(s,d) has intrinsic centre s+6d.
# Collect directed differences between every pair of blocks having the
# same canonical unsigned step.
centre_differences: dict[int, set[int]] = defaultdict(set)
repeated_step_sets: set[tuple[int, ...]] = set()

for cover in solutions:
    grouped: dict[int, list[int]] = defaultdict(list)
    for mask in cover:
        start, step = canonical_rep[mask]
        grouped[step].append((start + 6 * step) % L)

    repeated = tuple(sorted(step for step, values in grouped.items() if len(values) >= 2))
    require(repeated, "a normalized cover has no repeated step")
    repeated_step_sets.add(repeated)

    for step, centres in grouped.items():
        if len(centres) < 2:
            continue
        for first in centres:
            for second in centres:
                if first != second:
                    centre_differences[step].add((first - second) % L)

repeatable_steps = sorted(centre_differences)
require(
    repeatable_steps == [1, 2, 3, 4, 5, 44, 45],
    "wrong repeatable-step bank",
)
require(30 not in centre_differences, "step 30 unexpectedly repeats")

expected_difference_sizes = {
    1: 28,
    2: 24,
    3: 4,
    4: 6,
    5: 8,
    44: 8,
    45: 4,
}
positive_difference_generators = {
    1: (13, 17, 18, 19, 20, 26, 30, 31, 32, 33, 39, 43, 44, 45),
    2: (1, 3, 25, 26, 27, 29, 30, 31, 32, 38, 39, 40),
    3: (1, 2),
    4: (1, 2, 3),
    5: (1, 2, 3, 4),
    44: (1, 2, 43, 45),
    45: (13, 39),
}
require(
    {step: len(values) for step, values in centre_differences.items()}
    == expected_difference_sizes,
    "wrong directed centre-difference sizes",
)
require(
    centre_differences
    == {
        step: {
            signed % L
            for value in generators
            for signed in (value, -value)
        }
        for step, generators in positive_difference_generators.items()
    },
    "wrong directed centre-difference banks",
)

for values in centre_differences.values():
    require(
        all((-value) % L in values for value in values),
        "directed difference bank is not symmetric",
    )


# If two fixed labels have equal unsigned step d, their physical
# centre difference forces
#
#   round(91*n*Y) in d^{-1} S_d + {-1,0,1}  (mod 91).
#
# The directed symmetry absorbs both the guard reflection and the
# signed orientation of the repeated pair.
expanded_banks: dict[int, set[int]] = {}
for step, differences in centre_differences.items():
    inverse = pow(step, -1, L)
    expanded_banks[step] = {
        (inverse * difference + error) % L
        for difference in differences
        for error in (-1, 0, 1)
    }

expected_expanded_sizes = {
    1: 50,
    2: 44,
    3: 8,
    4: 12,
    5: 16,
    44: 12,
    45: 12,
}
require(
    {step: len(values) for step, values in expanded_banks.items()}
    == expected_expanded_sizes,
    "wrong thickened residue-bank sizes",
)

max_bank = max(map(len, expanded_banks.values()))
require(max_bank == 50, "wrong maximum parent residue bank")
parentwise_union = set().union(*expanded_banks.values())
require(
    len(parentwise_union) == 66,
    "wrong union for the invalid parentwise-pair shortcut",
)


# Exact finite controls for the nearest-integer sum and difference
# defects. Half-integer ties are skipped because the theorem removes
# endpoint-null phases.
observed_sum_errors: set[int] = set()
observed_difference_errors: set[int] = set()
for denominator in range(3, 32):
    values = [
        Fraction(numerator, denominator)
        for numerator in range(-2 * denominator, 2 * denominator + 1)
    ]
    for first in values:
        for second in values:
            observed_sum_errors.add(
                nearest(first) + nearest(second) - nearest(first + second)
            )
            observed_difference_errors.add(
                nearest(first) - nearest(second) - nearest(first - second)
            )

require(observed_sum_errors <= {-1, 0, 1}, "sum rounding defect escaped")
require(
    observed_difference_errors <= {-1, 0, 1},
    "difference rounding defect escaped",
)
require(observed_sum_errors == {-1, 0, 1}, "sum controls miss a boundary")
require(
    observed_difference_errors == {-1, 0, 1},
    "difference controls miss a boundary",
)


# For nonzero n, nY is Haar-uniform. The nearest integer to 91*nY,
# reduced modulo 91, therefore has 91 fibres of mass exactly 1/91:
# residues 1,...,90 are single intervals, while residue zero is the
# union of the two endpoint half-intervals.
uniform_fibre_masses = {
    residue: Fraction(1, L) for residue in range(1, L)
}
uniform_fibre_masses[0] = Fraction(1, 2 * L) + Fraction(1, 2 * L)
require(
    set(uniform_fibre_masses.values()) == {Fraction(1, L)},
    "nearest-residue fibres are not uniform",
)
require(sum(uniform_fibre_masses.values()) == 1, "fibre masses do not sum to one")

parent_cap = Fraction(max_bank, L)
blocker_safe_floor = Fraction(4, 7)
gap = blocker_safe_floor - parent_cap
require(parent_cap == Fraction(50, 91), "wrong parent cap")
require(blocker_safe_floor == Fraction(52, 91), "wrong blocker-safe floor")
require(gap == Fraction(2, 91) and gap > 0, "mass contradiction failed")


m_zero_before = {
    (0, 5, 0, 7),
    (1, 5, 1, 8),
    (2, 5, 2, 9),
}
m_positive_before = {
    (1, 5, 0, 7),
    (2, 0, 0, 2),
    (2, 5, 0, 7),
    (2, 5, 1, 8),
}
m_zero_after = sorted(row for row in m_zero_before if not (row[1] == 5 and row[2] == 0))
m_positive_after = sorted(
    row for row in m_positive_before if not (row[1] == 5 and row[2] == 0)
)
require(m_zero_after == [(1, 5, 1, 8), (2, 5, 2, 9)], "wrong M=0 residual")
require(
    m_positive_after == [(2, 0, 0, 2), (2, 5, 1, 8)],
    "wrong M>0 residual",
)


print("THM-2431 exact companion")
print(f"normalized_tilings={len(solutions)}")
print(f"independent_pair_triple_invoice={pair_triple_matches}")
print("repeatable_steps=" + ",".join(map(str, repeatable_steps)))
for step in repeatable_steps:
    print(
        f"step_{step}:"
        f"centre_differences={len(centre_differences[step])},"
        f"expanded_parent_residues={len(expanded_banks[step])}"
    )
print(f"maximum_parent_cap={parent_cap}")
print(f"invalid_parentwise_union_cap={len(parentwise_union)}/{L}")
print(f"blocker_safe_image_floor={blocker_safe_floor}")
print(f"strict_mass_gap={gap}")
print("M0_remaining=" + ",".join(map(str, m_zero_after)))
print("Mpositive_remaining=" + ",".join(map(str, m_positive_after)))
print("ALL CHECKS PASSED")
