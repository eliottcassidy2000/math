#!/usr/bin/env python3
"""Exact referee for THM-2197's scalar chord tournament no-go.

The audit checks:

* the genuine y=3/17 scalar masks in the covered and uncovered controls;
* rooted reflection symmetry of the covered carrier;
* absence of any tournament on the five mask vertices invariant under that
  reflection;
* equality of the two length profiles and inequality of their deficiencies;
* the Boolean deficiency product and singleton-test separation on all 512
  safe-sheet subsets.

Integer and Fraction arithmetic only; runtime checks survive ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


P = 13
Y = Fraction(3, 17)
ROOTS = tuple((Y + k) / P for k in range(P))
FORBIDDEN = frozenset({11, 12, 0, 1})
SAFE = frozenset(range(P)) - FORBIDDEN


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def danger_mask(coefficient: int, radius: Fraction) -> frozenset[int]:
    return frozenset(
        k
        for k, root in enumerate(ROOTS)
        if circle_norm(coefficient * root) < radius
    )


def chord_length(mask: frozenset[int]) -> int:
    require(len(mask) == 2, "chord length requested for a non-chord")
    left, right = sorted(mask)
    step = (right - left) % P
    return min(step, P - step)


def deficiency(masks: tuple[frozenset[int], ...]) -> frozenset[int]:
    covered = frozenset().union(*masks)
    return SAFE - covered


def reflection(sheet: int) -> int:
    return (12 - sheet) % P


def reflected(mask: frozenset[int]) -> frozenset[int]:
    return frozenset(reflection(sheet) for sheet in mask)


def invariant_tournament_count(permutation: tuple[int, ...]) -> int:
    """Enumerate labelled tournaments fixed by ``permutation``."""
    pairs = tuple(combinations(range(len(permutation)), 2))
    fixed = 0
    for bits in range(1 << len(pairs)):
        arcs = set()
        for index, (left, right) in enumerate(pairs):
            arcs.add((left, right) if bits & (1 << index) else (right, left))
        image = {
            (permutation[left], permutation[right])
            for left, right in arcs
        }
        if image == arcs:
            fixed += 1
    return fixed


def all_subsets(universe: tuple[int, ...]):
    for bits in range(1 << len(universe)):
        yield frozenset(
            universe[index]
            for index in range(len(universe))
            if bits & (1 << index)
        )


def main() -> None:
    guard_bad = danger_mask(1, Fraction(1, 7))
    require(guard_bad == FORBIDDEN, "guard block changed")
    require(
        all(
            circle_norm(multiplier * Y) > Fraction(1, 14)
            for multiplier in (1, 2, 3)
        ),
        "divided-deep hostile control became active",
    )

    covered_coefficients = (40, 12, 25, 92, 105)
    covered_masks = tuple(
        danger_mask(coefficient, Fraction(1, 14))
        for coefficient in covered_coefficients
    )
    expected_covered = (
        frozenset({6}),
        frozenset({2, 3}),
        frozenset({4, 5}),
        frozenset({9, 10}),
        frozenset({7, 8}),
    )
    require(covered_masks == expected_covered, "covered scalar masks changed")
    require(deficiency(covered_masks) == frozenset(), "covered control misses")

    reflection_action = tuple(
        covered_masks.index(reflected(mask))
        for mask in covered_masks
    )
    require(
        reflection_action == (0, 3, 4, 1, 2),
        "rooted reflection action changed",
    )
    require(
        invariant_tournament_count(reflection_action) == 0,
        "an invariant mask tournament unexpectedly exists",
    )

    uncovered_coefficients = (40, 12, 25, 105, 27)
    uncovered_masks = tuple(
        danger_mask(coefficient, Fraction(1, 14))
        for coefficient in uncovered_coefficients
    )
    expected_uncovered = (
        frozenset({6}),
        frozenset({2, 3}),
        frozenset({4, 5}),
        frozenset({7, 8}),
        frozenset({8, 9}),
    )
    require(
        uncovered_masks == expected_uncovered,
        "uncovered scalar masks changed",
    )
    require(
        deficiency(uncovered_masks) == frozenset({10}),
        "uncovered deficiency changed",
    )

    covered_profile = tuple(
        sorted([1] + [chord_length(mask) for mask in covered_masks[1:]])
    )
    uncovered_profile = tuple(
        sorted([1] + [chord_length(mask) for mask in uncovered_masks[1:]])
    )
    require(
        covered_profile == uncovered_profile == (1, 1, 1, 1, 1),
        "same-profile hostile control changed",
    )

    universe = tuple(sorted(SAFE))
    subsets = tuple(all_subsets(universe))
    product_checks = 0
    separation_checks = 0
    for left in subsets:
        left_deficiency = SAFE - left
        for right in subsets:
            right_deficiency = SAFE - right
            require(
                SAFE - (left | right)
                == left_deficiency & right_deficiency,
                "Boolean deficiency product failed",
            )
            product_checks += 1
            if left == right:
                continue
            witness = next(iter(left ^ right))
            context = SAFE - {witness}
            left_covers = (left | context) == SAFE
            right_covers = (right | context) == SAFE
            require(
                left_covers != right_covers,
                "singleton-test continuation did not separate states",
            )
            separation_checks += 1

    require(len(subsets) == 512, "deficiency state count changed")
    require(product_checks == 512 * 512, "product check count changed")
    require(separation_checks == 512 * 511, "separation check count changed")

    print("THM-2197 scalar chord deficiency/tournament audit")
    print(f"guard_forbidden={tuple(sorted(FORBIDDEN))}")
    print(f"covered_masks={tuple(tuple(sorted(mask)) for mask in covered_masks)}")
    print(f"reflection_action={reflection_action}")
    print("reflection_invariant_tournaments=0")
    print(f"uncovered_masks={tuple(tuple(sorted(mask)) for mask in uncovered_masks)}")
    print(f"profiles={covered_profile}; deficiencies=((),(10,))")
    print(
        f"deficiency_states={len(subsets)}; "
        f"product_checks={product_checks}; "
        f"separation_checks={separation_checks}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
