#!/usr/bin/env python3
"""Exact companion for THM-2427.

Checks the finite thirteen-root capacity used to prune the guard-top,
deep-c3 branch, the complete (k,t,b,W) enumeration, its primitive
M>0 reduction, and positive abstract and physical common-phase t=5
root-cover controls.

Only standard-library integer and Fraction arithmetic is used.
"""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle_norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def root_mask(offset, width):
    """Danger roots for ||(offset+j)/13|| < width/14."""

    return frozenset(
        root
        for root in range(P)
        if circle_norm((offset + root) / P) < Fraction(width, 14)
    )


def physical_root_mask(speed, width, base):
    """Danger roots of one actual integer speed at a common base."""

    return frozenset(
        root
        for root in range(P)
        if circle_norm(Fraction(speed) * (base + root) / P)
        < Fraction(width, 14)
    )


def base_root_masks(width):
    """All masks between, and on, the exact offset breakpoints."""

    boundaries = {Fraction(0), Fraction(1)}
    for root in range(P):
        candidates = (
            Fraction(P * width, 14) - root,
            Fraction(P * (14 - width), 14) - root,
        )
        boundaries.update(
            candidate
            for candidate in candidates
            if 0 <= candidate <= 1
        )
    ordered = sorted(boundaries)
    samples = set(ordered)
    samples.update(
        (left + right) / 2
        for left, right in zip(ordered, ordered[1:])
    )
    return {
        root_mask(sample, width)
        for sample in samples
        if sample < 1
    }


def affine_closure(masks):
    """All unit-affine relabellings of Z/13Z."""

    return {
        frozenset((unit * root + shift) % P for root in mask)
        for mask in masks
        for unit in range(1, P)
        for shift in range(P)
    }


def main():
    ordinary_masks = affine_closure(base_root_masks(1))
    guard_masks = affine_closure(base_root_masks(2))

    ordinary_sizes = sorted({len(mask) for mask in ordinary_masks})
    guard_sizes = sorted({len(mask) for mask in guard_masks})
    require(ordinary_sizes == [1, 2], "ordinary root sizes changed")
    require(guard_sizes == [3, 4], "guard root sizes changed")
    require(max(map(len, ordinary_masks)) == 2, "ordinary cap is not two")
    require(max(map(len, guard_masks)) == 4, "guard cap is not four")

    for top_ordinary in range(5):
        require(
            4 + 2 * top_ordinary < P,
            f"t={top_ordinary} unexpectedly reaches thirteen roots",
        )
    require(4 + 2 * 5 >= P, "t=5 capacity boundary changed")

    # A positive abstract root-word control: the capacity inequality is
    # sharp as a pruning threshold and must not delete t=5.
    ordinary_pairs = sorted(
        (mask for mask in ordinary_masks if len(mask) == 2),
        key=lambda mask: tuple(sorted(mask)),
    )
    guard_fours = sorted(
        (mask for mask in guard_masks if len(mask) == 4),
        key=lambda mask: tuple(sorted(mask)),
    )
    cover_witness = None
    for guard in guard_fours:
        for choices in combinations_with_replacement(ordinary_pairs, 5):
            covered = set(guard)
            for choice in choices:
                covered.update(choice)
            if len(covered) == P:
                cover_witness = (guard, choices)
                break
        if cover_witness is not None:
            break
    require(cover_witness is not None, "abstract t=5 control disappeared")

    # A genuine common-phase root partition.  Thus even keeping one
    # common base and 91-unit integer speeds cannot eliminate t=5
    # without septimal top-bin/owner/valuation data.
    physical_base = Fraction(1, 2)
    physical_guard_speed = 1
    physical_ordinary_speeds = (2, 3, 5, 11, 19)
    require(
        all(
            gcd(speed, 91) == 1
            for speed in (physical_guard_speed,) + physical_ordinary_speeds
        ),
        "physical t=5 control lost a 91-unit speed",
    )
    physical_guard = physical_root_mask(
        physical_guard_speed,
        2,
        physical_base,
    )
    physical_ordinary = tuple(
        physical_root_mask(speed, 1, physical_base)
        for speed in physical_ordinary_speeds
    )
    require(
        physical_guard == frozenset((0, 1, 11, 12)),
        "physical t=5 guard mask changed",
    )
    require(
        physical_ordinary
        == (
            frozenset((6,)),
            frozenset((4, 8)),
            frozenset((2, 10)),
            frozenset((3, 9)),
            frozenset((5, 7)),
        ),
        "physical t=5 ordinary masks changed",
    )
    physical_parts = (physical_guard,) + physical_ordinary
    require(
        all(
            physical_parts[left].isdisjoint(physical_parts[right])
            for left in range(len(physical_parts))
            for right in range(left)
        ),
        "physical t=5 masks stopped being disjoint",
    )
    require(
        set().union(*physical_parts) == set(range(P)),
        "physical t=5 masks stopped covering all roots",
    )

    # THM-2367 gives W<=k or W>=7.  The new theorem adds t=5 whenever
    # W>k.  Here H is top, so W=2+t+b.
    thm2367_types = []
    reduced_types = []
    for low_blockers in range(3):
        for top_ordinary in range(6):
            for top_low_blockers in range(low_blockers + 1):
                top_weight = 2 + top_ordinary + top_low_blockers
                if not (
                    top_weight <= low_blockers or top_weight >= 7
                ):
                    continue
                row = (
                    low_blockers,
                    top_ordinary,
                    top_low_blockers,
                    top_weight,
                )
                thm2367_types.append(row)
                if top_weight <= low_blockers or top_ordinary == 5:
                    reduced_types.append(row)

    expected_reduced = [
        (0, 5, 0, 7),
        (1, 5, 0, 7),
        (1, 5, 1, 8),
        (2, 0, 0, 2),
        (2, 5, 0, 7),
        (2, 5, 1, 8),
        (2, 5, 2, 9),
    ]
    require(
        reduced_types == expected_reduced,
        f"seven-type residual changed: {reduced_types}",
    )

    zero_depth_types = [
        row for row in reduced_types if row[2] == row[0]
    ]
    expected_zero_depth = [
        (0, 5, 0, 7),
        (1, 5, 1, 8),
        (2, 5, 2, 9),
    ]
    require(
        zero_depth_types == expected_zero_depth,
        "M=0 residual changed",
    )

    primitive_positive_depth = [
        row
        for row in reduced_types
        if not (row[1] == 5 and row[2] == row[0])
    ]
    expected_positive_depth = [
        (1, 5, 0, 7),
        (2, 0, 0, 2),
        (2, 5, 0, 7),
        (2, 5, 1, 8),
    ]
    require(
        primitive_positive_depth == expected_positive_depth,
        "primitive positive-depth residual changed",
    )

    guard, choices = cover_witness
    print("theorem=THM-2427")
    print("arithmetic=Fraction-and-finite-Z13")
    print(
        "ordinary-root-masks="
        f"{len(ordinary_masks)},sizes={ordinary_sizes},maximum=2"
    )
    print(
        "guard-root-masks="
        f"{len(guard_masks)},sizes={guard_sizes},maximum=4"
    )
    print("thirteen-root-capacity=13<=4+2t forces t=5")
    print(f"abstract-t5-guard={sorted(guard)}")
    print(
        "abstract-t5-ordinary="
        + ";".join(str(sorted(choice)) for choice in choices)
    )
    print("physical-t5-base=1/2,H=1,q=2,3,5,11,19")
    print(f"physical-t5-guard={sorted(physical_guard)}")
    print(
        "physical-t5-ordinary="
        + ";".join(str(sorted(mask)) for mask in physical_ordinary)
    )
    print("physical-t5-partition=all-13-roots")
    print(f"raw-THM2367-guard-top-types={len(thm2367_types)}")
    print(f"deep-c3-guard-top-residual-types={len(reduced_types)}")
    print(
        "residual-types="
        + ";".join(str(row) for row in reduced_types)
    )
    print(
        "M-zero-types="
        + ";".join(str(row) for row in zero_depth_types)
    )
    print(
        "primitive-M-positive-types="
        + ";".join(str(row) for row in primitive_positive_depth)
    )
    print("scope=no scalar-row decrement and no LRC(14) conclusion")
    print("status=PASS")


if __name__ == "__main__":
    main()
