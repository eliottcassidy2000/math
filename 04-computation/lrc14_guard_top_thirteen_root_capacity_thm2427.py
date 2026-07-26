#!/usr/bin/env python3
"""Exact companion for THM-2427.

Checks the finite thirteen-root capacity used to prune the guard-top,
deep-c3 branch, the complete (k,t,b,W) enumeration, its primitive
M>0 reduction, positive abstract and physical common-phase t=5
root-cover controls, and a strict physical 13-by-7 stalk partition
with all three blockers absent.

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


def physical_root_mask(modulus, speed, width, base):
    """Danger roots of one actual integer speed at a common base."""

    return frozenset(
        root
        for root in range(modulus)
        if circle_norm(Fraction(speed) * (base + root) / modulus)
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
        P,
        physical_guard_speed,
        2,
        physical_base,
    )
    physical_ordinary = tuple(
        physical_root_mask(P, speed, 1, physical_base)
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

    # A stronger positive control realizes the whole compositional
    # thirteen-root by seven-bin stalk, not merely its two marginals.
    # The quotient blockers have strict thirteen-adic depths 0,1,2
    # after division by thirteen and are safe everywhere on the stalk.
    joint_base = Fraction(11, 581)
    joint_guard_unit = 1
    joint_ordinary_units = (547, 1821, 3095, 4369, 5643)
    joint_units = (joint_guard_unit,) + joint_ordinary_units
    joint_widths = (2, 1, 1, 1, 1, 1)
    joint_blockers = (13 * 7, 13**2 * 7, 13**3 * 7)
    require(
        all(gcd(speed, 91) == 1 for speed in joint_units),
        "joint stalk lost a 91-unit top speed",
    )

    def joint_value(unit, depth, root13, root7):
        """Top phase on the depth-M compositional stalk."""

        top_speed = 7**depth * unit
        base = joint_base / 7**depth
        orbit_order = 7 ** (depth + 1)
        point = (
            (base + root13) / 13
            + Fraction(root7, orbit_order)
        )
        return Fraction(top_speed) * point

    # At depth zero, CRT identifies the product stalk with Z/91Z.
    # In that gauge the six masks are six consecutive blocks.
    joint_masks = []
    for unit, width in zip(joint_units, joint_widths):
        mask = frozenset(
            (7 * root13 + 13 * root7) % 91
            for root13 in range(13)
            for root7 in range(7)
            if circle_norm(joint_value(unit, 0, root13, root7))
            < Fraction(width, 14)
        )
        joint_masks.append(mask)
    expected_joint_masks = (
        frozenset((*range(13), *range(78, 91))),
        frozenset(range(13, 26)),
        frozenset(range(26, 39)),
        frozenset(range(39, 52)),
        frozenset(range(52, 65)),
        frozenset(range(65, 78)),
    )
    require(
        tuple(joint_masks) == expected_joint_masks,
        "joint Z/91 stalk masks changed",
    )
    require(
        all(
            joint_masks[left].isdisjoint(joint_masks[right])
            for left in range(len(joint_masks))
            for right in range(left)
        ),
        "joint Z/91 masks stopped being disjoint",
    )
    require(
        set().union(*joint_masks) == set(range(91)),
        "joint Z/91 masks stopped partitioning the stalk",
    )
    direct_joint_masks = tuple(
        physical_root_mask(
            91,
            unit,
            width,
            7 * joint_base,
        )
        for unit, width in zip(joint_units, joint_widths)
    )
    require(
        direct_joint_masks == tuple(joint_masks),
        "direct Z/91 path disagrees with product-stalk path",
    )

    # The thirteen-root section is itself an exact partition.
    joint_root13_masks = tuple(
        frozenset(
            root13
            for root13 in range(13)
            if circle_norm(joint_value(unit, 0, root13, 0))
            < Fraction(width, 14)
        )
        for unit, width in zip(joint_units, joint_widths)
    )
    expected_joint_root13 = (
        frozenset((0, 1, 12)),
        frozenset((2, 3)),
        frozenset((4, 5)),
        frozenset((6, 7)),
        frozenset((8, 9)),
        frozenset((10, 11)),
    )
    require(
        joint_root13_masks == expected_joint_root13,
        "joint thirteen-root section changed",
    )
    require(
        joint_root13_masks
        == tuple(
            physical_root_mask(13, unit, width, joint_base)
            for unit, width in zip(joint_units, joint_widths)
        ),
        "direct thirteen-root path disagrees with product section",
    )
    for blocker in joint_blockers:
        for root13 in range(13):
            for root7 in range(7):
                point = (
                    (joint_base + root13) / 13
                    + Fraction(root7, 7)
                )
                require(
                    circle_norm(Fraction(blocker) * point)
                    >= Fraction(1, 14),
                    "direct product path found a dangerous blocker",
                )

    # Exact self-similarity: after multiplying all top speeds by 7^M,
    # dividing the quotient base by 7^M, and enlarging the orbit to
    # 7^(M+1), the depth-M word is the same 13-by-7 word repeated
    # 7^M times.  The strict blockers remain safe at every point.
    for depth in range(9):
        orbit_order = 7 ** (depth + 1)
        for root13 in range(13):
            reduced_root13 = (7**depth * root13) % 13
            for root7 in range(7):
                occupancies = [
                    circle_norm(
                        joint_value(unit, depth, root13, root7)
                    )
                    < Fraction(width, 14)
                    for unit, width in zip(joint_units, joint_widths)
                ]
                expected_occupancies = [
                    circle_norm(
                        joint_value(unit, 0, reduced_root13, root7)
                    )
                    < Fraction(width, 14)
                    for unit, width in zip(joint_units, joint_widths)
                ]
                require(
                    occupancies == expected_occupancies,
                    f"depth-{depth} stalk reduction failed",
                )
                require(
                    sum(occupancies) == 1,
                    f"depth-{depth} stalk stopped being one-fold",
                )
            require(
                orbit_order // 7 == 7**depth,
                "septimal repetition count changed",
            )

        base = joint_base / 7**depth
        for blocker_index in range(1, 4):
            blocker = 13**blocker_index * 7 ** (depth + 1)
            quotient_blocker = blocker // 13
            blocker_phase = circle_norm(
                Fraction(quotient_blocker) * base
            )
            require(
                blocker_phase >= Fraction(1, 14),
                f"depth-{depth} blocker {blocker_index} became dangerous",
            )

    # The depth-zero packet is a positive chamber, not an endpoint
    # coincidence.  This is the exact minimum base displacement to a
    # top-word or blocker boundary.
    chamber_radii = []
    for unit, width in zip(joint_units, joint_widths):
        threshold = Fraction(width, 14)
        for root13 in range(13):
            for root7 in range(7):
                slack = abs(
                    circle_norm(
                        joint_value(unit, 0, root13, root7)
                    )
                    - threshold
                )
                chamber_radii.append(slack * 13 / unit)
    for blocker in joint_blockers:
        quotient_blocker = blocker // 13
        slack = abs(
            circle_norm(Fraction(quotient_blocker) * joint_base)
            - Fraction(1, 14)
        )
        chamber_radii.append(slack / quotient_blocker)
    joint_chamber_radius = min(chamber_radii)
    require(
        joint_chamber_radius == Fraction(1, 635614),
        "joint stalk positive-chamber radius changed",
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
    print("arithmetic=Fraction-and-finite-Z7xZ13")
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
    print("joint-stalk-quotient-base=11/581")
    print("joint-stalk-top=H:1,q:547,1821,3095,4369,5643")
    print("joint-stalk-blockers=91,1183,15379")
    print(
        "joint-Z91-guard="
        + ",".join(str(root) for root in sorted(joint_masks[0]))
    )
    print(
        "joint-Z91-ordinary="
        + ";".join(
            ",".join(str(root) for root in sorted(mask))
            for mask in joint_masks[1:]
        )
    )
    print(
        "joint-Z13-section="
        + ";".join(
            ",".join(str(root) for root in sorted(mask))
            for mask in joint_root13_masks
        )
    )
    print("joint-product-stalk=exact-one-fold-on-13x7")
    print("joint-depth-lift=verified-M0-through-M8")
    print(f"joint-positive-chamber-radius={joint_chamber_radius}")
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
