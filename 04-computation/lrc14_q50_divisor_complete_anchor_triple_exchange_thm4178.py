#!/usr/bin/env python3
"""Exact q=50 failure-atom certificate for THM-4178 anchor-triple exchange.

This tracked primary audit deliberately retains the full
30-label failure mask on every fixed-pool wall cell, including cells on which
one or more of the original anchors 120, 126, 143 fail.  The q=50 safe mass of
each labelled atom is then subset-zeta summed over deletion masks.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, lcm
import sys


sys.stdout.reconfigure(newline="\n")


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
ORIGINAL_ANCHORS = (120, 126, 143)
EXPECTED_ANCHORS = (
    (40, 143, 252),
    (80, 143, 252),
    (120, 126, 143),
    (120, 143, 252),
    (126, 143, 240),
    (143, 240, 252),
)
COMMON = 18_241_159_416_480
Q = 50


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def safe_prefix_tick(tick):
    """Integral of the q-safe comb through tick/COMMON.

    The returned numerator has denominator 14*q*COMMON.
    """
    product = Q * tick
    whole, rem = divmod(product, COMMON)
    scaled = 14 * rem
    if scaled <= COMMON:
        partial = 0
    elif scaled >= 13 * COMMON:
        partial = 12 * COMMON
    else:
        partial = scaled - COMMON
    return 12 * whole * COMMON + partial


def mask_of(labels):
    indices = {value: index for index, value in enumerate(POOL)}
    return sum(1 << indices[value] for value in labels)


def labels_of(mask):
    return tuple(value for index, value in enumerate(POOL) if mask >> index & 1)


def is_divisor_complete(anchor):
    return all(any(value % modulus == 0 for value in anchor) for modulus in range(2, 15))


def anchor_gcd(anchor):
    value = 0
    for label in anchor:
        value = gcd(value, label)
    return value


def find_cover_exact(edges, budget):
    """Find a hitting set of size at most budget, or prove none exists.

    The recursion uses a greedy disjoint-edge packing as a valid lower bound.
    Vertices and edges remain full 30-label bit masks.
    """
    failed = set()

    def search(chosen, remaining):
        key = (chosen, remaining)
        if key in failed:
            return None

        uncovered = 0
        matching_used = 0
        matching_size = 0
        for edge in edges:
            if edge & chosen:
                continue
            if not uncovered:
                uncovered = edge
            if not edge & matching_used:
                matching_used |= edge
                matching_size += 1
                if matching_size > remaining:
                    failed.add(key)
                    return None

        if not uncovered:
            return chosen
        if remaining == 0:
            failed.add(key)
            return None

        branch = uncovered
        while branch:
            bit = branch & -branch
            answer = search(chosen | bit, remaining - 1)
            if answer is not None:
                return answer
            branch ^= bit

        failed.add(key)
        return None

    answer = search(0, budget)
    return answer, len(failed)


def minimum_cover_through_eight(edges):
    searches = []
    for budget in range(9):
        cover, states = find_cover_exact(edges, budget)
        searches.append((budget, states))
        if cover is not None:
            require(cover.bit_count() == budget, ("nonminimal returned cover", budget, labels_of(cover)))
            require(all(edge & cover for edge in edges), ("invalid cover", budget, labels_of(cover)))
            return budget, cover, tuple(searches)
    return None, None, tuple(searches)


def build_full_failure_atoms():
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)
    require(common == COMMON, "common wall lattice changed")

    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    ticks = tuple(int(wall * COMMON) for wall in walls)
    require(len(walls) == 7_134, "wall count changed")
    require(all(F(tick, COMMON) == wall for tick, wall in zip(ticks, walls)), "wall embedding")

    atom_mass = defaultdict(int)
    cell_mask_count = Counter()
    original_optional_hist = Counter()
    original_anchor_mask = mask_of(ORIGINAL_ANCHORS)
    original_anchor_failing_cells = 0

    previous = safe_prefix_tick(ticks[0])
    for index, (left, right) in enumerate(zip(walls, walls[1:])):
        current = safe_prefix_tick(ticks[index + 1])
        midpoint = (left + right) / 2
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe_at(midpoint, speed):
                failure |= 1 << bit
        contribution = current - previous
        require(contribution >= 0, ("negative q-safe contribution", index))
        atom_mass[failure] += contribution
        cell_mask_count[failure] += 1

        if failure & original_anchor_mask:
            original_anchor_failing_cells += 1
        else:
            original_optional_hist[failure.bit_count()] += 1
        previous = current

    require(previous == 12 * Q * COMMON, "q-safe prefix endpoint")
    require(sum(atom_mass.values()) == 12 * Q * COMMON, "full atom mass")
    require(original_anchor_failing_cells == 3_177, "original-anchor failing-cell count")
    require(
        tuple(sorted(original_optional_hist.items()))
        == ((0, 150), (1, 328), (2, 518), (3, 678), (4, 728),
            (5, 666), (6, 472), (7, 242), (8, 102), (9, 38),
            (10, 20), (11, 6), (12, 2), (13, 2), (14, 2), (15, 2)),
        "original-anchor cell histogram",
    )
    return walls, atom_mass, cell_mask_count


def deletion_numerator(mask, atom_mass):
    numerator = 0
    subset = mask
    while True:
        numerator += atom_mass.get(subset, 0)
        if subset == 0:
            return numerator
        subset = (subset - 1) & mask


def main():
    require(len(POOL) == 30 and len(set(POOL)) == 30, "pool cardinality")

    valid_anchors = tuple(
        anchor
        for anchor in combinations(POOL, 3)
        if is_divisor_complete(anchor) and anchor_gcd(anchor) == 1
    )
    require(valid_anchors == EXPECTED_ANCHORS, ("primitive divisor-complete anchors changed", valid_anchors))

    type_hostile = (120, 126, 286)
    require(is_divisor_complete(type_hostile), "type hostile must be divisor-complete")
    require(anchor_gcd(type_hostile) == 2, "type hostile gcd")

    walls, atom_mass, cell_mask_count = build_full_failure_atoms()
    original_anchor_mask = mask_of(ORIGINAL_ANCHORS)
    original_anchor_failing_mass = sum(
        mass for failure, mass in atom_mass.items() if failure & original_anchor_mask
    )
    require(original_anchor_failing_mass > 0, "old-anchor-failing atoms were discarded")

    print("Q", Q)
    print("POOL", POOL)
    print("VALID_PRIMITIVE_DIVISOR_COMPLETE_ANCHORS", valid_anchors)
    print("TYPE_HOSTILE", type_hostile, "DIVISOR_COMPLETE", True, "GCD", anchor_gcd(type_hostile))
    print("WALLS", len(walls), "CELLS", len(walls) - 1, "COMMON", COMMON)
    print(
        "FULL_FAILURE_ATOMS",
        len(cell_mask_count),
        "NONZERO_Q_MASS_ATOMS",
        sum(mass > 0 for mass in atom_mass.values()),
        "TOTAL_Q_SAFE_NUMERATOR",
        sum(atom_mass.values()),
        "DENOMINATOR",
        14 * Q * COMMON,
    )
    print(
        "ORIGINAL_ANCHOR_FAILING_CELLS",
        sum(count for failure, count in cell_mask_count.items() if failure & original_anchor_mask),
        "ORIGINAL_ANCHOR_FAILING_Q_SAFE_NUMERATOR",
        original_anchor_failing_mass,
    )

    rows = []
    for deletion_arity in (5, 6):
        positive_masks = []
        equality_masks = []
        for vertices in combinations(range(len(POOL)), deletion_arity):
            mask = sum(1 << vertex for vertex in vertices)
            numerator = deletion_numerator(mask, atom_mass)
            difference = 9 * numerator - 8 * Q * COMMON
            if difference >= 0:
                positive_masks.append(mask)
            if difference == 0:
                equality_masks.append(mask)

        positive_masks = tuple(positive_masks)
        equality_masks = tuple(equality_masks)
        print(
            "GLOBAL_ARITY",
            deletion_arity,
            "TOTAL",
            comb(len(POOL), deletion_arity),
            "REPAIRS",
            len(positive_masks),
            "EQUALITIES",
            len(equality_masks),
        )

        for anchor in valid_anchors:
            anchor_mask = mask_of(anchor)
            active = tuple(mask for mask in positive_masks if not mask & anchor_mask)
            equalities = sum(not mask & anchor_mask for mask in equality_masks)
            require(len(active) > 0, ("empty repair hypergraph", anchor, deletion_arity))
            require(all(not edge & anchor_mask for edge in active), ("anchor deletion leaked", anchor, deletion_arity))
            require(all(edge.bit_count() == deletion_arity for edge in active), "edge arity")

            tau, cover, searches = minimum_cover_through_eight(active)
            cover_labels = labels_of(cover) if cover is not None else None
            status = f"TAU={tau}" if tau is not None else "TAU>8"
            if tau is not None:
                require(not cover & anchor_mask, ("cover used an anchor", anchor, deletion_arity, cover_labels))

            if anchor == ORIGINAL_ANCHORS and deletion_arity == 5:
                require((len(active), equalities, tau) == (2_063, 0, 5), "original d=5 control")
            if anchor == ORIGINAL_ANCHORS and deletion_arity == 6:
                require((len(active), equalities, tau) == (46_261, 0, 8), "original d=6 control")

            row = (
                anchor,
                deletion_arity,
                len(active),
                equalities,
                status,
                cover_labels,
                searches,
            )
            rows.append(row)
            print(
                "ROW",
                "ANCHOR",
                anchor,
                "ARITY",
                deletion_arity,
                "TOTAL",
                comb(27, deletion_arity),
                "EDGES",
                len(active),
                "EQUALITIES",
                equalities,
                status,
                "COVER",
                cover_labels,
                "SEARCH_STATES",
                searches,
                flush=True,
            )

    print("ROW_COUNT", len(rows))
    print("NONORIGINAL_ROW_COUNT", sum(anchor != ORIGINAL_ANCHORS for anchor, *_ in rows))


if __name__ == "__main__":
    main()
