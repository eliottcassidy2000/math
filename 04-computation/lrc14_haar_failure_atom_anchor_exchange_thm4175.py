#!/usr/bin/env python3
"""Exact audit for THM-4175 Haar failure atoms and q=50 anchor exchange.

The first half realizes the Boolean-zeta deletion transform on the THM-4156
pool while retaining every anchor-failure atom.  The second half reconstructs
both projected mixed-deletion and forced-anchor-deletion hypergraphs at the
hostile newcomer q=50, checks their exact failure/success boundaries, and
counts the resulting uniform primitive divisor-complete eleven-body
subfamilies.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, lcm


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(v for v in POOL if v not in ANCHORS)
Q = 50
COMMON = 18_241_159_416_480
THRESHOLD = F(4, 63)
TOMOGRAPHY_LABELS = (8, 15, 42, 63, 85, 88, 120, 126, 143, 252)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def safe_at(point, speed):
    residue = (speed * point) % 1
    return F(1, 14) <= residue <= F(13, 14)


def safe_prefix_numerator(tick, q):
    """Integral of G_q through tick/COMMON, denominator 14*q*COMMON."""
    whole, remainder = divmod(q * tick, COMMON)
    scaled = 14 * remainder
    if scaled <= COMMON:
        partial = 0
    elif scaled >= 13 * COMMON:
        partial = 12 * COMMON
    else:
        partial = scaled - COMMON
    return 12 * whole * COMMON + partial


def find_cover(edges, budget):
    """Exact bounded transversal search with a disjoint-edge lower bound."""
    failed = set()

    def search(chosen, remaining):
        key = (chosen, remaining)
        if key in failed:
            return None
        uncovered = 0
        matching_union = 0
        matching_size = 0
        for edge in edges:
            if edge & chosen:
                continue
            if not uncovered:
                uncovered = edge
            if not edge & matching_union:
                matching_union |= edge
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
            witness = search(chosen | bit, remaining - 1)
            if witness is not None:
                return witness
            branch ^= bit
        failed.add(key)
        return None

    return search(0, budget), len(failed)


def minimum_cover(edges, maximum):
    states = 0
    for budget in range(maximum + 1):
        witness, used = find_cover(edges, budget)
        states += used
        if witness is not None:
            require(witness.bit_count() <= budget, "cover budget mismatch")
            require(all(edge & witness for edge in edges), "invalid cover")
            return budget, witness, states
    return None, None, states


def mask_labels(mask, labels):
    return tuple(value for index, value in enumerate(labels) if mask >> index & 1)


def build_wall_atoms():
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)
    require(common == COMMON, "common lattice changed")

    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    ticks = tuple(int(point * COMMON) for point in walls)
    require(len(walls) == 7_134, "wall count changed")
    require(all(F(tick, COMMON) == point for tick, point in zip(ticks, walls)),
            "wall embedding failed")

    buckets = defaultdict(int)
    cell_rows = []
    previous = safe_prefix_numerator(ticks[0], Q)
    for index, (left, right) in enumerate(zip(walls, walls[1:])):
        current = safe_prefix_numerator(ticks[index + 1], Q)
        contribution = current - previous
        previous = current
        require(contribution >= 0, "negative q-safe cell contribution")
        midpoint = (left + right) / 2
        failure = sum(
            1 << label_index
            for label_index, speed in enumerate(POOL)
            if not safe_at(midpoint, speed)
        )
        buckets[failure] += contribution
        cell_rows.append((failure, contribution))
    require(len(cell_rows) == 7_133, "cell count changed")
    require(sum(buckets.values()) == 12 * Q * COMMON, "q-safe mass changed")
    return walls, ticks, dict(buckets), tuple(cell_rows)


def zeta_sum(buckets, deletion_mask):
    total = 0
    subset = deletion_mask
    while True:
        total += buckets.get(subset, 0)
        if subset == 0:
            return total
        subset = (subset - 1) & deletion_mask


def tomography_audit(buckets, cell_rows):
    pool_index = {speed: index for index, speed in enumerate(POOL)}
    selected_bits = tuple(1 << pool_index[speed] for speed in TOMOGRAPHY_LABELS)
    n = len(selected_bits)
    local_to_full = []
    for mask in range(1 << n):
        local_to_full.append(sum(selected_bits[i] for i in range(n) if mask >> i & 1))

    atom = [buckets.get(full, 0) for full in local_to_full]
    deletion = [zeta_sum(buckets, full) for full in local_to_full]

    # Direct retained-label cell integration is independent of subset summing.
    for local_mask, full_mask in enumerate(local_to_full):
        direct = sum(
            contribution
            for failure, contribution in cell_rows
            if not failure & ~full_mask
        )
        require(direct == deletion[local_mask], "direct/zeta deletion mismatch")

    # Boolean inversion recovers every exact atom supported in the selected set.
    recovered = []
    for mask in range(1 << n):
        value = 0
        subset = mask
        while True:
            value += (-1) ** (mask.bit_count() - subset.bit_count()) * deletion[subset]
            if subset == 0:
                break
            subset = (subset - 1) & mask
        recovered.append(value)
    require(recovered == atom, "Boolean atom inversion failed")
    require(all(value >= 0 for value in recovered), "negative recovered atom")

    # All mixed forward differences are sums of nonnegative failure atoms.
    difference_checks = 0
    for assignment in range(3 ** n):
        code = assignment
        base = 0
        delta = 0
        for index in range(n):
            digit = code % 3
            code //= 3
            if digit == 1:
                base |= 1 << index
            elif digit == 2:
                delta |= 1 << index
        value = 0
        subset = delta
        while True:
            sign = (-1) ** (delta.bit_count() - subset.bit_count())
            value += sign * deletion[base | subset]
            if subset == 0:
                break
            subset = (subset - 1) & delta
        atom_sum = sum(
            atom[delta | subset]
            for subset in range(1 << n)
            if not subset & ~base
        )
        require(value == atom_sum >= 0, "complete-monotonicity identity failed")
        difference_checks += 1

    # Symmetric deletion layers are the dual binomial transform of atom sizes.
    atom_layers = [0] * (n + 1)
    deletion_layers = [0] * (n + 1)
    for mask in range(1 << n):
        atom_layers[mask.bit_count()] += atom[mask]
        deletion_layers[mask.bit_count()] += deletion[mask]
    for r in range(n + 1):
        expected = sum(comb(n - s, r - s) * atom_layers[s] for s in range(r + 1))
        require(deletion_layers[r] == expected, "symmetric zeta layer failed")
    for s in range(n + 1):
        recovered_layer = sum(
            (-1) ** (s - r) * comb(n - r, s - r) * deletion_layers[r]
            for r in range(s + 1)
        )
        require(recovered_layer == atom_layers[s], "symmetric layer inversion failed")

    return {
        "labels": TOMOGRAPHY_LABELS,
        "subsets": 1 << n,
        "mixed_differences": difference_checks,
        "nonzero_atoms": sum(value > 0 for value in atom),
        "atom_layers": tuple(atom_layers),
        "deletion_layers": tuple(deletion_layers),
    }


def repair_summary(edges, raw_count, equalities, minimum_positive,
                   maximum_negative):
    tau, cover, states = minimum_cover(tuple(edges), 8)
    return {
        "raw_repairs": raw_count,
        "minimal_repairs": len(edges),
        "tau": tau,
        "cover": mask_labels(cover, OPTIONAL) if cover is not None else None,
        "states": states,
        "equalities": equalities,
        "minimum_positive": minimum_positive,
        "closest_negative": maximum_negative,
    }


def forced_anchor_hypergraph(buckets, removed, depth):
    pool_index = {speed: index for index, speed in enumerate(POOL)}
    ground_index = {speed: index for index, speed in enumerate(OPTIONAL)}
    anchor_bit = 1 << pool_index[removed]
    edges = []
    equalities = 0
    minimum_positive = None
    maximum_negative = None
    for deletion_tuple in combinations(OPTIONAL, depth):
        full_mask = anchor_bit | sum(1 << pool_index[speed] for speed in deletion_tuple)
        numerator = zeta_sum(buckets, full_mask)
        difference = 9 * numerator - 8 * Q * COMMON
        if difference >= 0:
            edges.append(sum(1 << ground_index[speed] for speed in deletion_tuple))
            if difference == 0:
                equalities += 1
            elif minimum_positive is None or difference < minimum_positive:
                minimum_positive = difference
        elif maximum_negative is None or difference > maximum_negative:
            maximum_negative = difference
    # Every edge has the same optional-deletion size, so no proper edge can
    # contain another and the combination iterator creates no duplicates.
    reduced = tuple(edges)
    return repair_summary(
        reduced, len(edges), equalities, minimum_positive, maximum_negative
    )


def mixed_anchor_hypergraph(buckets, removed, total_depth):
    """Project total-depth repairs from {removed} union OPTIONAL onto OPTIONAL."""
    pool_index = {speed: index for index, speed in enumerate(POOL)}
    ground_index = {speed: index for index, speed in enumerate(OPTIONAL)}
    deletion_ground = (removed,) + OPTIONAL
    edges = []
    equalities = 0
    minimum_positive = None
    maximum_negative = None
    for deletion_tuple in combinations(deletion_ground, total_depth):
        full_mask = sum(1 << pool_index[speed] for speed in deletion_tuple)
        numerator = zeta_sum(buckets, full_mask)
        difference = 9 * numerator - 8 * Q * COMMON
        if difference >= 0:
            edges.append(sum(
                1 << ground_index[speed]
                for speed in deletion_tuple
                if speed != removed
            ))
            if difference == 0:
                equalities += 1
            elif minimum_positive is None or difference < minimum_positive:
                minimum_positive = difference
        elif maximum_negative is None or difference > maximum_negative:
            maximum_negative = difference

    # Projecting away the omitted anchor creates edges of sizes t-1 and t.
    # The t-1 edges are minimal.  A t-edge is minimal exactly when none of
    # its t immediate submasks is an active t-1 edge.
    smaller = {edge for edge in edges if edge.bit_count() == total_depth - 1}
    larger = {edge for edge in edges if edge.bit_count() == total_depth}
    reduced_larger = tuple(sorted(
        edge for edge in larger
        if all((edge ^ bit) not in smaller for bit in iter_bits(edge))
    ))
    reduced = tuple(sorted(smaller)) + reduced_larger
    require(len(reduced) == len(set(reduced)), "duplicate projected repair")
    require(all(
        all((edge ^ bit) not in smaller for bit in iter_bits(edge))
        for edge in reduced_larger
    ), "nonminimal projected repair")
    return repair_summary(
        reduced, len(edges), equalities, minimum_positive, maximum_negative
    )


def iter_bits(mask):
    while mask:
        bit = mask & -mask
        yield bit
        mask ^= bit


def divisor_complete(body):
    return all(any(value % modulus == 0 for value in body) for modulus in range(2, 15))


def primitive(body):
    value = 0
    for speed in body:
        value = gcd(value, speed)
    return value == 1


def body_census():
    count_by_missing = Counter()
    divisor_by_missing = Counter()
    hostile_nonprimitive = 0
    for optional_body in combinations(OPTIONAL, 8):
        for removed in ANCHORS:
            # Omit q here: these are the subfamilies divisor-complete and
            # primitive uniformly for every newcomer q outside the pool.
            body = tuple(anchor for anchor in ANCHORS if anchor != removed) + optional_body
            if divisor_complete(body):
                divisor_by_missing[removed] += 1
                if primitive(body):
                    count_by_missing[removed] += 1
                elif removed == 143:
                    hostile_nonprimitive += 1

    missing_120_formula = (
        comb(27, 8) - comb(18, 8) - comb(17, 8) - comb(20, 8)
        + comb(11, 8) + comb(14, 8) + comb(12, 8) - comb(8, 8)
    )
    require(missing_120_formula == 2_029_699, "missing-120 formula")
    require(count_by_missing[120] == missing_120_formula, "missing-120 count")
    require(
        count_by_missing[126]
        == comb(27, 8) - comb(23, 8) - comb(25, 8) + comb(22, 8),
        "missing-126 inclusion-exclusion count",
    )
    require(count_by_missing[143] == comb(26, 7) - comb(20, 7),
            "missing-143 primitive count")
    require(hostile_nonprimitive == comb(20, 7), "even hostile count")
    require(divisor_by_missing[143] == comb(26, 7), "missing-143 divisor count")
    require(sum(count_by_missing.values()) == 3_577_935, "uniform total core count")
    require(sum(divisor_by_missing.values()) == 3_655_455,
            "odd-newcomer divisor/primitive total")
    return (
        tuple(sorted(divisor_by_missing.items())),
        tuple(sorted(count_by_missing.items())),
        hostile_nonprimitive,
    )


def main():
    require(len(POOL) == 30 and len(OPTIONAL) == 27, "pool sizes")
    require(Q not in POOL, "q must be a newcomer")
    walls, _ticks, buckets, cell_rows = build_wall_atoms()
    tomography = tomography_audit(buckets, cell_rows)

    mixed_rows = []
    mixed_expected = {
        120: ((2_467, 535, 5), (59_320, 13_158)),
        126: ((2_126, 1_417, 5), (51_810, 11_051)),
        143: ((2_476, 654, 7), (60_509, 14_321)),
    }
    for removed in ANCHORS:
        hostile = mixed_anchor_hypergraph(buckets, removed, 5)
        positive = mixed_anchor_hypergraph(buckets, removed, 6)
        (raw5, minimal5, tau5), (raw6, minimal6) = mixed_expected[removed]
        require(
            (hostile["raw_repairs"], hostile["minimal_repairs"], hostile["tau"])
            == (raw5, minimal5, tau5),
            "mixed hostile ledger changed",
        )
        require(
            (positive["raw_repairs"], positive["minimal_repairs"])
            == (raw6, minimal6),
            "mixed positive ledger changed",
        )
        require(positive["tau"] is None, "mixed depth-six repair regressed")
        require(hostile["equalities"] == positive["equalities"] == 0,
                "mixed threshold equality appeared")
        mixed_rows.append((removed, hostile, positive))

    boundary_depths = {120: (4, 5), 126: (5, 6), 143: (4, 5)}
    forced_rows = []
    for removed, (hostile_depth, positive_depth) in boundary_depths.items():
        hostile = forced_anchor_hypergraph(buckets, removed, hostile_depth)
        positive = forced_anchor_hypergraph(buckets, removed, positive_depth)
        require(hostile["tau"] is not None and hostile["tau"] <= 8,
                "hostile anchor-exchange level disappeared")
        require(positive["tau"] is None, "positive anchor-exchange level regressed")
        require(hostile["equalities"] == positive["equalities"] == 0,
                "threshold equality appeared")
        forced_rows.append((removed, hostile_depth, hostile, positive_depth, positive))

    divisor_counts, uniform_primitive_counts, hostile_nonprimitive = body_census()

    print("LRC14_HAAR_FAILURE_ATOM_ANCHOR_EXCHANGE")
    print("q", Q, "pool", len(POOL), "optional", len(OPTIONAL))
    print("walls", len(walls), "cells", len(cell_rows), "failure_atoms", len(buckets))
    print("threshold", f"{THRESHOLD.numerator}/{THRESHOLD.denominator}")
    print("tomography", tomography)
    print("mixed_total_depth_rows")
    for removed, hostile, positive in mixed_rows:
        print("removed", removed, "t5", hostile, "t6", positive)
    print("forced_anchor_optional_depth_rows")
    for removed, hostile_depth, hostile, positive_depth, positive in forced_rows:
        print("removed", removed, "d" + str(hostile_depth), hostile,
              "d" + str(positive_depth), positive)
    print("uniform_divisor_complete_counts", divisor_counts)
    print("uniform_primitive_divisor_complete_counts", uniform_primitive_counts)
    print("missing_143_even_nonprimitive_hostile", hostile_nonprimitive)
    print("uniform_even_q_primitive_total",
          sum(value for _, value in uniform_primitive_counts))
    print("uniform_odd_q_primitive_total", sum(value for _, value in divisor_counts))
    print("PASS")


if __name__ == "__main__":
    main()
