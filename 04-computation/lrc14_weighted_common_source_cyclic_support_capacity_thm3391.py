#!/usr/bin/env python3
"""Cheap exact audit for weighted common-source cyclic band covers.

This is the canonical exact companion for THM-3391.  It uses exact integers
and fractions, contains no optimization-dependent ``assert``, and exercises
the strict radius-1/14 boundary.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import ceil, gcd


RHO = Fraction(1, 14)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def lambda_sorted(n, source, numerator):
    """Exact maximum band load by circular sort/double/two-pointer.

    ``source`` is a tuple of source occurrences.  Repeated source residues and
    collisions under a nonunit numerator are deliberately retained.
    """
    images = sorted((numerator * residue) % n for residue in source)
    if not images:
        return 0, None
    size = len(images)
    doubled = images + [residue + n for residue in images]
    best = 0
    witness = None
    right = 0
    for left in range(size):
        if right < left:
            right = left
        while (
            right + 1 < left + size
            and 7 * (doubled[right + 1] - doubled[left]) < n
        ):
            right += 1
        load = right - left + 1
        candidate = (doubled[left] % n, doubled[right] % n, load)
        if load > best or (load == best and (witness is None or candidate < witness)):
            best = load
            witness = candidate
    return best, witness


def lambda_bitmap(n, source, numerator):
    """Exact maximum band load by pushforward bitmap and cyclic window."""
    pushforward = [0] * n
    for residue in source:
        pushforward[(numerator * residue) % n] += 1
    width = (n + 6) // 7
    extended = pushforward + pushforward[: width - 1]
    load = sum(extended[:width])
    best = load
    start = 0
    for candidate_start in range(1, n):
        load += extended[candidate_start + width - 1]
        load -= extended[candidate_start - 1]
        if load > best:
            best = load
            start = candidate_start
    return best, (start, width)


def full_fibre_capacity(n, numerator):
    common = gcd(numerator, n)
    orbit = n // common
    return common * ceil(Fraction(orbit, 7))


def window_edges(n, support, numerator):
    """Maximal pullback-band traces as masks on a distinct source support."""
    width = (n + 6) // 7
    edges = set()
    for start in range(n):
        positions = {(start + offset) % n for offset in range(width)}
        mask = 0
        for index, residue in enumerate(support):
            if (numerator * residue) % n in positions:
                mask |= 1 << index
        edges.add(mask)
    return tuple(sorted(edges))


def independent_cover_exists(n, support, numerators):
    target = (1 << len(support)) - 1
    edge_families = tuple(window_edges(n, support, a) for a in numerators)
    for edges in product(*edge_families):
        covered = 0
        for edge in edges:
            covered |= edge
        if covered == target:
            return True, edges
    return False, None


def q2_blocked(odd_speed, sheet, base):
    phase = Fraction(odd_speed, 2) * (base + sheet)
    residue = phase - phase.numerator // phase.denominator
    distance = min(residue, 1 - residue)
    return distance < RHO


def q2_cover_measure(first, second):
    """Exact measure of the correlated two-sheet full-cover locus on [0,1]."""
    boundaries = {Fraction(0), Fraction(1)}
    for speed in (first, second):
        radius = Fraction(1, 7 * speed)
        for sheet in (0, 1):
            for integer in range(-speed, 2 * speed + 1):
                centre = Fraction(2 * integer, speed) - sheet
                for boundary in (centre - radius, centre + radius):
                    if 0 < boundary < 1:
                        boundaries.add(boundary)
    ordered = sorted(boundaries)
    measure = Fraction(0)
    slabs = 0
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        covered = all(
            q2_blocked(first, sheet, midpoint)
            or q2_blocked(second, sheet, midpoint)
            for sheet in (0, 1)
        )
        if covered:
            measure += right - left
            slabs += 1
    return measure, slabs


def main():
    # Reported exact lambda universe: every subset through size three in
    # X_n, every multiplier (including zero and nonunits), 2 <= n <= 14.
    formula_packet = []
    for n in range(2, 15):
        for size in range(min(3, n) + 1):
            for support in combinations(range(n), size):
                for numerator in range(n):
                    sorted_load, sorted_witness = lambda_sorted(n, support, numerator)
                    bitmap_load, bitmap_witness = lambda_bitmap(n, support, numerator)
                    require(
                        sorted_load == bitmap_load,
                        ("lambda formula", n, support, numerator, sorted_load, bitmap_load),
                    )
                    formula_packet.append(
                        (n, support, numerator, sorted_load, sorted_witness, bitmap_witness)
                    )
    require(len(formula_packet) == 22_230, len(formula_packet))

    # Full-fibre formula and located image-order boundary through modulus 40.
    full_packet = []
    pair_packet = []
    for n in range(2, 41):
        full_source = tuple(range(n))
        for numerator in range(n):
            exact = lambda_sorted(n, full_source, numerator)[0]
            predicted = full_fibre_capacity(n, numerator)
            require(exact == predicted, ("full fibre", n, numerator, exact, predicted))
            full_packet.append((n, numerator, exact))
            for delta in range(1, n):
                source_order = n // gcd(n, delta)
                image_order = source_order // gcd(source_order, numerator)
                if 2 <= image_order <= 7:
                    load = lambda_sorted(n, (0, delta), numerator)[0]
                    require(load == 1, ("image-order pair", n, numerator, delta, image_order, load))
                    pair_packet.append((n, numerator, delta, image_order))
    require(len(full_packet) == 819, len(full_packet))
    require(len(pair_packet) == 4_574, len(pair_packet))

    # Weighted/multiplicity controls.  Raw occurrence count cannot be compared
    # with the unweighted ambient capacity after quotienting.
    repeated_single = (0,) * 100
    require(lambda_bitmap(14, repeated_single, 1)[0] == 100, "multiplicity collapse")
    weighted_antipodes = (0, 0, 0, 7, 7)
    require(lambda_bitmap(14, weighted_antipodes, 1)[0] == 3, "weighted antipodes")

    # Sharp one-blocker and unit-quantifier controls.
    require(lambda_sorted(14, (2, 3), 1)[0] == 2, "ambient equality cover")
    require(lambda_sorted(14, (0, 7), 1)[0] == 1, "order-two survivor")
    require(lambda_sorted(8, (0, 3), 1)[0] == 1, "fixed-unit survivor")
    require(lambda_sorted(8, (0, 3), 3)[0] == 2, "universal-unit hostile")
    require(lambda_sorted(4, (0, 2), 2)[0] == 2, "nonunit collapse")

    # Sum equality is decided by the Boolean cover, not the scalar sum.
    equality_support = (0, 1, 4, 5)
    equality_lambdas = tuple(lambda_sorted(8, equality_support, 1)[0] for _ in range(2))
    cover, cover_edges = independent_cover_exists(8, equality_support, (1, 1))
    require(sum(equality_lambdas) == len(equality_support), equality_lambdas)
    require(cover and cover_edges is not None, ("equality partition", cover_edges))

    # Losing common source identity is invalid: separate blockers can leave
    # different survivors while their union covers the whole source.
    source_masks = (0b01, 0b10)
    require(all(mask != 0b11 for mask in source_masks), source_masks)
    require(source_masks[0] | source_masks[1] == 0b11, source_masks)

    # Exact q=2 correlated graph and measure boundary for all odd speeds <=39.
    q2_packet = []
    for first in range(1, 40, 2):
        for second in range(1, 40, 2):
            measure, slabs = q2_cover_measure(first, second)
            predicted_edge = first + second > 7 * gcd(first, second)
            require((measure > 0) == predicted_edge, ("q2 graph", first, second, measure))
            require(measure <= Fraction(4, 7), ("q2 measure", first, second, measure))
            q2_packet.append((first, second, measure, slabs, predicted_edge))
    require(len(q2_packet) == 400, len(q2_packet))
    q2_nonedge = q2_cover_measure(1, 3)
    q2_edge = q2_cover_measure(1, 9)
    require(q2_nonedge[0] == 0, q2_nonedge)
    require(0 < q2_edge[0] <= Fraction(4, 7), q2_edge)
    require(Fraction(4, 7) == Fraction(52, 91) < Fraction(55, 91), "measure target")

    # Complete located order-r cosets: units block ceil(r/7), so m blockers
    # are defeated for 2 <= r <= 7 exactly when m < r.
    coset_packet = []
    for order in range(2, 8):
        capacity = ceil(Fraction(order, 7))
        require(capacity == 1, (order, capacity))
        for blockers in range(0, 9):
            exact_condition = blockers * capacity < order
            require(exact_condition == (blockers < order), (order, blockers))
            coset_packet.append((order, blockers, exact_condition))

    semantic_packet = (
        tuple(formula_packet),
        tuple(full_packet),
        tuple(pair_packet),
        (lambda_bitmap(14, repeated_single, 1), lambda_bitmap(14, weighted_antipodes, 1)),
        equality_lambdas,
        cover_edges,
        source_masks,
        tuple(q2_packet),
        tuple(coset_packet),
    )
    semantic_hash = digest(semantic_packet)

    print("WEIGHTED COMMON-SOURCE CYCLIC COVER EXACT PROBE")
    print("status=FINITE-EXACT;theorem=THM-3391;rho=1/14")
    print(
        "universe="
        f"lambda_formula:{len(formula_packet)};"
        f"full_fibre:{len(full_packet)};"
        f"image_order_pairs:{len(pair_packet)};"
        f"q2_ordered_pairs:{len(q2_packet)};"
        f"coset_thresholds:{len(coset_packet)}"
    )
    print(
        "formula=sort_double_strict_span_equals_pushforward_bitmap_window;"
        f"formula_sha256={digest(tuple(formula_packet))}"
    )
    print(
        "weighted_controls="
        f"single_residue_load:{lambda_bitmap(14, repeated_single, 1)[0]};"
        f"antipodal_load:{lambda_bitmap(14, weighted_antipodes, 1)[0]}"
    )
    print(
        "one_blocker_hostiles="
        "X14_adjacent_cover;X14_antipodal_survive;"
        "X8_fixed_unit_survive_but_unit3_covers;X4_nonunit_collapse"
    )
    print(
        "equality_partition="
        f"support:{equality_support};lambdas:{equality_lambdas};edges:{cover_edges};cover:true"
    )
    print("common_source_loss=two_individual_survivors_but_union_covers_source")
    print(
        "q2_boundary="
        f"nonedge_1_3:{q2_nonedge[0]};edge_1_9:{q2_edge[0]};"
        "all_cover_measures_le_4/7_lt_55/91"
    )
    print("unit_coset_rule=for_2_le_r_le_7_m_blockers_survive_iff_m_lt_r")
    print(f"semantic_sha256={semantic_hash}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
