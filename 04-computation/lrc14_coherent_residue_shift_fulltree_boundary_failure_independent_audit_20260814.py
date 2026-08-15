#!/usr/bin/env python3
"""Canonical self-contained audit of the first full-tree/literal boundaries.

This verifier does not import either census engine.  It constructs the six
periodic interval families literally, enumerates all 81 cross-K3,3 and all
1,296 K6 spanning trees, and computes union mass by an endpoint sweep.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations


BODY = (1, 2, 3, 4, 6, 12)
L = 168
CELL = 90
CASES = (
    (390, (393, 394, 402), (391, 392, 396),
     Q(-717_819, 1_739_713_360), Q(-717_819, 1_739_713_360),
     Q(699_359_657, 869_856_680)),
    (390, (393, 396, 402), (391, 392, 394),
     Q(-2_185_509_467, 63_483_722_064), Q(-100_664_755, 7_935_465_258),
     Q(1_435_246_319, 1_715_776_272)),
    (392, (393, 394, 398), (395, 396, 404),
     Q(-5_691_856, 138_190_745), Q(-5_691_856, 138_190_745),
     Q(1_553_722, 1_794_685)),
    (392, (394, 398, 404), (393, 395, 396),
     Q(-310_940_537, 4_550_843_675), Q(-133_244_522, 4_550_843_675),
     Q(113_423, 137_825)),
    (392, (395, 396, 404), (393, 394, 398),
     Q(-80_938_162, 2_313_097_605), Q(-80_938_162, 2_313_097_605),
     Q(58_693_679, 66_088_503)),
    (392, (395, 398, 404), (393, 394, 396),
     Q(-464_507_574, 11_402_675_725), Q(-39_048_713, 2_280_535_145),
     Q(4_497, 5_513)),
)
EXPECTED_FAILURE_SHIFTS = (390, 390, 392, 392, 392, 392)
EXPECTED_UNION_FAILURE_SHIFTS = (392, 392)
EXPECTED_SEMANTIC = "6eb12a82411220aaffdc415cecdda50c6ac9124b4048d506ec8b045720f3be7d"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def clip(left, right):
    left = max(Q(), left)
    right = min(Q(1), right)
    return (left, right) if left < right else None


def literal_arcs(multiplier):
    """Pull back ||multiplier*(CELL+u)/L||<1/14 literally."""
    require(multiplier > 0, multiplier)
    lo = multiplier * CELL // L - 1
    hi = multiplier * (CELL + 1) // L + 1
    rows = []
    for integer in range(lo, hi + 1):
        left = Q(L * (14 * integer - 1), 14 * multiplier) - CELL
        right = Q(L * (14 * integer + 1), 14 * multiplier) - CELL
        piece = clip(left, right)
        if piece is not None:
            rows.append(piece)
    return tuple(sorted(set(rows)))


def sweep_mass(families):
    rows = tuple(row for family in families for row in family)
    endpoints = sorted({endpoint for row in rows for endpoint in row})
    total = Q()
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        if any(a < midpoint < b for a, b in rows):
            total += right - left
    return total


def pair_mass(first, second):
    total = Q()
    for left1, right1 in first:
        for left2, right2 in second:
            total += max(Q(), min(right1, right2) - max(left1, left2))
    return total


def brute_maximum_tree(vertices, edges):
    maximum = None
    count = 0
    for chosen in combinations(edges, len(vertices) - 1):
        parent = {vertex: vertex for vertex in vertices}

        def root(vertex):
            while parent[vertex] != vertex:
                vertex = parent[vertex]
            return vertex

        for _weight, first, second in chosen:
            x, y = root(first), root(second)
            if x == y:
                break
            parent[x] = y
        else:
            count += 1
            credit = sum((edge[0] for edge in chosen), Q())
            maximum = credit if maximum is None else max(maximum, credit)
    require(maximum is not None, vertices)
    return maximum, count


def main():
    require(
        all(
            14 * ((label * CELL) % L) >= L
            and 14 * (((label * CELL) % L) + label) <= 13 * L
            for label in BODY
        ),
        "fixed cell is not body-safe",
    )
    digest = sha256()
    rows = []
    for shift, left, right, expected_cross, expected_full, expected_union in CASES:
        require(tuple(sorted(left + right)) == tuple(label + shift for label in BODY),
                (shift, left, right))
        vertices = tuple((3, residue) for residue in left) + tuple(
            (5, residue) for residue in right
        )
        events = {
            vertex: literal_arcs(vertex[0] * L - vertex[1])
            for vertex in vertices
        }
        singles = {vertex: sweep_mass((event,)) for vertex, event in events.items()}
        debt = sum(singles.values(), Q()) - Q(6, 7)
        cross_edges = tuple(
            (pair_mass(events[3, a], events[5, b]), (3, a), (5, b))
            for a in left for b in right
        )
        full_edges = tuple(
            (pair_mass(events[first], events[second]), first, second)
            for first, second in combinations(vertices, 2)
        )
        cross_credit, cross_count = brute_maximum_tree(vertices, cross_edges)
        full_credit, full_count = brute_maximum_tree(vertices, full_edges)
        union = sweep_mass(tuple(events.values()))
        cross_margin = cross_credit - debt
        full_margin = full_credit - debt
        require(cross_count == 81, (shift, cross_count))
        require(full_count == 1_296, (shift, full_count))
        require(cross_margin == expected_cross, (shift, cross_margin))
        require(full_margin == expected_full, (shift, full_margin))
        require(union == expected_union, (shift, union))
        row = (
            shift, left, right, cross_count, full_count, debt,
            cross_margin, full_margin, union, union - Q(6, 7),
        )
        rows.append(row)
        digest.update((repr(row) + "\n").encode())
    failure_shifts = tuple(row[0] for row in rows)
    union_failure_shifts = tuple(row[0] for row in rows if row[9] >= 0)
    require(failure_shifts == EXPECTED_FAILURE_SHIFTS, failure_shifts)
    require(union_failure_shifts == EXPECTED_UNION_FAILURE_SHIFTS,
            union_failure_shifts)
    semantic = digest.hexdigest()
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, semantic)
    print("LRC14 COHERENT-SHIFT FULL-TREE/LITERAL BOUNDARY FAILURE AUDIT")
    print("body", BODY, "L", L, "cell", CELL, "cases", len(rows))
    for row in rows:
        print("case", row)
    print("full_failure_shifts", failure_shifts)
    print("union_failure_shifts", union_failure_shifts)
    print("semantic_sha256", semantic)
    print("status=SELF-CONTAINED VERIFIED-EXACT; literal arcs, exhaustive trees, endpoint sweep")


if __name__ == "__main__":
    main()
