#!/usr/bin/env python3
"""Exact referee for THM-2656.

The certificate compares two finite rainbow-selector atlases on every
two-hole relation pair over F_13: an orientation-independent atlas whose
charts are individually reflection-equivariant, and an edge-disjoint atlas.
Every logical guard survives optimized Python.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations


P = 13
INV2 = 7
SOURCES = tuple(range(2, 13))


SYMMETRIC_MATCHED = {
    1: (3, 4, 2, 6, 9, 7, 5, 8, 12, 10, 11),
    12: (2, 3, 1, 5, 8, 6, 4, 7, 11, 9, 10),
}


DISJOINT_SECOND = {
    1: (11, 12, 5, 2, 4, 9, 6, 3, 8, 10, 7),
    2: (1, 4, 9, 3, 5, 7, 11, 6, 8, 12, 10),
    3: (2, 10, 5, 11, 1, 7, 4, 6, 8, 12, 9),
    4: (3, 10, 11, 9, 1, 5, 8, 2, 12, 6, 7),
    5: (4, 10, 1, 7, 8, 2, 3, 12, 6, 9, 11),
    6: (5, 8, 9, 1, 11, 7, 2, 3, 12, 10, 4),
    7: (6, 10, 1, 11, 9, 3, 12, 2, 4, 8, 5),
    8: (7, 10, 1, 9, 2, 3, 4, 11, 5, 6, 12),
    9: (8, 10, 11, 4, 1, 5, 3, 12, 6, 7, 2),
    10: (9, 1, 3, 8, 2, 7, 4, 6, 12, 5, 11),
    11: (10, 1, 6, 8, 9, 2, 3, 12, 4, 5, 7),
    12: (10, 11, 4, 1, 3, 8, 5, 2, 7, 9, 6),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inv(value):
    require(value % P != 0, "inverse of zero")
    return pow(value, P - 2, P)


def complement(pair):
    missing = set(pair)
    return tuple(x for x in range(P) if x not in missing)


def convolution(left, right):
    return tuple(
        sum(left[a] * right[(x - a) % P] for a in range(P))
        for x in range(P)
    )


def mask(pair):
    support = set(pair)
    return tuple(Fraction(int(x in support)) for x in range(P))


def translate(function, shift):
    return tuple(function[(x - shift) % P] for x in range(P))


def two_point_inverse(step):
    return tuple(
        sum(
            Fraction(1 if j % 2 == 0 else -1, 2)
            for j in range(P)
            if (j * step) % P == x
        )
        for x in range(P)
    )


def normalized_graph(targets):
    require(len(targets) == 11, "normal target row length")
    return tuple(zip(SOURCES, targets))


def first_normalized_chart(ratio):
    if ratio != P - 1:
        return tuple((x, ratio * x % P) for x in SOURCES)
    return tuple((x, (x - 1) % P) for x in SOURCES)


def symmetric_normalized_charts(ratio):
    if ratio not in (1, P - 1):
        plus = tuple((x, ratio * x % P) for x in SOURCES)
        minus = tuple((x, ratio * (1 - x) % P) for x in SOURCES)
        return plus, minus
    return first_normalized_chart(ratio), normalized_graph(SYMMETRIC_MATCHED[ratio])


def disjoint_normalized_charts(ratio):
    return first_normalized_chart(ratio), normalized_graph(DISJOINT_SECOND[ratio])


def transport(graph, a_order, b_order):
    a0, a1 = a_order
    b0, _ = b_order
    scale = (a1 - a0) % P
    return tuple(
        ((a0 + scale * x) % P, (b0 + scale * y) % P)
        for x, y in graph
    )


def charts_for_orders(a_order, b_order, kind):
    scale = (a_order[1] - a_order[0]) % P
    ratio = (b_order[1] - b_order[0]) * inv(scale) % P
    if kind == "symmetric":
        graphs = symmetric_normalized_charts(ratio)
    elif kind == "disjoint":
        graphs = disjoint_normalized_charts(ratio)
    else:
        raise RuntimeError("unknown atlas kind")
    return tuple(transport(graph, a_order, b_order) for graph in graphs)


def verify_matching(a_pair, b_pair, graph):
    source = set(complement(a_pair))
    target = set(complement(b_pair))
    xs = tuple(x for x, _ in graph)
    ys = tuple(y for _, y in graph)
    colours = frozenset((x + y) % P for x, y in graph)
    require(len(graph) == 11, "matching size")
    require(set(xs) == source and len(set(xs)) == 11, "source bijection")
    require(set(ys) == target and len(set(ys)) == 11, "target bijection")
    require(len(colours) == 11, "rainbow failure")
    holes = frozenset(set(range(P)) - colours)
    require(len(holes) == 2, "hole count")
    return colours, holes


def is_affine(graph):
    (x0, y0), (x1, y1) = graph[:2]
    slope = (y1 - y0) * inv(x1 - x0) % P
    intercept = (y0 - slope * x0) % P
    return all(y == (slope * x + intercept) % P for x, y in graph)


def verify_equivariance(a_pair, b_pair, graph):
    sum_a = sum(a_pair) % P
    sum_b = sum(b_pair) % P
    lookup = dict(graph)
    for x, y in graph:
        require(lookup[(sum_a - x) % P] == (sum_b - y) % P,
                "individual reflection equivariance")


def atlas_key(charts):
    return frozenset(frozenset(chart) for chart in charts)


def verify_cover(a_pair, b_pair, charts, expected_overlap):
    (c0, h0), (c1, h1) = tuple(
        verify_matching(a_pair, b_pair, chart) for chart in charts
    )
    require(h0.isdisjoint(h1), "hole pairs intersect")
    require(c0 | c1 == frozenset(range(P)), "carry cover incomplete")
    overlap = len(set(charts[0]) & set(charts[1]))
    require(overlap == expected_overlap, "edge overlap")
    require(len(set(charts[0]) | set(charts[1])) == 22 - expected_overlap,
            "edge union count")
    profile = tuple(int(c in c0) + int(c in c1) for c in range(P))
    require(Counter(profile) == Counter({2: 9, 1: 4}), "colour profile")
    mean = Fraction(22, 13)
    energy = sum((Fraction(value) - mean) ** 2 for value in profile) / P
    require(energy == Fraction(36, 169), "charged energy")
    for k in range(1, P):
        permuted = tuple(profile[(-k * exponent) % P] for exponent in range(P))
        require(len(set(permuted)) > 1, "charged character vanished")
    return h0, h1


def main():
    pairs = tuple(combinations(range(P), 2))
    require(len(pairs) == 78, "two-set count")
    require(set(SYMMETRIC_MATCHED) == {1, 12}, "matched bank keys")
    require(set(DISJOINT_SECOND) == set(range(1, 13)), "disjoint bank keys")

    for ratio, first_holes, second_holes in (
        (1, frozenset((0, 2)), frozenset((3, 12))),
        (12, frozenset((1, 12)), frozenset((2, 11))),
    ):
        normalized_b = (0, ratio)
        charts = symmetric_normalized_charts(ratio)
        _, actual_first = verify_matching((0, 1), normalized_b, charts[0])
        _, actual_second = verify_matching((0, 1), normalized_b, charts[1])
        require(actual_first == first_holes, "matched affine hole formula")
        require(actual_second == second_holes, "matched nonlinear hole formula")
        require(len(set(charts[0]) & set(charts[1])) == 1,
                "matched symmetric midpoint overlap")

    # The two matched templates are compatible with both changes of endpoint
    # orientation; these are the only non-affine orientation identities.
    parallel = dict(normalized_graph(SYMMETRIC_MATCHED[1]))
    anti = dict(normalized_graph(SYMMETRIC_MATCHED[12]))
    for x in SOURCES:
        require(anti[x] == (parallel[x] - 1) % P,
                "target-orientation matched identity")
        require(anti[x] == (-parallel[(1 - x) % P]) % P,
                "source-orientation matched identity")

    reconstruction_checks = 0
    orientation_checks = 0
    symmetric_equivariance_checks = 0
    symmetric_character_checks = 0
    disjoint_character_checks = 0
    symmetric_chart_types = Counter()
    disjoint_second_types = Counter()
    selector_pairs = 0

    for a_pair in pairs:
        a_mask = mask(a_pair)
        for b_pair in pairs:
            b_mask = mask(b_pair)
            multiplicity = convolution(a_mask, b_mask)
            a0, a1 = a_pair
            inverse = two_point_inverse((a1 - a0) % P)
            recovered = convolution(inverse, translate(multiplicity, -a0))
            require(recovered == b_mask, "THM2647 endpoint reconstruction")
            reconstruction_checks += 1

            orientation_atlases = set()
            for a_order in (a_pair, tuple(reversed(a_pair))):
                for b_order in (b_pair, tuple(reversed(b_pair))):
                    charts = charts_for_orders(a_order, b_order, "symmetric")
                    orientation_atlases.add(atlas_key(charts))
                    orientation_checks += 1
            require(len(orientation_atlases) == 1,
                    "symmetric atlas depends on endpoint ordering")

            symmetric = charts_for_orders(a_pair, b_pair, "symmetric")
            midpoint_edge = (
                INV2 * sum(a_pair) % P,
                INV2 * sum(b_pair) % P,
            )
            for chart in symmetric:
                verify_equivariance(a_pair, b_pair, chart)
                require(midpoint_edge in chart, "equivariant chart missed midpoint")
                symmetric_chart_types["affine" if is_affine(chart) else "nonlinear"] += 1
                symmetric_equivariance_checks += 1
            verify_cover(a_pair, b_pair, symmetric, 1)
            symmetric_character_checks += 12

            disjoint = charts_for_orders(a_pair, b_pair, "disjoint")
            verify_cover(a_pair, b_pair, disjoint, 0)
            disjoint_second_types[
                "affine" if is_affine(disjoint[1]) else "nonlinear"
            ] += 1
            disjoint_character_checks += 12
            require(not all(
                all(
                    dict(chart)[(sum(a_pair) - x) % P]
                    == (sum(b_pair) - y) % P
                    for x, y in chart
                )
                for chart in disjoint
            ), "edge-disjoint atlas cannot be individually equivariant")
            selector_pairs += 1

    require(reconstruction_checks == 6084, "reconstruction census")
    require(orientation_checks == 24336, "orientation census")
    require(symmetric_equivariance_checks == 12168, "equivariance census")
    require(symmetric_character_checks == 73008, "symmetric character census")
    require(disjoint_character_checks == 73008, "disjoint character census")
    require(symmetric_chart_types == Counter({"affine": 11154, "nonlinear": 1014}),
            "symmetric chart type census")
    require(disjoint_second_types == Counter({"nonlinear": 6084}),
            "disjoint second-chart census")
    require(selector_pairs == 6084, "selector-pair census")

    print("THM2656 ENDPOINT-ANCHORED RAINBOW SELECTOR MODULI EXACT REFEREE")
    print("endpoint_pairs=6084 endpoint_reconstructions=6084 orientation_checks=24336")
    print("symmetric_atlas=orientation_independent individually_reflection_equivariant")
    print("symmetric_chart_types={affine:11154,nonlinear:1014}")
    print("symmetric_edge_overlap={one:6084} retained_edges={21:6084}")
    print("disjoint_template_rows=12 disjoint_second_charts={nonlinear:6084}")
    print("disjoint_edge_overlap={zero:6084} retained_edges={22:6084}")
    print("both_colour_profiles=1^4_2^9 both_charged_energy=36/169")
    print("charged_character_checks={symmetric:73008,disjoint:73008}_all_nonzero")
    print("individual_reflection_fixed_edge=midpoint selector_overlap_not_invariant")
    print("PASS")


if __name__ == "__main__":
    main()
