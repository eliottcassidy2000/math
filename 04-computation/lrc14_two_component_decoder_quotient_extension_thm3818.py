#!/usr/bin/env python3
"""Exact companion for the rank-11, two-component THM-3818 branch.

This file is deliberately assertion-free so that normal and ``python -O``
executions have identical semantics.  It verifies two hostile families and a
small exhaustive control for the scale-slope formula; the general statements
and their proofs are recorded in THM-3818.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import product
from json import dumps
from math import gcd, lcm


Q = 91**6
B = 91**12
ATLAS_SUM_BOUND = 356

requirements = 0


def require(condition, message):
    global requirements
    requirements += 1
    if not condition:
        raise RuntimeError(message)


def factor_exponents(value):
    """Return the exact prime-exponent dictionary by deterministic trial division."""
    require(value >= 1, "factorization input must be positive")
    n = value
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def decoder_edge(x, y):
    """Literal THM-3818 full-decoder predicate for a positive speed pair."""
    require(x > 0 and y > 0 and x != y, "decoder inputs must be distinct positive speeds")
    g = gcd(x, y)
    p, q = sorted((x // g, y // g))
    primitive_sum = p + q
    if primitive_sum > ATLAS_SUM_BOUND:
        return False
    primitive_factors = factor_exponents(primitive_sum)
    if any(exponent > 2 for exponent in primitive_factors.values()):
        return False
    full_sum_factors = factor_exponents(g * primitive_sum)
    return all(prime % 3 == 2 for prime in full_sum_factors)


def decoder_edges(row):
    return tuple(
        (i, j)
        for i in range(len(row))
        for j in range(i + 1, len(row))
        if decoder_edge(row[i], row[j])
    )


def components(vertex_count, edges):
    adjacency = [set() for _ in range(vertex_count)]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    unseen = set(range(vertex_count))
    answer = []
    while unseen:
        root = min(unseen)
        stack = [root]
        unseen.remove(root)
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in sorted(adjacency[vertex], reverse=True):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer, key=lambda component: (len(component), component)))


def reduced_pair_height(x, y):
    g = gcd(x, y)
    return max(x // g, y // g)


def safe_at_one_fourteenth(row):
    return all(min(value % 14, 14 - value % 14) >= 1 for value in row)


def packet_signature(row, edges):
    """Full labelled (M,a,D,all residues) collection on the supplied edges."""
    packets = []
    for i, j in edges:
        x, y = row[i], row[j]
        require(x < y, "hostile rows are expected in increasing label order")
        g = gcd(x, y)
        p, q = x // g, y // g
        covector = tuple(q if k == i else -p if k == j else 0 for k in range(len(row)))
        modulus = x + y
        packets.append(
            (
                i,
                j,
                x**3 + y**3,
                covector,
                modulus,
                tuple(value % modulus for value in row),
            )
        )
    return tuple(packets)


def check_two_nontrivial_hostile():
    # A compact 2+11 witness.  Each component is connected through its
    # speed-one vertex before scaling; no shape is divisible by seven, so the
    # final row has the transparent safe time 1/14.
    shape_a = (1, 3)
    shape_b = (1, 3, 4, 9, 10, 16, 19, 22, 24, 33, 40)
    scale_a = 1
    scale_b = 2**42
    row = tuple(scale_a * value for value in shape_a) + tuple(
        scale_b * value for value in shape_b
    )
    edges = decoder_edges(row)
    comps = components(13, edges)

    require(tuple(sorted(map(len, comps))) == (2, 11), "hostile must have component sizes 2 and 11")
    require(13 - len(comps) == 11, "two-component decoder rank must be eleven")
    require(all((i < 2) == (j < 2) for i, j in edges), "no decoder edge may cross the cut")
    require(len(edges) == 25, "compact star components must have 25 full-decoder edges")

    height_a = max(reduced_pair_height(shape_a[i], shape_a[j]) for i in range(2) for j in range(i + 1, 2))
    height_b = max(reduced_pair_height(shape_b[i], shape_b[j]) for i in range(11) for j in range(i + 1, 11))
    require(height_a == 3 and height_b == 40, "same-component pair heights changed")
    require(max(height_a, height_b) < Q, "every same-component pair must lie below height Q")

    dominance_bound = 2 * Q * max(shape_a)
    require(scale_b > dominance_bound, "dominance must exclude every crossing support-at-most-three row")
    require(sum(row) < B, "hostile must remain inside the THM-763 sum box")
    require(gcd(*row) == 1, "hostile speed row must be primitive")
    require(safe_at_one_fourteenth(row), "hostile must be explicitly lonely at time 1/14")

    edge_a = (0, 1)
    edge_b = (2, 3)
    require(edge_a in edges and edge_b in edges, "one packet edge is required in each component")
    recovered_a = (row[edge_a[0]] + row[edge_a[1]]) // (
        shape_a[edge_a[0]] + shape_a[edge_a[1]]
    )
    recovered_b = (row[edge_b[0]] + row[edge_b[1]]) // (shape_b[0] + shape_b[1])
    require(recovered_a == scale_a and recovered_b == scale_b, "packet moduli must recover both scales")

    return {
        "component_sizes": tuple(sorted(map(len, comps))),
        "decoder_edges": len(edges),
        "decoder_rank": 13 - len(comps),
        "dominance_bound": dominance_bound,
        "max_pair_heights": (height_a, height_b),
        "recovered_scales": (recovered_a, recovered_b),
        "row": row,
        "row_sum": sum(row),
        "safe_time": "1/14",
    }


def check_singleton_packet_fibre():
    main_component = tuple(4**k for k in range(12))
    main_edges = decoder_edges(main_component)
    require(len(main_edges) == 29, "the twelve-vertex power-four component must have 29 edges")
    require(len(components(12, main_edges)) == 1, "the twelve-vertex component must be connected")

    packet_lcm = 1
    for i, j in main_edges:
        packet_lcm = lcm(packet_lcm, main_component[i] + main_component[j])
    require(packet_lcm == 22906142720, "packet-modulus lcm changed")

    singleton_0 = 2 * Q * max(main_component) + 1
    singleton_1 = singleton_0 + packet_lcm
    row_0 = main_component + (singleton_0,)
    row_1 = main_component + (singleton_1,)
    edges_0 = decoder_edges(row_0)
    edges_1 = decoder_edges(row_1)

    require(edges_0 == main_edges and edges_1 == main_edges, "singleton must remain decoder-isolated")
    require(tuple(sorted(map(len, components(13, edges_0)))) == (1, 12), "component sizes must be 1 and 12")
    require(13 - len(components(13, edges_0)) == 11, "singleton hostile decoder rank must be eleven")
    require(singleton_0 > 2 * Q * max(main_component), "dominance threshold must be strict")
    require(sum(row_1) < B, "both singleton hostiles must lie in the THM-763 box")
    require(packet_signature(row_0, edges_0) == packet_signature(row_1, edges_1), "full packet collections must coincide")
    require(safe_at_one_fourteenth(row_0) and safe_at_one_fourteenth(row_1), "both singleton rows must be safe at 1/14")

    dominance_tail_count = (B - sum(main_component) - singleton_0) // packet_lcm + 1
    require(dominance_tail_count == 14077914720208, "dominance-preserving packet fibre count changed")

    return {
        "component_sizes": (1, 12),
        "decoder_edges": len(edges_0),
        "decoder_rank": 11,
        "dominance_tail_count": dominance_tail_count,
        "main_component_sum": sum(main_component),
        "packet_lcm": packet_lcm,
        "singleton_pair": (singleton_0, singleton_1),
        "singleton_residue": singleton_0 % packet_lcm,
    }


def check_small_slope_atlas():
    """Exhaust all height-three crossing rows on a 2+2 toy quotient."""
    shape_a = (1, 4)
    shape_b = (1, 4)
    selected_scales = (2, 3)
    slopes = set()
    realized_rows = 0
    crossing_rows = 0
    for coefficients in product(range(-3, 4), repeat=4):
        support = tuple(i for i, value in enumerate(coefficients) if value)
        if not 2 <= len(support) <= 3:
            continue
        if not any(i < 2 for i in support) or not any(i >= 2 for i in support):
            continue
        crossing_rows += 1
        partial_a = sum(coefficients[i] * shape_a[i] for i in range(2))
        partial_b = sum(coefficients[i + 2] * shape_b[i] for i in range(2))
        require(partial_a != 0 and partial_b != 0, "a support-at-most-three crossing row cannot vanish on both sides")
        if partial_a * partial_b < 0:
            slopes.add(Fraction(-partial_b, partial_a))  # scale_a / scale_b
        if selected_scales[0] * partial_a + selected_scales[1] * partial_b == 0:
            realized_rows += 1
            require(
                Fraction(selected_scales[0], selected_scales[1])
                == Fraction(-partial_b, partial_a),
                "a crossing row must pin the quotient slope",
            )
    require(crossing_rows == 1008, "toy crossing-row universe changed")
    require(len(slopes) == 95, "toy positive slope atlas changed")
    require(realized_rows == 12, "toy selected-slope fibre changed")
    return {
        "coefficient_height": 3,
        "crossing_rows": crossing_rows,
        "positive_reduced_slopes": len(slopes),
        "realized_rows_at_scale_2_3": realized_rows,
    }


def main():
    two_nontrivial = check_two_nontrivial_hostile()
    singleton = check_singleton_packet_fibre()
    slope_control = check_small_slope_atlas()
    semantic = {
        "B": B,
        "Q": Q,
        "singleton": singleton,
        "slope_control": slope_control,
        "two_nontrivial": two_nontrivial,
    }
    semantic_hash = sha256(dumps(semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    print("THM-3818 TWO-COMPONENT QUOTIENT / PACKET-FIBRE EXACT PROBE")
    print(f"Q={Q}")
    print(f"B={B}")
    print("two_nontrivial=" + dumps(two_nontrivial, sort_keys=True, separators=(",", ":")))
    print("singleton=" + dumps(singleton, sort_keys=True, separators=(",", ":")))
    print("slope_control=" + dumps(slope_control, sort_keys=True, separators=(",", ":")))
    print(f"requirements={requirements}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
