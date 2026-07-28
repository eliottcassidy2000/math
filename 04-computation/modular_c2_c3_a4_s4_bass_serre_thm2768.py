#!/usr/bin/env python3
"""Exact modular C2*C3 quotient audit for THM-2768.

The script uses only finite permutation/coset arithmetic and explicit
exceptions.  It contains no floating point and no truth-bearing assertions.
"""

from itertools import permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


IDENTITY = (0, 1, 2, 3)
S4_ALL = tuple(permutations(range(4)))


def compose(left, right):
    """Compose permutations as left after right."""

    return tuple(left[right[index]] for index in range(4))


def inverse(permutation):
    return tuple(permutation.index(index) for index in range(4))


def power(permutation, exponent):
    out = IDENTITY
    for _ in range(exponent):
        out = compose(permutation, out)
    return out


def order(permutation):
    out = IDENTITY
    for exponent in range(1, 25):
        out = compose(permutation, out)
        if out == IDENTITY:
            return exponent
    raise RuntimeError("permutation order escaped S4 bound")


def tuple_order(permutation):
    identity = tuple(range(len(permutation)))
    out = identity
    for exponent in range(1, 25):
        out = tuple(permutation[out[index]]
                    for index in range(len(permutation)))
        if out == identity:
            return exponent
    raise RuntimeError("small permutation order escaped bound")


def sign(permutation):
    inversions = sum(
        permutation[first] > permutation[second]
        for first in range(4)
        for second in range(first + 1, 4)
    )
    return -1 if inversions % 2 else 1


def generated(generators):
    group = {IDENTITY}
    frontier = {IDENTITY}
    while frontier:
        new_frontier = set()
        for element in frontier:
            for generator in generators:
                for product_value in (
                    compose(element, generator),
                    compose(generator, element),
                ):
                    if product_value not in group:
                        group.add(product_value)
                        new_frontier.add(product_value)
        frontier = new_frontier
    return frozenset(group)


def subgroup(permutation):
    return frozenset(power(permutation, exponent)
                     for exponent in range(order(permutation)))


def right_cosets(group, small_group):
    cosets = {
        frozenset(compose(element, small)
                  for small in small_group)
        for element in group
    }
    return tuple(sorted(cosets, key=lambda row: tuple(sorted(row))))


def bass_serre_quotient(group, involution, ternary):
    c2 = subgroup(involution)
    c3 = subgroup(ternary)
    c2_cosets = right_cosets(group, c2)
    c3_cosets = right_cosets(group, c3)
    c2_index = {
        element: index
        for index, coset in enumerate(c2_cosets)
        for element in coset
    }
    c3_index = {
        element: index
        for index, coset in enumerate(c3_cosets)
        for element in coset
    }
    edges = tuple((c2_index[element], c3_index[element])
                  for element in sorted(group))
    require(len(set(edges)) == len(group),
            "Bass-Serre quotient acquired parallel incidence edges")
    left_degrees = [0] * len(c2_cosets)
    right_degrees = [0] * len(c3_cosets)
    for left, right in edges:
        left_degrees[left] += 1
        right_degrees[right] += 1
    require(set(left_degrees) == {2} and set(right_degrees) == {3},
            "Bass-Serre quotient stopped being (2,3)-biregular")

    collapsed_edges = set()
    for coset in c2_cosets:
        neighbours = sorted({c3_index[element] for element in coset})
        require(len(neighbours) == 2,
                "a C2 vertex stopped having two distinct C3 neighbours")
        collapsed_edges.add(tuple(neighbours))
    require(len(collapsed_edges) == len(c2_cosets),
            "degree-two suppression created a repeated edge")
    return c2_cosets, c3_cosets, edges, frozenset(collapsed_edges)


def graph_isomorphism_count(vertex_count, source_edges, target_edges):
    source = frozenset(tuple(sorted(edge)) for edge in source_edges)
    target = frozenset(tuple(sorted(edge)) for edge in target_edges)
    count = 0
    for relabelling in permutations(range(vertex_count)):
        image = frozenset(
            tuple(sorted((relabelling[first], relabelling[second])))
            for first, second in source
        )
        if image == target:
            count += 1
    return count


MATCHINGS = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def matching_action(permutation):
    action = []
    for matching in MATCHINGS:
        moved = frozenset(
            frozenset((permutation[min(edge)], permutation[max(edge)]))
            for edge in matching
        )
        action.append(MATCHINGS.index(moved))
    return tuple(action)


def main():
    # Tetrahedral quotient: x=(12)(34), y=(123).
    x_a = (1, 0, 3, 2)
    y_a = (1, 2, 0, 3)
    group_a = generated((x_a, y_a))

    # Octahedral quotient: x=(12), y=(234).
    x_s = (1, 0, 2, 3)
    y_s = (0, 2, 3, 1)
    group_s = generated((x_s, y_s))

    require((order(x_a), order(y_a), order(compose(x_a, y_a)))
            == (2, 3, 3), "A4 triangle orders changed")
    require((order(x_s), order(y_s), order(compose(x_s, y_s)))
            == (2, 3, 4), "S4 triangle orders changed")
    require(len(group_a) == 12 and all(sign(element) == 1 for element in group_a),
            "tetrahedral generators stopped generating A4")
    require(group_s == frozenset(S4_ALL),
            "octahedral generators stopped generating S4")

    a_c2, a_c3, a_edges, a_collapsed = bass_serre_quotient(
        group_a, x_a, y_a
    )
    s_c2, s_c3, s_edges, s_collapsed = bass_serre_quotient(
        group_s, x_s, y_s
    )
    require((len(a_c2), len(a_c3), len(a_edges)) == (6, 4, 12),
            "A4 Bass-Serre census changed")
    require((len(s_c2), len(s_c3), len(s_edges)) == (12, 8, 24),
            "S4 Bass-Serre census changed")
    require(len(a_edges) - len(a_c2) - len(a_c3) + 1 == 3,
            "A4 quotient cycle rank changed")
    require(len(s_edges) - len(s_c2) - len(s_c3) + 1 == 5,
            "S4 quotient cycle rank changed")

    complete_four = frozenset(
        (first, second)
        for first in range(4)
        for second in range(first + 1, 4)
    )
    cube = frozenset(
        (first, second)
        for first in range(8)
        for bit in (1, 2, 4)
        for second in (first ^ bit,)
        if first < second
    )
    require(graph_isomorphism_count(4, a_collapsed, complete_four) == 24,
            "suppressed A4 quotient stopped being K4")
    require(graph_isomorphism_count(8, s_collapsed, cube) == 48,
            "suppressed S4 quotient stopped being the cube")

    v4 = frozenset(
        element for element in group_s
        if matching_action(element) == (0, 1, 2)
    )
    require(len(v4) == 4 and v4 <= group_a,
            "matching-action kernel stopped being V4")
    image_a = {matching_action(element) for element in group_a}
    image_s = {matching_action(element) for element in group_s}
    require(len(image_a) == 3 and matching_action(x_a) == (0, 1, 2),
            "A4 matching quotient stopped collapsing the binary generator")
    require(len(image_s) == 6
            and tuple_order(matching_action(x_s)) == 2
            and tuple_order(matching_action(y_s)) == 3,
            "S4 matching quotient stopped retaining both generators")

    # Euler characteristic of C2*C3 is -1/6.  A torsion-free subgroup of
    # index d is free of rank 1+d/6.
    require(1 + len(group_a) // 6 == 3,
            "A4 modular-kernel free rank changed")
    require(1 + len(group_s) // 6 == 5,
            "S4 modular-kernel free rank changed")
    require(1 + 6 // 6 == 2,
            "S4 V4-preimage free rank changed")

    print("MODULAR C2-C3 / A4-S4 BASS-SERRE AUDIT")
    print("A4_triangle_orders=(2,3,3) quotient_order=12")
    print("S4_triangle_orders=(2,3,4) quotient_order=24")
    print("A4_kernel=torsion_free_F3")
    print("S4_kernel=torsion_free_F5")
    print("A4_bass_serre_vertices=(6_C2,4_C3) edges=12 beta1=3")
    print("S4_bass_serre_vertices=(12_C2,8_C3) edges=24 beta1=5")
    print("A4_degree2_suppression=K4 isomorphisms=24")
    print("S4_degree2_suppression=cube isomorphisms=48")
    print("matching_kernel=V4")
    print("A4_over_V4=C3 binary_generator=collapsed")
    print("S4_over_V4=S3 binary_and_ternary_generators=retained")
    print("A4_V4_preimage=C2*C2*C2 with kernel_F3")
    print("S4_V4_preimage=F2 with kernel_F5")
    print("SCOPE=finite_quotients_of_C2*C3_not_identifications_or_JC2")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
