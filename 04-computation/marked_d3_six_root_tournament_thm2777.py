#!/usr/bin/env python3
"""Exact marked-D3 determinant tournament audit for THM-2777.

Only integer matrix, permutation, tournament, and gcd arithmetic is used.
Truth checks use explicit exceptions; Python assertions and floating point
are absent.
"""

from itertools import combinations, permutations, product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(rows):
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def matvec(matrix, vector):
    return tuple(sum(row[index] * vector[index]
                     for index in range(3)) for row in matrix)


def permutation_sign(permutation):
    inversions = sum(
        permutation[first] > permutation[second]
        for first in range(3)
        for second in range(first + 1, 3)
    )
    return -1 if inversions % 2 else 1


def weyl_d3():
    matrices = set()
    for coordinate_permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1:
                continue
            rows = []
            for row_index in range(3):
                row = [0, 0, 0]
                row[coordinate_permutation[row_index]] = signs[row_index]
                rows.append(tuple(row))
            matrices.add(tuple(rows))
    return frozenset(matrices)


def two_minor_gcd(rows):
    values = []
    for row_first, row_second in combinations(range(3), 2):
        for col_first, col_second in combinations(range(3), 2):
            values.append(abs(
                rows[row_first][col_first] * rows[row_second][col_second]
                - rows[row_first][col_second] * rows[row_second][col_first]
            ))
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def tournament_data(roots, h):
    adjacency = [[False] * len(roots) for _ in roots]
    weights = {}
    determinants = {}
    for first, second in combinations(range(len(roots)), 2):
        value = determinant((roots[first], roots[second], h))
        require(value != 0, "marked root contraction acquired a tie")
        determinants[(first, second)] = value
        weights[(first, second)] = abs(value)
        if value > 0:
            adjacency[first][second] = True
        else:
            adjacency[second][first] = True
    return adjacency, weights, determinants


def directed_triangle_count(adjacency):
    count = 0
    for triple in combinations(range(len(adjacency)), 3):
        scores = [
            sum(adjacency[vertex][other]
                for other in triple if other != vertex)
            for vertex in triple
        ]
        if sorted(scores) == [1, 1, 1]:
            count += 1
    return count


def hamilton_path_count(adjacency):
    return sum(
        all(adjacency[path[index]][path[index + 1]]
            for index in range(len(path) - 1))
        for path in permutations(range(len(adjacency)))
    )


def automorphism_count(adjacency):
    size = len(adjacency)
    return sum(
        all(adjacency[first][second]
            == adjacency[permutation[first]][permutation[second]]
            for first in range(size) for second in range(size))
        for permutation in permutations(range(size))
    )


def switched_key(adjacency, mask):
    size = len(adjacency)
    bits = []
    for first, second in combinations(range(size), 2):
        value = adjacency[first][second]
        if ((mask >> first) & 1) != ((mask >> second) & 1):
            value = not value
        bits.append(int(value))
    return tuple(bits)


def unit(index):
    row = [0, 0, 0]
    row[index] = 1
    return tuple(row)


def add(left, right, right_scale=1):
    return tuple(left[index] + right_scale * right[index]
                 for index in range(3))


def chamber_roots(order):
    first, second, third = order
    return (
        add(unit(first), unit(second), -1),
        add(unit(first), unit(second)),
        add(unit(first), unit(third), -1),
        add(unit(first), unit(third)),
        add(unit(second), unit(third), -1),
        add(unit(second), unit(third)),
    )


def main():
    names = ("a12", "b12", "a13", "b13", "a23", "b23")
    roots = chamber_roots((0, 1, 2))
    h = (1, 1, 1)
    a_indices = {0, 2, 4}
    binary_pairs = {frozenset((0, 1)), frozenset((2, 3)), frozenset((4, 5))}
    ternary_pairs = {frozenset(pair) for pair in combinations(a_indices, 2)}

    group = weyl_d3()
    orbit = {matvec(matrix, h) for matrix in group}
    stabilizer = {matrix for matrix in group if matvec(matrix, h) == h}
    expected_orbit = {
        (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)
    }
    require(len(group) == 24 and orbit == expected_orbit,
            "W(D3) body-diagonal orbit changed")
    require(len(stabilizer) == 6
            and all(all(value >= 0 for row in matrix for value in row)
                    for matrix in stabilizer),
            "marked-state stabilizer stopped being coordinate S3")
    diagonal_v4 = {
        ((signs[0], 0, 0), (0, signs[1], 0), (0, 0, signs[2]))
        for signs in product((-1, 1), repeat=3)
        if signs[0] * signs[1] * signs[2] == 1
    }
    require({matvec(matrix, h) for matrix in diagonal_v4} == orbit,
            "diagonal V4 stopped acting simply transitively on body diagonals")

    require({index for index, root in enumerate(roots)
             if sum(root) == 0} == a_indices,
            "orthogonal A2 subsystem changed")

    adjacency, weights, determinants = tournament_data(roots, h)
    weight_one = {frozenset(pair) for pair, value in weights.items() if value == 1}
    weight_two = {frozenset(pair) for pair, value in weights.items() if value == 2}
    weight_three = {frozenset(pair) for pair, value in weights.items() if value == 3}
    require(len(weight_one) == 9 and weight_two == binary_pairs
            and weight_three == ternary_pairs,
            "binary/ternary determinant edge classification changed")

    # The half-Hadamard line dictionary identifies the six roots with K4 edges.
    k4_edges = (
        frozenset((2, 3)), frozenset((1, 4)),
        frozenset((2, 4)), frozenset((1, 3)),
        frozenset((3, 4)), frozenset((1, 2)),
    )
    for first, second in combinations(range(6), 2):
        pair = frozenset((first, second))
        disjoint = k4_edges[first].isdisjoint(k4_edges[second])
        both_opposite_triangle = (1 not in k4_edges[first]
                                  and 1 not in k4_edges[second])
        expected = 2 if disjoint else (3 if both_opposite_triangle else 1)
        require(weights[(first, second)] == expected,
                "K4 edge geometry stopped predicting determinant weight")

    # A unit two-minor makes each full frame cokernel cyclic of determinant
    # order.  Binary pairs row-reduce using h-b_ij=e_k; ternary pairs have
    # exact coordinate-sum quotient modulo three.
    for first, second in combinations(range(6), 2):
        frame = (roots[first], roots[second], h)
        require(two_minor_gcd(frame) == 1,
                "a marked full frame lost its cyclic Smith quotient")
        pair = frozenset((first, second))
        if pair in binary_pairs:
            b_index = first if first % 2 else second
            support = {index for index, value in enumerate(roots[b_index]) if value}
            complement = ({0, 1, 2} - support).pop()
            require(add(h, roots[b_index], -1) == unit(complement),
                    "binary frame stopped row-reducing to C2 plus a coordinate")
        if pair in ternary_pairs:
            require(all(sum(row) % 3 == 0 for row in frame),
                    "ternary frame left the coordinate-sum kernel")

    scores = tuple(sum(row) for row in adjacency)
    require(scores == (3, 2, 3, 4, 1, 2),
            "displayed tournament score sequence changed")
    triangles = directed_triangle_count(adjacency)
    paths = hamilton_path_count(adjacency)
    automorphisms = automorphism_count(adjacency)
    require(triangles == 6 and paths == 23 and automorphisms == 1,
            "displayed tournament census changed")

    # Reordering the A2 chamber by an even coordinate permutation preserves
    # the labelled positional tournament; an odd permutation reverses it.
    for order in permutations(range(3)):
        moved_adjacency, moved_weights, _ = tournament_data(chamber_roots(order), h)
        require(moved_weights == weights,
                "coordinate chamber permutation changed determinant magnitudes")
        sign = permutation_sign(order)
        require(all(moved_adjacency[first][second]
                    == (adjacency[first][second] if sign == 1
                        else adjacency[second][first])
                    for first in range(6) for second in range(6)
                    if first != second),
                "chamber parity stopped giving same/converse tournaments")

    switching_keys = {switched_key(adjacency, mask) for mask in range(1 << 6)}
    require(len(switching_keys) == 32,
            "unoriented root-line switching class changed")

    binary_text = ",".join(
        f"{names[first]}-{names[second]}"
        for first, second in ((0, 1), (2, 3), (4, 5))
    )
    ternary_text = ",".join(
        f"{names[first]}-{names[second]}"
        for first, second in ((0, 2), (0, 4), (2, 4))
    )

    print("MARKED D3 SIX-ROOT DETERMINANT TOURNAMENT AUDIT")
    print("W(D3)_order=24 h_orbit=4 stabilizer=W(A2)=S3_order6")
    print("half_Hadamard_edge_map=a12:23,b12:14,a13:24,b13:13,a23:34,b23:12")
    print("pair_determinants=abs1:9,abs2:3,abs3:3 ties=0")
    print("binary_opposite_pairs=" + binary_text)
    print("ternary_A2_pairs=" + ternary_text)
    print("Smith_cokernels=trivial:9,Z2:3,Z3:3")
    print("tournament_scores=(3,2,3,4,1,2) sorted=(1,2,2,3,3,4)")
    print(f"directed_3cycles={triangles} Hamilton_paths={paths} automorphisms={automorphisms}")
    print("chamber_parity=even:same,odd:converse")
    print(f"unoriented_root_line_switching_class={len(switching_keys)}")
    print("SCOPE=marked_oriented_frame_not_canonical_unmarked_tournament")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
