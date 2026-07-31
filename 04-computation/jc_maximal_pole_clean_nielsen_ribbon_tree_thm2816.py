#!/usr/bin/env python3
"""Exact permutation and Prüfer controls for THM-2816."""

import ast
from collections import defaultdict
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(ast.parse(Path(path).read_text(encoding="utf-8")))
    )


def positive_compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(1, total - length + 2):
        for tail in positive_compositions(total - first, length - 1):
            yield (first,) + tail


def noncrossing_matchings(points):
    points = tuple(points)
    if not points:
        yield ()
        return
    first = points[0]
    for partner_index in range(1, len(points), 2):
        partner = points[partner_index]
        inside = points[1:partner_index]
        outside = points[partner_index + 1 :]
        for inside_matching in noncrossing_matchings(inside):
            for outside_matching in noncrossing_matchings(outside):
                yield ((first, partner),) + inside_matching + outside_matching


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def cycles(permutation):
    seen = set()
    answer = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point)
            point = permutation[point]
        answer.append(tuple(cycle))
    return tuple(answer)


def rotated_matching(matching, shift, degree):
    return tuple(
        sorted(
            tuple(sorted(((left + shift) % degree, (right + shift) % degree)))
            for left, right in matching
        )
    )


def labelled_orbit_key(matching, labelled_cycles, degree):
    candidates = []
    for shift in range(degree):
        shifted_matching = rotated_matching(matching, shift, degree)
        shifted_cycles = tuple(
            tuple(sorted((point + shift) % degree for point in cycle))
            for cycle in labelled_cycles
        )
        candidates.append((shifted_matching, shifted_cycles))
    return min(candidates)


def direct_nielsen_atlas(e, degree):
    long_cycle = tuple((index + 1) % degree for index in range(degree))
    atlas = defaultdict(set)
    raw = 0
    for endpoints in combinations(range(degree), 2 * e):
        for matching in noncrossing_matchings(endpoints):
            zero_inertia = list(range(degree))
            for left, right in matching:
                zero_inertia[left], zero_inertia[right] = (
                    zero_inertia[right],
                    zero_inertia[left],
                )
            pole_cycles = cycles(compose(tuple(zero_inertia), long_cycle))
            require(len(pole_cycles) == e + 1, "noncrossing maximum-pole cycle count")
            raw += 1
            for labelled_cycles in permutations(pole_cycles):
                passport = tuple(len(cycle) for cycle in labelled_cycles)
                atlas[passport].add(
                    labelled_orbit_key(matching, labelled_cycles, degree)
                )

    catalan = comb(2 * e, e) // (e + 1)
    require(raw == catalan * comb(degree, 2 * e), "raw Catalan matching count")
    expected = factorial(e - 1) * comb(degree - e - 1, e - 1)
    compositions = tuple(positive_compositions(degree, e + 1))
    require(len(atlas) == len(compositions) == comb(degree - 1, e), "passport count")
    for passport in compositions:
        require(len(atlas[passport]) == expected, "direct ordered Nielsen count")
        require(prufer_ribbon_count(passport) == expected, "Pruefer/ribbon count")
    return raw, len(atlas), expected


def prufer_ribbon_count(passport):
    vertex_count = len(passport)
    e = vertex_count - 1
    direct = 0
    for excesses in weak_compositions(e - 1, vertex_count):
        direct += factorial(e - 1) * product_of(
            comb(part - 1, excess)
            for part, excess in zip(passport, excesses)
        )
    require(direct == factorial(e - 1) * comb(sum(passport) - e - 1, e - 1),
            "Vandermonde convolution")

    # Explicit tree/ribbon enumeration: each Pruefer word is one labelled
    # tree and contributes prod(k_i! binom(p_i-1,k_i)).
    tree_total = 0
    for word in product(range(vertex_count), repeat=max(0, vertex_count - 2)):
        excess_degrees = [0] * vertex_count
        for vertex in word:
            excess_degrees[vertex] += 1
        tree_total += product_of(
            factorial(excess) * comb(part - 1, excess)
            for part, excess in zip(passport, excess_degrees)
        )
    require(tree_total == direct, "Pruefer degree/ribbon cancellation")
    return tree_total


def weak_compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first,) + tail


def product_of(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    print("THM-2816 MAXIMAL-POLE CLEAN NIELSEN/RIBBON-TREE AUDIT")
    print("e,N | raw noncrossing | ordered passports | classes per passport")

    cases = (
        *((1, degree) for degree in range(2, 8)),
        *((2, degree) for degree in range(4, 10)),
        *((3, degree) for degree in range(6, 10)),
        *((4, degree) for degree in range(8, 10)),
        (5, 10),
    )
    for e, degree in cases:
        raw, passports, classes = direct_nielsen_atlas(e, degree)
        print(
            f"{e:1d},{degree:2d} | {raw:15d} |"
            f" {passports:17d} | {classes:20d}"
        )

    for e in range(1, 9):
        for degree in range(2 * e, 4 * e + 5):
            expected = factorial(e - 1) * comb(degree - e - 1, e - 1)
            passport = next(positive_compositions(degree, e + 1))
            require(
                expected
                == factorial(e - 1)
                * sum(
                    product_of(
                        comb(part - 1, excess)
                        for part, excess in zip(
                            passport,
                            excesses,
                        )
                    )
                    for excesses in weak_compositions(e - 1, e + 1)
                ),
                "fixed-passport Vandermonde control",
            )

    print("general formula: (e-1)!*binom(N-e-1,e-1)")
    print("e=3 ternary sequence N=6..10: 2,6,12,20,30")
    print("noncrossing, labelled-ribbon-tree, Pruefer, and Vandermonde counts: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
