#!/usr/bin/env python3
"""Exact small-universe audit for THM-2183."""

from functools import lru_cache
from itertools import permutations


def pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


def edge(mask: int, n: int, u: int, v: int) -> bool:
    """Whether u beats v; bit 0 means the smaller label beats the larger."""
    if u == v:
        raise ValueError("loops are not tournament edges")
    if u > v:
        return not edge(mask, n, v, u)
    index = pairs(n).index((u, v))
    return not bool((mask >> index) & 1)


def encode(n: int, predicate) -> int:
    mask = 0
    for index, (u, v) in enumerate(pairs(n)):
        if not predicate(u, v):
            mask |= 1 << index
    return mask


def relabel(mask: int, n: int, permutation: tuple[int, ...]) -> int:
    """Encode the target-labelled copy under source->target permutation."""
    inverse = [0] * n
    for source, target in enumerate(permutation):
        inverse[target] = source
    return encode(
        n,
        lambda u, v: edge(mask, n, inverse[u], inverse[v]),
    )


@lru_cache(maxsize=None)
def orbit(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                relabel(mask, n, permutation)
                for permutation in permutations(range(n))
            }
        )
    )


@lru_cache(maxsize=None)
def d_iso(mask_1: int, mask_2: int, n: int) -> int:
    return min(
        (mask_1 ^ copy).bit_count()
        for copy in orbit(mask_2, n)
    )


def order_join(left: int, a: int, right: int, b: int) -> int:
    return encode(
        a + b,
        lambda u, v: (
            edge(left, a, u, v)
            if v < a
            else (
                True
                if u < a
                else edge(right, b, u - a, v - a)
            )
        ),
    )


def theorem_audit() -> int:
    quadruples = 0
    for a in range(1, 4):
        masks_a = range(1 << len(pairs(a)))
        for b in range(1, 4):
            masks_b = range(1 << len(pairs(b)))
            for left_1 in masks_a:
                for left_2 in masks_a:
                    left_distance = d_iso(left_1, left_2, a)
                    for right_1 in masks_b:
                        joined_1 = order_join(left_1, a, right_1, b)
                        for right_2 in masks_b:
                            joined_2 = order_join(left_2, a, right_2, b)
                            expected = (
                                left_distance
                                + d_iso(right_1, right_2, b)
                            )
                            actual = d_iso(joined_1, joined_2, a + b)
                            if actual != expected:
                                raise AssertionError(
                                    (
                                        a,
                                        b,
                                        left_1,
                                        left_2,
                                        right_1,
                                        right_2,
                                        actual,
                                        expected,
                                    )
                                )
                            quadruples += 1
    return quadruples


def quotient_graph(n: int) -> tuple[tuple[int, ...], set[tuple[int, int]]]:
    representatives = tuple(
        sorted(
            {
                min(orbit(mask, n))
                for mask in range(1 << len(pairs(n)))
            }
        )
    )
    edges = {
        (first, second)
        for index, first in enumerate(representatives)
        for second in representatives[index + 1 :]
        if d_iso(first, second, n) == 1
    }
    return representatives, edges


def main() -> None:
    quadruples = theorem_audit()
    representatives, graph_edges = quotient_graph(4)
    triangle = (0, 2, 4)

    if len(representatives) != 4:
        raise AssertionError(representatives)
    if len(graph_edges) != 5:
        raise AssertionError(graph_edges)
    if not all(
        d_iso(first, second, 4) == 1
        for index, first in enumerate(triangle)
        for second in triangle[index + 1 :]
    ):
        raise AssertionError("claimed masks do not form a triangle")

    print("THM-2183 ORDER-JOIN METRIC PRODUCT -- exact audit")
    print("factor_orders=1..3 by 1..3")
    print(f"labelled_quadruples={quadruples}")
    print("metric_product_failures=0")
    print(f"n4_quotient_vertices={len(representatives)}")
    print(f"n4_quotient_edges={len(graph_edges)}")
    print(f"n4_triangle={triangle}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
