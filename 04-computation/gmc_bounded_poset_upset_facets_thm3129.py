#!/usr/bin/env python3
"""Exact partition-upset controls for THM-3129."""

from hashlib import sha256


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def partitions(n, cap=None):
    if n == 0:
        yield ()
        return
    if cap is None or cap > n:
        cap = n
    for first in range(cap, 0, -1):
        for tail in partitions(n - first, first):
            yield (first,) + tail


def hasse_edges(types):
    index = {mu: i for i, mu in enumerate(types)}
    edges = set()
    for i, mu in enumerate(types):
        for a in range(len(mu)):
            for b in range(a + 1, len(mu)):
                nu = [mu[j] for j in range(len(mu)) if j not in (a, b)]
                nu.append(mu[a] + mu[b])
                nu = tuple(sorted(nu, reverse=True))
                edges.add((i, index[nu]))
    return tuple(sorted(edges))


def upper_closures(q, edges):
    out = [[] for _ in range(q)]
    for fine, coarse in edges:
        out[fine].append(coarse)
    answer = []
    for root in range(q):
        seen = {root}
        todo = [root]
        while todo:
            current = todo.pop()
            for nxt in out[current]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        answer.append(frozenset(seen))
    return tuple(answer)


def enumerate_upsets(types, closures):
    """Enumerate upsets in a coarsest-first linear extension."""
    q = len(types)
    order = tuple(sorted(range(q), key=lambda i: (len(types[i]), types[i])))
    position = {vertex: i for i, vertex in enumerate(order)}
    require(
        all(position[coarse] <= position[fine]
            for fine in range(q) for coarse in closures[fine]),
        "declared order is not a coarsest-first linear extension",
    )
    answer = []

    def recurse(k, mask):
        if k == q:
            answer.append(mask)
            return
        vertex = order[k]
        recurse(k + 1, mask)
        strict_upper = (u for u in closures[vertex] if u != vertex)
        if all((mask >> u) & 1 for u in strict_upper):
            recurse(k + 1, mask | (1 << vertex))

    recurse(0, 0)
    return tuple(answer)


def connected(mask, adjacency):
    if mask == 0:
        return False
    start = (mask & -mask).bit_length() - 1
    seen = 1 << start
    todo = [start]
    while todo:
        vertex = todo.pop()
        for nxt in adjacency[vertex]:
            bit = 1 << nxt
            if mask & bit and not seen & bit:
                seen |= bit
                todo.append(nxt)
    return seen == mask


EXPECTED = {
    2: (2, 1, 3),
    3: (3, 2, 4),
    4: (5, 5, 7),
    5: (7, 9, 10),
    6: (11, 17, 27),
    7: (15, 28, 47),
    8: (22, 47, 168),
    9: (30, 73, 573),
    10: (42, 114, 3588),
    11: (56, 170, 19542),
}


def partition_census():
    lines = []
    for n, expected in EXPECTED.items():
        types = tuple(partitions(n))
        edges = hasse_edges(types)
        closures = upper_closures(len(types), edges)
        upsets = enumerate_upsets(types, closures)
        require(
            (len(types), len(edges), len(upsets)) == expected,
            f"partition-upset census drift at N={n}",
        )
        adjacency = [set() for _ in types]
        for fine, coarse in edges:
            adjacency[fine].add(coarse)
            adjacency[coarse].add(fine)
        full = (1 << len(types)) - 1
        digest = sha256()
        facets = 0
        for mask in sorted(upsets):
            digest.update(mask.to_bytes((len(types) + 7) // 8, "little"))
            if mask in (0, full):
                continue
            require(connected(mask, adjacency),
                    f"nonempty partition upset disconnected at N={n}")
            require(connected(full ^ mask, adjacency),
                    f"proper partition upset has disconnected complement at N={n}")
            facets += 1
        require(facets == len(upsets) - 2, "facet count is not upset count minus two")
        lines.append(
            f"N={n}:types={len(types)}:edges={len(edges)}:"
            f"upsets={len(upsets)}:nontrivial_facets={facets}:"
            f"digest={digest.hexdigest()}"
        )
    return lines


def disconnected_controls():
    # A bottom with two incomparable maxima.  Their union is a disconnected
    # upset and its indicator is the sum of the two component indicators.
    union = (0, 1, 1)
    left = (0, 1, 0)
    right = (0, 0, 1)
    require(tuple(left[i] + right[i] for i in range(3)) == union,
            "disconnected-upset decomposition failed")

    # Two incomparable minima below one top.  For U={top}, the complement is
    # disconnected and 1_(P\{a})+1_(P\{b})-1_U is the constant one vector.
    u = (0, 0, 1)
    not_a = (0, 1, 1)
    not_b = (1, 0, 1)
    difference = tuple(not_a[i] + not_b[i] - u[i] for i in range(3))
    require(difference == (1, 1, 1),
            "disconnected-complement quotient decomposition failed")
    return (
        "disconnected_upset=two_maxima:indicator=sum_components:PASS",
        "disconnected_complement=two_minima:quotient_difference=constant_one:PASS",
    )


def main():
    for line in partition_census():
        print(line)
    for line in disconnected_controls():
        print(line)
    print("bounded_partition_poset=every_nontrivial_upset_is_a_bond_facet")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
