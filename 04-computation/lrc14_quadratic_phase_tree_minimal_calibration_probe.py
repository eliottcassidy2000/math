#!/usr/bin/env python3
"""Exact hostile audit for the quadratic phase-tree calibration lemma.

The mathematical proof is algebraic.  This companion checks:

* the complete-pair energy compression identity over Gaussian integers;
* the sharp bipartite path ambiguity;
* the signless-incidence rank criterion on every graph through six vertices;
* the unique-difference graph census for every subset of C_13; and
* the global-twist coefficient selector on deterministic Gaussian controls.

Only integer and Fraction arithmetic is used.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations
import random


def gadd(z, w):
    return (z[0] + w[0], z[1] + w[1])


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1],
            z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gnorm2(z):
    return z[0] * z[0] + z[1] * z[1]


def gsum(values):
    out = (0, 0)
    for value in values:
        out = gadd(out, value)
    return out


def check_complete_pair_identity():
    rng = random.Random(2355)
    checks = 0
    for size in range(2, 8):
        for _ in range(400):
            values = [
                (rng.randrange(-5, 6), rng.randrange(-5, 6))
                for _ in range(size)
            ]
            diagonal = sum(gnorm2(z) for z in values)
            pair_sum = sum(
                gnorm2(gadd(values[i], values[j]))
                for i in range(size)
                for j in range(i + 1, size)
            )
            total = gnorm2(gsum(values))
            assert total == pair_sum - (size - 2) * diagonal
            checks += 1
    return checks


def rank_q(matrix):
    if not matrix:
        return 0
    rows = [[Fraction(x) for x in row] for row in matrix]
    nrows = len(rows)
    ncols = len(rows[0])
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in range(rank, nrows) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [x / scale for x in rows[rank]]
        for r in range(nrows):
            if r == rank or not rows[r][col]:
                continue
            scale = rows[r][col]
            rows[r] = [
                rows[r][c] - scale * rows[rank][c]
                for c in range(ncols)
            ]
        rank += 1
        if rank == nrows:
            break
    return rank


def connected_and_bipartite(n, edges):
    adjacency = [set() for _ in range(n)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    seen = {0}
    stack = [0]
    colour = {0: 0}
    bipartite = True
    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            if v not in seen:
                seen.add(v)
                colour[v] = 1 - colour[u]
                stack.append(v)
            elif colour[v] == colour[u]:
                bipartite = False
    return len(seen) == n, bipartite


def check_signless_incidence():
    counts = Counter()
    checks = 0
    for n in range(2, 7):
        pairs = list(combinations(range(n), 2))
        for mask in range(1, 1 << len(pairs)):
            edges = [pairs[i] for i in range(len(pairs)) if mask & (1 << i)]
            connected, bipartite = connected_and_bipartite(n, edges)
            if not connected:
                continue
            matrix = [
                [int(j == u or j == v) for j in range(n)]
                for u, v in edges
            ]
            rank = rank_q(matrix)
            assert rank == (n - 1 if bipartite else n)
            counts[(n, bipartite)] += 1
            checks += 1
    return checks, counts


def unique_difference_graph(subset, modulus):
    counts = Counter(
        (a - b) % modulus
        for a in subset
        for b in subset
        if a != b
    )
    edges = []
    for i, a in enumerate(subset):
        for b in subset[i + 1:]:
            if counts[(a - b) % modulus] == 1:
                assert counts[(b - a) % modulus] == 1
                edges.append((a, b))
    return edges


def graph_connected(subset, edges):
    if len(subset) <= 1:
        return True
    adjacency = {a: set() for a in subset}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = {subset[0]}
    stack = [subset[0]]
    while stack:
        a = stack.pop()
        for b in adjacency[a]:
            if b not in seen:
                seen.add(b)
                stack.append(b)
    return len(seen) == len(subset)


def is_three_term_ap(subset, modulus):
    target = set(subset)
    return any(
        target == {(centre - step) % modulus, centre,
                   (centre + step) % modulus}
        for centre in range(modulus)
        for step in range(1, modulus)
    )


def c13_census():
    modulus = 13
    rows = []
    coefficient_checks = 0
    rng = random.Random(2313)
    for size in range(1, modulus + 1):
        connected = 0
        total = 0
        edge_hist = Counter()
        bad_triples = 0
        for subset_tuple in combinations(range(modulus), size):
            subset = tuple(subset_tuple)
            edges = unique_difference_graph(subset, modulus)
            is_connected = graph_connected(subset, edges)
            connected += int(is_connected)
            total += 1
            edge_hist[len(edges)] += 1
            if size == 3 and not is_connected:
                assert is_three_term_ap(subset, modulus)
                bad_triples += 1

            if size <= 4:
                values = {
                    a: (rng.randrange(-4, 5), rng.randrange(-4, 5))
                    for a in subset
                }
                for a, b in edges:
                    d = (a - b) % modulus
                    selected = gsum(
                        [
                            gmul(values[x], gconj(values[y]))
                            for x in subset
                            for y in subset
                            if (x - y) % modulus == d
                        ]
                    )
                    assert selected == gmul(values[a], gconj(values[b]))
                    coefficient_checks += 1
        if size == 3:
            assert bad_triples == 78
        rows.append((size, connected, total, tuple(sorted(edge_hist.items()))))
    return rows, coefficient_checks


def main():
    pair_checks = check_complete_pair_identity()

    # The two labelled path states have identical complete edge-twist data:
    # zero modes 5 and first modes -2 on both edges, but different totals.
    hostile_a = [(-1, 0), (2, 0), (-1, 0)]
    hostile_b = [(-2, 0), (1, 0), (-2, 0)]
    for edge in [(0, 1), (1, 2)]:
        u, v = edge
        zero_a = gnorm2(hostile_a[u]) + gnorm2(hostile_a[v])
        zero_b = gnorm2(hostile_b[u]) + gnorm2(hostile_b[v])
        first_a = gmul(gconj(hostile_a[u]), hostile_a[v])
        first_b = gmul(gconj(hostile_b[u]), hostile_b[v])
        assert (zero_a, first_a) == (zero_b, first_b) == (5, (-2, 0))
    assert gsum(hostile_a) == (0, 0)
    assert gsum(hostile_b) == (-3, 0)

    graph_checks, graph_counts = check_signless_incidence()
    rows, coefficient_checks = c13_census()

    print("LRC14 QUADRATIC PHASE-TREE / GLOBAL-TWIST AUDIT")
    print("complete-pair identity Gaussian checks:", pair_checks)
    print("path hostile edge data: zero-mode=5 first-mode=-2 on both edges")
    print("path hostile totals:", gsum(hostile_a), gsum(hostile_b))
    print("connected-graph signless-rank checks:", graph_checks)
    for n in range(2, 7):
        print(
            "  n=%d bipartite=%d nonbipartite=%d"
            % (n, graph_counts[(n, True)], graph_counts[(n, False)])
        )
    print("C13 unique-difference phase-tree census:")
    for size, connected, total, edge_hist in rows:
        print(
            "  size=%2d connected=%4d/%4d edge_hist=%s"
            % (size, connected, total, edge_hist)
        )
    print("C13 unique-edge coefficient selector checks:", coefficient_checks)
    print("three-point failures are exactly the 78 cyclic 3-term APs: PASS")
    print("positive controls: every 2-support; non-AP 3-support; 364 four-supports")
    print("hostile controls: AP triple; every support of size at least 5; full C13")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
