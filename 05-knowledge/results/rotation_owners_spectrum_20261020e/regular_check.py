"""Exact endpoint audit of every regular tournament class at orders 7 and 9.

The order-nine bank is inherited and SHA-256 authenticated. Independent
degree-sequence counting, graph isomorphism, automorphism enumeration, and
orbit-mass equality additionally certify completeness of the regular subset.
Endpoint subset DP is cross-checked by direct path DFS. No new large census.
"""
from collections import Counter
import csv
from functools import cache
import gzip
import hashlib
from itertools import combinations
from math import factorial
from pathlib import Path
import networkx as nx


def need(value, message):
    if not value:
        raise RuntimeError(message)


def decode(code, n):
    result = [0] * n
    bit = n * (n - 1) // 2 - 1
    for i in range(n):
        for j in range(i + 1, n):
            if code >> bit & 1:
                result[i] |= 1 << j
            else:
                result[j] |= 1 << i
            bit -= 1
    return tuple(result)


def circulant(n, steps):
    return tuple(sum(1 << ((i + d) % n) for d in steps) for i in range(n))


def graph(a):
    g = nx.DiGraph()
    g.add_nodes_from(range(len(a)))
    g.add_edges_from((i, j) for i, row in enumerate(a)
                     for j in range(len(a)) if row >> j & 1)
    return g


@cache
def degree_count(degrees):
    """Number of labelled tournaments with this outdegree demand vector.

    Sorting only identifies equal count subproblems; every labelled choice
    of the first vertex's outneighbours is still summed separately.
    """
    n = len(degrees)
    if any(d < 0 or d >= n for d in degrees):
        return 0
    if sum(degrees) != n * (n - 1) // 2:
        return 0
    if n <= 1:
        return 1
    total = 0
    for wins in combinations(range(n - 1), degrees[0]):
        wins = set(wins)
        rest = tuple(sorted(degrees[i + 1] - (i not in wins)
                            for i in range(n - 1)))
        total += degree_count(rest)
    return total


def endpoint_dp(a):
    n = len(a)
    answer = []
    for source in range(n):
        dp = [[0] * n for _ in range(1 << n)]
        dp[1 << source][source] = 1
        for mask in range(1, 1 << n):
            for end, count in enumerate(dp[mask]):
                left = a[end] & ~mask if count else 0
                while left:
                    bit = left & -left
                    left -= bit
                    dp[mask | bit][bit.bit_length() - 1] += count
        answer.append(tuple(dp[-1]))
    return tuple(answer)


def endpoint_dfs(a):
    n = len(a)
    answer = [[0] * n for _ in range(n)]

    def visit(source, end, mask):
        if mask == (1 << n) - 1:
            answer[source][end] += 1
            return
        for nxt in range(n):
            if not mask >> nxt & 1 and a[end] >> nxt & 1:
                visit(source, nxt, mask | (1 << nxt))

    for source in range(n):
        visit(source, source, 1 << source)
    return tuple(tuple(row) for row in answer)


def audit(n, bank):
    d = (n - 1) // 2
    graphs = []
    mass, minimum = 0, None
    print("order", n)
    for code, a, expected_aut in bank:
        need(len(a) == n and all(row.bit_count() == d for row in a), "regular degree")
        need(all(not (a[i] >> i & 1) for i in range(n)), "no loops")
        need(all(bool(a[i] >> j & 1) != bool(a[j] >> i & 1)
                 for i, j in combinations(range(n), 2)), "exactly one arc per pair")
        g = graph(a)
        need(not any(nx.is_isomorphic(g, old) for old in graphs), "distinct isomorphism classes")
        graphs.append(g)
        aut = sum(1 for _ in nx.algorithms.isomorphism.DiGraphMatcher(g, g).isomorphisms_iter())
        need(aut == expected_aut, "independent automorphism order")
        mass += factorial(n) // aut
        endpoints = endpoint_dp(a)
        need(endpoints == endpoint_dfs(a), "independent endpoint DFS")
        off = [endpoints[i][j] for i in range(n) for j in range(n) if i != j]
        need(min(off) > 0, "Hamiltonian connectivity")
        minimum = min(off) if minimum is None else min(minimum, min(off))
        print("code=", code, "rows=", a, "Aut=", aut,
              "H=", sum(map(sum, endpoints)), "min_endpoint=", min(off),
              "endpoint_hist=", sorted(Counter(off).items()))
    exact = degree_count((d,) * n)
    need(mass == exact, "orbit sum exhausts independent labelled degree count")
    print("complete regular classes=", len(bank), "orbit_mass=", mass,
          "independent_labelled_count=", exact, "global_min_endpoint=", minimum)
    return minimum


def main():
    # Three pairwise nonisomorphic regular7 witnesses; orbit mass certifies coverage.
    bank7 = [
        ("carousel", circulant(7, (1, 2, 3)), 7),
        ("Paley", circulant(7, (1, 2, 4)), 21),
        ("triangle-reversal", (26, 28, 41, 112, 100, 67, 7), 3),
    ]
    minimum7 = audit(7, bank7)
    source = Path(__file__).resolve().parent.parent / "blue_quotient_20261020b/n9_classes.tsv.gz"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    need(digest == "959c2276c57108d55b56d4906edfbd761fb9737e708b14114d66dd4ad618db0c",
         "authenticated inherited order-nine bank")
    bank9, total = [], 0
    with gzip.open(source, "rt", encoding="utf-8") as stream:
        for record in csv.DictReader(stream, delimiter="\t"):
            total += 1
            code, aut = int(record["code"]), int(record["Aut"])
            a = decode(code, 9)
            if all(row.bit_count() == 4 for row in a):
                bank9.append((code, a, aut))
    need(total == 191536 and len(bank9) == 15, "inherited universe and regular subset")
    minimum9 = audit(9, bank9)
    need((minimum7, minimum9) == (1, 21), "minimum endpoint controls")
    for n, steps in ((3, (1,)), (5, (1, 2))):
        e = endpoint_dp(circulant(n, steps))
        zeros = [(i, j) for i in range(n) for j in range(n) if i != j and e[i][j] == 0]
        need(zeros, "small regular non-Hamiltonian-connected hostile")
        print("regular small hostile order=", n, "missing ordered endpoint pairs=", zeros)
    print("PASS: every regular tournament of order7 or9 is Hamiltonian-connected.")
    print("All-order inference additionally uses the proved connectivity floor and CITED Thomassen theorem.")


if __name__ == "__main__":
    main()
