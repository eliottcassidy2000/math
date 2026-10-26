#!/usr/bin/env python3
"""Independent brute-permutation audit of the order-five lumpability hostile.

No imports from verify_recursion.py. It checks every attachment and all 6!
permutations of every child; no degree pruning or inherited binary is used.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows_to_matrix(rows):
    n = len(rows)
    return [[(rows[i] >> j) & 1 for j in range(n)] for i in range(n)]


def chain_sizes(matrix):
    n = len(matrix)
    unseen = set(range(n))
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)
             if all(matrix[i][k] == matrix[j][k] for k in range(n) if k not in (i, j))]
    components = []
    while unseen:
        component = {min(unseen)}
        while True:
            enlarged = component | {j for i, j in edges if i in component} | {i for i, j in edges if j in component}
            if enlarged == component:
                break
            component = enlarged
        unseen -= component
        components.append(component)
    return sorted(map(len, components))


def exhaustive_aut_and_canonical(matrix):
    n = len(matrix)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    original = tuple(matrix[i][j] for i, j in pairs)
    best = None
    aut = 0
    for p in permutations(range(n)):
        word = tuple(matrix[p[i]][p[j]] for i, j in pairs)
        aut += word == original
        if best is None or word < best:
            best = word
    return aut, best


def main():
    results = []
    parents = [(16, 1, 3, 7, 14), (2, 24, 3, 5, 13)]
    expected = [{3: "4/3", 4: "1/3", 5: "7/3", 6: "2"},
                {3: "2", 4: "1/3", 5: "3", 6: "4/3"}]
    for index, rows in enumerate(parents):
        matrix = rows_to_matrix(rows)
        parent_aut, _ = exhaustive_aut_and_canonical(matrix)
        need(parent_aut == 1, "parents are rigid")
        sizes = chain_sizes(matrix)
        need(sizes == ([1, 1, 3] if index == 0 else [1, 2, 2]), "clone weights")
        raw_hist = Counter()
        aut_hist = Counter()
        weighted = Counter()
        codes = set()
        attachments = []
        for mask in range(32):
            child = [row[:] + [1 - (mask >> i & 1)] for i, row in enumerate(matrix)]
            child.append([(mask >> i) & 1 for i in range(5)] + [0])
            aut, canon = exhaustive_aut_and_canonical(child)
            b = len(chain_sizes(child))
            codes.add(canon)
            raw_hist[b] += 1
            aut_hist[aut] += 1
            weighted[b] += aut
            attachments.append([mask, b, aut])
        weights = {b: str(Fraction(value, 6)) for b, value in sorted(weighted.items())}
        need(weights == expected[index], "exact unlabelled transfer mass")
        need(raw_hist == {3: 8, 4: 2, 5: 14, 6: 8}, "same labelled core kernel")
        need(len(codes) == (24 if index == 0 else 23), "distinct children")
        results.append({"rows": rows, "chain_sizes": sizes,
                        "raw_core_histogram": dict(sorted(raw_hist.items())),
                        "child_aut_histogram": dict(sorted(aut_hist.items())),
                        "weighted_core_transition": weights,
                        "distinct_children": len(codes), "attachment_certificate": attachments})
    print(json.dumps({"status": "PASS", "child_permutations_checked": 46080,
                      "parents": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
