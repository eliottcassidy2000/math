#!/usr/bin/env python3
"""Exact six/seven-vertex comparisons and exhaustive D5/D6 bases.

Universe: every labelled switching class, represented by all-positive edges
at vertex zero. Primary path: integer Walsh transform of simple-cycle masks.
The independent C++ companion uses direct parity, without a transform.
"""

from collections import Counter
from itertools import combinations, permutations
from math import comb, factorial
import hashlib
import json

import numpy as np
import sympy as sp


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def cycles(n, k, index):
    result = []
    for vertices in combinations(range(n), k):
        root, *tail = vertices
        for order in permutations(tail):
            if order[0] > order[-1]:
                continue
            word = (root,) + order
            mask = 0
            for a, b in zip(word, word[1:] + word[:1]):
                if a and b:
                    mask ^= 1 << index[tuple(sorted((a, b)))]
            result.append(mask)
    require(len(result) == factorial(n) // (2*k*factorial(n-k)), "cycle count")
    require(len(set(result)) == len(result), "cycle mask uniqueness")
    return result


def transform(values):
    result = values.copy()
    width = 1
    while width < len(result):
        blocks = result.reshape(-1, 2, width)
        a, b = blocks[:, 0].copy(), blocks[:, 1].copy()
        blocks[:, 0], blocks[:, 1] = a+b, a-b
        width *= 2
    return result


def single_edge_masks(n, index):
    result = {1 << j for j in index.values()}
    for v in range(1, n):
        result.add(sum(1 << j for edge, j in index.items() if v in edge))
    require(len(result) == comb(n, 2), "labelled edge classes")
    return result


def main():
    rows = []
    for n in range(6, 9):
        index = {edge: j for j, edge in enumerate(combinations(range(1, n), 2))}
        size = 1 << len(index)
        counts = []
        for k in range(3, min(n, 7)+1):
            masks = cycles(n, k, index)
            frequency = np.zeros(size, dtype=np.int64)
            frequency[masks] = 1
            difference = len(masks) - transform(frequency)
            require(bool(np.all(difference % 2 == 0)), "integral negative counts")
            counts.append(difference // 2)
        require(all(int(c[0]) == 0 for c in counts), "balanced hostile")
        total = sum(counts[:4])
        minimum = int(total[1:].min())
        minimizers = set(map(int, np.flatnonzero(total == minimum)))
        edge_minimum = sum(factorial(n-2)//factorial(n-k) for k in range(3, 7))
        require(minimum == edge_minimum, "D5 base minimum")
        require(minimizers == single_edge_masks(n, index), "D5 equality classes")
        require(bool(np.all(3*counts[3] >= (n-4)*(n-5)*counts[1])), "lifted comparison")
        anti = size-1
        require(int(counts[1][anti]) == int(counts[3][anti]) == 0, "antibalanced even hostile")
        require(int(counts[0][anti]) == comb(n, 3), "antibalanced triangles")
        row = {"n": n, "classes": size, "minimum": minimum,
               "minimizers": len(minimizers), "gap": 2*minimum}
        if n >= 7:
            d6 = sum(counts)
            d6min = int(d6[1:].min())
            d6masks = set(map(int, np.flatnonzero(d6 == d6min)))
            require(d6min == sum(factorial(n-2)//factorial(n-k) for k in range(3, 8)), "D6 base minimum")
            require(d6masks == single_edge_masks(n, index), "D6 equality classes")
            comparison = 2*counts[4]-comb(n-4, 3)*counts[1]-comb(n-5, 2)*counts[2]
            require(bool(np.all(comparison >= 0)), "seven-vertex lifted comparison")
            row["d6"] = {"minimum": d6min, "minimizers": len(d6masks), "gap": 2*d6min}
            if n == 7:
                row["seven_vertex_comparison_minimum_nontrivial"] = int(comparison[1:].min())
        if n == 6:
            profiles = Counter(tuple(int(c[m]) for c in counts) for m in range(size))
            row["profiles"] = [[*profile, multiplicity] for profile, multiplicity in sorted(profiles.items())]
            sharp = [m for m in range(size) if counts[1][m] and 3*counts[3][m] == 2*counts[1][m]]
            require(bool(sharp), "nonzero sharp comparison witness")
            witness = sharp[0]
            row["sharp_witness"] = {"mask": witness,
                "negative_edges": [list(e) for e, j in index.items() if witness & (1 << j)]}
        rows.append(row)
    n = sp.Symbol("n")
    A = lambda x: (x-2)*(x**3-11*x*x+41*x-50)
    gap = sp.expand(n*A(n-1)-(n-5)*A(n)-n*(n-1)*(n-2)/3)
    require(sp.expand(gap-(n-5)*(n-4)*(n-3)*(3*n-25)/3) == 0, "induction identity")
    A6 = lambda x: A(x)+(x-2)*(x-3)*(x-4)*(x-5)*(x-6)
    gap6 = sp.expand(n*A6(n-1)-(n-6)*A6(n)-n*(n-1)*(n-2)/2)
    require(sp.expand(gap6-(n-4)*(n-3)*(2*n**3-42*n**2+285*n-620)/2) == 0, "D6 induction identity")
    result = {"status": "PASS", "method": "integer Walsh", "rows": rows,
              "induction_gap": str(sp.factor(gap)),
              "d6_induction_gap": str(sp.factor(gap6)),
              "scope": "D5 for every n>=6 and D6 for every n>=7; D>=7 remains open"}
    print(json.dumps(result, indent=2, sort_keys=True))
    print("semantic_sha256=" + hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest())


if __name__ == "__main__":
    main()
