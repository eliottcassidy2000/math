#!/usr/bin/env python3
from itertools import combinations

def count(rows, root):
    n = len(rows)
    full = (1 << n) - 1
    def build(starts):
        dp = [[0] * n for _ in range(1 << n)]
        for v in starts:
            dp[1 << v][v] = 1
        for s in range(1, full + 1):
            for u, ways in enumerate(dp[s]):
                if ways:
                    moves = rows[u] & (full ^ s)
                    while moves:
                        bit = moves & -moves
                        dp[s | bit][bit.bit_length() - 1] += ways
                        moves ^= bit
        return [1] + [sum(row) for row in dp[1:]]
    h = build(range(n))
    marked = build([root])
    covers = sum(marked[s] * h[full ^ s] for s in range(1, full + 1))
    return h[full], covers - 3 * h[full]

def expand(rows, root, slots):
    labels = [j for j in range(len(rows)) for _ in range(slots.get(j, 1))]
    matrix = [sum(1 << k for k, y in enumerate(labels)
                  if (j < k if x == y else rows[x] >> y & 1))
              for j, x in enumerate(labels)]
    return count(matrix, labels.index(root))

def verify():
    a = (88,73,107,16,6,91,24)
    b = (104,121,99,100,45,64,16)
    av = [expand(a, 0, sizes) for sizes in ({}, {2:2}, {5:2}, {2:2,5:2})]
    bv = [expand(b, 2, sizes) for sizes in ({}, {1:2})]
    if av != [(39,-28),(97,-82),(69,-54),(187,-160)]:
        raise RuntimeError(av)
    if bv != [(55,-28),(105,-10)]:
        raise RuntimeError(bv)
    peak = 0
    for n in range(1, 5):
        arcs = list(combinations(range(n), 2))
        for orient in range(1 << len(arcs)):
            rows = [0] * n
            for bit, (u,v) in enumerate(arcs):
                if orient >> bit & 1:
                    u,v = v,u
                rows[u] |= 1 << v
            peak = max(peak, count(rows, 0)[0])
    if peak != 5:
        raise RuntimeError(peak)
    print(av, bv, peak)

if __name__ == "__main__":
    verify()
