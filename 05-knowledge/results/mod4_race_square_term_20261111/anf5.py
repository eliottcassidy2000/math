#!/usr/bin/env python3
"""Algebraic normal form over F_2 of y(T) = c_odd(T) mod 2 = (H(T)-1)/2 mod 2 for tournaments on k vertices
(k = 3, 4, 5), in the arc variables x_ij (i<j, x_ij = 1 iff i -> j).  Monomials are classified by degree and by the
number of vertices they span, and by isomorphism type of the underlying graph."""
from itertools import combinations, permutations
from collections import Counter
import sys
def ham(n, out):
    full = (1 << n) - 1
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, full + 1):
        for v in range(n):
            c = f[S][v]
            if c:
                o = out[v] & ~S
                while o:
                    w = (o & -o).bit_length() - 1; o &= o - 1
                    f[S | 1 << w][w] += c
    return sum(f[full])
def anf(k):
    prs = list(combinations(range(k), 2)); N = len(prs)
    vals = [0] * (1 << N)
    for code in range(1 << N):
        out = [0] * k
        for t, (i, j) in enumerate(prs):
            if code >> t & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        vals[code] = ((ham(k, out) - 1) // 2) & 1
    a = vals[:]
    for t in range(N):
        for x in range(1 << N):
            if x >> t & 1: a[x] ^= a[x ^ (1 << t)]
    return prs, a
def graph_type(edges, k):
    # canonical form: sorted degree sequence + number of components, + edge multiset under all relabelings (small k)
    verts = sorted(set(v for e in edges for v in e))
    best = None
    for p in permutations(range(len(verts))):
        mp = {v: p[i] for i, v in enumerate(verts)}
        key = tuple(sorted(tuple(sorted((mp[a], mp[b]))) for a, b in edges))
        if best is None or key < best: best = key
    return (len(verts), best)
k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
prs, a = anf(k)
N = len(prs)
by = Counter(); types = Counter()
for x in range(1 << N):
    if a[x]:
        edges = [prs[t] for t in range(N) if x >> t & 1]
        deg = len(edges); span = len(set(v for e in edges for v in e))
        by[(deg, span)] += 1; types[(deg,) + graph_type(edges, k)] += 1
print("k = %d: ANF of c_odd mod 2 has %d monomials; by (degree, vertices spanned): %s" % (k, sum(by.values()), dict(sorted(by.items()))))
for key in sorted(types):
    print("   degree %d, spans %d vertices, type %s : %d monomials" % (key[0], key[1], key[2], types[key]))
