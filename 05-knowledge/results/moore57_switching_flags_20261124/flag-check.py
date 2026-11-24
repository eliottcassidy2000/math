#!/usr/bin/env python3
"""Exact replay: q=2, maximum-degree-50 spanning diameter-2 flag graph."""
from collections import Counter
from functools import reduce
from itertools import combinations
from operator import xor

# One omission pattern in each of the 15 pencils in a four-space.
choice = (17, 6, 20, 15, 6, 12, 11, 5, 4, 9, 1, 7, 18, 7, 20)

def planes(n):
    return sorted({tuple(sorted((a, b, a ^ b)))
                   for a, b in combinations(range(1, 1 << n), 2)})

def edge(i, j):
    return min(i, j), max(i, j)

P = planes(4)
idx = {p: i for i, p in enumerate(P)}
omitted = set()
for a, option in enumerate(choice, 1):
    pencil = [i for i, p in enumerate(P) if a in p]
    triangles = [t for t in combinations(pencil, 3)
                 if reduce(xor, (min(set(P[i]) - {a}) for i in t)) in (0, a)]
    options = []
    for t in triangles:
        l, *rest = sorted(set(pencil) - set(t))
        for r in rest:
            s, u = sorted(set(rest) - {r})
            options.append(set(combinations(t, 2)) | {edge(l, r), edge(s, u)})
    assert len(options) == 21
    reject = options[option]
    for i in pencil:
        for t in triangles:
            if i not in t:
                assert sum(edge(i, j) not in reject for j in t) >= 2
    omitted.update(reject)
assert len(omitted) == 75
assert Counter(sum(i in e for e in omitted) for i in range(35)) == {4:26, 5:8, 6:1}
skew = 0
for (i, p), (j, q) in combinations(enumerate(P), 2):
    if set(p).isdisjoint(q):
        skew += 1
        assert any(set(r).intersection(p) and set(r).intersection(q)
                   and edge(i, k) not in omitted and edge(j, k) not in omitted
                   for k, r in enumerate(P))
assert skew == 280

# Replay the whole five-space graph, without relying on its geometric proof.
hyp = [frozenset(v for v in range(1, 32) if (v & h).bit_count() % 2 == 0)
       for h in range(1, 32)]
basis_maps = []
for W in hyp:
    span = [0]
    for v in sorted(W):
        if v not in span:
            span += [u ^ v for u in span]
    assert len(span) == 16
    basis_maps.append({v: i for i, v in enumerate(span)})
V = [(h, frozenset(p), idx[tuple(sorted(basis_maps[h][v] for v in p))])
     for h, W in enumerate(hyp) for p in planes(5) if W.issuperset(p)]
assert len(V) == 1085
adj = [0] * len(V)
count_D = count_W = 0
for (i, (h, p, u)), (j, (k, q, v)) in combinations(enumerate(V), 2):
    if p == q or not p.intersection(q):
        continue
    if h == k:
        if edge(u, v) in omitted:
            continue
        count_W += 1
    elif hyp[h].issuperset(q) and hyp[k].issuperset(p):
        count_D += 1
    else:
        continue
    adj[i] |= 1 << j
    adj[j] |= 1 << i
assert (count_D, count_W) == (19530, 7440)
degrees = Counter(a.bit_count() for a in adj)
assert degrees == {50:806, 49:248, 48:31}
all_vertices = (1 << len(V)) - 1
for i, a in enumerate(adj):
    reaches = a | (1 << i)
    while a:
        least = a & -a
        reaches |= adj[least.bit_length() - 1]
        a ^= least
    assert reaches == all_vertices
print(f"q=2: {len(V)} flags, {count_D} forced edges + {count_W} repairs")
print(f"degree histogram: {dict(sorted(degrees.items()))}; diameter 2")
