#!/usr/bin/env python3
"""Classes behind the blue self-lines (grid-symmetric t with T_{t^c} isomorphic to T_t) at n = 4, 6, 8, 10:
for each self-line class: (H, |Aut|, HC, tilings t, symmetric s, closed c, score sequence, strong), and whether the
closed end of each self-line is the 'closing' or the 'opening' end; the n = 10 count 24 reproduces klein-S173."""
import subprocess, itertools, sys
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
def tiles_of(n): return [(a, b) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
def to_d6(n, arcs):
    A = [[0] * n for _ in range(n)]
    for u, v in arcs: A[u][v] = 1
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6): s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s
def labelg(lines, chunk=250000):
    out = []
    for s0 in range(0, len(lines), chunk):
        part = lines[s0:s0 + chunk]
        r = subprocess.run(['labelg', '-q'], input='\n'.join(part) + '\n', capture_output=True, text=True).stdout.split()
        assert len(r) == len(part); out.extend(r)
    return out
def ham_paths(n, adj):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in adj[v]:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full])
def ham_cycles(n, adj, adjm):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]; dp[1][0] = 1
    for mask in range(1, 1 << n):
        if not mask & 1: continue
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in adj[v]:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(n) if adjm[v][0])
for n in (4, 6, 8, 10):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    mirror = [idx[(n + 1 - b, n + 1 - a)] for a, b in tiles]; kclose = idx[(n, 1)]
    orbits = []; seen = set()
    for i, j in enumerate(mirror):
        if i in seen: continue
        seen |= {i, j}; orbits.append((i, j))
    allmask = (1 << m) - 1
    def arcs_of(mask):
        arcs = [(i, i - 1) for i in range(1, n)]
        for k, (a, b) in enumerate(tiles): arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    syms = []
    for bits in itertools.product([0, 1], repeat=len(orbits)):
        mk = 0
        for b, (i, j) in zip(bits, orbits):
            if b: mk |= (1 << i) | (1 << j)
        syms.append(mk)
    can = labelg([to_d6(n, arcs_of(mk)) for mk in syms]); canc = labelg([to_d6(n, arcs_of(allmask ^ mk)) for mk in syms])
    hits = [(mk, c) for mk, c1, c2 in zip(syms, can, canc) if c1 == c2 for c in [c1]]
    from collections import Counter, defaultdict
    per_class = defaultdict(list)
    for mk, c in hits: per_class[c].append(mk)
    symcount = Counter(can)
    print(f"\nn = {n}: symmetric tilings {len(syms)}, self-line tilings {len(hits)} = {len(hits)//2} blue self-lines in {len(per_class)} classes")
    for c, mks in sorted(per_class.items(), key=lambda kv: -len(kv[1])):
        arcs = arcs_of(mks[0]); adjm = [[0] * n for _ in range(n)]
        for u, v in arcs: adjm[u][v] = 1
        adj = [[v for v in range(n) if adjm[u][v]] for u in range(n)]
        H = ham_paths(n, adj); HC = ham_cycles(n, adj, adjm)
        D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
        A = sum(1 for _ in DiGraphMatcher(D, D).isomorphisms_iter())
        scores = sorted(sum(adjm[u]) for u in range(n))
        t = H // A; closed_hits = sum(1 for mk in mks if (mk >> kclose) & 1)
        print(f"  class: H = {H}, |Aut| = {A}, HC = {HC}, t = {t}, s = {symcount[c]}, closed tilings = {n*HC//A}, self-line tilings {len(mks)} (closed {closed_hits}), scores {scores}")
