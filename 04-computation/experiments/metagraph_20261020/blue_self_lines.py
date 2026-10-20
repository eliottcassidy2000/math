#!/usr/bin/env python3
"""Blue self-lines: symmetric tilings t with T_{t^c} isomorphic to T_t (the complement line stays inside the node).
Counts 1, 0, 2, 0, 4 for n = 4..8 suggest 2^{n/2-2} for even n and 0 for odd n. List them for n = 4, 6, 8 with the
class data (H, |Aut|, score sequence) and the tiling bit pattern on the pin grid."""
import subprocess, networkx as nx
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
def labelg(lines):
    r = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(r) == len(lines); return r
def ham(n, adj):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            for w in range(n):
                if not (mask >> w) & 1 and adj[v][w]: dp[mask | 1 << w][w] += c
    return sum(dp[full])
for n in (4, 5, 6, 7, 8):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    mirror = [idx[(n + 1 - b, n + 1 - a)] for a, b in tiles]
    orbits = []; seen = set()
    for i, j in enumerate(mirror):
        if i in seen: continue
        seen |= {i, j}; orbits.append((i, j))
    # enumerate symmetric tilings only
    import itertools
    syms = []
    for bits in itertools.product([0, 1], repeat=len(orbits)):
        mk = 0
        for b, (i, j) in zip(bits, orbits):
            if b: mk |= (1 << i) | (1 << j)
        syms.append(mk)
    allmask = (1 << m) - 1
    def arcs_of(mask):
        arcs = [(i, i - 1) for i in range(1, n)]
        for k, (a, b) in enumerate(tiles): arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    can = labelg([to_d6(n, arcs_of(mk)) for mk in syms]); canc = labelg([to_d6(n, arcs_of(allmask ^ mk)) for mk in syms])
    hits = [mk for mk, c1, c2 in zip(syms, can, canc) if c1 == c2]
    print(f"n = {n}: symmetric tilings {len(syms)}, blue self-lines (symmetric t with [T_t] = [T_t^c]): {len(hits) // 2} (tilings {len(hits)})")
    for mk in hits:
        arcs = arcs_of(mk); adj = [[0] * n for _ in range(n)]
        for u, v in arcs: adj[u][v] = 1
        D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
        aut = sum(1 for _ in DiGraphMatcher(D, D).isomorphisms_iter())
        scores = sorted(sum(r) for r in adj)
        grid = {}
        for k, (a, b) in enumerate(tiles): grid[(a - b - 1, b)] = (mk >> k) & 1
        rows = []
        for r in range(1, n - 1):
            rows.append(''.join(str(grid[(r, c)]) for c in range(1, n - r)))
        print(f"   tiling {mk:0{m}b}: H = {ham(n, adj)}, |Aut| = {aut}, scores {scores}, grid rows r=1..: {rows}")
