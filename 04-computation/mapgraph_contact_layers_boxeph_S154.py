#!/usr/bin/env python3
"""
mapgraph_contact_layers_boxeph_S154.py  (HYP-8235)

OWNER S154: map graphs (faces adjacent iff they share ANY boundary point, not
just an edge — Chen-Grigni-Papadimitriou; map graphs = half-squares of
incidence bipartite graphs) — extend abstractly to the repo.

TRANSLATION.  The tiling hypercube Q_m is the surface; the iso classes are the
FACES (fibers of the tiling fibration); the d=1 (wiggly) metagraph is the DUAL
GRAPH (edge contact); the d<=2 layer is the MAP-GRAPH closure (corner
contact).  The CGP half-square theorem suggests the identification test:

  HS(B) := half-square of the bridge incidence  B[tournament-class, even-class]
           (two tournament classes adjacent iff some EVEN class occurs in both
            fibers)  --  "faces meeting at an even-graph vertex."

MEASURED HERE (n = 5, 6, exact):
  (1) contact layers d<=k of the merged... (unmerged) class metagraph;
  (2) HS(B) vs each layer: containments, equalities, differences;
  (3) midpoint analysis of pure-2 contacts (d=2 but not d=1): do the two
      Q_m-midpoints of a connecting 4-cycle share an even class? tournament
      class?  (the 'corner' the two faces meet at);
  (4) the even-side mirror: even-class contact layers + HS(B^T).

Conventions: vertices 1..n; base path arcs k+1 -> k fixed; tiles (x,y), x-y>=2,
bit=1 flips to y -> x.  Even graph of a tiling = XOR over FLIPPED tiles of the
fundamental cycles {xy} u {(k,k+1): y<=k<x} (undirected).  boxeph-2026-07-20-S154.
"""

from itertools import permutations, combinations

def build(n):
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    m = len(tiles)
    perms = list(permutations(range(1, n + 1)))
    # tournament of a tiling: adj[i][j]=1 iff i->j ; vertices 1..n
    def tourn(bits):
        adj = {}
        for k in range(n, 1, -1): adj[(k, k - 1)] = 1   # base arcs k -> k-1
        for idx, (x, y) in enumerate(tiles):
            if (bits >> idx) & 1: adj[(y, x)] = 1
            else: adj[(x, y)] = 1
        return adj
    def tkey(adj):
        return frozenset(adj)
    def tcanon(adj):
        best = None
        for p in perms:
            pm = {i + 1: p[i] for i in range(n)}
            k = frozenset((pm[a], pm[b]) for (a, b) in adj)
            if best is None or sorted(k) < sorted(best): best = k
        return frozenset(best)
    def evengraph(bits):
        E = set()
        for idx, (x, y) in enumerate(tiles):
            if (bits >> idx) & 1:
                cyc = {frozenset((x, y))} | {frozenset((k, k + 1)) for k in range(y, x)}
                E ^= cyc
        return frozenset(E)
    def ecanon(E):
        best = None
        for p in perms:
            pm = {i + 1: p[i] for i in range(n)}
            k = sorted(tuple(sorted((pm[a], pm[b]))) for fs in E for (a, b) in [tuple(fs)])
            if best is None or k < best: best = k
        return tuple(best)
    tclass, eclass = {}, {}
    tmap, emap = {}, {}
    for bits in range(1 << m):
        tc = tcanon(tourn(bits))
        ec = ecanon(evengraph(bits))
        tclass.setdefault(tc, len(tclass))
        eclass.setdefault(ec, len(eclass))
        tmap[bits] = tclass[tc]
        emap[bits] = eclass[ec]
    return m, tmap, emap, len(tclass), len(eclass)

for n in (5, 6):
    m, tmap, emap, TC, EC = build(n)
    print("=" * 88)
    print("n=%d: m=%d, tilings=%d, tournament classes=%d (A000568), even classes=%d (A002854)"
          % (n, m, 1 << m, TC, EC))
    # contact layers
    def layer_edges(d):
        E = set()
        for b in range(1 << m):
            cb = tmap[b]
            for flip in combinations(range(m), d):
                b2 = b
                for f in flip: b2 ^= (1 << f)
                if b2 < b: continue
                c2 = tmap[b2]
                if c2 != cb: E.add((min(cb, c2), max(cb, c2)))
        return E
    L1 = layer_edges(1); L2 = layer_edges(2); L3 = layer_edges(3)
    C12 = L1 | L2; C123 = C12 | L3
    allp = set((i, j) for i in range(TC) for j in range(i + 1, TC))
    # bridge incidence and half-square
    B = {}
    for b in range(1 << m):
        B.setdefault(tmap[b], set()).add(emap[b])
    HS = set()
    for i in range(TC):
        for j in range(i + 1, TC):
            if B[i] & B[j]: HS.add((i, j))
    print("  dual (d=1): %d edges ; d<=2: %d ; d<=3: %d ; complete: %d"
          % (len(L1), len(C12), len(C123), len(allp)))
    print("  HS(B) ('meet at an even vertex'): %d edges" % len(HS))
    for name, S in (("d=1", L1), ("d<=2", C12), ("d<=3", C123), ("complete", allp)):
        rel = "==" if HS == S else ("SUPERSET of" if HS > S else ("SUBSET of" if HS < S else "INCOMPARABLE to"))
        print("    HS(B) %s %s  (HS-only: %d, %s-only: %d)"
              % (rel, name, len(HS - S), name, len(S - HS)))
    # pure-2 contacts and their corners
    P2 = C12 - L1
    print("  pure-2 contacts (map-graph edges absent from the dual): %d" % len(P2))
    same_e_mid = diff_e_mid = 0
    same_t_mid = 0
    shared_end_e = 0
    for b in range(1 << m):
        cb = tmap[b]
        for (i, j) in combinations(range(m), 2):
            b2 = b ^ (1 << i) ^ (1 << j)
            if b2 < b: continue
            c2 = tmap[b2]
            if c2 == cb: continue
            pr = (min(cb, c2), max(cb, c2))
            if pr not in P2: continue
            m1, m2 = b ^ (1 << i), b ^ (1 << j)
            if emap[m1] == emap[m2]: same_e_mid += 1
            else: diff_e_mid += 1
            if tmap[m1] == tmap[m2]: same_t_mid += 1
            if emap[b] == emap[b2]: shared_end_e += 1
    tot = same_e_mid + diff_e_mid
    print("  4-cycle corner stats over all pure-2 contact instances (%d):" % tot)
    print("    midpoints share EVEN class: %d/%d ; midpoints share TOURNAMENT class: %d/%d"
          % (same_e_mid, tot, same_t_mid, tot))
    print("    endpoints share EVEN class: %d/%d" % (shared_end_e, tot))
    # even-side mirror
    BT = {}
    for b in range(1 << m):
        BT.setdefault(emap[b], set()).add(tmap[b])
    HSe = set()
    for i in range(EC):
        for j in range(i + 1, EC):
            if BT[i] & BT[j]: HSe.add((i, j))
    Ee1 = set()
    for b in range(1 << m):
        for f in range(m):
            b2 = b ^ (1 << f)
            if b2 < b: continue
            if emap[b2] != emap[b]: Ee1.add((min(emap[b], emap[b2]), max(emap[b], emap[b2])))
    alle = set((i, j) for i in range(EC) for j in range(i + 1, EC))
    print("  EVEN side: dual (d=1) %d ; HS(B^T) %d ; complete %d ; HS==dual: %s ; HS==complete: %s"
          % (len(Ee1), len(HSe), len(alle), HSe == Ee1, HSe == alle))
print("DONE.")
