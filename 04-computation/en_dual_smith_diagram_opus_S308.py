#!/usr/bin/env python3
"""THE E_n DUAL SMITH DIAGRAM (opus-S308, HYP-6870).

The cycle-space bijection sends each tiling to an even graph (XOR of the
fundamental cycles of its bit-1 tiles over the base path). E_n = even-graph
iso classes (A002854: 2, 3, 7, 16, 54 for n = 3..7). The dual network: nodes =
E_n classes, edges = wiggly tile-flips (conductance = multiplicity), SOURCE =
the empty even graph's class (fiber 1 -- the dual of the transitive), SINK
rail = the maximal even graphs (K_n at odd n; K_n minus structure at even n),
axis = edge count |E|. Exact Fractions throughout (<= 54 nodes).

Outputs: R^E_n, the |E|-axis concordance (exact), spread-vs-gap, and the
CUT/CYCLE RECIPROCITY PROBE: R^G_n x R^E_n (BSST planar duals satisfy R R' = 1;
the Ihara/Bass cut+cycle duality suggests testing the analogue).
"""
import sys, os, itertools
from collections import defaultdict
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

A002854 = {3: 2, 4: 3, 5: 7, 6: 16, 7: 54}
RG = {4: F(3, 7), 5: F(112209451, 437566662)}   # exact from S307; R6 exact stored in .out
RG_float = {4: 3/7, 5: 0.256440, 6: 0.109041, 7: 0.069902}

def run(n):
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
    m = len(tiles)
    edge_idx = {}
    k = 0
    for u in range(1, n+1):
        for v in range(u+1, n+1):
            edge_idx[(u, v)] = k; k += 1
    E = k
    cyc = []
    for (x, y) in tiles:
        mask = 1 << edge_idx[(y, x)]
        for w in range(y, x):
            mask ^= 1 << edge_idx[(w, w+1)]
        cyc.append(mask)

    def graph_of(bits):
        g = 0
        for i in range(m):
            if bits >> i & 1: g ^= cyc[i]
        return g

    def ginv(g):
        adj = {v: set() for v in range(1, n+1)}
        for (u, v), i in edge_idx.items():
            if g >> i & 1: adj[u].add(v); adj[v].add(u)
        deg = {v: len(adj[v]) for v in adj}
        prof = {}
        for v in adj:
            nd = sorted(deg[u] for u in adj[v])
            tri = sum(1 for u in adj[v] for w in adj[v] if u < w and w in adj[u])
            prof[v] = (deg[v], tuple(nd), tri)
        prof2 = {}
        for v in adj:
            pn = sorted(prof[u] for u in adj[v])
            prof2[v] = (prof[v], tuple(pn))
        # global: sorted edge profile
        ep = sorted(tuple(sorted((prof[u], prof[v]))) for (u, v), i in edge_idx.items() if g >> i & 1)
        return (tuple(sorted(prof2.values())), tuple(ep))

    classes = {}; cls_of = {}; rep = {}
    for b in range(1 << m):
        g = graph_of(b)
        inv = ginv(g)
        if inv not in classes:
            classes[inv] = len(classes); rep[classes[inv]] = g
        cls_of[b] = classes[inv]
    C = len(classes)
    assert C == A002854[n], (n, C)
    ecount = {c: bin(rep[c]).count('1') for c in range(C)}
    # network
    W = defaultdict(int)
    for b in range(1 << m):
        ca = cls_of[b]
        for i in range(m):
            cb = cls_of[b ^ (1 << i)]
            key = (ca, cb) if ca <= cb else (cb, ca)
            W[key] += 1
    W = {kk: v // 2 for kk, v in W.items()}
    cond = {kk: v for kk, v in W.items() if kk[0] != kk[1]}
    src = cls_of[0]                       # empty even graph
    emax = max(ecount.values())
    sinks = [c for c in range(C) if ecount[c] == emax]
    nodes = [c for c in range(C) if c not in sinks]
    pos = {c: i for i, c in enumerate(nodes)}
    N = len(nodes)
    A = [[F(0)]*N for _ in range(N)]
    rhs = [F(0)]*N
    deg = defaultdict(lambda: F(0))
    for (a, bb), w in cond.items():
        deg[a] += w; deg[bb] += w
        if a in pos and bb in pos:
            A[pos[a]][pos[bb]] -= w; A[pos[bb]][pos[a]] -= w
    for c in nodes: A[pos[c]][pos[c]] += deg[c]
    rhs[pos[src]] = F(1)
    for i in range(N):
        piv = next(r for r in range(i, N) if A[r][i] != 0)
        A[i], A[piv] = A[piv], A[i]; rhs[i], rhs[piv] = rhs[piv], rhs[i]
        inv_ = F(1)/A[i][i]
        for r in range(i+1, N):
            if A[r][i] != 0:
                f = A[r][i]*inv_
                for cc in range(i, N): A[r][cc] -= f*A[i][cc]
                rhs[r] -= f*rhs[i]
    phi = [F(0)]*N
    for i in range(N-1, -1, -1):
        ssum = rhs[i]
        for cc in range(i+1, N): ssum -= A[i][cc]*phi[cc]
        phi[i] = ssum/A[i][i]
    pot = {c: phi[pos[c]] for c in nodes}
    for c in sinks: pot[c] = F(0)
    Reff = pot[src]
    # concordance vs the |E| axis (note: source has MIN |E|, so potential should DECREASE in...
    # orient: axis a(c) = emax - ecount[c] so source = max, sinks = 0, like the G-side)
    conc = disc = 0; discs = []
    for a, bb in itertools.combinations(range(C), 2):
        da = ecount[bb] - ecount[a]     # a-axis = -|E|
        if da == 0: continue
        dp = pot[a] - pot[bb]
        if dp == 0: continue
        if (da > 0) == (dp > 0): conc += 1
        else: disc += 1; discs.append((a, bb, ecount[a], ecount[bb]))
    # spread vs gap by level
    lv = defaultdict(list)
    for c in range(C): lv[ecount[c]].append(pot[c])
    levels = sorted(lv)
    overlaps = 0
    for l1, l2 in zip(levels, levels[1:]):
        if max(lv[l2]) > min(lv[l1]) and min(lv[l2]) < max(lv[l1]):
            pass
        if max(lv[l2]) > min(lv[l1]):   # higher-|E| level should have LOWER potential
            overlaps += 1 if min(lv[l2]) < max(lv[l1]) and max(lv[l2]) > min(lv[l1]) else 0
    print(f"===== E_{n}: {C} classes (A002854 OK); source = empty class (fiber 1), "
          f"sink rail = {len(sinks)} class(es) at |E| = {emax}")
    print(f"  DUAL RESISTANCE R^E_{n} = {Reff} = {float(Reff):.6f}")
    print(f"  |E|-axis concordance: {conc}/{conc+disc}" +
          (f"  DISCORDANT: {discs[:6]}" if disc else "  (perfect)"))
    if n in RG_float:
        print(f"  reciprocity probe: R^G_{n} x R^E_{n} = {RG_float[n]*float(Reff):.6f}")
    return Reff

if __name__ == '__main__':
    for n in (4, 5, 6, 7):
        run(n)
