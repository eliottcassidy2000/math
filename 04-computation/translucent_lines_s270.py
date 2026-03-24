#!/usr/bin/env python3
"""
translucent_lines_s270.py -- opus-2026-03-23-S270

TRANSLUCENT LINES: THE INTRA-CLASS TILING CONNECTIONS

Each node in the merged meta-graph G_n/Z_2 is an isomorphism class
containing many labeled tournaments (= tilings of the staircase).

Each tiling has m = C(n,2) cells. Flipping one cell either:
  (a) STAYS in the same class → TRANSLUCENT line (self-loop arc)
  (b) MOVES to a different class → OPAQUE line (meta-graph edge)

So each tiling has:
  - neutral_arcs TRANSLUCENT connections (to other tilings in same class)
  - (m - neutral_arcs) OPAQUE connections (to tilings in other classes)

The TRANSLUCENT SUBGRAPH within a class = the graph on labeled
tournaments where two are connected if they differ by one arc
AND are isomorphic. This is the FIBER GRAPH of the iso class.

THE KEY: each tiling connects to EXACTLY m other tilings (one per cell).
Of these m neighbors:
  - neutral_arcs are in the SAME class (translucent)
  - m - neutral_arcs are in DIFFERENT classes (opaque)

But the user says "T(n-2) other tilings" — let me understand this.
Actually, the user means: each tiling connects to m = C(n,2) others
by flipping one tile. But wait — they say T(n-2). Let me check:
at n=5: m = C(5,2) = 10. But T(n-2) = T(3) = number of arc-orbits
at n-2 = 3, which is T_3 = 4. That doesn't match m = 10.

Perhaps the user means something different: flipping one of the
"overlap tiles" (the C(n-2,2) cells in the overlap region of the
recursive decomposition). At n=5: C(3,2) = 3 overlap cells.

OR: the user means that each tiling connects to 2^{n-2} tilings
by flipping the n-2 boundary cells. That gives 2^3 = 8 at n=5.

Let me just compute the actual fiber graph structure.

Author: opus-2026-03-23-S270
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

for n in [4, 5]:
    banner("TRANSLUCENT LINES AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build iso classes
    cmap = {}; sizes = {}; tc = {}; reps = {}
    class_members = defaultdict(list)  # cid -> list of bits
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
            reps[cid] = [row[:] for row in A]
        cid = cmap[c]
        sizes[cid] += 1; tc[bits] = cid
        class_members[cid].append(bits)

    nc = len(cmap)
    info = {}
    for cid in range(nc):
        A = reps[cid]
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                     'SC': cmap[comp_c] == cid, 'size': sizes[cid],
                     'aut': factorial(n) // sizes[cid]}

    # Merged
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    print("  V_merged = %d, V_unmerged = %d, m = %d\n" % (nm, nc, m))

    # For each merged node: compute the TRANSLUCENT SUBGRAPH
    # = graph on labeled tournaments within the node, connected by
    #   single-arc-flip that preserves the iso class

    for mi in range(nm):
        cids = merged[mi]
        # Collect all tilings in this merged class
        tilings = []
        for cid in cids:
            tilings.extend(class_members[cid])

        fiber_size = len(tilings)
        H_val = info[cids[0]]['H']
        aut_size = info[cids[0]]['aut']

        # Build the translucent subgraph
        tiling_idx = {t: i for i, t in enumerate(tilings)}
        translucent_edges = 0
        opaque_edges = 0
        neutral_per_tiling = []

        for bits in tilings:
            cid_self = tc[bits]
            n_neutral = 0
            for k in range(m):
                flipped = bits ^ (1 << k)
                cid_flip = tc[flipped]
                # Is the flipped tiling in the SAME merged class?
                if mid_map[cid_flip] == mi:
                    n_neutral += 1
                    if flipped in tiling_idx and flipped > bits:
                        translucent_edges += 1
                else:
                    if flipped > bits:
                        opaque_edges += 1

            neutral_per_tiling.append(n_neutral)

        avg_neutral = np.mean(neutral_per_tiling) if neutral_per_tiling else 0
        avg_opaque = m - avg_neutral

        # The translucent graph: V = fiber_size, E = translucent_edges
        # Degree of each tiling = neutral arcs

        if fiber_size <= 120 or n <= 4:
            # Compute connected components of translucent subgraph
            adj_trans = defaultdict(set)
            for bits in tilings:
                for k in range(m):
                    flipped = bits ^ (1 << k)
                    if flipped in tiling_idx and mid_map[tc[flipped]] == mi:
                        adj_trans[bits].add(flipped)
                        adj_trans[flipped].add(bits)

            # BFS for components
            visited = set()
            components = 0
            for start in tilings:
                if start in visited: continue
                components += 1
                stack = [start]
                while stack:
                    v = stack.pop()
                    if v in visited: continue
                    visited.add(v)
                    for w in adj_trans[v]:
                        if w not in visited:
                            stack.append(w)

            # Diameter of translucent subgraph
            from collections import deque
            max_dist = 0
            for start in tilings[:min(20, len(tilings))]:
                dist = {start: 0}
                q = deque([start])
                while q:
                    v = q.popleft()
                    for w in adj_trans[v]:
                        if w not in dist:
                            dist[w] = dist[v] + 1
                            q.append(w)
                            max_dist = max(max_dist, dist[w])

            print("  Merged node %d (H=%d, |Aut|=%d, fiber=%d):" %
                  (mi, H_val, aut_size, fiber_size))
            print("    Translucent edges: %d" % translucent_edges)
            print("    Avg neutral arcs/tiling: %.2f / %d" % (avg_neutral, m))
            print("    Translucent components: %d" % components)
            print("    Translucent diameter: %d" % max_dist)
            print("    Neutral fraction: %.4f" % (avg_neutral / m))
            print()

    # ============================================================
    # SUMMARY STATISTICS
    # ============================================================
    banner("SUMMARY: TRANSLUCENT STRUCTURE AT n=%d" % n)

    # Total translucent edges across all nodes
    total_trans = 0
    total_opaque = 0
    for bits in range(2**m):
        for k in range(m):
            flipped = bits ^ (1 << k)
            if flipped > bits:
                if mid_map[tc[bits]] == mid_map[tc[flipped]]:
                    total_trans += 1
                else:
                    total_opaque += 1

    total_edges = total_trans + total_opaque
    print("  Total labeled edges (one per pair of adjacent tilings): %d" % total_edges)
    print("  Translucent (same class): %d (%.1f%%)" %
          (total_trans, 100*total_trans/total_edges))
    print("  Opaque (different class): %d (%.1f%%)" %
          (total_opaque, 100*total_opaque/total_edges))
    print("  Total = m × 2^{m-1} = %d × %d = %d ✓" %
          (m, 2**(m-1), m * 2**(m-1)))

    # The translucent fraction = Tr(F) / (2^m × m) = arc neutrality
    print("\n  Translucent fraction = %d/%d = %.6f = arc neutrality" %
          (2*total_trans, 2**m * m, 2*total_trans / (2**m * m)))

    # The translucent subgraph of the FULL cube Q_m:
    # This is the graph where vertices = all 2^m tournaments,
    # edges = pairs differing in one bit AND in the same iso class.
    # The connected components of this graph = the "neutral fibers."

    print("\n  THE TWO-LEVEL STRUCTURE:")
    print("  LEVEL 1 (coarse): G_n/Z_2 with %d nodes, connected by opaque lines" % nm)
    print("  LEVEL 2 (fine):   Each node has a translucent subgraph on its tilings")
    print("  The full Q_m = translucent subgraph ∪ opaque subgraph")
    print("  (disjoint edge partition of the m-cube)")

banner("THE TRANSLUCENT GRAPH AS CAYLEY GRAPH OF Aut(T)")

print("""
  For a class C with automorphism group Aut(C):

  The TRANSLUCENT GRAPH within C has:
  - Vertices: n!/|Aut(C)| labeled tournaments
  - Edges: neutral arc flips (stay in class C)

  The structure of this graph is determined by Aut(C):
  - If |Aut| = 1: each tournament has its own "neutral neighborhood"
    with no automorphism-induced equivalences.
    Translucent degree = # neutral arcs per tournament.

  - If |Aut| > 1: automorphisms map neutral arcs to neutral arcs,
    creating orbits within the translucent graph.
    The translucent graph has quotient structure.

  CONJECTURE: The translucent graph of class C is a CAYLEY GRAPH
  of S_n / Aut(C) with generators = the neutral arc types.

  For the TRANSITIVE tournament (|Aut| = 1 at n ≥ 4):
  - Fiber = n! labeled copies
  - Neutral arcs: those whose reversal gives another transitive
  - From S244: at n=5, neutral arcs of transitive = 4 (of 10)

  For the REGULAR tournament (|Aut| = n at odd n):
  - Fiber = n!/n = (n-1)! labeled copies
  - Neutral arcs: often 0 (regular tournaments are "rigid")
""")

# Check: is the translucent graph connected within each class?
# If connected: every pair of tilings in the same class is reachable
# by a sequence of neutral arc flips.

# If NOT connected: some tilings in the same class require "detours"
# through other classes to reach each other via single flips.

print("  TRANSLUCENT CONNECTIVITY:")
for n_val in [4, 5]:
    m_val = comb(n_val, 2)
    pairs_val = [(i,j) for i in range(n_val) for j in range(i+1,n_val)]
    cmap_val = {}; tc_val = {}; sizes_val = {}
    class_mem = defaultdict(list)
    for bits in range(2**m_val):
        A = [[0]*n_val for _ in range(n_val)]
        for k,(ii,jj) in enumerate(pairs_val):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n_val)
        if c not in cmap_val:
            cid = len(cmap_val); cmap_val[c] = cid; sizes_val[cid] = 0
        cid = cmap_val[c]; sizes_val[cid] += 1; tc_val[bits] = cid
        class_mem[cid].append(bits)

    # Check connectivity for each class
    all_connected = True
    for cid in range(len(cmap_val)):
        tilings = class_mem[cid]
        if len(tilings) <= 1: continue
        tset = set(tilings)
        visited = set()
        stack = [tilings[0]]
        while stack:
            v = stack.pop()
            if v in visited: continue
            visited.add(v)
            for k in range(m_val):
                w = v ^ (1 << k)
                if w in tset and w not in visited:
                    stack.append(w)
        if len(visited) < len(tilings):
            aut = factorial(n_val) // sizes_val[cid]
            print("  n=%d: Class %d (size=%d, |Aut|=%d) is DISCONNECTED translucently! "
                  "(%d/%d reachable)" %
                  (n_val, cid, len(tilings), aut, len(visited), len(tilings)))
            all_connected = False

    if all_connected:
        print("  n=%d: ALL classes are translucently CONNECTED!" % n_val)

print("\nDone. opus-2026-03-23-S270")
