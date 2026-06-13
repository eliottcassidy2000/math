#!/usr/bin/env python3
"""
tournament_graph_bijection_s319.py -- opus-2026-03-24-S319

THE TOURNAMENT-GRAPH BIJECTION VIA TILINGS

KEY INSIGHT: A tournament tiling IS a graph.

Fix the base path P = (n-1, n-2, ..., 1, 0).
A tournament T with HP P has:
  - n-1 BASE arcs: all "downhill" (i+1 → i for consecutive vertices)
  - m = C(n-1,2) TILE arcs: each either "downhill" (natural) or "uphill" (flipped)

REINTERPRETATION: view "downhill" = edge PRESENT, "uphill" = edge ABSENT.
Then:
  - The n-1 base arcs become n-1 edges ALWAYS present (= the HP as a path graph)
  - Tile bit 0 (natural/downhill) = edge PRESENT
  - Tile bit 1 (flipped/uphill) = edge ABSENT

So: tiling bits → graph on n vertices, where:
  G(tiling) = path P ∪ {tile edges with bit 0}

This graph ALWAYS contains the Hamiltonian path P.
The m tile bits determine which of the C(n-1,2) non-path edges are present.

TOTAL GRAPHS containing path P: exactly 2^m (one per tiling).
TOTAL GRAPHS on n vertices: 2^{C(n,2)} (one per edge subset of K_n).

Since C(n,2) = (n-1) + C(n-1,2) = (n-1) + m:
  2^{C(n,2)} = 2^{n-1} × 2^m

So: 2^{n-1} copies of the tiling space tile the full graph space.
Each copy fixes a different HP as the "base path."

THE BIJECTION:
  Tournaments ↔ Graphs-with-marked-HP (modulo choosing which HP to mark)

But more precisely:
  {tilings} = {graphs containing path P₀}

And by Rédei: every TOURNAMENT has at least one HP.
So: every tournament, when viewed through ANY of its HPs, gives a graph
that contains that HP. Different HPs give different graphs.

The number of (tournament, HP) pairs = Σ_T H(T) = the total HP count.
The number of (graph, marked HP) pairs = Σ_G #{HPs in G}.

These are NOT equal in general. But the TILING space gives
a BIJECTION between tournaments-with-fixed-HP and graphs-containing-fixed-HP.

Author: opus-2026-03-24-S319
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
banner("THE TILING-GRAPH MAP")
# ============================================================

n = 5
m = comb(n-1, 2)  # 6 tiles
base_set = set((i, i+1) for i in range(n-1))
all_pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
tile_pairs = [(i,j) for i,j in all_pairs if (i,j) not in base_set]
tile_pairs.sort()

print("  n=%d: %d total pairs, %d base-path edges, %d tile positions\n" %
      (n, comb(n,2), n-1, m))
print("  Base path edges: %s" % sorted(base_set))
print("  Tile positions: %s" % tile_pairs)
print()

# For each tiling: construct the corresponding GRAPH
# bit=0 (natural direction, "downhill") → edge PRESENT
# bit=1 (flipped, "uphill") → edge ABSENT

def tiling_to_graph(bits, n, tile_pairs, base_set):
    """Convert tiling bits to a graph adjacency matrix.
    Base-path edges always present. Tile bit 0 = present, 1 = absent."""
    adj = [[0]*n for _ in range(n)]
    # Base path edges (always present)
    for i, j in base_set:
        adj[i][j] = 1; adj[j][i] = 1
    # Tile edges (present iff bit = 0)
    for k, (i, j) in enumerate(tile_pairs):
        if not ((bits >> k) & 1):  # bit 0 = present
            adj[i][j] = 1; adj[j][i] = 1
    return adj

# What graphs does this produce?
print("  TILING → GRAPH MAP (first 16 tilings):\n")
print("  %-8s %-10s %-15s %-8s %-10s" %
      ("tiling", "#edges", "degree seq", "has HP?", "even?"))
print("  " + "-"*55)

graph_types = Counter()
even_count = 0; has_hp_count = 0

for bits in range(min(2**m, 64)):
    adj = tiling_to_graph(bits, n, tile_pairs, base_set)
    n_edges = sum(adj[i][j] for i in range(n) for j in range(i+1, n))
    deg_seq = tuple(sorted(sum(adj[v][w] for w in range(n) if w != v) for v in range(n)))

    # Check: does this graph contain the base path as HP?
    # By construction: YES (the base path edges are always present)
    has_hp = True  # Always true by construction
    has_hp_count += 1

    # Is it an even graph?
    is_even = all(d % 2 == 0 for d in deg_seq)
    if is_even: even_count += 1

    graph_types[deg_seq] += 1

    if bits < 16:
        print("  %-8s %-10d %-15s %-8s %-10s" %
              (format(bits, '06b'), n_edges, deg_seq, "✓", "E" if is_even else ""))

print("\n  Total: %d graphs, all contain HP P₀ ✓" % min(2**m, 64))
print("  Even graphs: %d / %d = %.1f%%" %
      (even_count, min(2**m, 64), 100*even_count/min(2**m, 64)))

# ============================================================
banner("THE EDGE COUNT DISTRIBUTION")
# ============================================================

# A graph from tiling(bits) has:
#   n_edges = (n-1) + #{tile bits that are 0}
#           = (n-1) + (m - hamming_weight(bits))
#           = (n-1) + m - |bits| = C(n,2) - |bits|

print("  EDGE COUNT = C(%d,2) - hamming_weight = %d - weight\n" % (n, comb(n,2)))

edge_dist = Counter()
for bits in range(2**m):
    wt = bin(bits).count('1')
    n_edges = comb(n, 2) - wt
    edge_dist[n_edges] += 1

for ne in sorted(edge_dist.keys()):
    bar = "█" * (edge_dist[ne] // 2)
    print("  %2d edges: %3d graphs  %s" % (ne, edge_dist[ne], bar))

print("\n  Edge range: [%d, %d] = [C(n,2)-m, C(n,2)] = [%d, %d]" %
      (comb(n,2)-m, comb(n,2), n-1, comb(n,2)))
print("  Minimum edges = n-1 = path P₀ only (all tiles absent)")
print("  Maximum edges = C(n,2) = complete graph K_n (all tiles present)")

# ============================================================
banner("EVEN GRAPHS IN THE TILING SPACE")
# ============================================================

# When is the graph from tiling(bits) an EVEN graph?
# Degree of vertex v = #{base edges incident to v} + #{tile edges incident to v with bit 0}
# Base edges incident to v: 1 if v ∈ {0, n-1}, else 2 (endpoints have 1 base neighbor)
# Wait: base path is (n-1, n-2, ..., 1, 0). Vertex 0 has base neighbor 1.
# Vertex n-1 has base neighbor n-2. All others have 2 base neighbors.

# Base degree:
base_deg = [0]*n
for i, j in base_set:
    base_deg[i] += 1; base_deg[j] += 1
print("  BASE DEGREES: %s (endpoints have 1, interior have 2)\n" % base_deg)

# Tile degree of v from tiling bits:
# tile_deg(v, bits) = #{tile pairs (i,j) involving v with bit=0}
# = #{tiles involving v} - #{tiles involving v with bit=1}

# For the graph to be EVEN: base_deg(v) + tile_deg(v, bits) must be even for all v.
# Since base_deg is fixed: tile_deg(v, bits) must have same parity as base_deg(v).
# For interior vertices (base_deg=2): tile_deg must be even.
# For endpoint vertices (base_deg=1): tile_deg must be odd.

print("  EVEN GRAPH CONDITION:")
print("  For each vertex v:")
print("    tile_deg(v) ≡ base_deg(v) mod 2")
print("    Interior vertices: tile_deg even")
print("    Endpoint vertices (0 and %d): tile_deg odd\n" % (n-1))

# Count even graphs
even_count_exact = 0
for bits in range(2**m):
    adj = tiling_to_graph(bits, n, tile_pairs, base_set)
    is_even = all(sum(adj[v][w] for w in range(n) if w != v) % 2 == 0 for v in range(n))
    if is_even:
        even_count_exact += 1

print("  Even graphs in tiling space: %d / %d = %.1f%%" %
      (even_count_exact, 2**m, 100*even_count_exact/2**m))
print("  Expected from cycle-space theory: 2^{C(n,2) - n + 1} / 2^{n-1}")
print("  = 2^{m - n + 1} / 1 = 2^{%d} = %d" % (m - n + 1, 2**(m-n+1) if m >= n-1 else 0))
print("  Wait: even graphs in K_n = 2^{C(n,2)-n+1}. In the tiling subspace (fixed HP):")
print("  The even-graph constraint removes n-1 degrees of freedom? No...")

# The constraint is: for each vertex, tile_deg ≡ base_deg mod 2.
# There are n constraints (one per vertex) but the sum is automatic (sum of degrees = 2|E|).
# So n-1 independent parity constraints on m bits.
# Number of solutions: 2^{m-(n-1)} = 2^{m-n+1} = 2^{C(n-1,2)-(n-1)} = 2^{C(n-2,2)}.

print("  Independent parity constraints: n-1 = %d" % (n-1))
print("  Solutions: 2^{m-(n-1)} = 2^{%d} = %d" % (m-n+1, 2**(m-n+1)))
print("  Actual: %d  Match: %s" % (even_count_exact, even_count_exact == 2**(m-n+1)))

# ============================================================
banner("THE BIJECTION: TOURNAMENT CLASSES ↔ EVEN GRAPH CLASSES")
# ============================================================

# In the tiling space:
# - Total tilings: 2^m
# - Tilings giving even graphs: 2^{m-n+1} = 2^{C(n-2,2)}
# - Tilings giving non-even graphs: 2^m - 2^{m-n+1}

# The EVEN GRAPH SUBSPACE of the tiling hypercube is a SUBCODE:
# it's the set of binary words satisfying n-1 linear parity constraints.
# This subcode has dimension m - (n-1) = C(n-2,2).

# The subcode IS a sub-hypercube Q_{m-n+1}.
# Its distance structure is inherited from Q_m.

# KEY: the tournament isomorphism classes partition Q_m.
# The even graph subspace cuts across these classes.
# Some tournament classes contain many even-graph tilings,
# others contain few or none.

print("  THE EVEN-GRAPH SUBCODE:\n")
print("  Full tiling space: Q_%d (%d tilings)" % (m, 2**m))
print("  Even-graph subcode: Q_%d (%d tilings)" % (m-n+1, 2**(m-n+1)))
print("  Dimension reduction: %d → %d (lost %d = n-1 dimensions)\n" % (m, m-n+1, n-1))

# How many tournament classes intersect the even subspace?
# For each tournament class: count how many of its tilings are even graphs.

def canon_tournament(A, nn):
    sc = [sum(A[i]) for i in range(nn)]
    sg = defaultdict(list)
    for v in range(nn): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs2):
        if not gs2: yield []; return
        for p in permutations(gs2[0]):
            for r in gp(gs2[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(nn) for j in range(nn) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, nn):
    dp = {}
    for v in range(nn): dp[(1<<v,v)] = 1
    full = (1<<nn)-1
    for S in range(1,1<<nn):
        for v in range(nn):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(nn):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(nn))

# Build tournament classes
tourn_cmap = {}; tourn_members = defaultdict(list)
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (bits >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    c = canon_tournament(A, n)
    if c not in tourn_cmap: tourn_cmap[c] = len(tourn_cmap)
    tourn_members[tourn_cmap[c]].append(bits)

nc = len(tourn_cmap)

# For each class: how many tilings are even graphs?
print("  TOURNAMENT CLASSES AND THEIR EVEN-GRAPH CONTENT:\n")
print("  %-6s %-4s %-6s %-8s %-8s %-12s" %
      ("class", "H", "fiber", "#even", "%even", "even_deg_seqs"))
print("  " + "-"*50)

for cid in sorted(range(nc)):
    tilings = tourn_members[cid]
    fiber = len(tilings)
    # Get H
    bits0 = tilings[0]
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (bits0 >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    H = Hcount(A, n)

    # Count even-graph tilings
    n_even = 0
    even_degs = set()
    for bits in tilings:
        adj = tiling_to_graph(bits, n, tile_pairs, base_set)
        is_ev = all(sum(adj[v][w] for w in range(n) if w != v) % 2 == 0 for v in range(n))
        if is_ev:
            n_even += 1
            deg = tuple(sorted(sum(adj[v][w] for w in range(n) if w != v) for v in range(n)))
            even_degs.add(deg)

    pct = 100 * n_even / fiber if fiber > 0 else 0
    deg_str = str(even_degs)[:12] if even_degs else "-"
    print("  %-6d %-4d %-6d %-8d %-8.1f %-12s" %
          (cid, H, fiber, n_even, pct, deg_str))

# ============================================================
banner("THE KEY IDENTITY")
# ============================================================

print("""
  THE TOURNAMENT-GRAPH BRIDGE:

  Fix a Hamiltonian path P₀ = (n-1, n-2, ..., 1, 0).

  TILING SPACE Q_m has:
    2^m tilings = 2^{C(n-1,2)} binary words

  Each tiling is SIMULTANEOUSLY:
  1. A TOURNAMENT (orientation of K_n with P₀ as HP)
  2. A GRAPH (edge subset of K_n containing P₀)

  The MAP: bit 0 (natural arc direction) = edge PRESENT
           bit 1 (flipped arc direction) = edge ABSENT

  EVEN GRAPH SUBSPACE:
  The n-1 parity constraints (one per vertex, one redundant) define
  a subcode of dimension m - (n-1) = C(n-2,2).
  This subcode has 2^{C(n-2,2)} elements.

  THE "TRIVIAL VERTEX":
  The tournament has n vertices but only m = C(n-1,2) tile bits.
  A graph on n vertices has C(n,2) edge decisions.
  The difference: C(n,2) - m = n-1 = the base-path edges.

  These n-1 edges are ALWAYS present in the graph interpretation.
  They form the HP P₀. The HP is "free" — it costs no bits.

  So: the tournament effectively operates on n-1 "free" dimensions
  (the tiles) while the graph needs all C(n,2) dimensions.
  The n-1 base-path edges are the "trivial" layer.

  In the Lie algebra language:
  - Tile bits = compound roots of A_{n-1} (the free part)
  - Base-path bits = simple roots (the fixed/trivial part)
  - The tournament "quotients out" the simple roots

  The "trivial vertex" is the one that contributes the LAST row
  of the staircase — vertex 0 (or n-1). Deleting it gives Mode A:
  n → n-1, removing one layer of the hypotenuse.

  THE REDÉI CONNECTION:
  Every tournament has ≥ 1 HP (Rédei).
  In the graph interpretation: every tournament tiling IS a graph
  containing the base path. But the CONVERSE is also true:
  every graph containing a specific HP can be viewed as a tournament
  tiling (with that HP as base path).

  So: {graphs with a marked HP} ↔ {tournament tilings}
  This is a BIJECTION (not just a surjection).

  The number of graphs containing HP P₀ = 2^m = 2^{C(n-1,2)}.
  Each such graph corresponds to exactly one tiling.
  And each tiling corresponds to exactly one labeled tournament with HP P₀.
""")

print("Done. opus-2026-03-24-S319")
