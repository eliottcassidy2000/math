#!/usr/bin/env python3
"""
dual_bijection_s321.py -- opus-2026-03-24-S321

THE DUAL BIJECTION: HP AS MISSING EDGES

Two interpretations of a tiling:

INTERPRETATION A (S319): "HP-present"
  Base-path edges: PRESENT
  Tile bit 0 (natural): PRESENT
  Tile bit 1 (flipped): ABSENT
  → Graph G_A has edges at base + natural tiles. Edge count = C(n,2) - weight.

INTERPRETATION B: "HP-absent" (anti-HP)
  Base-path edges: ABSENT
  Tile bit 0 (natural): ABSENT
  Tile bit 1 (flipped): PRESENT
  → Graph G_B has edges ONLY at flipped tiles. Edge count = weight.

G_B = complement of G_A within K_n:
  G_A ∪ G_B = K_n, G_A ∩ G_B = ∅.

INTERPRETATION C: "HP-present, tiles-flipped"
  Base-path edges: PRESENT
  Tile bit 0: ABSENT
  Tile bit 1: PRESENT
  → Graph G_C = base_path ∪ {flipped tiles}. Edge count = (n-1) + weight.

INTERPRETATION D: "HP-absent, tiles-natural"
  Base-path edges: ABSENT
  Tile bit 0: PRESENT
  Tile bit 1: ABSENT
  → Graph G_D = {natural tiles only}. Edge count = m - weight.

So we have FOUR graph interpretations from each tiling!
And: G_A + G_B = K_n, G_C + G_D = K_n \ base_path.

THE KEY QUESTION: For the classes that CAN'T reach PV=0 in Interpretation A,
can they reach PV=0 in Interpretation B, C, or D?

Author: opus-2026-03-24-S321
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

def find_all_HPs(A, n):
    paths = []
    def dfs(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        v = path[-1]
        for w in range(n):
            if w not in visited and A[v][w]:
                dfs(path + [w], visited | {w})
    for s in range(n):
        dfs([s], {s})
    return paths

def tourn_canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs2):
        if not gs2: yield []; return
        for p in permutations(gs2[0]):
            for r in gp(gs2[1:]): yield list(p)+r
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

def is_even_graph(adj, n):
    return all(sum(adj[v][w] for w in range(n) if w != v) % 2 == 0 for v in range(n))

def parity_violations(adj, n):
    return sum(sum(adj[v][w] for w in range(n) if w != v) % 2 for v in range(n))

def hp_to_four_graphs(A, hp, n):
    """Given tournament A and HP, produce all four graph interpretations."""
    pos = {hp[i]: i for i in range(n)}
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Classify each pair
    base_edges = set()
    for k in range(n-1):
        i, j = min(hp[k], hp[k+1]), max(hp[k], hp[k+1])
        base_edges.add((i,j))

    tile_natural = set()  # non-base, arc goes HP-earlier → HP-later
    tile_flipped = set()  # non-base, arc goes HP-later → HP-earlier

    for i, j in all_p:
        if (i,j) in base_edges:
            continue
        # Who is earlier in HP?
        if pos[i] < pos[j]:
            # i earlier. Arc i→j = i beats j = "natural"
            if A[i][j]:
                tile_natural.add((i,j))
            else:
                tile_flipped.add((i,j))
        else:
            # j earlier. Arc j→i = j beats i = "natural"
            if A[j][i]:
                tile_natural.add((i,j))
            else:
                tile_flipped.add((i,j))

    # Four interpretations
    def make_adj(edges):
        adj = [[0]*n for _ in range(n)]
        for i,j in edges:
            adj[i][j] = 1; adj[j][i] = 1
        return adj

    GA = make_adj(base_edges | tile_natural)           # HP-present, natural=present
    GB = make_adj(tile_flipped)                         # HP-absent, flipped=present
    GC = make_adj(base_edges | tile_flipped)           # HP-present, flipped=present
    GD = make_adj(tile_natural)                         # HP-absent, natural=present

    return GA, GB, GC, GD

# ============================================================
banner("THE FOUR GRAPH INTERPRETATIONS AT n=5")
# ============================================================

n = 5
m_arc = comb(n, 2)
all_pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

# Build tournament classes
tourn_classes = defaultdict(list)
for bits in range(2**m_arc):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(all_pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    c = tourn_canon(A, n)
    tourn_classes[c].append([row[:] for row in A])

nc = len(tourn_classes)
print("  n=%d: %d tournament iso classes\n" % (n, nc))

# For each class: check all 4 interpretations across all HPs
print("  %-4s %-6s %-8s %-8s %-8s %-8s %-8s" %
      ("H", "#HPs", "A_pv0", "B_pv0", "C_pv0", "D_pv0", "ANY_pv0"))
print("  " + "-"*55)

total_any = 0
for tc, members in sorted(tourn_classes.items(), key=lambda x: Hcount(x[1][0], n)):
    A0 = members[0]
    H = Hcount(A0, n)
    hps = find_all_HPs(A0, n)

    best = {'A': n+1, 'B': n+1, 'C': n+1, 'D': n+1}
    count_pv0 = {'A': 0, 'B': 0, 'C': 0, 'D': 0}

    for hp in hps:
        GA, GB, GC, GD = hp_to_four_graphs(A0, hp, n)
        for label, G in [('A', GA), ('B', GB), ('C', GC), ('D', GD)]:
            pv = parity_violations(G, n)
            best[label] = min(best[label], pv)
            if pv == 0:
                count_pv0[label] += 1

    any_pv0 = any(count_pv0[x] > 0 for x in 'ABCD')
    if any_pv0: total_any += 1

    print("  %-4d %-6d %-8d %-8d %-8d %-8d %-8s" %
          (H, len(hps), count_pv0['A'], count_pv0['B'],
           count_pv0['C'], count_pv0['D'],
           "✓" if any_pv0 else "✗"))

print("\n  Classes with ANY interpretation giving even graph: %d / %d (%.1f%%)" %
      (total_any, nc, 100*total_any/nc))

# ============================================================
banner("THE UNION: A ∪ B ∪ C ∪ D")
# ============================================================

print("""
  COMBINING ALL FOUR INTERPRETATIONS:

  Interpretation A: HP edges present, natural tiles present
  Interpretation B: HP edges absent, flipped tiles present (= complement of A)
  Interpretation C: HP edges present, flipped tiles present
  Interpretation D: HP edges absent, natural tiles present (= complement of C)

  PARITY STRUCTURE:
  In interpretation A: endpoint degrees = 1 + tile_deg (odd + ?)
  In interpretation B: endpoint degrees = 0 + tile_deg_flipped (0 + ?)
  In interpretation C: endpoint degrees = 1 + tile_deg_flipped (odd + ?)
  In interpretation D: endpoint degrees = 0 + tile_deg (0 + ?)

  For EVEN graph: all degrees even.
  A: endpoints need odd tile_deg → possible but constrained
  B: endpoints need even tile_deg (of flipped tiles) → different constraint!
  C: endpoints need odd tile_deg (of flipped tiles)
  D: endpoints need even tile_deg (of natural tiles)

  B and D remove the base-path parity issue entirely!
  In B: no base edges → endpoint parity = tile_deg_flipped parity
  In D: no base edges → endpoint parity = tile_deg parity

  The constraint for B to be even: every vertex has even "flipped-tile-degree"
  The constraint for D to be even: every vertex has even "natural-tile-degree"

  These are DIFFERENT constraints from A's. And they may be satisfiable
  for classes where A isn't!
""")

# ============================================================
banner("WHICH INTERPRETATION REACHES EACH CLASS?")
# ============================================================

# Detailed breakdown
print("  PER-CLASS BEST PV BY INTERPRETATION:\n")
print("  %-4s %-6s %-6s %-6s %-6s %-6s %-6s" %
      ("H", "#HPs", "A_best", "B_best", "C_best", "D_best", "BEST"))
print("  " + "-"*42)

all_best = Counter()
for tc, members in sorted(tourn_classes.items(), key=lambda x: Hcount(x[1][0], n)):
    A0 = members[0]
    H = Hcount(A0, n)
    hps = find_all_HPs(A0, n)

    best = {'A': n+1, 'B': n+1, 'C': n+1, 'D': n+1}
    for hp in hps:
        GA, GB, GC, GD = hp_to_four_graphs(A0, hp, n)
        for label, G in [('A', GA), ('B', GB), ('C', GC), ('D', GD)]:
            pv = parity_violations(G, n)
            best[label] = min(best[label], pv)

    overall_best = min(best.values())
    best_labels = [k for k, v in best.items() if v == overall_best]
    all_best[overall_best] += 1

    print("  %-4d %-6d %-6d %-6d %-6d %-6d %-6d (%s)" %
          (H, len(hps), best['A'], best['B'], best['C'], best['D'],
           overall_best, "+".join(best_labels)))

print("\n  OVERALL BEST PV DISTRIBUTION:")
for pv in sorted(all_best.keys()):
    print("    PV=%d: %d classes (%.1f%%)" % (pv, all_best[pv], 100*all_best[pv]/nc))

pv0_total = all_best.get(0, 0)
print("\n  TOTAL CLASSES REACHABLE TO EVEN (any interpretation): %d / %d (%.1f%%)" %
      (pv0_total, nc, 100*pv0_total/nc))

# ============================================================
banner("THE RESULT")
# ============================================================

print("""
  COMBINING A + B + C + D:

  With FOUR interpretations and adaptive HP selection:
  Can we reach PV=0 for EVERY tournament class?

  If YES: the union of the four interpretations gives a complete
  bijection from tournaments to even graphs (via choosing the right
  interpretation and HP for each class).

  If NO: some classes are "truly dark" — they can't be represented
  as even graphs in ANY interpretation with ANY HP.
""")

print("Done. opus-2026-03-24-S321")
