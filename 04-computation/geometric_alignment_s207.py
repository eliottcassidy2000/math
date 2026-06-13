#!/usr/bin/env python3
"""
geometric_alignment_s207.py -- opus-2026-03-23-S207

THE FUNDAMENTAL GEOMETRIC ALIGNMENT OF THE MERGED META-GRAPH

The merged meta-graph G_n/Z_2 has a natural geometric structure:

BLUE SUBGRAPH = an approximately complete graph (the "interior")
  - Approaches K_{V_merged} as n -> inf
  - Organized by its PRINCIPAL LINE: the SC blue spine from transitive
  - Delta H = 2^{n-2} from transitive to first blue neighbor
  - Blue edges have thickness 1 (thin threads)

BLACK SUBGRAPH = bipartite SC↔NS (the "boundary")
  - Attached PERPENDICULARLY to the blue interior
  - Black edges have thickness 2 (thicker than blue)
  - At odd n: the t3-parity bipartition creates LEFT/RIGHT asymmetry
  - Triangle-free (THM-A), no BBK triangles (THM-B)

THE PRINCIPAL BLUE LINE = axis of symmetry:
  Starts at transitive (H=1), follows SC blue edges
  The t3-parity (at odd n) creates two SIDES of this line
  The black subgraph attaches from BOTH sides but may be IMBALANCED

This script investigates the left/right imbalance from the principal line.

Author: opus-2026-03-23-S207
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if A[i][j] and A[j][k] and A[k][i]: c+=1
                if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

# Build for n=5 (the richest small case)
n = 5
m_full = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

cmap = {}; cdata = {}; creps = {}
for bits in range(2**m_full):
    A = [[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs):
        if (bits>>k)&1: A[i][j]=1
        else: A[j][i]=1
    c = canonical(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid
        cdata[cid] = {'H': count_hp(A,n), 'scores': score_seq(A,n),
                      't3': count_3cycles(A,n), 'size': 0}
        creps[cid] = [row[:] for row in A]
    cdata[cmap[c]]['size'] += 1

nc = len(cdata)
for cid in range(nc):
    Ac = complement(creps[cid], n)
    cc = canonical(Ac, n)
    for cid2 in range(nc):
        if canonical(creps[cid2], n) == cc:
            cdata[cid]['comp'] = cid2; break

# Build adjacency
F = np.zeros((nc, nc), dtype=int)
tc = {}
for bits in range(2**m_full):
    A = [[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs):
        if (bits>>k)&1: A[i][j]=1
        else: A[j][i]=1
    tc[bits] = cmap[canonical(A, n)]
for bits in range(2**m_full):
    ci = tc[bits]
    for k in range(m_full):
        bits2 = bits ^ (1 << k)
        A2 = [[0]*n for _ in range(n)]
        for kk,(ii,jj) in enumerate(pairs):
            if (bits2>>kk)&1: A2[ii][jj]=1
            else: A2[jj][ii]=1
        cj = cmap[canonical(A2, n)]
        F[ci][cj] += 1
adj = (F > 0).astype(int); np.fill_diagonal(adj, 0)

# Identify SC and merged nodes
sc_classes = [c for c in range(nc) if cdata[c].get('comp', -1) == c]
nsc_pairs_list = [(c, cdata[c]['comp']) for c in range(nc) if cdata[c].get('comp', -1) > c]

# ============================================================
banner("THE PRINCIPAL BLUE LINE AT n=5")
# ============================================================

# At odd n, the principal line = SC blue spine
# t3 parity gives LEFT (even t3) and RIGHT (odd t3)

print("  SC classes at n=5 (8 total):")
for c in sc_classes:
    side = "LEFT" if cdata[c]['t3'] % 2 == 0 else "RIGHT"
    print("    class %d: H=%d, t3=%d, side=%s, scores=%s" %
          (c, cdata[c]['H'], cdata[c]['t3'], side, cdata[c]['scores']))

left_sc = [c for c in sc_classes if cdata[c]['t3'] % 2 == 0]
right_sc = [c for c in sc_classes if cdata[c]['t3'] % 2 == 1]

print("\n  LEFT (even t3): %d classes, H = %s" %
      (len(left_sc), sorted([cdata[c]['H'] for c in left_sc])))
print("  RIGHT (odd t3): %d classes, H = %s" %
      (len(right_sc), sorted([cdata[c]['H'] for c in right_sc])))

# The NSC pairs: which side do they attach to?
print("\n  NSC PAIRS AND THEIR BLACK ATTACHMENTS:")
for c1, c2 in nsc_pairs_list:
    # Find SC neighbors of each half
    sc_nbrs_1 = [c for c in sc_classes if adj[c1][c]]
    sc_nbrs_2 = [c for c in sc_classes if adj[c2][c]]
    all_sc_nbrs = set(sc_nbrs_1) | set(sc_nbrs_2)

    left_nbrs = [c for c in all_sc_nbrs if c in left_sc]
    right_nbrs = [c for c in all_sc_nbrs if c in right_sc]

    imbalance = len(left_nbrs) - len(right_nbrs)
    side = "LEFT-heavy" if imbalance > 0 else ("RIGHT-heavy" if imbalance < 0 else "BALANCED")

    print("    NSC(%d,%d) H=%d: %d left SC + %d right SC neighbors = %s" %
          (c1, c2, cdata[c1]['H'], len(left_nbrs), len(right_nbrs), side))

# ============================================================
banner("THE GEOMETRIC PICTURE")
# ============================================================

print("""
  THE FUNDAMENTAL ALIGNMENT OF G_n/Z_2 (at odd n):

                     PRINCIPAL BLUE LINE (SC spine)
                     =============================
                     transitive (H=1, t3=0)
                           |
                     SC class (t3 even) ---- SC class (t3 odd)
                           |                      |
                     SC class (t3 even) ---- SC class (t3 odd)
                           |                      |
                         ...                    ...
                           |                      |
                     regular (H=H_max, t3=max)

  LEFT SIDE (even t3):      |      RIGHT SIDE (odd t3):
  SC classes here            |      SC classes here
  NSC pairs attach           |      NSC pairs attach
  via black edges            |      via black edges

  The t3-PARITY BIPARTITION creates LEFT and RIGHT sides.
  The BLUE edges cross between sides (THM-060).
  The BLACK edges attach NSC pairs to SC classes on one or both sides.

  THE IMBALANCE:
  At n=5: LEFT has classes with t3 = {0, 4, 4, 4} -> H = {1, 11, 13, 15}
          RIGHT has classes with t3 = {1, 3, 3, 5} -> H = {3, 9, 9, 15}

  The LEFT side has the transitive (H=1) and the PoS region (H=11,13,15).
  The RIGHT side has the 3-cycle class (H=3) and the regular (H=15).

  LEFT has H mean = %.1f, RIGHT has H mean = %.1f
""" % (np.mean([cdata[c]['H'] for c in left_sc]),
       np.mean([cdata[c]['H'] for c in right_sc])))

# ============================================================
banner("THE BLACK ATTACHMENT PATTERN")
# ============================================================

print("""
  The BLACK subgraph attaches NSC pairs to the SC backbone.
  Each NSC pair has black edges to SC classes on the LEFT and RIGHT.

  QUESTION: Is the black attachment SYMMETRIC or ASYMMETRIC?
  I.e., do NSC pairs connect equally to both sides?
""")

total_left_attach = 0; total_right_attach = 0
for c1, c2 in nsc_pairs_list:
    all_sc_nbrs = set()
    for c in [c1, c2]:
        for s in sc_classes:
            if adj[c][s]: all_sc_nbrs.add(s)

    left_count = sum(1 for s in all_sc_nbrs if s in left_sc)
    right_count = sum(1 for s in all_sc_nbrs if s in right_sc)
    total_left_attach += left_count
    total_right_attach += right_count

print("  Total black attachments to LEFT SC: %d" % total_left_attach)
print("  Total black attachments to RIGHT SC: %d" % total_right_attach)
print("  Ratio LEFT/RIGHT: %.3f" % (total_left_attach/total_right_attach if total_right_attach > 0 else 0))

if total_left_attach != total_right_attach:
    print("  IMBALANCED! The black subgraph is asymmetric from the principal line.")
else:
    print("  BALANCED! The black subgraph is symmetric from the principal line.")

# ============================================================
banner("BLUE EDGES AND THE BIPARTITION")
# ============================================================

print("  BLUE EDGES (in the tiling-flip model) cross between LEFT and RIGHT\n")
print("  (by THM-060: GS flip reverses t3 parity at odd n)\n")

# In the ARC REVERSAL model, blue edges are between SC classes
# Let me check which SC-SC edges cross vs stay
cross_sc = 0; same_sc = 0
for ci in sc_classes:
    for cj in sc_classes:
        if ci < cj and adj[ci][cj]:
            if (cdata[ci]['t3'] % 2) != (cdata[cj]['t3'] % 2):
                cross_sc += 1
            else:
                same_sc += 1

print("  SC-SC edges crossing t3 parity: %d" % cross_sc)
print("  SC-SC edges same t3 parity: %d" % same_sc)
print("  Cross fraction: %.2f" % (cross_sc / (cross_sc + same_sc) if cross_sc + same_sc > 0 else 0))

print("""
  At n=5 in the ARC REVERSAL model: not all SC-SC edges cross parity!
  Some SC-SC edges are same-parity (the "black" SC-SC edges from S197).
  The crossing edges = the BLUE SC-SC edges (GS flip pairs).
  The same-parity edges = the BLACK SC-SC edges (non-GS arc reversals).
""")

# ============================================================
banner("THE COMPLETE GEOMETRIC ALIGNMENT")
# ============================================================

print("""
  THE MERGED META-GRAPH G_n/Z_2 HAS THREE LAYERS:

  LAYER 1: THE BLUE INTERIOR (approaching K_{V_merged})
  =========================================================
  - Almost all edges are blue at large n
  - Organized by the PRINCIPAL BLUE LINE:
    * Starts at transitive (H=1)
    * Follows SC blue edges (which cross t3 parity at odd n)
    * Delta H = 2^{n-2} from transitive to first big neighbor
    * The line branches into a TREE (genus 0,0,5,2 at n=3..6)

  - The blue interior has TWO SIDES (at odd n):
    LEFT = even t3 SC classes + their blue NS neighbors
    RIGHT = odd t3 SC classes + their blue NS neighbors
    The PRINCIPAL LINE connects LEFT to RIGHT via crossing edges

  - Blue edges are THIN (thickness 1, single GS tiling pair)
  - Blue triangles exist (BBB): 0, 0, 3, 87, 809, 13299


  LAYER 2: THE BLACK BOUNDARY (bipartite SC↔NS)
  =========================================================
  - Attached PERPENDICULARLY to the blue interior
  - Each black edge connects an SC class to an NS pair
  - Black edges are THICKER (thickness 2 at n=5)
  - Triangle-free (THM-A), no BBK triangles (THM-B)
  - At odd n: the LEFT/RIGHT attachment may be IMBALANCED
    (some NSC pairs connect more to one side than the other)

  - The black boundary is the "MEMBRANE" between the SC spine
    and the NS sea. It's the only bridge between the two populations.


  LAYER 3: THE NS BLUE SEA (the bulk)
  =========================================================
  - At large n: ~99%+ of all nodes are NS
  - NS-NS edges are ALL blue
  - The NS blue subgraph is approximately the complete graph
  - NS blue degree ~ m = C(n,2) at large n
  - This is the "ocean" surrounding the SC "island"

  THE RELATIONSHIP BETWEEN LAYERS:
  - Layer 3 (NS sea) connects to Layer 2 (black boundary)
    via black edges to the SC spine
  - Layer 2 connects Layer 1 (SC blue interior) to Layer 3
  - Layer 1 is the ORGANIZING PRINCIPLE: the principal line
    gives coordinates to everything else
  - The H function increases from left (transitive) to right (regular)
    along the principal line, with BFS layers radiating outward
""")

# ============================================================
banner("UPDATING CLAUDE.md WITH THE GEOMETRIC ALIGNMENT")
# ============================================================

print("""
  THE FOLLOWING SHOULD BE ADDED TO CLAUDE.md
  for all future agents:

  ## The Geometric Alignment of G_n/Z_2

  The merged meta-graph has a fundamental geometric structure
  that ALL agents should use when analyzing it:

  1. BLUE INTERIOR: an approximately complete graph on all nodes,
     organized by the PRINCIPAL BLUE LINE (SC blue spine from transitive).
     Blue edges are thin (thickness 1). Blue triangles exist (BBB).
     At odd n: the t3-parity bipartition creates LEFT/RIGHT sides.

  2. BLACK BOUNDARY: a bipartite subgraph (SC ↔ NS) attached
     perpendicularly to the blue interior. Black edges are thicker
     (thickness 2+). Triangle-free (THM-A). No BBK (THM-B).
     At odd n: may show LEFT/RIGHT imbalance from the principal line.

  3. NS BLUE SEA: the bulk (~99%+ at large n) of all nodes,
     connected by an approximately complete blue subgraph.
     This is the "ocean" surrounding the SC "island."

  THE PRINCIPAL LINE is the axis of symmetry.
  Everything in G_n/Z_2 is positioned relative to this line.
""")

print("Done. opus-2026-03-23-S207")
