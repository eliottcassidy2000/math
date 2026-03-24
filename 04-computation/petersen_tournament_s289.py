#!/usr/bin/env python3
"""
petersen_tournament_s289.py -- opus-2026-03-24-S289

THE PETERSEN GRAPH AS A TOURNAMENT ON THE META-GRAPH NODES

The n=5 merged meta-graph has 10 nodes (= vertices of the Petersen graph).
The Petersen graph K(5,2) has 15 edges and 30 non-edges.

CONSTRUCTION:
  Orient the 15 Petersen edges in the "transitive direction"
  (from lower H to higher H, or by the natural ordering).
  Orient the 30 non-edges in the OPPOSITE direction.

This creates a tournament T_P on 10 vertices.

QUESTION: What is H(T_P)? What are its properties?
Does this tournament live at a special point in the meta-graph?

Author: opus-2026-03-24-S289
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

# The 10 nodes of the n=5 merged meta-graph, ordered by H
nodes = [
    {'mid': 0, 'H': 1,  'score': (0,1,2,3,4), 'label': 'trans'},
    {'mid': 1, 'H': 3,  'score': (1,1,1,3,4), 'label': 'H3_NS'},
    {'mid': 3, 'H': 3,  'score': (0,2,2,2,4), 'label': 'H3_SC'},
    {'mid': 2, 'H': 5,  'score': (1,1,2,2,4), 'label': 'H5_NS'},
    {'mid': 4, 'H': 9,  'score': (1,1,2,3,3), 'label': 'H9a'},
    {'mid': 5, 'H': 9,  'score': (1,1,2,3,3), 'label': 'H9b'},
    {'mid': 6, 'H': 11, 'score': (1,2,2,2,3), 'label': 'H11'},
    {'mid': 8, 'H': 13, 'score': (1,2,2,2,3), 'label': 'H13'},
    {'mid': 7, 'H': 15, 'score': (1,2,2,2,3), 'label': 'H15a'},
    {'mid': 9, 'H': 15, 'score': (2,2,2,2,2), 'label': 'reg'},
]

# Meta-graph adjacency
meta_adj_raw = {
    (0,1): 1, (0,2): 1, (0,3): 1, (0,4): 1,
    (1,2): 1, (1,5): 1,
    (2,3): 1, (2,4): 1, (2,5): 1, (2,6): 1, (2,8): 1,
    (3,7): 1,
    (4,5): 1, (4,6): 1, (4,7): 1,
    (5,6): 1, (5,8): 1,
    (6,8): 1, (6,9): 1,
    (7,8): 1,
    (8,9): 1,
}

# Build proper adjacency matrix using mids
mid_to_idx = {nodes[i]['mid']: i for i in range(10)}

meta_adj = np.zeros((10, 10), dtype=int)
for (mi, mj), v in meta_adj_raw.items():
    i = mid_to_idx[mi]; j = mid_to_idx[mj]
    meta_adj[i][j] = 1; meta_adj[j][i] = 1

ne_meta = meta_adj.sum() // 2
print("Meta-graph: V=10, E=%d\n" % ne_meta)

# ============================================================
banner("THE PETERSEN GRAPH ON THESE 10 NODES")
# ============================================================

# Petersen K(5,2): 2-subsets of [5], adjacent iff disjoint.
# But we need to map Petersen vertices to meta-graph nodes.
# There's no canonical mapping. Let me try: use the H-ordering.

# The meta-graph has 21 edges and 24 non-edges.
# Petersen has 15 edges and 30 non-edges.
# The meta-graph has MORE edges than Petersen.

# Instead of mapping Petersen to meta-graph nodes, let me just
# BUILD the tournament the user describes:
# The 10 nodes in H-order: indices 0..9 as above.
# We need 15 "Petersen edges" among these 10 nodes.
# The meta-graph has 21 edges — too many for Petersen.
# The meta-graph's NON-edges have 24 — too many for Petersen complement's 30.

# So let me just USE the meta-graph adjacency AS the Petersen-like structure:
# 21 edges oriented "transitive" (low H → high H)
# 24 non-edges oriented "anti-transitive" (high H → low H)

# TOURNAMENT FROM META-GRAPH:
# For each pair (i,j) with i < j (in H-ordering):
# If meta_adj[i][j] = 1: orient i → j (transitive direction, lower H wins)
# If meta_adj[i][j] = 0: orient j → i (anti-transitive, higher H wins)

print("  CONSTRUCTING TOURNAMENT T_meta FROM THE META-GRAPH:\n")
print("  Rule: meta-edge → orient LOW H to HIGH H (transitive)")
print("  Rule: non-edge → orient HIGH H to LOW H (anti-transitive)\n")

T = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        # i has lower or equal H to j (since sorted by H)
        if meta_adj[i][j]:
            # Meta-edge: orient i → j (low H beats high H = transitive)
            T[i][j] = 1
        else:
            # Non-edge: orient j → i (high H beats low H = anti-transitive)
            T[j][i] = 1

# Verify tournament
for i in range(10):
    for j in range(10):
        if i != j:
            assert T[i][j] + T[j][i] == 1, "Not a tournament!"

# Score sequence
scores = [sum(T[i]) for i in range(10)]
score_sorted = sorted(scores)
print("  Score sequence: %s" % score_sorted)

# 3-cycle count
c3 = 0
for i in range(10):
    for j in range(i+1, 10):
        for k in range(j+1, 10):
            if T[i][j] and T[j][k] and T[k][i]: c3 += 1
            if T[i][k] and T[k][j] and T[j][i]: c3 += 1
print("  3-cycles: %d" % c3)

# H count
def H10(A):
    n = 10
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    full = (1<<n) - 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            c = dp.get((S,v), 0)
            if c == 0: continue
            for w in range(n):
                if S & (1<<w): continue
                if A[v][w]:
                    dp[(S|(1<<w), w)] = dp.get((S|(1<<w), w), 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

H_T = H10(T)
print("  H(T_meta) = %d" % H_T)

# Szele bound for n=10: 10!/2^9 = 3628800/512 = 7087.5
szele = factorial(10) // (2**9)
print("  Szele bound: %d" % szele)
print("  H/Szele = %.4f" % (H_T / szele))

# Is it transitive?
print("  Is transitive: %s" % (score_sorted == list(range(10))))

# The adjacency: which edges are "transitive" (meta-edge) vs "anti" (non-edge)?
print("\n  ORIENTATION TABLE:")
print("  %-15s | %-15s | direction | meta-edge?" % ("from (H)", "to (H)"))
print("  " + "-"*60)
for i in range(10):
    for j in range(i+1, 10):
        Hi = nodes[i]['H']; Hj = nodes[j]['H']
        if T[i][j]:
            direction = "→ (low→high)"
            is_meta = meta_adj[i][j]
        else:
            direction = "← (high→low)"
            is_meta = not meta_adj[i][j]
        # Only print first few
        if i < 3 or (i == 0 and j < 5):
            print("  %d(H=%d)%s | %d(H=%d)%s | %s | %s" %
                  (i, Hi, " "*(5-len(str(i))),
                   j, Hj, " "*(5-len(str(j))),
                   direction, "yes" if meta_adj[i][j] else "NO"))

# ============================================================
banner("THE PURE PETERSEN TOURNAMENT")
# ============================================================

# Now construct a tournament DIRECTLY from the Petersen graph.
# Petersen K(5,2): vertices = 2-subsets of {0,1,2,3,4}.
# Edges = disjoint subsets.

subsets = list(combinations(range(5), 2))
print("  Petersen vertices (2-subsets of [5]):")
for i, s in enumerate(subsets):
    print("    %d: %s" % (i, s))

# Build Petersen adjacency
pet = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        if len(set(subsets[i]) & set(subsets[j])) == 0:
            pet[i][j] = 1; pet[j][i] = 1

ne_pet = pet.sum() // 2
print("\n  Petersen: V=10, E=%d" % ne_pet)

# Orient Petersen edges: use lexicographic order as "transitive"
# subset i < subset j in lex order → i beats j on Petersen edges
# non-Petersen edges: j beats i
T_pet = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        if pet[i][j]:
            T_pet[i][j] = 1  # Petersen edge: lex-lower wins
        else:
            T_pet[j][i] = 1  # non-edge: lex-higher wins

scores_pet = sorted([sum(T_pet[i]) for i in range(10)])
c3_pet = 0
for i in range(10):
    for j in range(i+1, 10):
        for k in range(j+1, 10):
            if T_pet[i][j] and T_pet[j][k] and T_pet[k][i]: c3_pet += 1
            if T_pet[i][k] and T_pet[k][j] and T_pet[j][i]: c3_pet += 1

H_pet = H10(T_pet)

print("  Score sequence: %s" % scores_pet)
print("  3-cycles: %d" % c3_pet)
print("  H(T_Petersen) = %d" % H_pet)
print("  H/Szele = %.4f" % (H_pet / szele))

# ============================================================
banner("COMPARING THE TWO TOURNAMENTS")
# ============================================================

print("  T_meta (from meta-graph adjacency):")
print("    Scores: %s" % sorted(scores))
print("    c3: %d" % c3)
print("    H: %d" % H_T)

print("\n  T_Petersen (from Petersen + lex orientation):")
print("    Scores: %s" % scores_pet)
print("    c3: %d" % c3_pet)
print("    H: %d" % H_pet)

# How many arcs differ?
diff = sum(1 for i in range(10) for j in range(10) if i != j and T[i][j] != T_pet[i][j])
print("\n  Arc difference: %d / %d arcs differ" % (diff // 2, 45))

# ============================================================
banner("THE STRUCTURE OF T_meta")
# ============================================================

print("  T_meta is the tournament where:")
print("  - Lower-H classes beat higher-H classes IF they're meta-graph neighbors")
print("  - Higher-H classes beat lower-H IF they're NOT meta-graph neighbors")
print()
print("  This inverts the usual 'transitive = low beats high' for NON-neighbors.")
print("  The non-neighbors are the LONG-RANGE connections (large dH).")
print("  So: short-range = transitive, long-range = anti-transitive.")
print()

# Which orientation gives the MOST Hamiltonian paths?
# Try: pure transitive (all low→high)
T_trans = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        T_trans[i][j] = 1  # always low→high

H_trans = H10(T_trans)
print("  Pure transitive (all low→high): H = %d" % H_trans)
print("  T_meta (edges transitive, non-edges anti): H = %d" % H_T)
print("  T_Petersen (Petersen edges transitive): H = %d" % H_pet)

# How does the meta-graph structure affect H?
# The meta-graph edges (21) go transitive (low H → high H).
# The non-edges (24) go anti-transitive (high H → low H).
# Since non-edges tend to connect DISTANT H classes,
# the anti-transitive edges create long backward arcs.

# These backward arcs create 3-cycles, increasing H beyond 1 (transitive).
print("\n  Score analysis:")
print("  Pure transitive: scores = %s" % sorted([sum(T_trans[i]) for i in range(10)]))
print("  T_meta: scores = %s" % sorted(scores))

print("\nDone. opus-2026-03-24-S289")
