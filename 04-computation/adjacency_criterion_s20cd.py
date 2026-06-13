#!/usr/bin/env python3
"""
adjacency_criterion_s20cd.py -- kind-pasteur-2026-03-22-S20cd

THE ADJACENCY CRITERION: When are two iso classes connected in G_n?

Two classes C1, C2 are adjacent iff there exists a tournament T in C1
and an arc (a,b) such that flip(T, a, b) is in C2.

QUESTION: Can adjacency be determined from CLASS INVARIANTS alone?
If yes: G_n is fully determined analytically.
If no: the fine structure requires tournament-level computation.

TEST: Do (H, score, |Aut|, c3) jointly determine adjacency?

Author: kind-pasteur-2026-03-22-S20cd
"""
import sys
import numpy as np
from math import comb, factorial
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

print("=" * 70)
print("  THE ADJACENCY CRITERION")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Build classes with full invariants
canon_map = defaultdict(list)
H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    H_map[bits] = H
    s = A.sum(axis=1).astype(int)
    S2 = int(sum(s*s))
    c3 = comb(n,3) - (S2 - comb(n,2)) // 2
    score = tuple(sorted(s))
    cf = canonical(A, n)
    canon_map[cf].append((bits, H, score, c3, S2))

classes = []
cf_to_id = {}
for cf, members in sorted(canon_map.items(), key=lambda x: x[1][0][1]):
    cid = len(classes)
    cf_to_id[cf] = cid
    bits0, H, score, c3, S2 = members[0]
    aut_size = factorial(n) // len(members)
    classes.append({
        'id': cid, 'H': H, 'score': score, 'c3': c3, 'S2': S2,
        'aut': aut_size, 'size': len(members), 'members': set(b for b, *_ in members)
    })

bits_to_class = {}
for c in classes:
    for b in c['members']:
        bits_to_class[b] = c['id']

N = len(classes)

# Build adjacency
adj = np.zeros((N, N), dtype=int)
for c in classes:
    T = list(c['members'])[0]
    for k in range(m):
        nb = T ^ (1 << k)
        nb_class = bits_to_class[nb]
        if nb_class != c['id']:
            adj[c['id']][nb_class] = 1
            adj[nb_class][c['id']] = 1

# ================================================================
# 1. NECESSARY CONDITIONS FOR ADJACENCY
# ================================================================
print(f"\n{'='*70}")
print(f"  1. NECESSARY CONDITIONS FOR ADJACENCY")
print(f"{'='*70}\n")

# An arc flip changes score[i] by -1 and score[j] by +1.
# So the SORTED score can change by at most swapping two adjacent values.
# This means: adjacent classes must have "close" score sequences.

# What is "close"? The L1 distance between sorted scores.
print(f"  SCORE DISTANCE between adjacent vs non-adjacent pairs:")

adj_dists = []
nonadj_dists = []
for i in range(N):
    for j in range(i+1, N):
        si = np.array(classes[i]['score'])
        sj = np.array(classes[j]['score'])
        d = int(np.sum(np.abs(si - sj)))
        if adj[i][j]:
            adj_dists.append(d)
        else:
            nonadj_dists.append(d)

from collections import Counter
print(f"  Adjacent pairs: L1 distances = {sorted(Counter(adj_dists).items())}")
print(f"  Non-adjacent pairs: L1 distances = {sorted(Counter(nonadj_dists).items())}")

# Is L1 distance <= 2 necessary for adjacency?
max_adj_dist = max(adj_dists)
print(f"\n  Max L1 distance for adjacent pairs: {max_adj_dist}")
print(f"  NECESSARY CONDITION: score L1 distance <= {max_adj_dist}")

# Is L1 distance <= 2 sufficient?
d2_nonadj = sum(1 for d in nonadj_dists if d <= max_adj_dist)
d2_adj = sum(1 for d in adj_dists if d <= max_adj_dist)
print(f"  Pairs with L1 <= {max_adj_dist}: {d2_adj} adjacent + {d2_nonadj} non-adjacent")
print(f"  SUFFICIENT? {'YES' if d2_nonadj == 0 else 'NO'}")

# ================================================================
# 2. |delta_H| CONSTRAINT
# ================================================================
print(f"\n{'='*70}")
print(f"  2. |delta_H| CONSTRAINT")
print(f"{'='*70}\n")

adj_dH = []
nonadj_dH = []
for i in range(N):
    for j in range(i+1, N):
        dH = abs(classes[i]['H'] - classes[j]['H'])
        if adj[i][j]:
            adj_dH.append(dH)
        else:
            nonadj_dH.append(dH)

print(f"  Adjacent: |dH| in {sorted(set(adj_dH))}")
print(f"  Non-adjacent: |dH| in {sorted(set(nonadj_dH))}")
print(f"  Max adjacent |dH|: {max(adj_dH)}")
print(f"  Min non-adjacent |dH|: {min(nonadj_dH) if nonadj_dH else 'N/A'}")

# ================================================================
# 3. WHICH INVARIANTS DETERMINE ADJACENCY?
# ================================================================
print(f"\n{'='*70}")
print(f"  3. ADJACENCY FROM INVARIANTS")
print(f"{'='*70}\n")

# For each pair (C1, C2), compute invariant-based features:
# - |dH|
# - L1 score distance
# - |dc3|
# - |d(|Aut|)|
# See if some combination perfectly predicts adjacency.

features_adj = []
features_nonadj = []

for i in range(N):
    for j in range(i+1, N):
        dH = abs(classes[i]['H'] - classes[j]['H'])
        dS = int(np.sum(np.abs(np.array(classes[i]['score']) - np.array(classes[j]['score']))))
        dc3 = abs(classes[i]['c3'] - classes[j]['c3'])
        same_score = (classes[i]['score'] == classes[j]['score'])
        both_sc = (classes[i]['aut'] > 1 and classes[j]['aut'] > 1)

        feat = (dH, dS, dc3, same_score)
        if adj[i][j]:
            features_adj.append(feat)
        else:
            features_nonadj.append(feat)

# Find features that appear in BOTH adjacent and non-adjacent
adj_feat_set = set(features_adj)
nonadj_feat_set = set(features_nonadj)
ambiguous = adj_feat_set & nonadj_feat_set

print(f"  Feature tuples (|dH|, L1_score, |dc3|, same_score):")
print(f"    Adjacent only: {len(adj_feat_set - nonadj_feat_set)} distinct")
print(f"    Non-adjacent only: {len(nonadj_feat_set - adj_feat_set)} distinct")
print(f"    AMBIGUOUS (in both): {len(ambiguous)} distinct")

if ambiguous:
    print(f"\n  AMBIGUOUS features (can't determine adjacency):")
    for feat in sorted(ambiguous):
        n_adj = features_adj.count(feat)
        n_nonadj = features_nonadj.count(feat)
        print(f"    {feat}: {n_adj} adj, {n_nonadj} non-adj")

    print(f"\n  CONCLUSION: Invariants (dH, L1_score, dc3, same_score)")
    print(f"  do NOT fully determine adjacency at n={n}.")
    print(f"  There are {len(ambiguous)} ambiguous feature tuples.")
    print(f"  The FINE STRUCTURE of G_n requires tournament-level information")
    print(f"  beyond these simple invariants.")
else:
    print(f"\n  CONCLUSION: These invariants DO determine adjacency!")

# ================================================================
# 4. THE FULL ADJACENCY MATRIX
# ================================================================
print(f"\n{'='*70}")
print(f"  4. THE FULL ADJACENCY MATRIX OF G_{n}")
print(f"{'='*70}\n")

print(f"  G_{n} adjacency (1=adjacent, .=not):")
print(f"  {'':>4s}", end="")
for j in range(N):
    print(f" {j:>2d}", end="")
print()
for i in range(N):
    print(f"  {i:>2d}:", end="")
    for j in range(N):
        if i == j:
            print(f"  .", end="")
        elif adj[i][j]:
            print(f"  1", end="")
        else:
            print(f"  .", end="")
    print(f"  H={classes[i]['H']}")

# ================================================================
# 5. GRAPH PROPERTIES OF THE ADJACENCY
# ================================================================
print(f"\n{'='*70}")
print(f"  5. GRAPH PROPERTIES")
print(f"{'='*70}\n")

# Clique number
max_clique = 0
from itertools import combinations
for size in range(N, 0, -1):
    found = False
    for subset in combinations(range(N), size):
        if all(adj[i][j] for i in subset for j in subset if i != j):
            max_clique = size
            found = True
            break
    if found: break

print(f"  Clique number (max complete subgraph): {max_clique}")

# Chromatic number (upper bound from greedy)
colors = [-1] * N
for v in range(N):
    used = set(colors[w] for w in range(N) if adj[v][w] and colors[w] >= 0)
    c = 0
    while c in used: c += 1
    colors[v] = c
chromatic_upper = max(colors) + 1
print(f"  Chromatic number (greedy upper bound): {chromatic_upper}")

# Independence number (from earlier: 5 at n=5)
print(f"  Independence number: {5} (from earlier computation)")

# Diameter
dist = np.full((N, N), N+1)
np.fill_diagonal(dist, 0)
for i in range(N):
    for j in range(N):
        if adj[i][j]: dist[i][j] = 1
for k in range(N):
    for i in range(N):
        for j in range(N):
            if dist[i][k] + dist[k][j] < dist[i][j]:
                dist[i][j] = dist[i][k] + dist[k][j]
diameter = int(dist.max())
print(f"  Diameter: {diameter}")

# Girth
girth = N + 1
for i in range(N):
    for j in range(i+1, N):
        if not adj[i][j]: continue
        for k in range(j+1, N):
            if adj[j][k] and adj[k][i]:
                girth = min(girth, 3)
print(f"  Girth: {girth}")

# Average clustering coefficient
cc_vals = []
for v in range(N):
    neighbors = [w for w in range(N) if adj[v][w]]
    if len(neighbors) < 2:
        cc_vals.append(0)
        continue
    n_edges = sum(1 for i in range(len(neighbors)) for j in range(i+1, len(neighbors)) if adj[neighbors[i]][neighbors[j]])
    max_edges = len(neighbors) * (len(neighbors) - 1) // 2
    cc_vals.append(n_edges / max_edges if max_edges > 0 else 0)

avg_cc = sum(cc_vals) / len(cc_vals)
print(f"  Average clustering coefficient: {avg_cc:.4f}")

# The complement graph
comp_edges = sum(1 for i in range(N) for j in range(i+1, N) if not adj[i][j])
print(f"  Complement edges: {comp_edges}")
print(f"  G_n + complement = K_{N}: {n_edges + comp_edges} = C({N},2) = {comb(N,2)}")

n_edges = sum(1 for i in range(N) for j in range(i+1, N) if adj[i][j])
print(f"\n  SUMMARY: G_{n} = ({N} vertices, {n_edges} edges)")
print(f"  Clique={max_clique}, Chromatic<={chromatic_upper}, Indep=5, Diam={diameter}, CC={avg_cc:.3f}")
