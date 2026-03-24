#!/usr/bin/env python3
"""
wiggly_subgraphs_s274.py -- opus-2026-03-23-S274

THE WIGGLY CLASS SUBGRAPHS

Each wiggly class X (= one tile position) defines a subgraph:
  Nodes: merged iso classes
  Edges: {C, D} is an edge iff ∃ labeled T in C such that
         flipping tile X in T gives a tournament in class D.

From local=global (S262): EVERY class X sees EVERY meta-graph edge.
So the subgraph for each class has the SAME edges as G_n/Z_2.

But the WEIGHTS differ: W_X(C, D) = #{labeled T in C : flip X gives class D}.
These weights encode the FINE STRUCTURE of each tile position.

ANALYSIS:
1. Weight distributions per wiggly class
2. Which classes have the most uniform vs skewed weights
3. Self-loop weights per class (= neutral arcs at position X)
4. The STAIRCASE GEOMETRY of the weight landscape
5. Do different wiggly classes see different H-gradients?

Author: opus-2026-03-23-S274
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

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

banner("WIGGLY CLASS SUBGRAPHS AT n=%d" % n)

# Build everything
cmap = {}; sizes = {}; tc = {}; reps = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    sizes[cmap[c]] += 1; tc[bits] = cmap[c]
nc = len(cmap)

info = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                 'SC': cmap[comp_c] == cid, 'size': sizes[cid]}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

# Build per-wiggly-class weight matrices
# W_k[mi][mj] = #{labeled T : tc[T] in merged[mi], tc[T^k] in merged[mj]}
W = [np.zeros((nm, nm), dtype=int) for _ in range(m)]

for bits in range(2**m):
    mi = mid_map[tc[bits]]
    for k in range(m):
        mj = mid_map[tc[bits ^ (1 << k)]]
        W[k][mi][mj] += 1

# Total weight matrix (sum over all classes)
W_total = sum(W)

# The FULL meta-graph adjacency (unweighted)
A_adj = np.zeros((nm, nm), dtype=int)
for mi in range(nm):
    for mj in range(nm):
        if mi != mj and W_total[mi][mj] > 0:
            A_adj[mi][mj] = 1

ne = A_adj.sum() // 2

print("  V_merged = %d, E_merged = %d, m = %d wiggly classes\n" % (nm, ne, m))

# ============================================================
# PER-CLASS WEIGHT ANALYSIS
# ============================================================

banner("WEIGHT MATRICES PER WIGGLY CLASS")

H_vec = [info[merged[mi][0]]['H'] for mi in range(nm)]

for k in range(m):
    arc = pairs[k]
    # Range in staircase
    r = arc[1] - arc[0]

    # Self-loop weight for this class
    self_weight = sum(W[k][mi][mi] for mi in range(nm))
    cross_weight = W[k].sum() - self_weight

    # How uniform is the cross-weight distribution?
    cross_vals = []
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if A_adj[mi][mj]:
                w = W[k][mi][mj] + W[k][mj][mi]
                cross_vals.append(w)

    avg_w = np.mean(cross_vals) if cross_vals else 0
    std_w = np.std(cross_vals) if cross_vals else 0
    cv = std_w / avg_w if avg_w > 0 else 0

    # Expected ΔH for this class
    dH_sum = 0; dH_count = 0
    for bits in range(2**m):
        mi = mid_map[tc[bits]]
        mj = mid_map[tc[bits ^ (1 << k)]]
        dH_sum += abs(H_vec[mj] - H_vec[mi])
        dH_count += 1
    avg_dH = dH_sum / dH_count

    print("  Class %d: arc %s (range %d)" % (k, arc, r))
    print("    Self-loop weight: %d (%.1f%% of %d total)" %
          (self_weight, 100*self_weight/W[k].sum(), W[k].sum()))
    print("    Cross-weight: avg=%.1f, std=%.1f, CV=%.3f" % (avg_w, std_w, cv))
    print("    Avg |ΔH| per flip: %.2f" % avg_dH)
    print()

# ============================================================
banner("STAIRCASE GEOMETRY OF WEIGHTS")
# ============================================================

# Group wiggly classes by their RANGE (= distance in staircase)
range_groups = defaultdict(list)
for k in range(m):
    r = pairs[k][1] - pairs[k][0]
    range_groups[r].append(k)

print("  WEIGHTS BY RANGE (staircase distance):\n")
print("  %-6s | %-6s | %-12s | %-12s | %-10s | %-10s" %
      ("Range", "#arcs", "Self-loop %", "Avg cross-w", "CV", "Avg |ΔH|"))
print("  " + "-"*65)

for r in sorted(range_groups.keys()):
    arcs = range_groups[r]
    # Aggregate stats
    total_self = sum(sum(W[k][mi][mi] for mi in range(nm)) for k in arcs)
    total_all = sum(W[k].sum() for k in arcs)
    self_pct = 100 * total_self / total_all if total_all > 0 else 0

    all_cross = []
    for k in arcs:
        for mi in range(nm):
            for mj in range(mi+1, nm):
                if A_adj[mi][mj]:
                    all_cross.append(W[k][mi][mj] + W[k][mj][mi])

    avg_cw = np.mean(all_cross) if all_cross else 0
    cv_cw = np.std(all_cross) / avg_cw if avg_cw > 0 else 0

    # Avg ΔH
    total_dH = 0; total_count = 0
    for k in arcs:
        for bits in range(2**m):
            mi = mid_map[tc[bits]]
            mj = mid_map[tc[bits ^ (1 << k)]]
            total_dH += abs(H_vec[mj] - H_vec[mi])
            total_count += 1
    avg_dH = total_dH / total_count if total_count > 0 else 0

    print("  %-6d | %-6d | %10.1f%% | %12.1f | %10.3f | %10.2f" %
          (r, len(arcs), self_pct, avg_cw, cv_cw, avg_dH))

print("""
  THE ΔH = 2^{range} FORMULA:
  Higher-range arcs should create LARGER H-changes.
  Range 1: adjacent vertices, smallest ΔH
  Range n-1: source-sink (apex), largest ΔH ≈ 2^{n-2}
""")

# ============================================================
banner("WHICH EDGES ARE WITNESSED BY WHICH MULTIPLICITIES?")
# ============================================================

# For each meta-graph edge: compute the weight per wiggly class
print("  EDGE WEIGHTS PER WIGGLY CLASS (sorted by total weight):\n")

edge_weights = {}
for mi in range(nm):
    for mj in range(mi+1, nm):
        if not A_adj[mi][mj]: continue
        weights = []
        for k in range(m):
            w = W[k][mi][mj] + W[k][mj][mi]
            weights.append(w)
        edge_weights[(mi, mj)] = weights

# Sort by total weight
for (mi, mj), weights in sorted(edge_weights.items(), key=lambda x: -sum(x[1])):
    Hi = H_vec[mi]; Hj = H_vec[mj]
    total_w = sum(weights)
    # Group by range
    range_weights = defaultdict(int)
    for k in range(m):
        r = pairs[k][1] - pairs[k][0]
        range_weights[r] += weights[k]

    rw_str = ", ".join("r%d:%d" % (r, w) for r, w in sorted(range_weights.items()) if w > 0)
    print("  {%d(H=%d),%d(H=%d)}: total=%d [%s]" %
          (mi, Hi, mj, Hj, total_w, rw_str))

# ============================================================
banner("SELF-LOOP WEIGHTS PER WIGGLY CLASS")
# ============================================================

# Self-loops = neutral arcs. Per wiggly class k:
# W_k[mi][mi] = #{T in class mi : flipping arc k gives T' in same class}

print("  SELF-LOOPS PER NODE PER WIGGLY CLASS:\n")
for mi in range(nm):
    H = H_vec[mi]
    sl_by_range = defaultdict(int)
    total_sl = 0
    for k in range(m):
        r = pairs[k][1] - pairs[k][0]
        sl_by_range[r] += W[k][mi][mi]
        total_sl += W[k][mi][mi]

    if total_sl > 0:
        rsl_str = ", ".join("r%d:%d" % (r, w) for r, w in sorted(sl_by_range.items()) if w > 0)
        fiber = sum(sizes[cid] for cid in merged[mi])
        print("  Node %d (H=%d, fiber=%d): SL=%d [%s]" %
              (mi, H, fiber, total_sl, rsl_str))

# ============================================================
banner("THE WIGGLY WEIGHT AS A TOURNAMENT INVARIANT")
# ============================================================

print("""
  For each labeled tournament T and each wiggly class k:
  Define w_k(T) = 1 if flipping arc k preserves the iso class, 0 otherwise.

  Then:
  - SL(T) = sum_k w_k(T) = total neutral arcs of T
  - W_k(class C) = sum_{T in C} w_k(T) = weight of self-loop at C via class k

  The WIGGLY PROFILE of T = (w_1(T), w_2(T), ..., w_m(T)) ∈ {0,1}^m.

  This is a BINARY VECTOR encoding which arcs of T are neutral.
  The Hamming weight of this vector = SL(T) = neutral arc count.

  THE WIGGLY PROFILE IS A SECONDARY TOURNAMENT INVARIANT:
  It's not determined by the iso class alone — different labeled
  copies of the same tournament have different wiggly profiles
  (because S_n permutes the arc positions).

  But the MULTISET of wiggly profiles over all labeled copies
  IS an invariant of the iso class. This multiset is exactly
  the weight distribution W_k(C) for each k.
""")

# Compute wiggly profiles for all tournaments at n=5
profile_dist = defaultdict(Counter)  # mi -> Counter of profiles
for bits in range(2**m):
    mi = mid_map[tc[bits]]
    profile = tuple(1 if mid_map[tc[bits ^ (1 << k)]] == mi else 0 for k in range(m))
    profile_dist[mi][profile] += 1

print("  WIGGLY PROFILE DISTRIBUTION PER CLASS:\n")
for mi in range(nm):
    H = H_vec[mi]
    profiles = profile_dist[mi]
    n_distinct = len(profiles)
    total_tilings = sum(profiles.values())
    total_weight = sum(sum(p) * c for p, c in profiles.items())
    avg_neutral = total_weight / total_tilings if total_tilings > 0 else 0

    print("  Node %d (H=%d): %d tilings, %d distinct profiles, avg neutral=%.2f" %
          (mi, H, total_tilings, n_distinct, avg_neutral))
    for profile, count in sorted(profiles.items(), key=lambda x: -x[1])[:3]:
        print("    %s: %d tilings (weight %d)" %
              (profile, count, sum(profile)))

# ============================================================
banner("THE WIGGLY SPECTRUM: EIGENVALUES PER CLASS")
# ============================================================

# Each wiggly class k defines a weight matrix W_k.
# Its eigenvalues reveal the structure of that arc position.

print("  TOP EIGENVALUE OF EACH WIGGLY CLASS WEIGHT MATRIX:\n")
for k in range(m):
    arc = pairs[k]; r = arc[1] - arc[0]
    W_k_sym = (W[k] + W[k].T).astype(float)
    np.fill_diagonal(W_k_sym, 0)
    evals = np.sort(np.linalg.eigvalsh(W_k_sym))[::-1]
    print("  Class %d (arc %s, range %d): λ_max = %.2f, λ_2 = %.2f" %
          (k, arc, r, evals[0], evals[1] if len(evals) > 1 else 0))

print("\nDone. opus-2026-03-23-S274")
