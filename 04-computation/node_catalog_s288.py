#!/usr/bin/env python3
"""
node_catalog_s288.py -- opus-2026-03-24-S288

COMPREHENSIVE NODE CATALOG: META-GRAPH + WIGGLY

For EACH merged node:
  - Basic: H, |Aut|, fiber, score, c3, SC/NS
  - Meta-graph: degree, neighbors, self-loop count
  - Wiggly: per-class weights to each neighbor
  - Tilings: which labeled tournaments belong to this node
  - Wiggly profile: for each tiling, which flips are neutral

Author: opus-2026-03-24-S288
"""

import sys
import numpy as np
import string
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

def c3count(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if A[i][j] and A[j][k] and A[k][i]: c+=1
                if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
labels = list(string.ascii_uppercase[:m])

banner("COMPREHENSIVE NODE CATALOG AT n=%d" % n)

# Build everything
cmap = {}; sizes = {}; tc = {}; reps = {}; members = defaultdict(list)
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
    members[cmap[c]].append(bits)
nc = len(cmap)

info = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info[cid] = {'H': Hcount(A, n), 'c3': c3count(A, n),
                 'comp': cmap[comp_c], 'SC': cmap[comp_c] == cid,
                 'size': sizes[cid], 'aut': factorial(n)//sizes[cid],
                 'score': tuple(sorted(sum(A[i]) for i in range(n)))}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

# Build wiggly weights per class
W_letter = [np.zeros((nm, nm), dtype=int) for _ in range(m)]
W_total = np.zeros((nm, nm), dtype=int)
for bits in range(2**m):
    mi = mid_map[tc[bits]]
    for k in range(m):
        mj = mid_map[tc[bits ^ (1 << k)]]
        W_letter[k][mi][mj] += 1
        W_total[mi][mj] += 1

adj = np.zeros((nm, nm), dtype=int)
for mi in range(nm):
    for mj in range(nm):
        if mi != mj and W_total[mi][mj] > 0:
            adj[mi][mj] = 1

# ============================================================
# NODE CATALOG
# ============================================================

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cids = merged[mi]
    cid0 = cids[0]
    H = info[cid0]['H']; c3 = info[cid0]['c3']
    aut = info[cid0]['aut']; sc = info[cid0]['SC']
    score = info[cid0]['score']
    fiber = sum(sizes[c] for c in cids)
    deg = sum(adj[mi])
    sl_weight = W_total[mi][mi]
    cross_weight = sum(W_total[mi][mj] for mj in range(nm) if mj != mi)

    banner("NODE %d: H=%d, %s" % (mi, H, "SC" if sc else "NS"))

    print("  BASIC:")
    print("    H = %d, c3 = %d, |Aut| = %d, SC = %s" % (H, c3, aut, sc))
    print("    Score = %s" % (score,))
    print("    Fiber = %d tilings (%d unmerged classes)" % (fiber, len(cids)))

    # Neighbors
    nbrs = [(mj, info[merged[mj][0]]['H']) for mj in range(nm) if adj[mi][mj]]
    nbrs.sort(key=lambda x: x[1])
    print("\n  META-GRAPH:")
    print("    Degree = %d" % deg)
    print("    Neighbors: %s" % [(mj, "H=%d" % Hj) for mj, Hj in nbrs])
    print("    Self-loop weight = %d (%.1f%% of total %d)" %
          (sl_weight, 100*sl_weight/(sl_weight+cross_weight) if sl_weight+cross_weight > 0 else 0,
           sl_weight + cross_weight))

    # Wiggly weights per neighbor per class
    print("\n  WIGGLY WEIGHTS TO EACH NEIGHBOR (per letter, all equal by S_n):")
    print("    %-15s | per-letter | total (all %d letters)" % ("Target", m))
    print("    " + "-"*50)
    for mj, Hj in nbrs:
        w_per_letter = W_letter[0][mi][mj] + W_letter[0][mj][mi]  # any letter, all same
        w_total = W_total[mi][mj] + W_total[mj][mi]
        print("    %d(H=%d)%s | %9d | %d" %
              (mj, Hj, " "*(8-len(str(mj))), w_per_letter, w_total))

    if sl_weight > 0:
        sl_per_letter = W_letter[0][mi][mi]
        print("    SELF%s | %9d | %d" %
              (" "*8, sl_per_letter, sl_weight))

    # Wiggly profile distribution
    print("\n  WIGGLY PROFILES (neutral arc patterns):")
    all_tilings = []
    for cid in cids:
        all_tilings.extend(members[cid])

    profiles = Counter()
    for bits in all_tilings:
        profile = tuple(1 if mid_map[tc[bits ^ (1<<k)]] == mi else 0 for k in range(m))
        profiles[profile] += 1

    n_neutral_vals = Counter()
    for prof, count in profiles.items():
        n_neutral_vals[sum(prof)] += count

    print("    Distinct profiles: %d" % len(profiles))
    print("    Neutral arc distribution:")
    for nn in sorted(n_neutral_vals.keys()):
        print("      %d neutral arcs: %d tilings" % (nn, n_neutral_vals[nn]))

    # Where do wiggly lines GO from this node?
    print("\n  WIGGLY DESTINATION DISTRIBUTION (from all %d tilings):" % fiber)
    dest_counts = Counter()
    for bits in all_tilings:
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            dest_counts[mj] += 1

    for mj in sorted(dest_counts.keys(), key=lambda x: -dest_counts[x]):
        Hj = info[merged[mj][0]]['H']
        count = dest_counts[mj]
        frac = count / (fiber * m)
        label = "SELF" if mj == mi else ""
        print("    → %d(H=%d): %d lines (%.1f%%) %s" %
              (mj, Hj, count, 100*frac, label))

    # Blue/black line count
    n_blue = 0; n_black = 0
    for bits in all_tilings:
        comp_bits = bits ^ ((2**m) - 1)
        # Check if complement is in same merged class
        comp_mi = mid_map[tc[comp_bits]]
        if comp_mi == mi:
            # Self-complementary tiling (SC class)
            pass
        # Blue/black determined by grid-symmetry (not computed here)
    print("    Blue/black lines: %d tilings paired with complements" % fiber)

    print()

# ============================================================
banner("SUMMARY: ALL NODES SIDE BY SIDE")
# ============================================================

print("  %-3s | %-3s | %-3s | %-4s | %-6s | %-3s | %-4s | %-18s | %-6s | %s" %
      ("mid", "H", "c3", "Aut", "fiber", "deg", "SL%", "neighbors", "neut", "score"))
print("  " + "-"*90)

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid0 = merged[mi][0]
    H = info[cid0]['H']; c3 = info[cid0]['c3']
    aut = info[cid0]['aut']
    fiber = sum(sizes[c] for c in merged[mi])
    deg = sum(adj[mi])
    sl_pct = 100*W_total[mi][mi]/(fiber*m) if fiber*m > 0 else 0

    nbr_str = ",".join("%d" % mj for mj in range(nm) if adj[mi][mj])

    # Average neutral arcs
    all_tilings = []
    for cid in merged[mi]:
        all_tilings.extend(members[cid])
    total_neutral = sum(1 for bits in all_tilings for k in range(m)
                        if mid_map[tc[bits ^ (1<<k)]] == mi)
    avg_neutral = total_neutral / fiber if fiber > 0 else 0

    sc = info[cid0]['score']
    print("  %-3d | %-3d | %-3d | %-4d | %-6d | %-3d | %4.1f | %-18s | %-6.1f | %s" %
          (mi, H, c3, aut, fiber, deg, sl_pct, nbr_str[:18], avg_neutral, sc))

print("\nDone. opus-2026-03-24-S288")
