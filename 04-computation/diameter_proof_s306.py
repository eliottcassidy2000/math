#!/usr/bin/env python3
"""
diameter_proof_s306.py -- opus-2026-03-24-S306

PROVING diam(G_n) = n-2

Verified: n=4: diam=1, n=5: diam=3, n=6: diam=4.
Pattern: 1, 3, 4 for n=4,5,6. This is n-3, n-2, n-2.
Wait — n=4: diam=1=n-3. n=5: diam=3=n-2. n=6: diam=4=n-2.

So the conjecture is diam = n-2 for n ≥ 5, and diam = 1 for n = 4.

WHY diam ≤ n-2:
  The wiggly metagraph has V_n nodes.
  From the completeness theorem: any two classes C, D can be
  connected by a path of wiggly edges of length ≤ diam.
  Each wiggly edge = 1 tile flip.
  So we need at most diam tile flips to get from any class to any other.

WHY diam ≥ n-2 (for n ≥ 5):
  Need to show there exist classes C, D with min Hamming distance = n-2.
  The TRANSITIVE class has 1 tiling (at bits=0: all arcs go high→low).
  The REGULAR class (if it exists at odd n) has few tilings.
  If the minimum Hamming distance between any transitive tiling
  and any regular tiling is n-2, then diam ≥ n-2.

APPROACH: compute the min Hamming distance from the transitive tiling
to each class, and identify which class is farthest.

Author: opus-2026-03-24-S306
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

# ============================================================
banner("DISTANCE FROM TRANSITIVE TO EACH CLASS")
# ============================================================

for n in [4, 5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    # Build classes
    cmap = {}; tc = {}; members = defaultdict(list); sizes = {}
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap); sizes[cmap[c]] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid; members[cid].append(bits)

    nc = len(cmap)
    info = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        info[cid] = {'comp': cmap.get(comp_c, -1), 'H': Hcount(A, n)}

    # Merge
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        if comp >= 0 and comp != cid:
            merged[mid] = (cid, comp); mid_map[cid] = mid; mid_map[comp] = mid
        else:
            merged[mid] = (cid,); mid_map[cid] = mid
        mid += 1
    nm = mid

    # The transitive tiling is bits=0 (all arcs go high→low = natural direction)
    trans_mid = mid_map[tc[0]]
    trans_H = info[merged[trans_mid][0]]['H']

    print("  n=%d: m=%d, %d merged classes, transitive at mid=%d (H=%d)\n" %
          (n, m, nm, trans_mid, trans_H))

    # For each merged class: find min Hamming distance from transitive tiling (bits=0)
    print("  %-4s %-4s %-6s %-8s %-40s" %
          ("mid", "H", "fiber", "min_d", "closest tiling"))
    print("  " + "-"*65)

    max_min_d = 0
    farthest_class = -1

    for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
        H = info[merged[mi][0]]['H']
        all_tilings = []
        for c in merged[mi]:
            all_tilings.extend(members[c])
        fiber = len(all_tilings)

        # Min Hamming distance from bits=0 to any tiling in this class
        min_d = m + 1
        closest = -1
        for t in all_tilings:
            d = bin(t).count('1')  # distance from 0 = Hamming weight
            if d < min_d:
                min_d = d
                closest = t

        if min_d > max_min_d:
            max_min_d = min_d
            farthest_class = mi

        closest_str = format(closest, '0%db' % m) if closest >= 0 else "-"
        print("  %-4d %-4d %-6d %-8d %-40s" %
              (mi, H, fiber, min_d, closest_str))

    print("\n  FARTHEST class from transitive: mid=%d (H=%d), min_dist=%d" %
          (farthest_class, info[merged[farthest_class][0]]['H'], max_min_d))
    print("  This gives diam ≥ %d" % max_min_d)

    # Now check: is the FULL diameter = max over ALL pairs?
    # We already know from S301 that diam = max min_dist over all pairs.
    # The transitive is often one endpoint. Let's verify.

    # Compute BFS distances in the wiggly metagraph
    adj = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = mid_map[tc[bits]]
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            if mi != mj:
                adj[mi][mj] = 1; adj[mj][mi] = 1

    from scipy.sparse.csgraph import shortest_path
    dist_matrix = shortest_path(adj.astype(float), directed=False, unweighted=True)
    graph_diam = int(dist_matrix.max())

    print("  Graph diameter (BFS) = %d" % graph_diam)
    print("  diam = %d = n - %d" % (graph_diam, n - graph_diam))

    # Which pair achieves the diameter?
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if dist_matrix[mi][mj] == graph_diam:
                Hi = info[merged[mi][0]]['H']
                Hj = info[merged[mj][0]]['H']
                print("  Diameter pair: {%d(H=%d), %d(H=%d)}" % (mi, Hi, mj, Hj))

    print()

# ============================================================
banner("THE TRANSITIVE DISTANCE PROFILE")
# ============================================================

# For each n: the distance from the transitive class to each other class
# This gives the "distance profile" — is it palindromic?

for n in [5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    cmap = {}; tc = {}; members = defaultdict(list); sizes = {}
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap); sizes[cmap[c]] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid; members[cid].append(bits)

    nc = len(cmap)
    info2 = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        info2[cid] = {'comp': cmap.get(comp_c, -1), 'H': Hcount(A, n)}

    merged2 = {}; mid_map2 = {}; mid2 = 0
    for cid in range(nc):
        if cid in mid_map2: continue
        comp = info2[cid]['comp']
        if comp >= 0 and comp != cid:
            merged2[mid2] = (cid, comp); mid_map2[cid] = mid2; mid_map2[comp] = mid2
        else:
            merged2[mid2] = (cid,); mid_map2[cid] = mid2
        mid2 += 1
    nm2 = mid2

    adj = np.zeros((nm2, nm2), dtype=int)
    for bits in range(total):
        mi = mid_map2[tc[bits]]
        for k in range(m):
            mj = mid_map2[tc[bits ^ (1 << k)]]
            if mi != mj: adj[mi][mj] = 1

    from scipy.sparse.csgraph import shortest_path
    dist = shortest_path(adj.astype(float), directed=False, unweighted=True)

    trans_mid = mid_map2[tc[0]]

    # Distance profile from transitive
    dist_profile = Counter()
    for mi in range(nm2):
        d = int(dist[trans_mid][mi])
        dist_profile[d] += 1

    print("  n=%d: Distance profile from transitive (mid=%d, H=%d):" %
          (n, trans_mid, info2[merged2[trans_mid][0]]['H']))
    for d in sorted(dist_profile.keys()):
        bar = "█" * dist_profile[d]
        print("    d=%d: %d classes  %s" % (d, dist_profile[d], bar))

    profile_list = [dist_profile.get(d, 0) for d in range(max(dist_profile.keys())+1)]
    is_palindrome = profile_list == profile_list[::-1]
    print("  Profile: %s" % profile_list)
    print("  Palindromic: %s" % is_palindrome)
    print()

# ============================================================
banner("WHY diam = n-2: THE ARGUMENT")
# ============================================================

print("""
  THE DIAMETER CONJECTURE: diam(G_n) = n - 2 for n ≥ 5.

  EVIDENCE:
    n=4: diam=1 (exception: only 3 classes, all adjacent)
    n=5: diam=3 = n-2 ✓
    n=6: diam=4 = n-2 ✓
    n=7: diam≤7 (from S300 sampling), expect diam=5=n-2

  THE UPPER BOUND: diam ≤ n-2.
    ARGUMENT: any two iso classes can be connected by changing
    at most n-2 tiles. This is because:
    - The staircase has n-2 rows.
    - Two tournaments in different classes differ in their
      treatment of at least one row.
    - Each row change can be accomplished by one tile flip per row.
    - So at most n-2 tile flips suffice (one per row).

    BUT WAIT — this gives an upper bound on the HAMMING distance
    between representatives, not on the GRAPH distance.
    The graph distance ≤ Hamming distance ≤ n-2.
    Actually no: the Hamming distance between optimal representatives
    could be less than n-2 even if graph distance is larger.

    THE CORRECT ARGUMENT for the upper bound:
    From the completeness theorem, graph distance = min Hamming distance.
    So we need: for any two classes C, D, there exist t ∈ C, t' ∈ D
    with Hamming distance ≤ n-2.

    This is NOT obvious. Some classes might be Hamming-far from all
    tilings in all other classes. Need a constructive argument.

  THE LOWER BOUND: diam ≥ n-2 for n ≥ 5.
    ARGUMENT: there exist classes C, D where the minimum Hamming
    distance between any tiling in C and any in D is exactly n-2.

    From the data:
    n=5: the pair {mid=0(H=9), mid=9(H=15)} has min_dist=3=n-2.
    n=6: the pair {mid=0(H=?), mid=?(H=?)} achieves min_dist=4=n-2.

    These are typically the transitive class vs a high-|Aut| class.
""")

# ============================================================
banner("THE HAMMING WEIGHT PROFILE OF THE TRANSITIVE CLASS")
# ============================================================

# The transitive tiling has bits=0 (weight 0).
# Every other tiling has weight ≥ 1.
# The distance from transitive to any other tiling = weight of that tiling.
# So: min_dist(transitive, class C) = min weight of any tiling in C.

for n in [5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    cmap = {}; tc = {}; members = defaultdict(list)
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]; members[cmap[c]].append(bits)

    nc = len(cmap)
    info3 = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        info3[cid] = {'H': Hcount(A, n)}

    # For each class: weight distribution of tilings
    print("  n=%d: Weight distribution per class:\n" % n)
    print("  %-6s %-4s %-30s %-8s" % ("class", "H", "weight distribution", "min_wt"))
    print("  " + "-"*55)

    for cid in sorted(range(nc), key=lambda x: info3[x]['H']):
        H = info3[cid]['H']
        wt_dist = Counter(bin(b).count('1') for b in members[cid])
        min_wt = min(wt_dist.keys())
        wt_str = ", ".join("%d:%d" % (w, c) for w, c in sorted(wt_dist.items()))
        if len(members[cid]) <= 20:
            print("  %-6d %-4d %-30s %-8d" % (cid, H, wt_str, min_wt))

    # The class with the HIGHEST minimum weight is the farthest from transitive
    max_min_wt = 0
    farthest = -1
    for cid in range(nc):
        min_wt = min(bin(b).count('1') for b in members[cid])
        if min_wt > max_min_wt:
            max_min_wt = min_wt
            farthest = cid

    print("\n  FARTHEST from transitive: class %d (H=%d), min_weight=%d" %
          (farthest, info3[farthest]['H'], max_min_wt))
    print("  This contributes diam ≥ %d" % max_min_wt)
    print()

print("Done. opus-2026-03-24-S306")
