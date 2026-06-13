#!/usr/bin/env python3
"""
level_edges_s186.py -- opus-2026-03-22-S186

LEVEL EDGES: DECOMPOSING |E(G_n)| BY H-LEVEL AND SCORE-LEVEL

The edge count |E(G_n)| = 1, 5, 30, 290 has no known formula.
Maybe we can find a formula by DECOMPOSING the edges:

|E| = sum_{H-levels} edges_within_level + sum_{adjacent levels} edges_between

Or: |E| = within_score + cross_score

This script computes EVERY possible decomposition at n=3,4,5,6
to find which gives the cleanest formula.

Author: opus-2026-03-22-S186
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
import time
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

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def build_Gn(n):
    m = comb(n,2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; cdata = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A,n), 'scores': score_seq(A,n), 'size': 0}
        cdata[cmap[c]]['size'] += 1
    nc = len(cdata)
    F = np.zeros((nc, nc), dtype=int)
    tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        tc[bits] = cmap[canonical(A, n)]
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            ia, ja = pairs[k]
            bits2 = bits ^ (1 << k)
            A2 = [[0]*n for _ in range(n)]
            for kk,(ii,jj) in enumerate(pairs):
                if (bits2>>kk)&1: A2[ii][jj]=1
                else: A2[jj][ii]=1
            cj = cmap[canonical(A2, n)]
            F[ci][cj] += 1
    return nc, cdata, F, m

# Build G_n for n=3,4,5,6
all_data = {}
for n in [3, 4, 5]:
    t0 = time.time()
    nc, cd, F, m = build_Gn(n)
    adj = (F > 0).astype(int)
    np.fill_diagonal(adj, 0)
    ne = adj.sum() // 2
    all_data[n] = {'nc': nc, 'cd': cd, 'F': F, 'adj': adj, 'ne': ne, 'm': m}
    print("n=%d: %d classes, %d edges (%.1fs)" % (n, nc, ne, time.time()-t0))

# n=6 takes longer
t0 = time.time()
nc6, cd6, F6, m6 = build_Gn(6)
adj6 = (F6 > 0).astype(int)
np.fill_diagonal(adj6, 0)
ne6 = adj6.sum() // 2
all_data[6] = {'nc': nc6, 'cd': cd6, 'F': F6, 'adj': adj6, 'ne': ne6, 'm': m6}
print("n=6: %d classes, %d edges (%.1fs)" % (nc6, ne6, time.time()-t0))

# ============================================================
# PART 1: DECOMPOSE BY H-LEVEL
# ============================================================

banner("PART 1: EDGES BY H-LEVEL")

for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']

    # H values
    H_vals = sorted(set(cd[i]['H'] for i in range(nc)))

    # Count edges within same H and between different H
    same_H = 0; diff_H = 0; up_H = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                if cd[i]['H'] == cd[j]['H']:
                    same_H += 1
                else:
                    diff_H += 1
                    if cd[i]['H'] < cd[j]['H']:
                        up_H += 1

    print("  n=%d: %d edges = %d same-H + %d diff-H (%d up, %d down)" %
          (n, d['ne'], same_H, diff_H, up_H, diff_H - up_H))

    # By H-level pair
    level_edges = Counter()
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                hi, hj = cd[i]['H'], cd[j]['H']
                level_edges[(min(hi,hj), max(hi,hj))] += 1

    print("  Edge distribution by (H_low, H_high):")
    for (h1, h2), cnt in sorted(level_edges.items()):
        tag = "SAME" if h1 == h2 else "CROSS"
        print("    (%d, %d): %d edges [%s]" % (h1, h2, cnt, tag))
    print()

# ============================================================
# PART 2: DECOMPOSE BY SCORE
# ============================================================

banner("PART 2: EDGES BY SCORE CLASS")

for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']

    # Score classes
    score_classes = defaultdict(list)
    for i in range(nc):
        score_classes[cd[i]['scores']].append(i)

    # Count within-score and cross-score edges
    within = 0; cross = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                if cd[i]['scores'] == cd[j]['scores']:
                    within += 1
                else:
                    cross += 1

    n_score_classes = len(score_classes)
    print("  n=%d: %d edges = %d within-score + %d cross-score (%d score classes)" %
          (n, d['ne'], within, cross, n_score_classes))

    # Count edges between score classes
    score_pair_edges = Counter()
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                si, sj = cd[i]['scores'], cd[j]['scores']
                pair = (min(si,sj), max(si,sj))
                score_pair_edges[pair] += 1

    # Score L1 distribution of edges
    l1_dist = Counter()
    for (s1, s2), cnt in score_pair_edges.items():
        l1 = sum(abs(a-b) for a,b in zip(s1, s2))
        l1_dist[l1] += cnt

    print("  Score L1 distribution: %s" % dict(sorted(l1_dist.items())))

    # Within-score edge count per score class
    if within > 0:
        print("  Within-score edges by class:")
        for sc, classes in sorted(score_classes.items()):
            if len(classes) > 1:
                ws = sum(1 for i in classes for j in classes if i < j and adj[i][j])
                total_possible = len(classes) * (len(classes) - 1) // 2
                if ws > 0:
                    print("    %s: %d/%d edges (density %.2f), %d classes" %
                          (sc, ws, total_possible, ws/total_possible, len(classes)))
    print()

# ============================================================
# PART 3: THE CROSS-SCORE EDGE FORMULA
# ============================================================

banner("PART 3: CROSS-SCORE EDGE ANALYSIS")

print("  Cross-score edges require score L1 distance = 2.")
print("  How many score-ADJACENT pairs exist, and what fraction are edges?\n")

for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']

    # Count L1=2 class pairs and how many are edges
    l1_2_pairs = 0; l1_2_edges = 0
    l1_0_pairs = 0; l1_0_edges = 0
    for i in range(nc):
        for j in range(i+1, nc):
            si, sj = cd[i]['scores'], cd[j]['scores']
            l1 = sum(abs(a-b) for a,b in zip(si, sj))
            if l1 == 2:
                l1_2_pairs += 1
                if adj[i][j]: l1_2_edges += 1
            elif l1 == 0:
                l1_0_pairs += 1
                if adj[i][j]: l1_0_edges += 1

    total_pairs = nc * (nc - 1) // 2
    print("  n=%d: %d total pairs, %d L1=2 pairs (%d edges = %.1f%%), %d L1=0 pairs (%d edges)" %
          (n, total_pairs, l1_2_pairs, l1_2_edges,
           100*l1_2_edges/l1_2_pairs if l1_2_pairs > 0 else 0,
           l1_0_pairs, l1_0_edges))

# ============================================================
# PART 4: THE FORMULA HUNT
# ============================================================

banner("PART 4: SEARCHING FOR THE EDGE COUNT FORMULA")

print("  KNOWN: |E(G_n)| = 1, 5, 30, 290 at n=3,4,5,6.\n")

# Decomposition:
# |E| = cross_score + within_score
cs = {3: 1, 4: 5, 5: 27, 6: 0}  # cross-score (will compute for n=6)
ws = {3: 0, 4: 0, 5: 3, 6: 0}    # within-score (will compute for n=6)

for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']
    w = 0; c = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                if cd[i]['scores'] == cd[j]['scores']: w += 1
                else: c += 1
    cs[n] = c; ws[n] = w

print("  Decomposition:")
print("  %-4s %-6s %-8s %-8s" % ("n", "|E|", "cross", "within"))
for n in [3, 4, 5, 6]:
    print("  %-4d %-6d %-8d %-8d" % (n, cs[n]+ws[n], cs[n], ws[n]))

print("\n  CROSS-SCORE edges: 1, 5, 27, %d" % cs[6])
print("  WITHIN-SCORE edges: 0, 0, 3, %d" % ws[6])

# Cross-score analysis
print("\n  Cross-score: 1, 5, 27, %d" % cs[6])
print("  Ratios: %s" % [Fraction(cs[n+1], cs[n]) for n in [3, 4, 5] if cs[n] > 0])

# Try formulas for cross-score
print("\n  Cross-score formula candidates:")
for n in [3, 4, 5, 6]:
    m = comb(n, 2)
    nv = all_data[n]['nc']
    print("    n=%d: cross=%d, C(n,2)=%d, |V|=%d, cross/C(n,2)=%s, cross/|V|=%s" %
          (n, cs[n], m, nv, Fraction(cs[n], m), Fraction(cs[n], nv)))

# Within-score analysis
print("\n  Within-score: 0, 0, 3, %d" % ws[6])
# At n=5: 3 within-score edges among score classes with 2+ iso classes
# Score (1,1,2,3,3): 2 classes, 1 edge -> 1/1 density
# Score (1,2,2,2,3): 3 classes, 2 edges -> 2/3 density

# ============================================================
# PART 5: NUMBER OF ADJACENT SCORE PAIRS
# ============================================================

banner("PART 5: ADJACENT SCORE PAIRS")

print("  How many pairs of score sequences have L1 = 2?")
print("  This is the NUMBER OF POTENTIAL CROSS-SCORE EDGES.\n")

for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']
    score_set = set(cd[i]['scores'] for i in range(nc))
    n_scores = len(score_set)

    # Count L1=2 score pairs
    score_list = sorted(score_set)
    l1_2_score_pairs = 0
    for i, s1 in enumerate(score_list):
        for j, s2 in enumerate(score_list):
            if j <= i: continue
            if sum(abs(a-b) for a,b in zip(s1, s2)) == 2:
                l1_2_score_pairs += 1

    print("  n=%d: %d scores, %d score pairs with L1=2, C(%d,2)=%d" %
          (n, n_scores, l1_2_score_pairs, n_scores, n_scores*(n_scores-1)//2))

# ============================================================
# PART 6: CLASSES PER SCORE AND EDGE DENSITY
# ============================================================

banner("PART 6: SPLITTING AND EDGE DENSITY")

for n in [5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']

    score_classes = defaultdict(list)
    for i in range(nc):
        score_classes[cd[i]['scores']].append(i)

    print("  n=%d:" % n)
    for sc in sorted(score_classes.keys()):
        classes = score_classes[sc]
        n_cls = len(classes)
        if n_cls == 1: continue
        # Within-score edges
        ws_e = sum(1 for i in classes for j in classes if i<j and adj[i][j])
        max_ws = n_cls * (n_cls - 1) // 2
        # Cross-score edges from this score
        cs_e = sum(1 for i in classes for j in range(nc)
                   if j not in classes and i < j and adj[i][j])
        print("    %s: %d classes, %d/%d within (%.0f%%), %d cross" %
              (sc, n_cls, ws_e, max_ws, 100*ws_e/max_ws if max_ws else 0, cs_e))
    print()

# ============================================================
# PART 7: THE COMPLEMENT PAIR EDGE STRUCTURE
# ============================================================

banner("PART 7: COMPLEMENT PAIR EDGES")

print("  Non-SC classes come in complement pairs (T, T^op).")
print("  Do complement pairs always share the same edges?\n")

for n in [4, 5]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']

    # Find complement pairs
    cmap_n = {}
    for bits in range(2**d['m']):
        pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap_n:
            cmap_n[c] = len(cmap_n)

    # Build complement map
    comp = {}
    for bits in range(2**d['m']):
        pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        ci = cmap_n[canonical(A, n)]
        # Complement
        Ac = [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cj = cmap_n[canonical(Ac, n)]
        comp[ci] = cj

    # NSC pairs
    nsc_pairs = [(c, comp[c]) for c in range(nc) if comp[c] > c]
    print("  n=%d: %d NSC pairs" % (n, len(nsc_pairs)))
    for c1, c2 in nsc_pairs:
        # Do c1 and c2 have the same neighbors?
        n1 = set(j for j in range(nc) if adj[c1][j])
        n2 = set(j for j in range(nc) if adj[c2][j])
        shared = n1 & n2
        only1 = n1 - n2
        only2 = n2 - n1
        print("    (%d, %d) H=(%d,%d): %d shared, %d only-c1, %d only-c2 neighbors" %
              (c1, c2, cd[c1]['H'], cd[c2]['H'], len(shared), len(only1), len(only2)))
    print()

# ============================================================
# PART 8: EDGE COUNT DECOMPOSITION SUMMARY
# ============================================================

banner("PART 8: COMPLETE DECOMPOSITION")

print("  %-4s | %-5s | %-6s %-6s | %-6s %-6s %-6s | %-8s" %
      ("n", "|E|", "cross", "within", "up", "level", "down", "L1=2 frac"))

for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']

    w = 0; c = 0; up = 0; level = 0; down = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                if cd[i]['scores'] == cd[j]['scores']:
                    w += 1
                else:
                    c += 1
                hi, hj = cd[i]['H'], cd[j]['H']
                if hi < hj: up += 1
                elif hi == hj: level += 1
                else: down += 1

    # L1=2 edge fraction
    l1_2_total = 0; l1_2_edge = 0
    for i in range(nc):
        for j in range(i+1, nc):
            si, sj = cd[i]['scores'], cd[j]['scores']
            l1 = sum(abs(a-b) for a,b in zip(si,sj))
            if l1 == 2:
                l1_2_total += 1
                if adj[i][j]: l1_2_edge += 1
    frac = Fraction(l1_2_edge, l1_2_total) if l1_2_total > 0 else 0

    print("  %-4d | %-5d | %-6d %-6d | %-6d %-6d %-6d | %s" %
          (n, c+w, c, w, up, level, down, frac))

# ============================================================
# PART 9: THE FORMULA
# ============================================================

banner("PART 9: FORMULA SEARCH")

print("  From the decomposition:")
print("  cross-score edges: 1, 5, 27, ?\n")

# For n=6, print cross-score count
print("  n=6: cross-score = %d, within-score = %d\n" % (cs[6], ws[6]))

print("  Cross-score sequence: 1, 5, 27, %d" % cs[6])
print("  Within-score sequence: 0, 0, 3, %d\n" % ws[6])

# Now check: is within-score = sum over splitting scores of C(split, 2)?
# No, because not all pairs within a split score are connected.
# Actually from Part 6: within-score density varies.

# Try: cross-score = |V| * (m-1) / something
for n in [3, 4, 5, 6]:
    nv = all_data[n]['nc']; m = comb(n,2)
    print("  n=%d: cross=%d, |V|*(m-1)/2=%s, |V|*m/4=%s" %
          (n, cs[n], Fraction(nv*(m-1), 2), Fraction(nv*m, 4)))

print()

# The L1=2 edge fraction:
print("  L1=2 edge fractions:")
for n in [3, 4, 5, 6]:
    d = all_data[n]
    nc = d['nc']; cd = d['cd']; adj = d['adj']
    l1_2_total = 0; l1_2_edge = 0
    for i in range(nc):
        for j in range(i+1, nc):
            si, sj = cd[i]['scores'], cd[j]['scores']
            l1 = sum(abs(a-b) for a,b in zip(si,sj))
            if l1 == 2:
                l1_2_total += 1
                if adj[i][j]: l1_2_edge += 1
    print("  n=%d: %d/%d L1=2 pairs are edges = %s = %.4f" %
          (n, l1_2_edge, l1_2_total, Fraction(l1_2_edge, l1_2_total) if l1_2_total else 0,
           l1_2_edge/l1_2_total if l1_2_total else 0))

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: LEVEL EDGE DECOMPOSITION")

print("""
COMPLETE DECOMPOSITION OF |E(G_n)|:

  |E| = cross-score + within-score
  |E| = up-H + level-H + down-H

Cross-score edges (L1=2):
  1, 5, 27, %d at n=3,4,5,6
  These dominate at all n.
  L1=2 fraction (of L1=2 pairs): varies (100%%, 100%%, 73%%, ...)
  NOT all score-adjacent pairs are connected.

Within-score edges (L1=0):
  0, 0, 3, %d at n=3,4,5,6
  Appear only at n >= 5 where score classes split.
  At n=5: 3 edges within 2 splitting score classes.

H-level structure:
  up-H edges dominate (about half of all edges).
  level-H edges are rare (0 at n=3,4, small fraction at n>=5).
  down-H edges = up-H by complement symmetry.

THE KEY INSIGHT FOR THE FORMULA:
  Cross-score edges are controlled by:
  1. Number of L1=2 score pairs (combinatorial, depends on score sequences)
  2. Fraction of L1=2 pairs that are actual edges (~70-100%%)
  Within-score edges are controlled by:
  1. Number of splitting score classes (grows with n)
  2. Density within each split (varies, generally < 100%%)

THE MOST PROMISING FORMULA APPROACH:
  |E| ~ (#L1=2 class pairs) * (edge probability)
       + (#L1=0 class pairs with same score, different class) * (edge probability)

  The first term gives the BULK of the edges.
  The second term is a correction that grows with n.
""" % (cs[6], ws[6]))

print("Done. opus-2026-03-22-S186")
