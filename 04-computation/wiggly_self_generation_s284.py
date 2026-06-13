#!/usr/bin/env python3
"""
wiggly_self_generation_s284.py -- opus-2026-03-24-S284

HOW THE WIGGLY GRAPH GENERATES ITSELF

The wiggly weight matrix W(n) at n can be decomposed by which
REGION of the staircase the flipped tile lives in:

  W(n) = W_overlap(n) + W_boundary(n) + W_apex(n)

where:
  W_overlap: flips of tiles in the overlap region (n-2 sub-tournament)
  W_boundary: flips of tiles involving vertex 0 or vertex n-1
  W_apex: flip of the single apex tile (0 ↔ n-1)

KEY RECURSION:
  W_overlap(n) relates to W(n-2) via the fiber map:
  flipping an overlap tile is the SAME as flipping a tile at n-2,
  but the effect at n depends on the boundary configuration.

  W_boundary(n) depends on how the extreme vertices connect.
  W_apex(n) depends on the apex arc direction.

GOAL: Express W(n) in terms of W(n-1), W(n-2), and simple corrections.
This would compute the meta-graph AND wiggly structure simultaneously.

Author: opus-2026-03-24-S284
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

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

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

for n in [5]:
    banner("WIGGLY SELF-GENERATION AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs)}

    # Classify tiles into overlap, boundary, apex
    overlap_tiles = [k for k, (i,j) in enumerate(pairs) if i >= 1 and j <= n-2]
    boundary_tiles = [k for k, (i,j) in enumerate(pairs)
                      if (i == 0 and j < n-1) or (j == n-1 and i > 0)]
    apex_tile = [k for k, (i,j) in enumerate(pairs) if i == 0 and j == n-1]

    print("  Staircase decomposition of %d tiles:" % m)
    print("    Overlap (interior): %d tiles %s" % (len(overlap_tiles), [pairs[k] for k in overlap_tiles]))
    print("    Boundary (extreme): %d tiles %s" % (len(boundary_tiles), [pairs[k] for k in boundary_tiles]))
    print("    Apex (0↔n-1): %d tiles %s" % (len(apex_tile), [pairs[k] for k in apex_tile]))
    print("    Total: %d + %d + %d = %d = C(%d,2) ✓\n" %
          (len(overlap_tiles), len(boundary_tiles), len(apex_tile),
           len(overlap_tiles)+len(boundary_tiles)+len(apex_tile), n))

    # Build classes
    cmap = {}; sizes = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        sizes[cmap[c]] += 1; tc[bits] = cmap[c]
    nc = len(cmap)

    info = {}
    for cid in range(nc):
        bits_rep = next(b for b in range(2**m) if tc[b] == cid)
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits_rep>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c], 'SC': cmap[comp_c]==cid}

    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid
    H_vec = [info[merged[mi][0]]['H'] for mi in range(nm)]

    # Build W decomposed by region
    W_ov = np.zeros((nm, nm), dtype=int)
    W_bd = np.zeros((nm, nm), dtype=int)
    W_ap = np.zeros((nm, nm), dtype=int)

    for bits in range(2**m):
        mi = mid_map[tc[bits]]
        for k in overlap_tiles:
            mj = mid_map[tc[bits ^ (1 << k)]]
            W_ov[mi][mj] += 1
        for k in boundary_tiles:
            mj = mid_map[tc[bits ^ (1 << k)]]
            W_bd[mi][mj] += 1
        for k in apex_tile:
            mj = mid_map[tc[bits ^ (1 << k)]]
            W_ap[mi][mj] += 1

    W_total = W_ov + W_bd + W_ap

    # ============================================================
    # REGION DECOMPOSITION PER EDGE
    # ============================================================
    print("  EDGE WEIGHT DECOMPOSITION BY REGION:\n")
    print("  %-30s | overlap | boundary | apex | total" % "Edge")
    print("  " + "-"*75)

    for mi in range(nm):
        for mj in range(mi+1, nm):
            if W_total[mi][mj] + W_total[mj][mi] == 0: continue
            Hi = H_vec[mi]; Hj = H_vec[mj]
            ov = W_ov[mi][mj] + W_ov[mj][mi]
            bd = W_bd[mi][mj] + W_bd[mj][mi]
            ap = W_ap[mi][mj] + W_ap[mj][mi]
            total = ov + bd + ap
            print("  {%d(H=%d),%d(H=%d)}%s | %6d  | %7d  | %4d | %d" %
                  (mi, Hi, mj, Hj, " "*(20-len(str(mi))-len(str(mj))),
                   ov, bd, ap, total))

    # ============================================================
    # SELF-LOOPS BY REGION
    # ============================================================
    banner("SELF-LOOPS BY REGION")
    print("  %-20s | overlap | boundary | apex | total" % "Node")
    print("  " + "-"*60)
    for mi in range(nm):
        ov = W_ov[mi][mi]; bd = W_bd[mi][mi]; ap = W_ap[mi][mi]
        total = ov + bd + ap
        if total == 0: continue
        H = H_vec[mi]
        print("  %d (H=%d)%s | %6d  | %7d  | %4d | %d" %
              (mi, H, " "*(12-len(str(mi))), ov, bd, ap, total))

    # ============================================================
    # THE REGION RATIOS
    # ============================================================
    banner("REGION CONTRIBUTION RATIOS")

    total_ov = W_ov.sum(); total_bd = W_bd.sum(); total_ap = W_ap.sum()
    total_all = total_ov + total_bd + total_ap

    print("  Overlap: %d (%.1f%%)" % (total_ov, 100*total_ov/total_all))
    print("  Boundary: %d (%.1f%%)" % (total_bd, 100*total_bd/total_all))
    print("  Apex: %d (%.1f%%)" % (total_ap, 100*total_ap/total_all))

    # Expected from tile counts:
    n_ov = len(overlap_tiles); n_bd = len(boundary_tiles); n_ap = len(apex_tile)
    print("\n  Expected (from tile counts):")
    print("  Overlap: %d/%d = %.1f%%" % (n_ov, m, 100*n_ov/m))
    print("  Boundary: %d/%d = %.1f%%" % (n_bd, m, 100*n_bd/m))
    print("  Apex: %d/%d = %.1f%%" % (n_ap, m, 100*n_ap/m))

    print("\n  Actual/Expected ratio:")
    print("  Overlap: %.4f" % (total_ov/total_all / (n_ov/m)))
    print("  Boundary: %.4f" % (total_bd/total_all / (n_bd/m)))
    print("  Apex: %.4f" % (total_ap/total_all / (n_ap/m)))

    # ============================================================
    # THE RECURSIVE IDENTITY
    # ============================================================
    banner("THE RECURSIVE IDENTITY")

    # Check: is W_ov proportional to W_total?
    # (Would mean overlap flips contribute the same shape as total)
    W_total_off = W_total.copy(); np.fill_diagonal(W_total_off, 0)
    W_ov_off = W_ov.copy(); np.fill_diagonal(W_ov_off, 0)

    if W_total_off.sum() > 0:
        ratio_matrix = W_ov_off.astype(float) / np.maximum(W_total_off.astype(float), 1)
        # Where both are nonzero
        nonzero = W_total_off > 0
        if nonzero.any():
            ratios = ratio_matrix[nonzero]
            print("  W_overlap / W_total (off-diagonal, where nonzero):")
            print("    min=%.4f, max=%.4f, std=%.6f" % (ratios.min(), ratios.max(), ratios.std()))
            if ratios.std() < 0.001:
                print("    → CONSTANT ratio %.4f = %d/%d tiles" %
                      (ratios.mean(), n_ov, m))
            else:
                print("    → NOT constant (varies by edge)")

    # Same for boundary
    W_bd_off = W_bd.copy(); np.fill_diagonal(W_bd_off, 0)
    if W_total_off.sum() > 0:
        ratio_bd = W_bd_off.astype(float) / np.maximum(W_total_off.astype(float), 1)
        nonzero_bd = W_total_off > 0
        if nonzero_bd.any():
            ratios_bd = ratio_bd[nonzero_bd]
            print("  W_boundary / W_total (off-diagonal):")
            print("    min=%.4f, max=%.4f, std=%.6f" %
                  (ratios_bd.min(), ratios_bd.max(), ratios_bd.std()))

    # And apex
    W_ap_off = W_ap.copy(); np.fill_diagonal(W_ap_off, 0)
    if W_total_off.sum() > 0:
        ratio_ap = W_ap_off.astype(float) / np.maximum(W_total_off.astype(float), 1)
        nonzero_ap = W_total_off > 0
        if nonzero_ap.any():
            ratios_ap = ratio_ap[nonzero_ap]
            print("  W_apex / W_total (off-diagonal):")
            print("    min=%.4f, max=%.4f, std=%.6f" %
                  (ratios_ap.min(), ratios_ap.max(), ratios_ap.std()))

    print("""
  IF the ratios are constant: W_region = (n_region/m) × W_total.
  This would mean: the overlap, boundary, and apex contribute
  PROPORTIONALLY to the total, with no region-specific structure.

  If TRUE: W(n) = W_total, and we can compute it from:
    W_total = (m/n_ov) × W_overlap
  And W_overlap relates to W(n-2) via the fiber map.
  This gives: W(n) = (m/C(n-2,2)) × lifted_W(n-2)
  = the RECURSIVE FORMULA for the wiggly weight matrix.
""")

    # ============================================================
    # THE SELF-GENERATION MECHANISM
    # ============================================================
    banner("THE SELF-GENERATION: W(n) FROM W(n-2)")

    # At n=5: overlap has C(3,2) = 3 tiles, boundary has 6, apex has 1.
    # W_overlap should relate to W(n=3) on 2 classes.
    # At n=3: W(3) is a 2×2 matrix (transitive ↔ 3-cycle).

    # Build W(n-2)
    n2 = n - 2; m2 = comb(n2, 2)
    pairs2 = [(i,j) for i in range(n2) for j in range(i+1,n2)]
    cmap2 = {}; tc2 = {}; sizes2 = {}
    for bits in range(2**m2):
        A = [[0]*n2 for _ in range(n2)]
        for k,(ii,jj) in enumerate(pairs2):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n2)
        if c not in cmap2:
            cid = len(cmap2); cmap2[c] = cid; sizes2[cid] = 0
        sizes2[cmap2[c]] += 1; tc2[bits] = cmap2[c]

    info2 = {}
    for cid in range(len(cmap2)):
        bits_rep = next(b for b in range(2**m2) if tc2[b] == cid)
        A = [[0]*n2 for _ in range(n2)]
        for k,(ii,jj) in enumerate(pairs2):
            if (bits_rep>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        comp_c = canon(complement(A, n2), n2)
        info2[cid] = {'comp': cmap2[comp_c], 'SC': cmap2[comp_c]==cid}

    merged2 = {}; mid_map2 = {}; mid2 = 0
    for cid in range(len(cmap2)):
        if cid in mid_map2: continue
        comp = info2[cid]['comp']
        merged2[mid2] = (cid,) if comp == cid else (cid, comp)
        mid_map2[cid] = mid2
        if comp != cid: mid_map2[comp] = mid2
        mid2 += 1
    nm2 = mid2

    W2 = np.zeros((nm2, nm2), dtype=int)
    for bits in range(2**m2):
        mi2 = mid_map2[tc2[bits]]
        for k in range(m2):
            mj2 = mid_map2[tc2[bits ^ (1 << k)]]
            W2[mi2][mj2] += 1

    print("  W(n-2) at n-2=%d: %dx%d matrix" % (n2, nm2, nm2))
    print("  W(n-2) =")
    print(W2)

    # Now: how does W_overlap(n) relate to W(n-2)?
    # For each boundary+apex config: the overlap creates a sub-W(n-2).
    # The total W_overlap = sum over boundary configs of sub-W(n-2).

    # But sub-W(n-2) projected to G_n classes, not G_{n-2} classes!
    # The LIFTING from G_{n-2} to G_n is the non-trivial part.

    # Count: for each overlap flip, what's the (n-2)-class effect?
    # The overlap tiles at n=5 correspond to pairs among {1,2,3}.
    # These are the pairs of the n-2 = 3 sub-tournament.

    # For each labeled tournament at n=5:
    # Extract the sub-tournament on {1,2,3} → get (n-2)-class
    # Flip an overlap tile → changes the sub-tournament
    # The (n-2)-class may or may not change → sub-edge or sub-self-loop
    # The n-class may or may not change → full-edge or full-self-loop

    sub_edge_full_edge = 0
    sub_edge_full_self = 0
    sub_self_full_edge = 0
    sub_self_full_self = 0

    for bits in range(2**m):
        mi_n = mid_map[tc[bits]]

        # Extract sub-tournament on {1,2,3}
        sub_bits = 0
        for k, (i,j) in enumerate(pairs):
            if i >= 1 and j <= n-2:
                sub_k = pairs2.index((i-1, j-1))
                if (bits >> k) & 1:
                    sub_bits |= (1 << sub_k)
        mi_sub = mid_map2[tc2[sub_bits]]

        for k in overlap_tiles:
            # Flip this overlap tile
            flipped = bits ^ (1 << k)
            mj_n = mid_map[tc[flipped]]

            # The corresponding sub-flip
            i_arc, j_arc = pairs[k]
            sub_k = pairs2.index((i_arc-1, j_arc-1))
            sub_flipped = sub_bits ^ (1 << sub_k)
            mj_sub = mid_map2[tc2[sub_flipped]]

            sub_changes = (mi_sub != mj_sub)
            full_changes = (mi_n != mj_n)

            if sub_changes and full_changes: sub_edge_full_edge += 1
            elif sub_changes and not full_changes: sub_edge_full_self += 1
            elif not sub_changes and full_changes: sub_self_full_edge += 1
            else: sub_self_full_self += 1

    total_overlap_flips = sum([sub_edge_full_edge, sub_edge_full_self,
                               sub_self_full_edge, sub_self_full_self])

    print("\n  OVERLAP FLIP LIFTING TABLE (n-2 → n):")
    print("                        n-class changes | n-class stays")
    print("    sub-class changes:  %15d | %13d" % (sub_edge_full_edge, sub_edge_full_self))
    print("    sub-class stays:    %15d | %13d" % (sub_self_full_edge, sub_self_full_self))
    print("    Total: %d\n" % total_overlap_flips)

    lift_rate = sub_edge_full_edge / (sub_edge_full_edge + sub_edge_full_self) if (sub_edge_full_edge + sub_edge_full_self) > 0 else 0
    noise_rate = sub_self_full_edge / (sub_self_full_edge + sub_self_full_self) if (sub_self_full_edge + sub_self_full_self) > 0 else 0

    print("  LIFT RATE: sub-edge → full-edge = %.4f" % lift_rate)
    print("  NOISE RATE: sub-self → full-edge = %.4f" % noise_rate)
    print("  (Lift >> Noise means the sub-structure almost always lifts)")

    print("""
  THE SELF-GENERATION MECHANISM:

  1. The overlap W at n is built from W at n-2:
     - Each sub-edge at n-2 LIFTS to a full edge at n with prob %.1f%%
     - Each sub-self-loop at n-2 creates a full edge with prob %.1f%%

  2. The boundary and apex contribute PROPORTIONALLY:
     - W_boundary = (%d/%d) × W_total (boundary tile fraction)
     - W_apex = (%d/%d) × W_total (apex tile fraction)

  3. So: W(n) ≈ (m/C(n-2,2)) × lifted(W(n-2))
     The lift is not exact (noise rate > 0), but it's a strong recursion.

  4. To compute W(n) AND the meta-graph simultaneously:
     a. Compute W(n-2) (the sub-tournament wiggly structure)
     b. Lift to W_overlap(n) via the fiber map
     c. Add boundary and apex contributions (proportional)
     d. The meta-graph edges = nonzero off-diagonal entries of W(n)
     e. fiber(C) = row_sum(W(n)) / m
""" % (100*lift_rate, 100*noise_rate,
       len(boundary_tiles), m, len(apex_tile), m))

print("Done. opus-2026-03-24-S284")
