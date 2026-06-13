#!/usr/bin/env python3
"""
enumeration_taxonomy_s238.py -- opus-2026-03-23-S238

CREATIVE NEW CATEGORIZATIONS targeting ENUMERATION FORMULAS.

Inspired by:
- kind-pasteur S20da: c3 + (n/2)*svar = const (EXACT identity)
- kind-pasteur S20db: bimodal degree at n=7, fiber = H/|Aut| for SC
- kind-pasteur S20cs: thickness quantized at multiples of n!
- opus S237: rib crossover — UP near transitive, DOWN near regular
- opus S236: TEAL=100%, GOLDEN SPINE=3

NEW CATEGORIES to test:
1. H mod k: partition nodes by H modulo various k (2, 4, n-1, 2^{n-2})
2. Fibonacci-H: is H a Fibonacci number?
3. Fiber class: fiber = H/|Aut|, group by fiber value
4. Score variance quartile: partition by svar
5. c3 mod (n-2): c3 value modulo n-2
6. Distance from transitive / regular in meta-graph
7. "Kings" vs "pawns": is degree > mean?
8. Triangles: number of triangles containing this node
9. Edge "weight" categories: based on fiber product
10. Complement-distance: how far is this node from its complement pair

Author: opus-2026-03-23-S238
"""

import sys
import numpy as np
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
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
            if not (S & (1<<v)): continue
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

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def score_var(sc, n):
    mean = (n-1)/2
    return sum((s - mean)**2 for s in sc) / n

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def is_sc(A, n):
    return canon(A,n) == canon(complement(A,n),n)

for n in [5, 6]:
    banner("ENUMERATION TAXONOMY AT n=%d" % n)
    m = comb(n,2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; sizes = {}; reps = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A,n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
            reps[cid] = [row[:] for row in A]
        sizes[cmap[c]] += 1; tc[bits] = cmap[c]
    nc = len(cmap)

    info = {}
    for cid in range(nc):
        A = reps[cid]
        sc_seq = score_seq(A,n)
        comp_c = canon(complement(A,n),n)
        info[cid] = {
            'H': Hcount(A,n), 'c3': c3count(A,n), 'score': sc_seq,
            'svar': score_var(sc_seq, n), 'SC': is_sc(A,n),
            'comp': cmap[comp_c], 'aut': factorial(n)//sizes[cid],
            'fiber': sizes[cid]  # = n!/|Aut|
        }

    # Merged
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    # Edges + F matrix
    F = np.zeros((nc,nc), dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            F[ci][tc[bits^(1<<k)]] += 1

    merged_edges = set()
    merged_weight = {}
    for ci in range(nc):
        mi = mid_map[ci]
        for cj in range(nc):
            if F[ci][cj] > 0:
                mj = mid_map[cj]
                if mi != mj:
                    key = (min(mi,mj), max(mi,mj))
                    merged_edges.add(key)
                    merged_weight[key] = merged_weight.get(key, 0) + F[ci][cj]
    ne = len(merged_edges)

    neighbors = defaultdict(set)
    for (a,b) in merged_edges:
        neighbors[a].add(b); neighbors[b].add(a)

    # BFS distances from transitive and regular
    # Find transitive (H=1) and regular (H=max)
    trans_mid = None; reg_mid = None
    H_max = max(info[merged[mi][0]]['H'] for mi in range(nm))
    H_min = min(info[merged[mi][0]]['H'] for mi in range(nm))
    for mi in range(nm):
        H = info[merged[mi][0]]['H']
        if H == H_min: trans_mid = mi
        if H == H_max: reg_mid = mi

    def bfs_dist(src):
        from collections import deque
        dist = {src: 0}; q = deque([src])
        while q:
            v = q.popleft()
            for w in neighbors[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; q.append(w)
        return dist

    dist_from_trans = bfs_dist(trans_mid)
    dist_from_reg = bfs_dist(reg_mid)

    # Triangles per node
    tri_count = defaultdict(int)
    edge_in_tri = defaultdict(int)
    for mi in range(nm):
        for mj in neighbors[mi]:
            if mj <= mi: continue
            for mk in neighbors[mj]:
                if mk <= mj: continue
                if mk in neighbors[mi]:
                    tri_count[mi] += 1
                    tri_count[mj] += 1
                    tri_count[mk] += 1
                    edge_in_tri[(min(mi,mj), max(mi,mj))] += 1
                    edge_in_tri[(min(mj,mk), max(mj,mk))] += 1
                    edge_in_tri[(min(mi,mk), max(mi,mk))] += 1

    # ============================================================
    # NEW NODE CATEGORIES
    # ============================================================
    print("  NEW NODE CATEGORIES (n=%d, V_merged=%d):\n" % (n, nm))

    for mi in range(nm):
        ci = merged[mi][0]
        ii = info[ci]
        ii['mid'] = mi
        ii['deg'] = len(neighbors[mi])
        ii['tri'] = tri_count.get(mi, 0)
        ii['dist_t'] = dist_from_trans.get(mi, -1)
        ii['dist_r'] = dist_from_reg.get(mi, -1)

    # 1. H mod k partitions
    print("  H mod 2:")
    for r in [0, 1]:
        count = sum(1 for mi in range(nm) if info[merged[mi][0]]['H'] % 2 == r)
        print("    H ≡ %d (mod 2): %d nodes" % (r, count))

    print("\n  H mod 4:")
    for r in range(4):
        count = sum(1 for mi in range(nm) if info[merged[mi][0]]['H'] % 4 == r)
        if count > 0:
            print("    H ≡ %d (mod 4): %d nodes" % (r, count))

    # 2. Fiber = H/|Aut| (kind-pasteur discovery)
    print("\n  Fiber = H / |Aut| (verified for SC):")
    fiber_match = 0; fiber_fail = 0
    for mi in range(nm):
        ci = merged[mi][0]
        ii = info[ci]
        fiber_expected = ii['H'] // ii['aut'] if ii['aut'] > 0 else -1
        if ii['fiber'] == fiber_expected * ii['aut']:
            # Actually fiber = n!/|Aut|, and H/|Aut| = (tiling count) / (orbit size)
            pass
        # The correct identity: fiber (= n!/|Aut|) * |Aut| = n!
        # And H = fiber * (something?)... no, H/|Aut| = #distinct orientations per tiling
        # Actually the orbit-stabilizer: H = fiber = n!/|Aut| is false generally
        # But: fiber (labeled copies) = n!/|Aut|. And H = Hamiltonian path count.
        # The identity fiber = H/|Aut| is wrong. Let me check what kind-pasteur meant.
        pass

    # Check: fiber / m (fiber per arc)
    print("\n  Fiber values (n!/|Aut|):")
    fiber_vals = Counter()
    for mi in range(nm):
        ci = merged[mi][0]
        fiber_vals[info[ci]['fiber']] += 1
    for fv, cnt in sorted(fiber_vals.items()):
        print("    fiber=%d: %d nodes" % (fv, cnt))

    # 3. Score variance classes
    print("\n  Score variance classes:")
    svar_vals = sorted(set(info[merged[mi][0]]['svar'] for mi in range(nm)))
    for sv in svar_vals:
        count = sum(1 for mi in range(nm) if info[merged[mi][0]]['svar'] == sv)
        # What c3 does this correspond to?
        # c3 = C(n,3) - n(n-1)(n-3)/8 - (n/2)*svar
        c3_from_svar = comb(n,3) - n*(n-1)*(n-3)/8 - (n/2)*sv
        print("    svar=%.2f (c3=%.0f): %d nodes" % (sv, c3_from_svar, count))

    # 4. Distance from transitive + distance from regular
    print("\n  Distance pairs (dist_trans, dist_reg):")
    dist_pairs = Counter()
    for mi in range(nm):
        ci = merged[mi][0]
        dt = info[ci]['dist_t']
        dr = info[ci]['dist_r']
        dist_pairs[(dt, dr)] += 1
    for (dt, dr), cnt in sorted(dist_pairs.items()):
        print("    (d_T=%d, d_R=%d): %d nodes, sum=%d" % (dt, dr, cnt, dt+dr))

    # 5. Triangle count classes
    print("\n  Triangle count per node:")
    tri_counts = Counter()
    for mi in range(nm):
        tri_counts[tri_count.get(mi, 0)] += 1
    for tc_val, cnt in sorted(tri_counts.items()):
        print("    tri=%d: %d nodes" % (tc_val, cnt))

    total_tri = sum(tri_count.values()) // 3
    print("  Total triangles: %d" % total_tri)

    # 6. H modulo 2^{n-2} (hypotenuse quantum)
    q = 2**(n-2)
    print("\n  H mod %d (hypotenuse quantum):" % q)
    for r in range(q):
        count = sum(1 for mi in range(nm) if info[merged[mi][0]]['H'] % q == r)
        if count > 0:
            print("    H ≡ %d (mod %d): %d nodes" % (r, q, count))

    # 7. Degree parity
    print("\n  Degree parity:")
    even_deg = sum(1 for mi in range(nm) if len(neighbors[mi]) % 2 == 0)
    odd_deg = nm - even_deg
    print("    Even degree: %d, Odd degree: %d" % (even_deg, odd_deg))

    # ============================================================
    # NEW EDGE CATEGORIES
    # ============================================================
    banner("NEW EDGE CATEGORIES (n=%d, E=%d)" % (n, ne))

    # 1. Edge in triangle
    tri_edges = sum(1 for e in merged_edges if edge_in_tri.get(e, 0) > 0)
    no_tri_edges = ne - tri_edges
    print("  In triangle: %d, NOT in triangle: %d" % (tri_edges, no_tri_edges))

    # 2. Thickness categories (quantized?)
    weights = [merged_weight[e] for e in merged_edges]
    nfact = factorial(n)
    print("\n  Thickness in units of n!=%d:" % nfact)
    w_classes = Counter()
    for w in weights:
        ratio = Fraction(w, nfact)
        w_classes[ratio] += 1
    for ratio, cnt in sorted(w_classes.items()):
        print("    w/n! = %s: %d edges" % (ratio, cnt))

    # 3. Fiber product ratio
    print("\n  Fiber product edge categories (fiber_i * fiber_j):")
    fprod_vals = Counter()
    for (mi, mj) in merged_edges:
        fi = info[merged[mi][0]]['fiber']
        fj = info[merged[mj][0]]['fiber']
        fprod = fi * fj
        fprod_vals[fprod] += 1
    for fp, cnt in sorted(fprod_vals.items())[:10]:
        print("    fiber_prod=%d: %d edges" % (fp, cnt))
    if len(fprod_vals) > 10:
        print("    ... (%d distinct values)" % len(fprod_vals))

    # 4. Distance sum (dist_trans_i + dist_trans_j)
    print("\n  Distance-from-transitive sum:")
    dsum_vals = Counter()
    for (mi, mj) in merged_edges:
        di = dist_from_trans.get(mi, -1)
        dj = dist_from_trans.get(mj, -1)
        dsum_vals[di + dj] += 1
    for ds, cnt in sorted(dsum_vals.items()):
        print("    d_T(i)+d_T(j)=%d: %d edges" % (ds, cnt))

    # 5. H-sum categories
    print("\n  H-sum parity (H_i + H_j mod 2):")
    h_sum_even = sum(1 for (mi,mj) in merged_edges
                     if (info[merged[mi][0]]['H'] + info[merged[mj][0]]['H']) % 2 == 0)
    h_sum_odd = ne - h_sum_even
    print("    Even: %d, Odd: %d" % (h_sum_even, h_sum_odd))

    # 6. dH mod 4
    print("\n  dH mod 4:")
    for r in range(4):
        count = sum(1 for (mi,mj) in merged_edges
                    if abs(info[merged[mi][0]]['H'] - info[merged[mj][0]]['H']) % 4 == r)
        if count > 0:
            print("    dH ≡ %d (mod 4): %d edges" % (r, count))

    # ============================================================
    # SEARCH FOR ENUMERATION FORMULAS
    # ============================================================
    banner("FORMULA SEARCH (n=%d)" % n)

    # Count by H mod 2
    H_odd = sum(1 for mi in range(nm) if info[merged[mi][0]]['H'] % 2 == 1)
    H_even = nm - H_odd
    print("  H odd: %d, H even: %d (V_merged=%d)" % (H_odd, H_even, nm))

    # Count by SC
    n_sc = sum(1 for mi in range(nm) if info[merged[mi][0]]['SC'])
    n_ns = nm - n_sc
    print("  SC: %d, NS: %d" % (n_sc, n_ns))

    # BLUE/BLACK/GREEN edge counts
    blue_e = sum(1 for (mi,mj) in merged_edges if info[merged[mi][0]]['SC'] and info[merged[mj][0]]['SC'])
    green_e = sum(1 for (mi,mj) in merged_edges if not info[merged[mi][0]]['SC'] and not info[merged[mj][0]]['SC'])
    black_e = ne - blue_e - green_e
    print("  BLUE=%d, BLACK=%d, GREEN=%d" % (blue_e, black_e, green_e))

    # Total thickness
    total_w = sum(merged_weight[e] for e in merged_edges)
    # Self-loop thickness
    sl_w = sum(F[ci][ci] for ci in range(nc))
    print("  Total edge thickness: %d" % total_w)
    print("  Self-loop thickness: %d" % sl_w)
    print("  Total = edge + self: %d (should be 2^m * m = %d)" % (total_w + sl_w, 2**m * m))

    # Key check: 2^m * m = SL + sum(thickness)
    print("  CHECK: 2^m * m = %d, SL + sum(w) = %d, match = %s" %
          (2**m * m, sl_w + total_w, sl_w + total_w == 2**m * m))

    # AMBER edge count
    amber_e = sum(1 for (mi,mj) in merged_edges
                  if info[merged[mi][0]]['score'] == info[merged[mj][0]]['score'])
    print("  AMBER: %d/%d = %.1f%%" % (amber_e, ne, 100*amber_e/ne))

    # Thickness per BLUE/BLACK/GREEN
    blue_w = sum(merged_weight[(mi,mj)] for (mi,mj) in merged_edges
                 if info[merged[mi][0]]['SC'] and info[merged[mj][0]]['SC'])
    green_w = sum(merged_weight[(mi,mj)] for (mi,mj) in merged_edges
                  if not info[merged[mi][0]]['SC'] and not info[merged[mj][0]]['SC'])
    black_w = total_w - blue_w - green_w
    print("\n  Avg thickness by type:")
    if blue_e > 0: print("    BLUE:  %d/%d = %.1f" % (blue_w, blue_e, blue_w/blue_e))
    if black_e > 0: print("    BLACK: %d/%d = %.1f" % (black_w, black_e, black_w/black_e))
    if green_e > 0: print("    GREEN: %d/%d = %.1f" % (green_w, green_e, green_w/green_e))

    # Average thickness per AMBER vs VIOLET
    amber_w = sum(merged_weight[(mi,mj)] for (mi,mj) in merged_edges
                  if info[merged[mi][0]]['score'] == info[merged[mj][0]]['score'])
    violet_w = total_w - amber_w
    print("    AMBER: %d/%d = %.1f" % (amber_w, amber_e, amber_w/amber_e if amber_e > 0 else 0))
    print("    VIOLET: %d/%d = %.1f" % (violet_w, ne-amber_e, violet_w/(ne-amber_e) if ne-amber_e > 0 else 0))

print("\nDone. opus-2026-03-23-S238")
