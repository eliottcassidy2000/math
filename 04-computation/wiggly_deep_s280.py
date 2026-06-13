#!/usr/bin/env python3
"""
wiggly_deep_s280.py -- opus-2026-03-24-S280

DEEP WIGGLY META-GRAPH: COMBINING ALL TECHNIQUES

Compute the wiggly weight matrix W and extract:
1. Fiber formula: fiber = W_row/m
2. Transition matrix P = W / (fiber × m)
3. H as eigenvector of W
4. Self-loop fraction per node = neutral arc probability
5. Wiggly Betti numbers (clique complex of weighted graph)
6. The base weight formula (from perfect factorization)
7. Wiggly profiles per class (which arcs are neutral)
8. Cross-tabulation: wiggly weight vs H, score, SC/NS

PUSH TO n=6 (34 merged nodes, feasible).

Author: opus-2026-03-24-S280
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

for n in [5, 6]:
    banner("DEEP WIGGLY ANALYSIS AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

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
                     'SC': cmap[comp_c] == cid, 'size': sizes[cid],
                     'aut': factorial(n)//sizes[cid],
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

    # Build W matrix
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(2**m):
        mi = mid_map[tc[bits]]
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            W[mi][mj] += 1

    # Node properties
    H_vec = np.array([info[merged[mi][0]]['H'] for mi in range(nm)], dtype=float)
    fiber_vec = np.array([sum(sizes[cid] for cid in merged[mi]) for mi in range(nm)], dtype=float)
    aut_vec = np.array([info[merged[mi][0]]['aut'] for mi in range(nm)], dtype=float)
    SC_vec = np.array([1.0 if info[merged[mi][0]]['SC'] else 0.0 for mi in range(nm)])

    # Unweighted adjacency
    A_adj = np.zeros((nm, nm));
    for mi in range(nm):
        for mj in range(nm):
            if mi != mj and W[mi][mj] > 0: A_adj[mi][mj] = 1.0
    deg = A_adj.sum(axis=1)
    ne = int(A_adj.sum()) // 2

    # ============================================================
    # 1. FIBER FROM W
    # ============================================================
    print("  1. FIBER = W_row / m:")
    fiber_check = np.array([sum(W[mi]) // m for mi in range(nm)])
    print("     All match: %s" % np.all(fiber_check == fiber_vec))

    # ============================================================
    # 2. W EIGENVECTORS + H CORRELATION
    # ============================================================
    # Use W (off-diagonal, symmetric)
    W_sym = W.astype(float).copy()
    np.fill_diagonal(W_sym, 0)

    w_evals, w_evecs = np.linalg.eigh(W_sym)
    idx = np.argsort(-w_evals)
    w_evals = w_evals[idx]; w_evecs = w_evecs[:, idx]

    print("\n  2. W EIGENVALUES (top 5): %s" % [round(e, 2) for e in w_evals[:5]])

    # H correlation
    H_c = H_vec - np.mean(H_vec)
    if np.linalg.norm(H_c) > 0:
        H_n = H_c / np.linalg.norm(H_c)
        projs = (w_evecs.T @ H_n)**2
        print("     H on eigvec 0: %.4f, eigvec 1: %.4f" % (projs[0], projs[1]))
        best_proj = max(range(min(nm, 10)), key=lambda i: projs[i])
        print("     Best H projection: eigvec %d (%.4f of variance)" % (best_proj, projs[best_proj]))

    # Fiber correlation
    F_c = fiber_vec - np.mean(fiber_vec)
    if np.linalg.norm(F_c) > 0:
        F_n = F_c / np.linalg.norm(F_c)
        f_projs = (w_evecs.T @ F_n)**2
        best_f = max(range(min(nm, 10)), key=lambda i: f_projs[i])
        print("     Best fiber projection: eigvec %d (%.4f)" % (best_f, f_projs[best_f]))

    # ============================================================
    # 3. TRANSITION MATRIX + MARKOV
    # ============================================================
    row_sums = W.sum(axis=1).astype(float)
    P = np.diag(1.0 / np.maximum(row_sums, 1)) @ W.astype(float)
    pi_stat = row_sums / row_sums.sum()

    m_evals = np.sort(np.real(np.linalg.eigvals(P)))[::-1]
    gap = 1 - m_evals[1] if len(m_evals) > 1 else 0

    print("\n  3. MARKOV CHAIN:")
    print("     Spectral gap: %.4f (cf 2/n = %.4f)" % (gap, 2.0/n))
    print("     π ∝ fiber: %s" % np.allclose(pi_stat, fiber_vec/fiber_vec.sum(), atol=1e-10))

    # ============================================================
    # 4. SELF-LOOP ANALYSIS
    # ============================================================
    sl_frac = np.array([W[mi][mi] / row_sums[mi] if row_sums[mi] > 0 else 0 for mi in range(nm)])
    avg_neutral = np.array([W[mi][mi] / fiber_vec[mi] if fiber_vec[mi] > 0 else 0 for mi in range(nm)])

    print("\n  4. SELF-LOOP STRUCTURE:")
    print("     Nodes with SL: %d/%d" % (sum(1 for mi in range(nm) if W[mi][mi] > 0), nm))
    print("     Avg SL fraction: %.4f" % np.mean(sl_frac))
    print("     Avg neutral arcs/tiling: %.2f / %d" % (np.mean(avg_neutral), m))
    print("     Corr(SL_frac, H): %.4f" % np.corrcoef(sl_frac, H_vec)[0,1])
    print("     Corr(SL_frac, fiber): %.4f" % np.corrcoef(sl_frac, fiber_vec)[0,1])

    # ============================================================
    # 5. BASE WEIGHT DECOMPOSITION
    # ============================================================
    # From S277: weight(edge, range r) = base_weight × (n-r)
    # base_weight = W[mi][mj] / (n-1) for range-1 arcs at this edge
    # Let's compute base weights for all edges

    W_range = [np.zeros((nm, nm), dtype=int) for _ in range(n)]
    for bits in range(2**m):
        mi = mid_map[tc[bits]]
        for k in range(m):
            r = pairs[k][1] - pairs[k][0]
            mj = mid_map[tc[bits ^ (1 << k)]]
            W_range[r][mi][mj] += 1

    base_weights = []
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if A_adj[mi][mj] == 0: continue
            # Base weight = weight at range 1 / (n-1)
            w_r1 = W_range[1][mi][mj] + W_range[1][mj][mi]
            base = w_r1 / (n - 1)
            base_weights.append(base)

    print("\n  5. BASE WEIGHTS:")
    bw_counter = Counter(int(b) for b in base_weights)
    print("     Distribution: %s" % dict(sorted(bw_counter.items())))
    print("     All multiples of 2^{n-2}=%d: %s" %
          (2**(n-2), all(int(b) % 2**(n-2) == 0 for b in base_weights)))

    # ============================================================
    # 6. THE H-WEIGHTED STRUCTURE
    # ============================================================
    # For each edge: compute dH, base_weight, and check correlations
    edges_data = []
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if A_adj[mi][mj] == 0: continue
            dH = abs(H_vec[mi] - H_vec[mj])
            total_w = W[mi][mj] + W[mj][mi]
            w_r1 = W_range[1][mi][mj] + W_range[1][mj][mi]
            base = w_r1 / (n-1) if n > 1 else 0
            edges_data.append({'mi': mi, 'mj': mj, 'dH': dH,
                              'total_w': total_w, 'base': base,
                              'fiber_prod': fiber_vec[mi] * fiber_vec[mj]})

    dH_arr = np.array([e['dH'] for e in edges_data])
    tw_arr = np.array([e['total_w'] for e in edges_data])
    bw_arr = np.array([e['base'] for e in edges_data])
    fp_arr = np.array([e['fiber_prod'] for e in edges_data])

    print("\n  6. EDGE CORRELATIONS:")
    print("     Corr(base_weight, dH): %.4f" % np.corrcoef(bw_arr, dH_arr)[0,1])
    print("     Corr(base_weight, fiber_prod): %.4f" % np.corrcoef(bw_arr, fp_arr)[0,1])
    print("     Corr(total_weight, fiber_prod): %.4f" % np.corrcoef(tw_arr, fp_arr)[0,1])

    # Is base_weight = fiber_prod / something?
    ratios = bw_arr / np.maximum(fp_arr, 1)
    print("     base/fiber_prod: min=%.6f, max=%.6f, std=%.6f" %
          (ratios.min(), ratios.max(), ratios.std()))
    if ratios.std() < 0.001:
        print("     → base_weight = %.6f × fiber_prod (CONSTANT!)" % ratios.mean())
    else:
        print("     → NOT a constant ratio")

    # ============================================================
    # 7. THE MAGIC FORMULA SEARCH
    # ============================================================
    print("\n  7. FORMULA SEARCH: base_weight = f(fiber_i, fiber_j, H_i, H_j, m)?")

    # Try: base = c × min(fiber_i, fiber_j)?
    min_fibers = np.array([min(fiber_vec[e['mi']], fiber_vec[e['mj']]) for e in edges_data])
    r_min = bw_arr / np.maximum(min_fibers, 1)
    print("     base/min(fiber): min=%.4f, max=%.4f, std=%.4f" %
          (r_min.min(), r_min.max(), r_min.std()))

    # Try: base = c × (fiber_i + fiber_j)?
    sum_fibers = np.array([fiber_vec[e['mi']] + fiber_vec[e['mj']] for e in edges_data])
    r_sum = bw_arr / np.maximum(sum_fibers, 1)
    print("     base/sum(fiber): min=%.4f, max=%.4f, std=%.4f" %
          (r_sum.min(), r_sum.max(), r_sum.std()))

    # Try: base = c × sqrt(fiber_i × fiber_j)?
    sqrt_fibers = np.sqrt(fp_arr)
    r_sqrt = bw_arr / np.maximum(sqrt_fibers, 1)
    print("     base/sqrt(fiber_prod): min=%.4f, max=%.4f, std=%.4f" %
          (r_sqrt.min(), r_sqrt.max(), r_sqrt.std()))

    # ============================================================
    # 8. SUMMARY TABLE
    # ============================================================
    banner("SUMMARY TABLE (n=%d)" % n)
    print("  %-3s | %-3s | %-4s | %-6s | %-6s | %-8s | %-8s | %-8s" %
          ("mid", "H", "Aut", "fiber", "deg", "SL_frac", "neutral", "score"))
    print("  " + "-"*65)
    for mi in range(nm):
        H = int(H_vec[mi]); aut = int(aut_vec[mi]); fib = int(fiber_vec[mi])
        d = int(deg[mi]); sl = sl_frac[mi]; neut = avg_neutral[mi]
        sc = info[merged[mi][0]]['score']
        print("  %-3d | %-3d | %-4d | %-6d | %-6d | %-8.4f | %-8.2f | %s" %
              (mi, H, aut, fib, d, sl, neut, sc))

    print()

print("Done. opus-2026-03-24-S280")
