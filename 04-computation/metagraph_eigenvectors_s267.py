#!/usr/bin/env python3
"""
metagraph_eigenvectors_s267.py -- opus-2026-03-23-S267

EIGENVECTORS OF THE MERGED META-GRAPH G_n/Z_2

The adjacency matrix of G_n/Z_2 encodes ALL connectivity information.
Its eigenvectors reveal the fundamental modes of tournament space.

From S167: at n=5, the Markov chain eigenvalues are ALL multiples of 1/5.
Spectral gap = 2/n at n=4,5.

This script computes:
1. Full spectrum of G_n/Z_2 adjacency matrix at n=5,6
2. Eigenvectors and what they correlate with (H, score, c3, SC/NS)
3. The Markov chain (transition matrix) spectrum
4. Whether eigenvectors align with the SPINE/HALO/SEA decomposition
5. Whether H is an eigenvector (or close to one)

Author: opus-2026-03-23-S267
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict
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

for n in [5, 6]:
    banner("EIGENVECTORS OF G_%d/Z_2" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build iso classes
    cmap = {}; sizes = {}; reps = {}; tc = {}
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

    # Class info
    info = {}
    for cid in range(nc):
        A = reps[cid]
        comp_c = canon(complement(A, n), n)
        info[cid] = {
            'H': Hcount(A, n), 'c3': c3count(A, n),
            'score': tuple(sorted(sum(A[i]) for i in range(n))),
            'SC': cmap[comp_c] == cid, 'comp': cmap[comp_c],
            'size': sizes[cid], 'aut': factorial(n) // sizes[cid]
        }

    # Merged graph
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    # F matrix
    F = np.zeros((nc, nc), dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m): F[ci][tc[bits^(1<<k)]] += 1

    # Merged adjacency (weighted)
    W = np.zeros((nm, nm))
    for ci in range(nc):
        mi = mid_map[ci]
        for cj in range(nc):
            if F[ci][cj] > 0:
                mj = mid_map[cj]
                W[mi][mj] += F[ci][cj]

    # Unweighted adjacency
    A_adj = (W > 0).astype(float)
    np.fill_diagonal(A_adj, 0)

    # Node properties for merged nodes
    H_vec = np.array([info[merged[mi][0]]['H'] for mi in range(nm)], dtype=float)
    c3_vec = np.array([info[merged[mi][0]]['c3'] for mi in range(nm)], dtype=float)
    SC_vec = np.array([1.0 if info[merged[mi][0]]['SC'] else 0.0 for mi in range(nm)])
    fiber_vec = np.array([info[merged[mi][0]]['size'] for mi in range(nm)], dtype=float)
    # Fiber for merged: sum of sizes
    for mi in range(nm):
        fiber_vec[mi] = sum(sizes[cid] for cid in merged[mi])

    deg_vec = A_adj.sum(axis=1)

    print("  V_merged = %d, E_merged = %d\n" % (nm, int(A_adj.sum()) // 2))

    # ============================================================
    # UNWEIGHTED ADJACENCY SPECTRUM
    # ============================================================
    eigenvalues, eigenvectors = np.linalg.eigh(A_adj)
    idx = np.argsort(-eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    print("  ADJACENCY SPECTRUM:")
    for i, ev in enumerate(eigenvalues):
        print("    λ_%d = %.6f" % (i, ev))

    # Spectral gap
    if len(eigenvalues) >= 2:
        gap = eigenvalues[0] - eigenvalues[1]
        print("\n  Spectral gap: λ_0 - λ_1 = %.6f" % gap)
        print("  λ_0 = %.6f (should be max degree)" % eigenvalues[0])

    # ============================================================
    # MARKOV CHAIN (RANDOM WALK) SPECTRUM
    # ============================================================
    D_inv = np.diag(1.0 / np.maximum(deg_vec, 1))
    M = D_inv @ A_adj  # transition matrix

    m_eigenvalues = np.linalg.eigvals(M)
    m_eigenvalues = np.sort(np.real(m_eigenvalues))[::-1]

    print("\n  MARKOV (RANDOM WALK) SPECTRUM:")
    for i, ev in enumerate(m_eigenvalues):
        print("    μ_%d = %.6f" % (i, ev))

    if len(m_eigenvalues) >= 2:
        m_gap = 1 - m_eigenvalues[1]
        print("\n  Markov spectral gap: 1 - μ_1 = %.6f" % m_gap)
        print("  Compare 2/n = %.6f" % (2.0/n))

    # ============================================================
    # EIGENVECTOR CORRELATIONS
    # ============================================================
    banner("EIGENVECTOR CORRELATIONS (n=%d)" % n)

    # Normalize H, c3, SC, fiber to zero mean
    features = {
        'H': H_vec, 'c3': c3_vec, 'SC': SC_vec,
        'fiber': fiber_vec, 'degree': deg_vec, '1': np.ones(nm)
    }

    print("  Correlation of eigenvectors with node properties:\n")
    print("  %-8s" % "λ_i", end="")
    for fname in ['H', 'c3', 'SC', 'fiber', 'degree', '1']:
        print(" | %-8s" % fname, end="")
    print()
    print("  " + "-"*70)

    for i in range(min(nm, 10)):
        ev = eigenvectors[:, i]
        print("  λ=%-6.3f" % eigenvalues[i], end="")
        for fname, fvec in features.items():
            # Pearson correlation
            if np.std(ev) > 1e-10 and np.std(fvec) > 1e-10:
                corr = np.corrcoef(ev, fvec)[0, 1]
            else:
                corr = 0
            print(" | %-+8.4f" % corr, end="")
        print()

    # ============================================================
    # IS H AN EIGENVECTOR?
    # ============================================================
    banner("IS H AN EIGENVECTOR? (n=%d)" % n)

    # Check: A_adj @ H = λ * H?
    AH = A_adj @ H_vec
    # If H is an eigenvector: AH / H should be constant
    ratios = AH / np.maximum(H_vec, 1)
    print("  A × H / H (should be constant if H is eigenvector):")
    for mi in range(nm):
        print("    mid=%d H=%d: (A×H)/H = %.4f" % (mi, int(H_vec[mi]), ratios[mi]))

    print("\n  Std of ratios: %.6f (0 = perfect eigenvector)" % np.std(ratios))

    # Project H onto eigenbasis
    H_centered = H_vec - np.mean(H_vec)
    H_norm = H_centered / np.linalg.norm(H_centered)
    projections = eigenvectors.T @ H_norm
    print("\n  H projection onto eigenbasis:")
    for i in range(min(nm, 10)):
        print("    |⟨v_%d, H⟩|² = %.6f (λ=%.3f)" %
              (i, projections[i]**2, eigenvalues[i]))

    # ============================================================
    # WEIGHTED SPECTRUM (using F matrix weights)
    # ============================================================
    banner("WEIGHTED SPECTRUM (n=%d)" % n)

    W_sym = (W + W.T) / 2
    np.fill_diagonal(W_sym, 0)

    w_eigenvalues = np.sort(np.linalg.eigvalsh(W_sym))[::-1]
    print("  Weighted adjacency spectrum:")
    for i, ev in enumerate(w_eigenvalues[:min(nm, 10)]):
        print("    λ_w_%d = %.2f" % (i, ev))

    # ============================================================
    # LAPLACIAN SPECTRUM
    # ============================================================
    banner("LAPLACIAN SPECTRUM (n=%d)" % n)

    D = np.diag(deg_vec)
    L = D - A_adj
    L_eigenvalues = np.sort(np.linalg.eigvalsh(L))

    print("  Laplacian spectrum:")
    for i, ev in enumerate(L_eigenvalues[:min(nm, 10)]):
        print("    μ_L_%d = %.6f" % (i, ev))

    # Algebraic connectivity = 2nd smallest Laplacian eigenvalue
    if len(L_eigenvalues) >= 2:
        alg_conn = L_eigenvalues[1]
        print("\n  Algebraic connectivity (Fiedler value): %.6f" % alg_conn)
        print("  Graph IS connected: %s" % (alg_conn > 1e-10))

    # Fiedler vector = eigenvector of 2nd smallest Laplacian eigenvalue
    L_vals, L_vecs = np.linalg.eigh(L)
    fiedler = L_vecs[:, 1]

    print("\n  Fiedler vector (determines graph bipartition):")
    for mi in range(nm):
        sc = "SC" if info[merged[mi][0]]['SC'] else "NS"
        print("    mid=%d H=%d %s: Fiedler=%.4f" %
              (mi, int(H_vec[mi]), sc, fiedler[mi]))

    # Does the Fiedler vector separate SC from NS?
    sc_fiedler = [fiedler[mi] for mi in range(nm) if info[merged[mi][0]]['SC']]
    ns_fiedler = [fiedler[mi] for mi in range(nm) if not info[merged[mi][0]]['SC']]
    if sc_fiedler and ns_fiedler:
        print("\n  SC Fiedler values: mean=%.4f, range=[%.4f, %.4f]" %
              (np.mean(sc_fiedler), min(sc_fiedler), max(sc_fiedler)))
        print("  NS Fiedler values: mean=%.4f, range=[%.4f, %.4f]" %
              (np.mean(ns_fiedler), min(ns_fiedler), max(ns_fiedler)))
        overlap = (min(max(sc_fiedler), max(ns_fiedler)) -
                   max(min(sc_fiedler), min(ns_fiedler)))
        print("  SC-NS overlap: %.4f (0 = perfect separation)" % max(0, overlap))

print("\nDone. opus-2026-03-23-S267")
