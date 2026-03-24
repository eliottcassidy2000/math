#!/usr/bin/env python3
"""
eigenvalue_deep_s268.py — Deep eigenvalue analysis of G_n/Z_2
opus-2026-03-23-S268

Building on S267 results (n=5,6), extend to n=7,8 and investigate:

1. Spectral gap 2/n conjecture: verified at n=5 (0.480 vs 0.400), n=6 (0.331 vs 0.333)
2. H ≈ 2nd eigenvector: 72% at n=5, 79% at n=6 — does it improve?
3. Spine/ribs/sea spectral structure
4. Lie-theoretic interpretation of eigenvalues (roots of A_{n-1}?)
5. Spectral density distribution — does it converge to Wigner semicircle?
6. Cheeger constant and expansion
7. Ramanujan-ness: is λ_1 ≤ 2√(d-1)?
"""

import sys, time
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict
import subprocess

sys.stdout.reconfigure(line_buffering=True)

def build_merged_metagraph(n):
    """Build G_n/Z_2 using nauty. Returns adjacency matrix + metadata."""
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1,n)]

    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]

    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m-1-k)): a[i][j] = 1
            else: a[j][i] = 1
        return a

    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

    def c3(a):
        c = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if a[i][j] and a[j][k] and a[k][i]: c += 1
                    if a[i][k] and a[k][j] and a[j][i]: c += 1
        return c

    def H(a):
        dp = {}
        for v in range(n): dp[(1<<v,v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: dp[(S|(1<<u),u)] = dp.get((S|(1<<u),u),0) + val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

    def cp(a):
        sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[sc[v]].append(v)
        gs = [sg[s] for s in sorted(sg.keys())]
        best = None; cnt = 0
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gp(gs2[1:]): yield list(p)+rest
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
            cnt += 1
            if cnt > 500000: break
        return best

    cls = []
    hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid': i, 'bits': bits, 'adj': adj, 'score': s,
                    'c3': cc, 'H': h, 'sc': sc, 'comp_canon': cc2, 'canon': canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)

    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)

    edges = set()
    for cl in cls:
        b = cl['bits']
        for arc in range(m):
            fb = b ^ (1 << (m-1-arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                edges.add((min(cl['cid'], nb), max(cl['cid'], nb)))

    # Merge complement pairs
    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1
    me = set()
    for (a,b) in edges:
        x, y = mn[a], mn[b]
        if x != y: me.add((min(x,y), max(x,y)))

    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'score': cl0['score'],
                  'c3': cl0['c3']}

    V = mid; E = len(me)

    # Build adjacency matrix
    A = np.zeros((V, V))
    for (a, b) in me:
        A[a, b] = 1; A[b, a] = 1

    return A, V, E, mi, me

def analyze_spectrum(n, A, V, E, mi, me):
    """Full spectral analysis."""
    print(f"\n{'='*72}")
    print(f"  EIGENVALUE ANALYSIS: G_{n}/Z_2  (V={V}, E={E})")
    print(f"{'='*72}")

    # Adjacency spectrum
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    idx = np.argsort(-eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    print(f"\n  ADJACENCY SPECTRUM ({V} eigenvalues):")
    for i in range(min(V, 20)):
        print(f"    λ_{i} = {eigenvalues[i]:.6f}")
    if V > 20:
        print(f"    ... ({V-20} more)")
        print(f"    λ_{V-1} = {eigenvalues[-1]:.6f}")

    # Degree sequence
    degrees = np.sum(A, axis=1)
    avg_deg = np.mean(degrees)
    max_deg = np.max(degrees)
    min_deg = np.min(degrees)

    print(f"\n  DEGREE: avg={avg_deg:.2f}, min={min_deg:.0f}, max={max_deg:.0f}")
    print(f"  λ_0 = {eigenvalues[0]:.4f} (≈ max_deg = {max_deg:.0f})")

    # Spectral gap
    adj_gap = eigenvalues[0] - eigenvalues[1]
    print(f"\n  ADJACENCY GAP: λ_0 - λ_1 = {adj_gap:.4f}")

    # Markov (random walk) spectrum
    D_inv = np.diag(1.0 / degrees)
    M = D_inv @ A  # transition matrix
    markov_evals = np.sort(np.linalg.eigvals(M).real)[::-1]

    markov_gap = 1 - markov_evals[1]
    two_over_n = 2.0 / n

    print(f"\n  MARKOV SPECTRUM:")
    for i in range(min(V, 10)):
        print(f"    μ_{i} = {markov_evals[i]:.6f}")
    print(f"\n  Markov gap: 1 - μ_1 = {markov_gap:.6f}")
    print(f"  2/n = {two_over_n:.6f}")
    print(f"  Ratio gap/(2/n) = {markov_gap/two_over_n:.4f}")

    # Laplacian spectrum
    L = np.diag(degrees) - A
    lap_evals = np.sort(np.linalg.eigvalsh(L))
    algebraic_conn = lap_evals[1]

    print(f"\n  LAPLACIAN:")
    print(f"    Algebraic connectivity (λ_1^L) = {algebraic_conn:.6f}")
    print(f"    Spectral radius = {lap_evals[-1]:.4f}")

    # Cheeger constant bound: h ≥ λ_1^L / (2 * d_max)
    cheeger_lower = algebraic_conn / (2 * max_deg)
    print(f"    Cheeger lower bound: h ≥ {cheeger_lower:.4f}")

    # Ramanujan check: is graph a good expander?
    # For d-regular: Ramanujan iff |λ_1| ≤ 2√(d-1)
    # For non-regular: check against spectral radius
    ramanujan_bound = 2 * np.sqrt(avg_deg - 1)
    is_ramanujan = abs(eigenvalues[1]) <= ramanujan_bound
    print(f"\n  RAMANUJAN CHECK:")
    print(f"    |λ_1| = {abs(eigenvalues[1]):.4f}")
    print(f"    2√(d̄-1) = {ramanujan_bound:.4f}")
    print(f"    Ramanujan-like: {is_ramanujan}")

    # H as eigenvector
    H_vals = np.array([mi[i]['H'] for i in range(V)], dtype=float)
    H_normalized = H_vals - np.mean(H_vals)
    H_normalized = H_normalized / np.linalg.norm(H_normalized)

    # Project H onto eigenbasis
    projections = np.abs(eigenvectors.T @ H_normalized) ** 2

    print(f"\n  H PROJECTION ONTO EIGENBASIS:")
    top_projs = sorted(enumerate(projections), key=lambda x: -x[1])[:5]
    total_top5 = sum(p for _, p in top_projs)
    for rank, (i, p) in enumerate(top_projs):
        print(f"    #{rank+1}: λ_{i} = {eigenvalues[i]:.4f}, |⟨v_{i}, H⟩|² = {p:.4f} ({100*p:.1f}%)")
    print(f"    Top 5 total: {100*total_top5:.1f}%")

    # Correlation with H
    best_corr_i = top_projs[0][0]
    best_corr = np.corrcoef(eigenvectors[:, best_corr_i], H_vals)[0, 1]
    print(f"\n  Best H correlate: eigenvector #{best_corr_i}, r = {best_corr:.4f}")

    # SC/NS spectral decomposition
    sc_mask = np.array([mi[i]['sc'] for i in range(V)])
    n_sc = np.sum(sc_mask)
    n_ns = V - n_sc

    if n_sc > 0 and n_ns > 0:
        # Spine/ribs/sea edge classification
        spine = sum(1 for (a,b) in me if mi[a]['sc'] and mi[b]['sc'])
        ribs = sum(1 for (a,b) in me if mi[a]['sc'] != mi[b]['sc'])
        sea = sum(1 for (a,b) in me if not mi[a]['sc'] and not mi[b]['sc'])

        print(f"\n  SC/NS STRUCTURE:")
        print(f"    SC nodes: {n_sc}, NS nodes: {n_ns}")
        print(f"    Spine: {spine}, Ribs: {ribs}, Sea: {sea}")

        # Fiedler vector (2nd Laplacian eigenvector) and SC separation
        lap_vecs = np.linalg.eigh(L)[1]
        fiedler = lap_vecs[:, 1]  # 2nd smallest eigenvalue's vector

        sc_fiedler = fiedler[sc_mask.astype(bool)]
        ns_fiedler = fiedler[~sc_mask.astype(bool)]

        print(f"\n  FIEDLER SEPARATION (SC vs NS):")
        print(f"    SC mean = {np.mean(sc_fiedler):.4f}, std = {np.std(sc_fiedler):.4f}")
        print(f"    NS mean = {np.mean(ns_fiedler):.4f}, std = {np.std(ns_fiedler):.4f}")
        if np.std(sc_fiedler) + np.std(ns_fiedler) > 0:
            separation = abs(np.mean(sc_fiedler) - np.mean(ns_fiedler)) / (np.std(sc_fiedler) + np.std(ns_fiedler))
            print(f"    Separation index: {separation:.4f}")

    # Spectral density
    print(f"\n  SPECTRAL DENSITY:")
    # Normalize eigenvalues to [-1, 1] range
    evals_norm = eigenvalues / np.max(np.abs(eigenvalues))
    hist, bin_edges = np.histogram(evals_norm, bins=10)
    for i in range(len(hist)):
        bar = '#' * hist[i]
        center = (bin_edges[i] + bin_edges[i+1]) / 2
        print(f"    [{center:+.2f}] {bar} ({hist[i]})")

    # Trace formulas
    trace_A2 = np.trace(A @ A)  # = 2E
    trace_A3 = np.trace(A @ A @ A)  # = 6 × triangles
    n_triangles = int(round(trace_A3 / 6))

    print(f"\n  TRACE FORMULAS:")
    print(f"    tr(A²) = {trace_A2:.0f} = 2E = {2*E} {'✓' if abs(trace_A2-2*E)<0.5 else '✗'}")
    print(f"    tr(A³) = {trace_A3:.0f} → triangles = {n_triangles}")
    print(f"    Σ λ_i² = {sum(eigenvalues**2):.4f}")
    print(f"    Σ λ_i³ = {sum(eigenvalues**3):.4f}")

    # Return key data
    return {
        'eigenvalues': eigenvalues,
        'markov_evals': markov_evals,
        'markov_gap': markov_gap,
        'adj_gap': adj_gap,
        'algebraic_conn': algebraic_conn,
        'H_proj_best': top_projs[0][1],
        'H_proj_idx': top_projs[0][0],
        'is_ramanujan': is_ramanujan,
        'n_triangles': n_triangles,
        'avg_deg': avg_deg,
    }

if __name__ == '__main__':
    print("=" * 72)
    print("  DEEP EIGENVALUE ANALYSIS OF G_n/Z_2")
    print("  opus-2026-03-23-S268")
    print("=" * 72)

    results = {}
    for n in range(3, 9):
        t0 = time.time()
        print(f"\n  Building G_{n}/Z_2...", end=' ', flush=True)
        A, V, E, mi, me = build_merged_metagraph(n)
        dt_build = time.time() - t0
        print(f"done ({dt_build:.1f}s)")

        r = analyze_spectrum(n, A, V, E, mi, me)
        results[n] = r

    # Summary table
    print(f"\n{'='*72}")
    print("  SUMMARY")
    print(f"{'='*72}")
    print(f"  {'n':>3} {'V':>6} {'E':>7} {'gap_M':>7} {'2/n':>6} {'ratio':>6} "
          f"{'H_proj%':>8} {'at_λ':>5} {'Raman':>6} {'algcon':>7} {'tri':>6}")
    for n in sorted(results.keys()):
        r = results[n]
        A, V, E, mi, me = build_merged_metagraph(n)
        gap_r = r['markov_gap'] / (2.0/n)
        print(f"  {n:>3} {V:>6} {E:>7} {r['markov_gap']:>7.4f} {2.0/n:>6.4f} {gap_r:>6.3f} "
              f"{100*r['H_proj_best']:>7.1f}% {r['H_proj_idx']:>5} "
              f"{'✓' if r['is_ramanujan'] else '✗':>6} {r['algebraic_conn']:>7.3f} "
              f"{r['n_triangles']:>6}")

    print(f"\n  KEY FINDINGS:")
    print(f"  1. Markov gap / (2/n) → ? as n grows")
    print(f"  2. H projection in 2nd eigenvector → ? as n grows")
    print(f"  3. Ramanujan-like: how does |λ_1| compare to 2√(d̄-1)?")

    print("\nDONE.")
