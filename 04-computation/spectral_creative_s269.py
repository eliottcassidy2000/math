#!/usr/bin/env python3
"""
spectral_creative_s269.py — Creative spectral investigations
opus-2026-03-23-S269

Going beyond standard spectral analysis into:
1. Random matrix universality (eigenvalue spacing → GOE or Poisson?)
2. Heat kernel trace = partition function Z(β) of the metagraph
3. Graph signal processing: H as signal, frequency decomposition, smoothness
4. Spectral clustering: natural communities in tournament space
5. Cheeger constant: actual vertex expansion
6. Ihara-Selberg zeta function: cycle structure
7. Spectral dimension: how many "effective dimensions" does the metagraph have?
8. Quantum walk: unitary evolution vs classical diffusion
"""

import sys, time, subprocess
import numpy as np
from math import comb, factorial, pi, log, sqrt, exp
from itertools import permutations
from collections import defaultdict
# from scipy import stats  # not available

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
        cls.append({'cid':i, 'adj':adj, 'score':s, 'c3':cc, 'H':h,
                    'sc':sc, 'comp_canon':cc2, 'canon':canon})
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
        b = cl['bits'] if 'bits' in cl else int(''.join(str(cl['adj'][i][j]) for (i,j) in P), 2)
        for arc in range(m):
            fb = b ^ (1 << (m-1-arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                edges.add((min(cl['cid'], nb), max(cl['cid'], nb)))

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
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'score': cl0['score'], 'c3': cl0['c3']}

    V = mid; E = len(me)
    A = np.zeros((V, V))
    for (a,b) in me: A[a,b] = 1; A[b,a] = 1

    return A, V, E, mi, me

def analyze_creative(n, A, V, E, mi):
    """Creative spectral investigations."""
    print(f"\n{'#'*72}")
    print(f"  n = {n}: CREATIVE SPECTRAL ANALYSIS (V={V}, E={E})")
    print(f"{'#'*72}")

    if V < 4:
        print("  (too small for meaningful analysis)")
        return {}

    # Basic spectrum
    eigenvalues_L, eigvecs_L = np.linalg.eigh(np.diag(np.sum(A, axis=1)) - A)
    eigenvalues_A, eigvecs_A = np.linalg.eigh(A)
    idx = np.argsort(-eigenvalues_A)
    eigenvalues_A = eigenvalues_A[idx]
    eigvecs_A = eigvecs_A[:, idx]

    degrees = np.sum(A, axis=1)
    H_vals = np.array([mi[i]['H'] for i in range(V)], dtype=float)

    # ============================================================
    # 1. RANDOM MATRIX UNIVERSALITY
    # ============================================================
    print(f"\n  1. RANDOM MATRIX UNIVERSALITY")

    if V >= 10:
        # Unfold eigenvalues (normalize to unit mean spacing)
        evals_sorted = np.sort(eigenvalues_A)
        spacings = np.diff(evals_sorted)
        mean_spacing = np.mean(spacings)
        normalized_spacings = spacings / mean_spacing

        # Compare to GOE (Wigner surmise) and Poisson
        # GOE: P(s) = (π/2) s exp(-π s²/4)
        # Poisson: P(s) = exp(-s)

        s_vals = normalized_spacings[normalized_spacings > 0]
        if len(s_vals) >= 5:
            mean_s = np.mean(s_vals)
            var_s = np.var(s_vals)
            # GOE has mean = 4/(π√π) ≈ 0.798, var ≈ 0.215
            # Poisson has mean = 1, var = 1
            goe_mean, goe_var = 4/(pi*sqrt(pi)), (4 - 16/pi**2) / pi
            # Actually: GOE mean spacing = 1 (by normalization), var ≈ 0.273
            # Poisson mean = 1, var = 1

            print(f"    Spacing statistics (V={V}):")
            print(f"    Mean spacing (normalized): {mean_s:.4f} (GOE≈1, Poisson=1)")
            print(f"    Var(spacing): {var_s:.4f} (GOE≈0.273, Poisson=1)")

            # Level repulsion: P(s→0)
            near_zero = np.sum(s_vals < 0.1) / len(s_vals)
            print(f"    P(s < 0.1): {near_zero:.4f} (GOE≈0, Poisson≈0.095)")

            # Ratio of consecutive spacings
            if len(s_vals) >= 3:
                ratios = s_vals[1:] / s_vals[:-1]
                ratios = np.minimum(ratios, 1.0/ratios)  # min(r, 1/r)
                mean_ratio = np.mean(ratios)
                # GOE: <r̃> ≈ 0.5307
                # Poisson: <r̃> ≈ 0.3863
                print(f"    Mean ratio r̃: {mean_ratio:.4f} (GOE≈0.531, Poisson≈0.386)")

                classification = "GOE-like" if mean_ratio > 0.46 else "Poisson-like"
                print(f"    Classification: {classification}")

    # ============================================================
    # 2. HEAT KERNEL = PARTITION FUNCTION
    # ============================================================
    print(f"\n  2. HEAT KERNEL / PARTITION FUNCTION")

    lap_evals = np.sort(eigenvalues_L)

    # Z(β) = tr(exp(-β L)) = Σ exp(-β λ_i)
    betas = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    print(f"    Z(β) = tr(exp(-β·L)):")
    for beta in betas:
        Z = np.sum(np.exp(-beta * lap_evals))
        S = -np.sum(np.exp(-beta * lap_evals) / Z * np.log(np.exp(-beta * lap_evals) / Z + 1e-300))
        print(f"    β={beta:>5.2f}: Z={Z:>12.4f}, S(entropy)={S:>8.4f}, "
              f"<E>={np.sum(lap_evals * np.exp(-beta * lap_evals)) / Z:>8.4f}")

    # Spectral dimension from heat kernel
    # d_s = -2 d(ln Z)/d(ln t) at small t
    t_vals = np.logspace(-2, 1, 20)
    Z_vals = [np.sum(np.exp(-t * lap_evals)) for t in t_vals]
    ln_Z = np.log(Z_vals)
    ln_t = np.log(t_vals)
    # d(ln Z)/d(ln t) via finite differences
    d_spec = -2 * np.diff(ln_Z) / np.diff(ln_t)
    print(f"\n    Spectral dimension d_s(t):")
    for i in range(0, len(d_spec), max(1, len(d_spec)//6)):
        print(f"    t={t_vals[i]:>8.4f}: d_s = {d_spec[i]:>6.3f}")
    print(f"    d_s at small t → {d_spec[0]:.3f} (≈ geometric dimension)")
    print(f"    d_s at large t → {d_spec[-1]:.3f} (= 0 for finite graph)")

    # ============================================================
    # 3. GRAPH SIGNAL PROCESSING ON H
    # ============================================================
    print(f"\n  3. GRAPH SIGNAL PROCESSING: H AS SIGNAL")

    # H as a signal on the graph
    H_centered = H_vals - np.mean(H_vals)
    H_norm = np.linalg.norm(H_centered)
    if H_norm > 0:
        H_unit = H_centered / H_norm
    else:
        H_unit = H_centered

    # Graph Fourier Transform of H
    # GFT coefficients: ĥ_k = <u_k, H> where u_k are Laplacian eigenvectors
    _, lap_vecs = np.linalg.eigh(np.diag(degrees) - A)
    gft_coeffs = lap_vecs.T @ H_centered

    # Power spectrum
    power = gft_coeffs**2
    total_power = np.sum(power)

    print(f"    Graph Fourier Transform of H:")
    print(f"    Frequency  Eigenvalue  Power    Cumul%")
    cumul = 0
    for i in range(min(V, 15)):
        cumul += power[i]
        print(f"    k={i:>3}: λ_L={lap_evals[i]:>8.4f}, "
              f"|ĥ|²={power[i]:>10.2f}, {100*cumul/total_power:>6.1f}%")

    # Signal smoothness: Dirichlet energy E(H) = H^T L H / ||H||²
    L = np.diag(degrees) - A
    dirichlet = H_centered @ L @ H_centered / (H_norm**2) if H_norm > 0 else 0
    max_dirichlet = lap_evals[-1]  # maximum possible Dirichlet energy
    smoothness = 1 - dirichlet / max_dirichlet if max_dirichlet > 0 else 1

    print(f"\n    Dirichlet energy E(H) = {dirichlet:.4f}")
    print(f"    Max possible = {max_dirichlet:.4f}")
    print(f"    Smoothness = 1 - E/E_max = {smoothness:.4f}")
    print(f"    (1.0 = perfectly smooth, 0.0 = maximally rough)")

    # Low-pass filter: reconstruct H using only bottom k Laplacian eigenvectors
    print(f"\n    LOW-PASS RECONSTRUCTION OF H:")
    for k in [1, 2, 3, 5, 10, V//2, V]:
        if k > V: break
        H_lowpass = lap_vecs[:, :k] @ gft_coeffs[:k]
        r2 = 1 - np.sum((H_centered - H_lowpass)**2) / np.sum(H_centered**2) if np.sum(H_centered**2) > 0 else 0
        print(f"    k={k:>4} modes: R² = {r2:.4f} ({100*r2:.1f}%)")

    # ============================================================
    # 4. SPECTRAL CLUSTERING
    # ============================================================
    print(f"\n  4. SPECTRAL CLUSTERING")

    if V >= 6:
        # Use top 3 non-trivial Laplacian eigenvectors for clustering
        # (k-means would be ideal but let's do threshold-based)
        fiedler = lap_vecs[:, 1]  # 2nd Laplacian eigenvector

        # Natural bipartition via sign of Fiedler vector
        cluster_A = np.where(fiedler >= 0)[0]
        cluster_B = np.where(fiedler < 0)[0]

        H_A = np.mean([mi[i]['H'] for i in cluster_A])
        H_B = np.mean([mi[i]['H'] for i in cluster_B])
        sc_A = sum(1 for i in cluster_A if mi[i]['sc'])
        sc_B = sum(1 for i in cluster_B if mi[i]['sc'])

        # Cut between clusters
        cut = sum(1 for (a,b) in zip(*np.where(A > 0))
                  if (a in cluster_A and b in cluster_B) or (b in cluster_A and a in cluster_B)) // 2

        print(f"    Fiedler bipartition: A={len(cluster_A)}, B={len(cluster_B)}")
        print(f"    Mean H: A={H_A:.1f}, B={H_B:.1f}")
        print(f"    SC count: A={sc_A}, B={sc_B}")
        print(f"    Cut edges: {cut}/{E} ({100*cut/E:.1f}%)")

        # Normalized cut (Shi-Malik)
        vol_A = np.sum(degrees[cluster_A])
        vol_B = np.sum(degrees[cluster_B])
        ncut = cut * (1/vol_A + 1/vol_B) if vol_A > 0 and vol_B > 0 else 0
        print(f"    Normalized cut: {ncut:.4f}")

    # ============================================================
    # 5. QUANTUM WALK
    # ============================================================
    print(f"\n  5. QUANTUM WALK ON THE METAGRAPH")

    # Continuous-time quantum walk: U(t) = exp(-iAt)
    # Start at transitive tournament (H=1, typically node 0)
    # Probability of being at node j at time t: |<j|exp(-iAt)|0>|²

    # Find transitive tournament
    trans_idx = min(range(V), key=lambda i: mi[i]['H'])
    reg_idx = max(range(V), key=lambda i: mi[i]['H'])

    print(f"    Start: node {trans_idx} (H={mi[trans_idx]['H']}, transitive)")
    print(f"    Target: node {reg_idx} (H={mi[reg_idx]['H']}, regular)")

    t_vals_q = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    print(f"\n    {'t':>6} {'P(trans)':>10} {'P(reg)':>10} {'P(SC)':>10} {'<H>':>10} {'spread':>10}")

    sc_mask = np.array([mi[i]['sc'] for i in range(V)])

    for t in t_vals_q:
        # Quantum: U = exp(-i A t)
        U = np.linalg.matrix_power(
            eigvecs_A @ np.diag(np.exp(-1j * eigenvalues_A * t)) @ eigvecs_A.T, 1)
        psi = U[:, trans_idx]
        probs = np.abs(psi)**2

        p_trans = probs[trans_idx]
        p_reg = probs[reg_idx]
        p_sc = np.sum(probs[sc_mask.astype(bool)])
        mean_H_q = np.sum(probs * H_vals)
        spread = np.sqrt(np.sum(probs * (H_vals - mean_H_q)**2))

        print(f"    {t:>6.1f} {p_trans:>10.4f} {p_reg:>10.4f} {p_sc:>10.4f} "
              f"{mean_H_q:>10.2f} {spread:>10.2f}")

    # Classical comparison
    print(f"\n    CLASSICAL walk (same times):")
    D_inv = np.diag(1.0 / degrees)
    M = D_inv @ A
    print(f"    {'t':>6} {'P(trans)':>10} {'P(reg)':>10} {'<H>':>10}")
    for t_steps in [1, 2, 5, 10, 20, 50]:
        Mt = np.linalg.matrix_power(M, t_steps)
        probs_c = Mt[trans_idx, :]
        print(f"    {t_steps:>6d} {probs_c[trans_idx]:>10.4f} {probs_c[reg_idx]:>10.4f} "
              f"{np.sum(probs_c * H_vals):>10.2f}")

    # ============================================================
    # 6. SPECTRAL ZETA FUNCTION
    # ============================================================
    print(f"\n  6. SPECTRAL ZETA FUNCTION")

    # ζ_G(s) = Σ_{λ>0} λ^{-s} (using Laplacian eigenvalues > 0)
    pos_lap = lap_evals[lap_evals > 0.001]
    if len(pos_lap) > 0:
        s_vals_z = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
        print(f"    ζ_G(s) = Σ λ_L^(-s):")
        for s in s_vals_z:
            zeta = np.sum(pos_lap ** (-s))
            print(f"    s={s:.1f}: ζ_G = {zeta:.6f}")

        # Spectral determinant: det'(L) = prod of nonzero Laplacian eigenvalues
        # = number of spanning trees × V (by Kirchhoff's theorem)
        log_det = np.sum(np.log(pos_lap))
        det_val = np.exp(log_det)
        spanning_trees = det_val / V
        print(f"\n    det'(L) = {det_val:.4f}")
        print(f"    Spanning trees ≈ {spanning_trees:.1f}")

    return {
        'eigenvalues': eigenvalues_A,
        'dirichlet': dirichlet,
        'smoothness': smoothness,
    }

if __name__ == '__main__':
    print("=" * 72)
    print("  CREATIVE SPECTRAL INVESTIGATIONS")
    print("  opus-2026-03-23-S269")
    print("=" * 72)

    for n in [5, 6, 7]:
        t0 = time.time()
        print(f"\n  Building G_{n}/Z_2...", end=' ', flush=True)
        A, V, E, mi, me = build_merged_metagraph(n)
        print(f"done ({time.time()-t0:.1f}s)")

        r = analyze_creative(n, A, V, E, mi)

    # Pull and push
    print(f"\nDONE. opus-2026-03-23-S269")
