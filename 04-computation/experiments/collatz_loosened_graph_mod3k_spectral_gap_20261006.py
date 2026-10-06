#!/usr/bin/env python3
"""
collatz_loosened_graph_mod3k_spectral_gap_20261006.py

KLS-inspired question (Song--Zhang 2026 prove a dimension-free Cheeger/Poincare constant for isotropic
log-concave measures): does the Le--Smith loosened Collatz graph E (arcs n -> 3n+1 for every n, n -> n/2 for
even n; the object of the repo's E-SCC conjecture HYP-9120 / Q2) have a DIMENSION-FREE spectral gap on its
3-adic quotients?  Modulo 3^k both generators descend: f(x) = 3x+1 and h(x) = x * 2^{-1} (2 is a primitive root
mod 3^k, so h is a single cycle of length 2*3^(k-1) on the units U_k).  f maps U_k 3-to-1 onto the class
1 (mod 3).  E_k := the Schreier digraph of {f, h} on U_k.  We compute, for k = 2..KMAX:
  * the second-largest singular/eigen quantities of the symmetrised lazy walk P = (1/2)(I + (M + M^T)/2)/..., i.e.
    the spectral gap 1 - lambda_2 of the undirected Schreier graph's normalised Laplacian (degrees vary: 2 out,
    1 + 3[x = 1 mod 3] in),
  * the Cheeger lower bound h >= (1 - lambda_2)/2 and the upper bound h <= sqrt(2(1 - lambda_2)),
  * the mixing of the directed walk: ||P^t 1_x - pi||_TV after t = c k steps, and the second eigenvalue modulus of P.
Also the same for the Syracuse-type graph with only h (a cycle: gap ~ 1/|U_k|^2, control) and for the random
"Chung--Diaconis--Graham" comparison x -> 2x, x -> 2x+1 on Z/3^k (an affine expander-like walk).
Reproduce: python3 collatz_loosened_graph_mod3k_spectral_gap_20261006.py [KMAX]
"""
import sys, math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9

def units(k):
    N = 3 ** k
    return [x for x in range(N) if x % 3 != 0]

def schreier(k, gens):
    N = 3 ** k
    U = units(k); idx = {x: i for i, x in enumerate(U)}; n = len(U)
    rows = []; cols = []
    for x in U:
        for g in gens:
            y = g(x, N)
            rows.append(idx[x]); cols.append(idx[y])
    M = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))   # out-arcs (with multiplicity)
    return U, M

def gens_E(k):
    inv2 = pow(2, -1, 3 ** k)
    return [lambda x, N: (3 * x + 1) % N, lambda x, N, inv2=inv2: (x * inv2) % N]

def gens_cycle(k):
    inv2 = pow(2, -1, 3 ** k)
    return [lambda x, N, inv2=inv2: (x * inv2) % N]

def gens_CDG(k):
    return [lambda x, N: (2 * x) % N, lambda x, N: (2 * x + 1) % N]   # not unit-preserving; use on all of Z/3^k below

def normalized_laplacian_gap(M):
    A = (M + M.T).tocsr()                      # undirected multigraph adjacency
    d = np.asarray(A.sum(axis=1)).ravel()
    Dm = sp.diags(1 / np.sqrt(d))
    L = sp.identity(A.shape[0]) - Dm @ A @ Dm
    k = min(6, A.shape[0] - 2)
    vals = spla.eigsh(L, k=k, which="SM", return_eigenvectors=False, tol=1e-10)
    vals = np.sort(vals)
    return vals[1], vals  # lambda_2 of the normalized Laplacian (gap), spectrum head

def directed_walk_second_eig(M, dense_limit=4500):
    d = np.asarray(M.sum(axis=1)).ravel()
    P = sp.diags(1 / d) @ M                    # row-stochastic (each out-arc with equal weight)
    n = M.shape[0]
    if n > dense_limit:
        return None, P, None
    vals = np.linalg.eigvals(P.toarray())
    mods = np.sort(np.abs(vals))[::-1]
    # spectrum structure: count moduli in bins
    structure = {"=1": int(np.sum(np.abs(mods - 1) < 1e-8)),
                 "=1/2": int(np.sum(np.abs(mods - 0.5) < 1e-6)),
                 "=0": int(np.sum(mods < 1e-8)),
                 "other": int(np.sum((np.abs(mods - 1) >= 1e-8) & (np.abs(mods - 0.5) >= 1e-6) & (mods >= 1e-8)))}
    return mods[1], P, structure

def tv_mixing(P, t_list, start=0):
    n = P.shape[0]
    # stationary distribution via power iteration
    pi = np.ones(n) / n
    for _ in range(5000):
        pi_new = P.T @ pi
        if np.abs(pi_new - pi).sum() < 1e-14:
            pi = pi_new; break
        pi = pi_new
    out = {}
    v = np.zeros(n); v[start] = 1.0
    t = 0
    for T in t_list:
        while t < T:
            v = P.T @ v; t += 1
        out[T] = 0.5 * np.abs(v - pi).sum()
    return out, pi

print("Loosened Collatz graph E_k on the units of Z/3^k: generators f(x) = 3x+1 (3-to-1 onto 1 mod 3), h(x) = x/2 (a single cycle)")
print(f"{'k':>2s} {'|U_k|':>7s} {'gap lam2(NormLap)':>18s} {'Cheeger in [g/2, sqrt(2g)]':>28s} {'|lam2| directed':>16s} {'TV after 4k':>12s} {'TV after 8k':>12s} {'cycle-only gap':>15s}")
rows = []
for k in range(2, KMAX + 1):
    U, M = schreier(k, gens_E(k))
    gap, head = normalized_laplacian_gap(M)
    lam2d, P, structure = directed_walk_second_eig(M)
    tv, pi = tv_mixing(P, [k, 2 * k, 4 * k, 8 * k])
    Lk = 2 * 3 ** (k - 1)
    gapc = 1 - math.cos(2 * math.pi / Lk)          # exact normalized-Laplacian gap of the h-cycle alone
    rows.append((k, len(U), gap, lam2d, tv[4 * k], tv[8 * k], gapc))
    l2 = f"{lam2d:16.6f}" if lam2d is not None else f"{'(n too large)':>16s}"
    print(f"{k:2d} {len(U):7d} {gap:18.6f} {'[%.4f, %.4f]' % (gap/2, math.sqrt(2*gap)):>28s} {l2} {tv[4*k]:12.2e} {tv[8*k]:12.2e} {gapc:15.2e}")
    if structure:
        print(f"      directed spectrum moduli: {structure};  TV after k, 2k steps: {tv[k]:.3e}, {tv[2*k]:.3e};  stationary max/min: {pi.max()/pi.min():.3f}")
        # predicted spectrum: {1, -1/2} u { e^{2 pi i j / L_m}/2 : 2 <= m <= k, 0 <= j < L_m, 3 does not divide j }
        pred = [1.0, -0.5]
        for mm in range(2, k + 1):
            Lm = 2 * 3 ** (mm - 1)
            pred += [np.exp(2j * np.pi * jj / Lm) / 2 for jj in range(Lm) if jj % 3 != 0]
        d = np.asarray(M.sum(axis=1)).ravel(); Pd = (sp.diags(1 / d) @ M).toarray()
        ev = np.linalg.eigvals(Pd)
        pred = np.array(sorted(pred, key=lambda z: (round(z.real, 9), round(z.imag, 9))))
        got = np.array(sorted(ev, key=lambda z: (round(z.real, 9), round(z.imag, 9))))
        ok = len(pred) == len(got) and np.max(np.abs(pred - got)) < 1e-6
        print(f"      predicted spectrum (1, -1/2, and e^(2 pi i j/L_m)/2 for 2<=m<=k, 3 does not divide j) matches: {ok}  (max dev {np.max(np.abs(pred-got)) if len(pred)==len(got) else 'size mismatch'})")
print("\nIf the gap stays bounded below as k grows, E_k is an expander family (a KLS-type dimension-free constant);")
print("if it decays like 1/k or 1/3^k it is not.  The cycle-only column is the control (gap ~ 1/|U_k|^2).")
# fit decay
ks = np.array([r[0] for r in rows[2:]]); gaps = np.array([r[2] for r in rows[2:]])
if len(ks) >= 3:
    slope = np.polyfit(np.log(3.0 ** ks), np.log(gaps), 1)[0]
    print(f"log-log slope of the undirected gap vs |Z/3^k| over k >= 4: {slope:.3f}  (0 = expander, -1 = cycle-like in 1/N, -2 = pure cycle; -log_3 2 = -0.631 = halving per level)")
    ratios = [rows[i][2] / rows[i + 1][2] for i in range(len(rows) - 1)]
    print("successive gap ratios gap_k/gap_(k+1):", ["%.4f" % r for r in ratios])

print("\nComparison: Chung--Diaconis--Graham affine walk x -> 2x, x -> 2x+1 on all of Z/3^k (classical fast mixer)")
for k in range(2, min(KMAX, 8) + 1):
    N = 3 ** k
    rows_i = []; cols_i = []
    for x in range(N):
        for y in ((2 * x) % N, (2 * x + 1) % N):
            rows_i.append(x); cols_i.append(y)
    M = sp.csr_matrix((np.ones(len(rows_i)), (rows_i, cols_i)), shape=(N, N))
    gap, _ = normalized_laplacian_gap(M)
    lam2d, P, structure = directed_walk_second_eig(M)
    print(f"  k={k}: N={N:6d}  NormLap gap {gap:.6f}   directed |lam2| {lam2d if lam2d is None else round(lam2d, 6)}   spectrum {structure}")
