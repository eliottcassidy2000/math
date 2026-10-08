#!/usr/bin/env python3
"""Optimal block-balance margin by the ellipsoid method (quasi-concave maximization in P = Q^-1).

For a block matrix A (PSD), whitening by Q gives W = Q^-1/2 A Q^-1/2, which has the same nonzero spectrum as
X_A(P) = A^1/2 P A^1/2, linear in P.  The margin of A under Q is 1 - 2 lambda_max(X_A(P)) / tr(A P); for fixed mu the
constraint margin >= mu reads lambda_max(X_A(P)) - (1 - mu)/2 tr(A P) <= 0, convex in P.  So the minimum margin over a
family of block matrices is quasi-concave in P, and the central-cut ellipsoid method (cut through the centre with the
subgradient of the worst constraint) converges to the global optimum.  The optimum is then rationalized and checked
exactly by block_balance2.exact_ok.  Usage: python3 block_lmi.py MAPNAME k [iters]"""
import math, sys, time
import numpy as np
from block_balance import int_coords
from block_balance2 import block_arrays, exact_ok, integer_forms, MAPS

def sym_basis(n):
    """Orthonormal basis of traceless symmetric n x n matrices (Frobenius)."""
    B = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = E[j, i] = 1 / math.sqrt(2); B.append(E)
    for k in range(1, n):
        E = np.zeros((n, n))
        for i in range(k): E[i, i] = 1
        E[k, k] = -k; E /= np.linalg.norm(E); B.append(E)
    return np.array(B)

def psd_sqrt(As):
    w, V = np.linalg.eigh(As)
    w = np.clip(w, 0, None)
    return (V * np.sqrt(w)[:, None, :]) @ np.transpose(V, (0, 2, 1))

def margin_all(P, Ah, As):
    X = Ah @ P @ Ah
    w, V = np.linalg.eigh(X)
    tr = np.einsum('nij,ij->n', As, P)
    return 1 - 2 * w[:, -1] / tr, V[:, :, -1], w[:, -1], tr

def ellipsoid_opt(As, iters=600, verbose=False):
    n = As.shape[1]
    Ah = psd_sqrt(As)
    B = sym_basis(n); dim = len(B)
    c = np.zeros(dim); E = np.eye(dim) * 4.0          # ellipsoid {y : (y-c)^T E^-1 (y-c) <= 1}
    best_m, best_P = -np.inf, None
    for it in range(iters):
        P = np.eye(n) / n + np.tensordot(c, B, axes=1)
        wP, VP = np.linalg.eigh(P)
        if wP[0] <= 1e-9:
            u = VP[:, 0]
            g = -np.array([u @ Bk @ u for Bk in B])     # feasible region has u^T P u > current value: cut g.(y - c) <= 0
        else:
            m, top, lam, tr = margin_all(P, Ah, As)
            k = int(np.argmin(m)); mk = m[k]
            if mk > best_m: best_m, best_P = mk, P.copy()
            # constraint h(P) = lam_max(Ah P Ah) - (1 - mk)/2 tr(A P); gradient wrt P: Ah w w^T Ah - (1-mk)/2 A
            wv = Ah[k] @ top[k]
            G = np.outer(wv, wv) - (1 - mk) / 2 * As[k]
            g = np.array([np.sum(G * Bk) for Bk in B])
        # central cut: keep {y : g.(y - c) <= 0}
        Eg = E @ g; gEg = g @ Eg
        if gEg <= 1e-30: break
        gt = Eg / math.sqrt(gEg)
        c = c - gt / (dim + 1)
        E = (dim * dim / (dim * dim - 1.0)) * (E - (2.0 / (dim + 1)) * np.outer(gt, gt))
        if verbose and it % 100 == 0: print(f"      it {it}: best margin {best_m:+.6f}", flush=True)
    return best_P, best_m

if __name__ == '__main__':
    name = sys.argv[1]; k = int(sys.argv[2]); iters = int(sys.argv[3]) if len(sys.argv) > 3 else 800
    d, m, r = MAPS[name]
    rho, v = int_coords(m)
    t0 = time.time()
    H, A = block_arrays(d, m, r, v, rho, k)
    flat = A.reshape(-1, rho, rho)
    nz = np.any(flat.reshape(len(flat), -1) != 0, axis=1)
    U = np.unique(flat[nz].reshape(-1, rho * rho), axis=0).reshape(-1, rho, rho)
    Uf = U.astype(float)
    Uf = Uf / np.trace(Uf, axis1=1, axis2=2)[:, None, None]
    P, mg = ellipsoid_opt(Uf, iters=iters)
    amax = 2 / (1 - mg) - 2 if mg > 0 else float('nan')
    Q = np.linalg.inv(P)
    print(f"{name} k = {k}: {len(U)} distinct nonzero blocks; optimal margin {mg:+.5f} (alpha_max {amax:.4f}); "
          f"Q ~ {np.round(Q / np.abs(Q).max(), 4).tolist()}  [{time.time() - t0:.1f}s]", flush=True)
    if mg > 0:
        for Qi in integer_forms(Q, rho):
            ok, bad = exact_ok(Qi, U)
            if ok:
                print(f"   EXACT: integer form Q = {Qi} balances all {len(U)} distinct block matrices  [{time.time() - t0:.1f}s]", flush=True)
                break
        else:
            print("   exact rationalization failed at the tried scales", flush=True)
