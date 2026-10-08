#!/usr/bin/env python3
"""Block balance, vectorized (numpy) with an exact integer verification.  See block_balance.py for the criterion:
one positive definite Q with tr(A Q^-1) > 2 lambda_max(A Q^-1) for every nonzero k-step block matrix
A_k(h) = d^k Sigma_k(h), h = (M mod d^k, e mod d^k) in H_k x Z/d^k, makes the debt walk transient (THM-4609 (3) applied
to the walk sampled every k steps).

Exact check, integer only: with adj(Q) and det(Q) > 0, the condition is that
        S' = tr(A adj Q) Q - 2 det(Q) A
is positive definite; we test its leading principal minors with Python integers on every distinct A.
Usage: python3 block_balance2.py MAPNAME k  (maps listed in MAPS)"""
import math, sys, time
import numpy as np
from block_balance import int_coords, ratio_group, coupling_D, std_r

def block_arrays(d, m, r, v, rho, k):
    """Return (H_k list, A array of shape (|H_k|, d^k, rho, rho), int64)."""
    D = np.zeros((d, d, rho, rho), dtype=np.int64)
    for a in range(1, d):
        if math.gcd(a, d) != 1: continue
        for b in range(d): D[a, b] = np.array(coupling_D(v, a, b, d, rho), dtype=np.int64)
    prevH = None; prevA = None
    for s in range(1, k + 1):
        mod = d ** s; H = np.array(ratio_group(m, mod), dtype=np.int64); E = np.arange(mod, dtype=np.int64)
        a = (H % d)[:, None] * np.ones((1, mod), dtype=np.int64)
        b = np.ones((len(H), 1), dtype=np.int64) * (E % d)[None, :]
        A = D[a, b] * (d ** (s - 1))
        if s > 1:
            pmod = d ** (s - 1)
            pos = -np.ones(pmod, dtype=np.int64); pos[prevH] = np.arange(len(prevH))
            Mg = H[:, None] * np.ones((1, mod), dtype=np.int64); Eg = np.ones((len(H), 1), dtype=np.int64) * E[None, :]
            mm = np.array(m, dtype=np.int64); rr = np.array(r, dtype=np.int64)
            inv = np.array([pow(x, -1, mod) for x in m], dtype=np.int64)
            for j in range(d):
                i = (a * j + b) % d
                # N = m_i e + r_i - m_i inv(m_j) M r_j  (mod d^s); all intermediate products < 2^63 for d^s <= 10^5
                t = (mm[i] * inv[j]) % mod
                t = (t * Mg) % mod
                t = (t * rr[j]) % mod
                N = (mm[i] * Eg + rr[i] - t) % mod
                assert np.all(N % d == 0)
                e2 = (N // d) % pmod
                M2 = (((mm[i] * inv[j]) % mod) * Mg) % pmod
                idx = pos[M2]
                assert np.all(idx >= 0)
                A = A + prevA[idx, e2]
        prevH, prevA = H, A
    return prevH, prevA

def margins(L, As):
    Li = np.linalg.inv(L)
    W = Li @ As @ Li.T
    ev = np.linalg.eigvalsh(W)
    tr = ev.sum(axis=1)
    return (tr - 2 * ev[:, -1]) / tr

def optimize(As, rho, rnd, Q0=None, iters=4000):
    """Maximize the minimum margin over As by random local search on the Cholesky factor, with an active set."""
    if Q0 is None: Q0 = As.mean(axis=0)
    L = np.linalg.cholesky(Q0 / np.trace(Q0))
    full = margins(L, As); best = full.min()
    for rnd_round in range(6):
        act = As[np.argsort(full)[:3000]]
        step = 0.1; bl = L; bm = margins(L, act).min()
        for it in range(iters):
            Ln = bl + np.tril(rnd.normal(0, step, (rho, rho)))
            if np.any(np.diag(Ln) <= 1e-4): continue
            mn = margins(Ln, act).min()
            if mn > bm: bm, bl = mn, Ln
            if it % 800 == 799: step *= 0.6
        full_n = margins(bl, As)
        if full_n.min() > best: best, L, full = full_n.min(), bl, full_n
        else: break
    return L @ L.T, best

def adj3(Q):
    n = len(Q)
    def minor(i, j): return det_int([[Q[x][y] for y in range(n) if y != j] for x in range(n) if x != i])
    return [[(-1) ** (i + j) * minor(j, i) for j in range(n)] for i in range(n)]

def det_int(M):
    n = len(M)
    if n == 0: return 1
    if n == 1: return M[0][0]
    if n == 2: return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return sum((-1) ** c * M[0][c] * det_int([row[:c] + row[c + 1:] for row in M[1:]]) for c in range(n))

def exact_ok(Qi, distinct):
    """Qi: integer symmetric matrix (list of lists). distinct: array (N, rho, rho) of int64. Exact integer test."""
    rho = len(Qi)
    if not all(det_int([row[:t] for row in Qi[:t]]) > 0 for t in range(1, rho + 1)): return False, None
    adj = adj3(Qi); dq = det_int(Qi)
    bad = None
    for A in distinct:
        A = [[int(x) for x in row] for row in A]
        t = sum(A[x][y] * adj[y][x] for x in range(rho) for y in range(rho))
        S = [[t * Qi[x][y] - 2 * dq * A[x][y] for y in range(rho)] for x in range(rho)]
        if not all(det_int([row[:u] for row in S[:u]]) > 0 for u in range(1, rho + 1)):
            bad = A; return False, bad
    return True, None

def integer_forms(Q, rho):
    Qn = Q / np.abs(Q).max()
    for S in (4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 50, 64, 100, 128, 200, 256, 500, 1000):
        Qi = [[int(round(Qn[i][j] * S)) for j in range(rho)] for i in range(rho)]
        Qi = [[(Qi[i][j] + Qi[j][i]) // 2 if (Qi[i][j] + Qi[j][i]) % 2 == 0 else Qi[i][j] for j in range(rho)] for i in range(rho)]
        Qi = [[Qi[min(i, j)][max(i, j)] for j in range(rho)] for i in range(rho)]
        g = 0
        for row in Qi:
            for x in row: g = math.gcd(g, x)
        if g > 1: Qi = [[x // g for x in row] for row in Qi]
        yield Qi

MAPS = {
    'Z5_12371': (5, [1, 2, 3, 7, 1], [0, 3, 4, 4, 1]),
    'Z5_11237': (5, [1, 1, 2, 3, 7], std_r([1, 1, 2, 3, 7])),
    'Z5_123711': (5, [1, 2, 3, 7, 11], std_r([1, 2, 3, 7, 11])),
    'Z4_1357': (4, [1, 3, 5, 7], std_r([1, 3, 5, 7])),
    'Z7_1111235': (7, [1, 1, 1, 1, 2, 3, 5], std_r([1, 1, 1, 1, 2, 3, 5])),
    'Z7_1112357': (7, [1, 1, 1, 2, 3, 5, 7], std_r([1, 1, 1, 2, 3, 5, 7])),
}

if __name__ == '__main__':
    name = sys.argv[1]; k = int(sys.argv[2])
    d, m, r = MAPS[name]
    for i in range(d): assert (m[i] * i + r[i]) % d == 0
    rho, v = int_coords(m)
    Lam = sum(math.log(x / d) for x in m) / d
    t0 = time.time()
    H, A = block_arrays(d, m, r, v, rho, k)
    flat = A.reshape(-1, rho, rho)
    nz = np.any(flat.reshape(len(flat), -1) != 0, axis=1)
    U = np.unique(flat[nz].reshape(-1, rho * rho), axis=0).reshape(-1, rho, rho)
    ranks = np.linalg.matrix_rank(U.astype(float))
    print(f"{name} m = {m} r = {r} rank {rho} Lambda = {Lam:+.4f}; k = {k}: {len(flat)} states, {int((~nz).sum())} frozen, "
          f"{len(U)} distinct nonzero blocks, min rank {ranks.min()}  [{time.time() - t0:.1f}s]", flush=True)
    if ranks.min() < 3:
        bad = U[np.argmin(ranks)]
        print(f"   impossible: a block matrix of rank {ranks.min()}: {bad.tolist()}", flush=True); sys.exit(0)
    rnd = np.random.default_rng(17)
    Q, mg = optimize(U.astype(float), rho, rnd)
    amax = 2 / (1 - mg) - 2 if mg > 0 else float('nan')
    print(f"   best float margin {mg:+.5f} (alpha_max {amax:.4f}); Q ~ {np.round(Q / np.abs(Q).max(), 4).tolist()}  [{time.time() - t0:.1f}s]", flush=True)
    if mg > 0:
        for Qi in integer_forms(Q, rho):
            ok, bad = exact_ok(Qi, U)
            if ok:
                print(f"   EXACT: integer form Q = {Qi} balances all {len(U)} distinct block matrices  [{time.time() - t0:.1f}s]", flush=True)
                break
        else:
            print("   exact rationalization failed at the tried scales", flush=True)
