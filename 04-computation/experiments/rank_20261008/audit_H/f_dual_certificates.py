#!/usr/bin/env python3
"""audit H: EXACT non-existence certificates for fixed-length block forms (THM-4611 3' "none has a fixed-length
certificate for k <= 4"; Reading "optimal margin -0.275, -0.091 at k = 2, 3" for Z_5 (1,2,3,7,1)).

A form P = Q^-1 balances block A iff tr(A P) > 2 lambda_max(A^1/2 P A^1/2).  For any rational z with z^T A z > 0,
    B(A, z) = A - 2 (A z)(A z)^T / (z^T A z)
satisfies <B, P> <= tr(AP) - 2 lambda_max(...) ... more precisely <B(A,z), P> = tr(AP) - 2 w^T A^1/2 P A^1/2 w for the unit
vector w = A^1/2 z/|A^1/2 z|, so balance of A by P implies <B(A, z), P> > 0.  Hence if  sum_l c_l B(A_l, z_l)  is negative
definite for some c_l >= 0, NO positive definite form balances all the A_l: an exact certificate of non-existence.
We find (c_l, z_l) by a cutting-plane LP (max_P min_l <B_l, P>, trace-1 P, PSD cuts), take the LP duals, rationalize
c_l and z_l, and verify negative definiteness of sum_l c_l B(A_l, z_l) EXACTLY (Fractions, principal-minor-sum test).
The primal side gives the optimal margin approximately (min over blocks of the margin at the LP optimum P).
Usage: python3 f_dual_certificates.py"""
import itertools, math, sys, time
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from hdp import Tables

def sym_vec(B):
    return np.array([B[0, 0], B[1, 1], B[2, 2], 2 * B[0, 1], 2 * B[0, 2], 2 * B[1, 2]])

def to_P(x):
    return np.array([[x[0], x[3], x[4]], [x[3], x[1], x[5]], [x[4], x[5], x[2]]])

def psd_sqrt(A):
    w, V = np.linalg.eigh(A)
    w = np.clip(w, 0, None)
    return (V * np.sqrt(w)[:, None, :]) @ np.transpose(V, (0, 2, 1)), (V * (1 / np.sqrt(np.where(w > 1e-12, w, np.inf)))[:, None, :]) @ np.transpose(V, (0, 2, 1))

def evaluate(P, Ah, An, mu=0.0):
    X = Ah @ P @ Ah
    lam, V = np.linalg.eigh(X)
    tr = np.einsum('nij,ij->n', An, P)
    g = (1 - mu) * tr - 2 * lam[:, -1]
    return g, V[:, :, -1], tr, lam[:, -1]

def cutting_plane(An, Ah, Aih, mu=0.0, iters=200, verbose=False):
    rnd = np.random.default_rng(0)
    cuts_A = []; cuts_w = []; rows = []
    psd_rows = []
    for y in list(np.eye(3)) + list(rnd.normal(size=(60, 3))):
        y = y / np.linalg.norm(y); psd_rows.append(sym_vec(np.outer(y, y)))
    P = np.eye(3) / 3
    best = None
    for it in range(iters):
        g, W, tr, lam = evaluate(P, Ah, An, mu)
        order = np.argsort(g)[:60]
        new = 0
        for i in order:
            B = (1 - mu) * An[i] - 2 * np.outer(Ah[i] @ W[i], Ah[i] @ W[i])
            rows.append(sym_vec(B)); cuts_A.append(i); cuts_w.append(W[i]); new += 1
        nb = len(rows)
        A_ub = np.vstack([np.hstack([-np.array(rows), np.ones((nb, 1))]),
                          np.hstack([-np.array(psd_rows), np.zeros((len(psd_rows), 1))])])
        b_ub = np.zeros(len(A_ub))
        A_eq = np.array([[1, 1, 1, 0, 0, 0, 0]], dtype=float)
        bounds = [(0, 1)] * 3 + [(-1, 1)] * 3 + [(None, None)]
        res = linprog(np.array([0, 0, 0, 0, 0, 0, -1.0]), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=bounds, method='highs')
        assert res.status == 0, res.message
        x = res.x; s = x[6]; Pn = to_P(x[:6])
        ev, EV = np.linalg.eigh(Pn)
        if ev[0] < -1e-12:
            y = EV[:, 0]; psd_rows.append(sym_vec(np.outer(y, y)))
        gt, _, _, _ = evaluate(Pn, Ah, An, mu) if ev[0] >= -1e-12 else (np.array([-np.inf]), None, None, None)
        true_min = gt.min()
        if best is None or true_min > best[0]:
            best = (true_min, Pn.copy())
        if verbose:
            print(f"      it {it}: LP s = {s:+.6f}, true min g = {true_min:+.6f}, cuts {nb}", flush=True)
        P = Pn if ev[0] >= -1e-12 else (EV * np.clip(ev, 1e-9, None)) @ EV.T / np.clip(ev, 1e-9, None).sum()
        if ev[0] >= -1e-12 and true_min >= s - 1e-10:
            break
    duals = -res.ineqlin.marginals[:nb]
    return s, best, duals, cuts_A, cuts_w

def frac_vec(v, den=10 ** 9):
    return [Fr(int(round(x * den)), den) for x in v]

def exact_cert(U, cuts_A, cuts_w, duals, mu=Fr(0)):
    """sum_l c_l [ (1-mu) A_l - 2 (A_l z_l)(A_l z_l)^T/(z_l^T A_l z_l) ] negative definite?  (A_l integer blocks)"""
    N = [[Fr(0)] * 3 for _ in range(3)]
    used = 0
    for i, w, c in zip(cuts_A, cuts_w, duals):
        if c <= 1e-12:
            continue
        A = U[i]
        Af = A.astype(float)
        # z with A^1/2 z parallel to w:  z = A^-1/2 w  (A full rank here)
        lam, V = np.linalg.eigh(Af)
        z = V @ ((V.T @ w) / np.sqrt(lam))
        z = frac_vec(z / np.abs(z).max())
        Ai = [[int(A[p][q]) for q in range(3)] for p in range(3)]
        Az = [sum(Ai[p][q] * z[q] for q in range(3)) for p in range(3)]
        zAz = sum(z[p] * Az[p] for p in range(3))
        assert zAz > 0
        cc = Fr(int(round(c * 10 ** 12)), 10 ** 12) / Fr(int(np.trace(A)))   # weight for the trace-normalised cut
        for p in range(3):
            for q in range(3):
                N[p][q] += cc * ((1 - mu) * Ai[p][q] - 2 * Az[p] * Az[q] / zAz)
        used += 1
    negN = [[-x for x in row] for row in N]
    e1 = sum(negN[p][p] for p in range(3))
    e2 = sum(negN[p][p] * negN[q][q] - negN[p][q] * negN[q][p] for p, q in itertools.combinations(range(3), 2))
    e3 = (negN[0][0] * (negN[1][1] * negN[2][2] - negN[1][2] * negN[2][1]) - negN[0][1] * (negN[1][0] * negN[2][2] - negN[1][2] * negN[2][0])
          + negN[0][2] * (negN[1][0] * negN[2][1] - negN[1][1] * negN[2][0]))
    return (e1 > 0 and e2 > 0 and e3 > 0), used, [float(x) for x in np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in N]))]

def blocks(m, k):
    tb = Tables(5, m)
    T = tb.table(k).reshape(-1, 3, 3)
    nz = np.any(T.reshape(len(T), -1) != 0, axis=1)
    U = np.unique(T[nz].reshape(-1, 9), axis=0).reshape(-1, 3, 3)
    return U

def margins_at(P, U):
    Uf = U.astype(float)
    Ah, _ = psd_sqrt(Uf)
    lam = np.linalg.eigvalsh(Ah @ P @ Ah)[:, -1]
    tr = np.einsum('nij,ij->n', Uf, P)
    return 1 - 2 * lam / tr

if __name__ == '__main__':
    jobs = [([1, 2, 3, 7, 1], 2), ([1, 2, 3, 7, 1], 3)]
    for m in ([1, 4, 1, 11, 34], [1, 6, 11, 11, 4], [1, 1, 6, 39, 11], [1, 14, 4, 1, 29], [1, 4, 31, 1, 6]):
        for k in (2, 3, 4):
            jobs.append((m, k))
    if len(sys.argv) > 1:
        jobs = jobs[int(sys.argv[1]):int(sys.argv[2])]
    for m, k in jobs:
        t0 = time.time()
        U = blocks(m, k)
        Uf = U.astype(float)
        An = Uf / np.trace(Uf, axis1=1, axis2=2)[:, None, None]
        Ah, Aih = psd_sqrt(An)
        rk = np.linalg.matrix_rank(Uf)
        if rk.min() < 3:
            print(f"Z_5 m = {m}, k = {k}: {len(U)} blocks; a block of rank {rk.min()} -> no form (exact by rank)", flush=True)
            continue
        s, best, duals, cuts_A, cuts_w = cutting_plane(An, Ah, Aih)
        Pbest = best[1]
        mg = margins_at(Pbest, U).min()
        ok, used, eig = exact_cert(U, cuts_A, cuts_w, duals)
        print(f"Z_5 m = {m}, k = {k}: {len(U)} blocks; LP max-min s = {s:+.6f}; margin at that P {mg:+.5f}; dual certificate with "
              f"{used} cuts: sum c B negative definite EXACTLY: {ok} (eigs {[round(x, 6) for x in eig]})  [{time.time() - t0:.1f}s]", flush=True)
