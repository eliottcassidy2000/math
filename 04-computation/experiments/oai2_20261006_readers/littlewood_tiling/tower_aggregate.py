#!/usr/bin/env python3
"""Exact aggregates for the tower rows, comparison with Walsh (Sylvester) rows,
Rudin-Shapiro, Legendre (Paley) rows; nearest-Walsh distance of tower rows."""
import numpy as np, sys
from tower_merit import tower, autocorr_int

def sylvester(k):
    S = np.array([[1]], dtype=np.int8)
    for _ in range(k):
        S = np.block([[S, S], [S, -S]]).astype(np.int8)
    return S

def fwht(X):
    X = X.astype(np.int64).copy()
    m, N = X.shape
    h = 1
    while h < N:
        X = X.reshape(m, N // (2 * h), 2, h)
        a = X[:, :, 0, :].copy(); b = X[:, :, 1, :].copy()
        X[:, :, 0, :] = a + b; X[:, :, 1, :] = a - b
        X = X.reshape(m, N)
        h *= 2
    return X

def rows_sumsq(X, bs):
    out = []
    for s in range(0, X.shape[0], bs):
        C = autocorr_int(X[s:s+bs])
        out.append((C[:, 1:] ** 2).sum(axis=1))
    return np.concatenate(out)

def rudin_shapiro(k):
    P = np.array([1], dtype=np.int8); Q = np.array([1], dtype=np.int8)
    for _ in range(k):
        P, Q = np.concatenate([P, Q]), np.concatenate([P, -Q])
    return P, Q

if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("k N  sum_rows sumC^2 (tower) | ratio to prev | Walsh: Fmax Fmean | tower: max |WHT|/N (closest Walsh agreement) min..max over rows | RS F")
    prev = None
    for k in range(1, kmax + 1):
        H = tower(k); N = H.shape[0]
        bs = max(1, (1 << 22) // N)
        s = rows_sumsq(H, bs)
        tot = int(s.sum())
        S = sylvester(k)
        sw = rows_sumsq(S, bs)
        Fw = N * N / (2.0 * sw)
        W = fwht(H)
        agree = np.abs(W).max(axis=1) / N
        P, Q = rudin_shapiro(k)
        sp = rows_sumsq(P[None, :], 1)[0]
        Frs = N * N / (2.0 * sp)
        print(f"{k:2d} {N:6d} {tot:22d} | {('%.5f' % (tot/prev)) if prev else '   -   '} | {Fw.max():.4f} {Fw.mean():.4f} | "
              f"{agree.min():.4f}..{agree.max():.4f} mean {agree.mean():.4f} | {Frs:.4f}", flush=True)
        prev = tot
