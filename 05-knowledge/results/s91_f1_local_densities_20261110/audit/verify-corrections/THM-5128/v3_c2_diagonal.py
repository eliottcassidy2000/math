#!/usr/bin/env python3
"""Verifier check 3 (THM-5128 (C), VERIFIED paragraph): the remainder of (C2) at p = 5, sheet 0, for X = 2^18,
2^20, 2^22, with D_sigma as defined in (C2) (triples with h1 != h2 whose smaller state is p-divisible) and with
the session's variant (first_order_exact.py also counts the diagonal h1 = h2 in F_small)."""
import numpy as np

LOGH = 22
H = 1 << LOGH
psi = np.zeros(H, dtype=np.int64)
for n in range(1, LOGH + 1):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    par = h >> 1
    psi[lo:hi] = np.where(h & 1 == 0, psi[par] + par, psi[par] + (par >> 1) + 1)
hh = np.arange(H, dtype=np.int64)
x = 2 * psi + 1
y = 2 * hh + 1
g = np.gcd(x, y)

for LOGX in (18, 20, 22):
    X = 1 << LOGX
    for p in (5,):
        cnt = np.zeros(p)
        Dstrict = 0.0
        Ddiag = 0.0
        lim = X + 1                               # sheet 0: M <= X + 1
        h1 = 0
        T = 0.0
        while True:
            y1 = 2 * h1 + 1
            if y1 * y1 > lim:
                break
            h2max = min(H, (lim // y1 + 1) // 2 + 1)
            Mv = int(x[h1]) * x[h1:h2max] + y1 * y[h1:h2max]
            sel = np.nonzero(Mv <= lim)[0]
            if len(sel):
                wt = np.where(sel == 0, 1.0, 2.0)
                cnt += np.bincount(Mv[sel] % p, weights=wt, minlength=p)
                T += wt.sum()
                if g[h1] % p == 0:
                    Dstrict += wt[sel != 0].sum()
                    Ddiag += wt.sum()
            h1 += 1
        out = []
        for D in (Dstrict, Ddiag):
            pred = np.full(p, (T - D) / p)
            pred[0] += D
            out.append(p * np.abs(cnt - pred).max() / T)
        print("X = 2^%d, p = %d, sheet 0: #T = %d; p max|R| = %.6f with D of (C2) (h1 != h2), %.6f with the diagonal counted"
              % (LOGX, p, T, out[0], out[1]))
