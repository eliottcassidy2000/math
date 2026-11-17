#!/usr/bin/env python3
"""profiles8.py: the profile (h, e, B, M) of every tournament class of orders 1-8 (gentourng order), saved to
profiles8.npz: for each class the count h, the end-counts e, the root responses B, the pair-cover matrix M(y, x) =
ordered pairs of nonempty paths covering T with the first ending at y and the second starting at x (THM-5266 (A)),
and the strong flag.  Orders 1 and 2 built in."""
import numpy as np
exec(open('nearjoin_check.py').read().split('# (1)')[0])
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def profile(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T)); full = N - 1
    e = [dp[full][x] for x in range(n)]; s = [dpr[full][x] for x in range(n)]
    M = [[0]*n for _ in range(n)]
    for m in range(1, N):
        c = full & ~m
        if c == 0: continue
        for y in range(n):
            a = dp[m][y]
            if not a or not (m >> y) & 1: continue
            for x in range(n):
                if (c >> x) & 1 and dpr[c][x]: M[y][x] += a * dpr[c][x]
    B = [s[x] + sum(M[y][x] for y in range(n)) for x in range(n)]
    return H[full], e, B, M
def classes1(n): return [[[0]]] if n == 1 else ([[[0, 1], [0, 0]]] if n == 2 else classes(n))
out = {}
for n in range(1, 9):
    cl = classes1(n); k = len(cl)
    A = np.zeros((k, n, n), dtype=np.int8); h = np.zeros(k, dtype=np.int64); e = np.zeros((k, n), dtype=np.int64)
    B = np.zeros((k, n), dtype=np.int64); M = np.zeros((k, n, n), dtype=np.int64); st = np.zeros(k, dtype=np.int8)
    for i, T in enumerate(cl):
        hh, ee, BB, MM = profile(T); A[i] = T; h[i] = hh; e[i] = ee; B[i] = BB; M[i] = MM; st[i] = strong(T)
    out[f'A{n}'] = A; out[f'h{n}'] = h; out[f'e{n}'] = e; out[f'B{n}'] = B; out[f'M{n}'] = M; out[f'strong{n}'] = st
    print(f"order {n}: {k} classes profiled; max h {h.max()}, max M entry {M.max()}")
np.savez_compressed('profiles8.npz', **out)
