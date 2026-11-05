#!/usr/bin/env python3
"""fourier.py -- exponential sums of the ribbon counts: S_n(theta) = sum_{|u| = n} exp(2 pi i theta H(u)),
normalised by the number of words 4^(n-1); the decay away from rationals with small denominators is the
analytic form of the collective invariant (HYP-9474).  Also the Lyapunov exponent of log H."""
import sys, math
from collections import Counter
import numpy as np
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
def load_intrep(path):
    with open(path) as f:
        r, nl = map(int, f.readline().split())
        y0 = list(map(int, f.readline().split())); hc = list(map(int, f.readline().split()))
        A = {}
        for _ in range(nl):
            mask = int(f.readline().strip()); A[mask] = [list(map(int, f.readline().split())) for _ in range(r)]
    return r, y0, hc, A
r, y0, hc, A = load_intrep('span2_intrep.txt')
Am = {m: np.array(M, dtype=np.int64) for m, M in A.items()}
masks = sorted(Am); hidx = hc.index(1)
cur = Counter({tuple(y0): 1})
thetas = [(1, 3), (1, 5), (1, 7), (1, 8), (1, 9), (1, 16), (2, 7), (3, 11), (5, 13)]
irr = [math.sqrt(2) - 1, (math.sqrt(5) - 1) / 2, math.pi - 3, 0.123456789]
for n in range(1, NMAX + 1):
    nxt = Counter()
    for y, c in cur.items():
        v = np.array(y, dtype=np.int64)
        for m in masks:
            if n >= 2 and (m >> (n - 1)): continue
            if n == 1 and m: continue
            nxt[tuple((Am[m] @ v).tolist())] += c
    cur = nxt
    if n >= 8:
        Hs = np.array([y[hidx] for y in cur], dtype=np.float64); ws = np.array([c for c in cur.values()], dtype=np.float64)
        W = ws.sum()
        lam = (ws * np.log(Hs)).sum() / W / n
        parts = []
        for a, q in thetas:
            S = (ws * np.exp(2j * math.pi * a / q * Hs)).sum()
            parts.append(f"{a}/{q}:{abs(S)/W:.3f}")
        parts2 = [f"{abs((ws * np.exp(2j * math.pi * t * Hs)).sum())/W:.4f}" for t in irr]
        print(f"n={n}: words {int(W)}, (1/n) E log H = {lam:.4f}; |S_n(theta)|/4^(n-1) at rationals {' '.join(parts)}; at irrationals {' '.join(parts2)}", flush=True)
