#!/usr/bin/env python3
"""Meet-depth anatomy of first-arrival pair sums (mac-mini-2026-10-18-sixcolour). Two sources x1, x2 first meet at the
node y = P^k1(x1) = P^k2(x2) (paths disjoint before y; k = 0 means the source is y itself). The pair contributes
(-1)^(t1 + t2) c_0(x1) c_0(x2) to the clock-lag Gram at lag |t1 - t2| of every stopping set entered at or after y,
with t the literal clock from the source to y. Prints the sums / N by (k1 <= k2, lag D) for k1, k2 <= K, restricted
to meets that happen at or before first entry into the stopping set V (file VFILE): both paths must avoid V strictly
before y; y itself may be the entry seed. With VFILE = none, all meets count (V = {1}).
Usage: python3 pair_series_meets.py N K VFILE"""
import sys
import numpy as np
from collections import defaultdict

def c0_array(N):
    M = 4 * N
    m = np.arange(1, M + 1, 2, dtype=np.int64)
    r = m % 3
    K = np.where(r == 2, 2, -1) + np.where(m < 2 * N, np.where(r == 0, 1, np.where(r == 1, -2, 1)), 0)
    s = m.copy()
    while True:
        msk = (s & 7) == 5
        if not msk.any(): break
        s[msk] = (s[msk] - 1) >> 2
    c = np.zeros((M + 1) // 2, dtype=np.int64)
    np.add.at(c, (s - 1) // 2, K)
    return c

def step(x):
    t = 3 * x + 1
    a = np.where(t % 4 == 0, 2, 1)
    z = t >> a
    strips = np.zeros_like(z)
    while True:
        msk = (z & 7) == 5
        if not msk.any(): break
        z[msk] = (z[msk] - 1) >> 2; strips[msk] += 1
    return z, a + 2 * strips

N, KMAX, VF = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
Vset = np.array(sorted(int(l.split()[0]) for l in open(VF) if l.strip() and not l.startswith('#')), dtype=np.int64) if VF != 'none' else np.array([1], dtype=np.int64)
c = c0_array(N); idx = np.flatnonzero(c); x = 2 * idx + 1; cv = c[idx]
ys = [x]; ts = [np.zeros_like(x)]
for k in range(KMAX):
    z, tau = step(ys[-1]); ys.append(z); ts.append(ts[-1] + tau)
rec_y, rec_b, rec_k, rec_t, rec_c = [], [], [], [], []
alive = np.ones_like(x, dtype=bool)                  # no node strictly before y_k lies in V
for k in range(KMAX + 1):
    keep = alive.copy()
    alive &= ~np.isin(ys[k], Vset)                  # after y_k, the path is stopped if y_k is in V
    rec_y.append(ys[k][keep]); rec_b.append((ys[k - 1] if k else -np.ones_like(x))[keep]); rec_k.append(np.full(keep.sum(), k)); rec_t.append(ts[k][keep]); rec_c.append(cv[keep])
Y = np.concatenate(rec_y); B = np.concatenate(rec_b); Kk = np.concatenate(rec_k); T = np.concatenate(rec_t); C = np.concatenate(rec_c)
o = np.lexsort((B, Y)); Y, B, Kk, T, C = Y[o], B[o], Kk[o], T[o], C[o]
brk = np.flatnonzero(np.diff(Y)) + 1; st = np.concatenate(([0], brk)); en = np.concatenate((brk, [len(Y)]))
acc = defaultdict(float)
for g in np.flatnonzero(en - st >= 2):
    s0, e0 = st[g], en[g]
    bb, kk, tt, cc = B[s0:e0], Kk[s0:e0], T[s0:e0], C[s0:e0]
    n = e0 - s0
    for i in range(n):
        for j in range(i + 1, n):
            if bb[i] == bb[j] and bb[i] != -1: continue          # same branch: they met earlier
            k1, k2 = (kk[i], kk[j]) if kk[i] <= kk[j] else (kk[j], kk[i])
            D = abs(int(tt[i]) - int(tt[j]))
            acc[(int(k1), int(k2), D)] += (-1) ** (int(tt[i]) + int(tt[j])) * int(cc[i]) * int(cc[j])
print('N = %d, K = %d, stopping set %s; pair sums / N by (k1, k2) and lag D (|value| >= 0.002 shown)' % (N, KMAX, VF))
tot = defaultdict(float)
for (k1, k2, D), v in sorted(acc.items()):
    tot[D] += v / N
for k1 in range(KMAX + 1):
    for k2 in range(k1, KMAX + 1):
        row = [(D, acc[(k1, k2, D)] / N) for D in range(0, 40) if abs(acc.get((k1, k2, D), 0)) / N >= 0.002]
        if row: print('  (k1,k2)=(%d,%d): ' % (k1, k2) + '  '.join('D%d %+.4f' % (D, v) for D, v in row))
print('  totals by lag D (all k1, k2 <= K): ' + '  '.join('D%d %+.4f' % (D, v) for D, v in sorted(tot.items()) if abs(v) >= 0.002))
