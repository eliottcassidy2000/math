#!/usr/bin/env python3
"""Depth-one meets of the signed ROOT transport (mac-mini-2026-10-18-sixcolour): pairs of distinct sources x, x' with the
same image y = P(x) = P(x'). Their step clocks tau, tau' are listed in THM-4690(2): 1 (the unpaired preimage when
y = 2 mod 3) and 2j+2, 2j+3 for B^j(y) = 1 (mod 3). The pair contributes sigma sigma' c_0 c_0' = (-1)^(tau+tau') c_0 c_0'
to the lag-|tau - tau'| Gram of every stopping set beyond y. Prints, per lag D, the total of these contributions / N."""
import sys
import numpy as np
from fractions import Fraction as Fr

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

def run(N, DMAX=20):
    c = c0_array(N)
    idx = np.flatnonzero(c)
    x = 2 * idx + 1
    t = 3 * x + 1
    a = np.where(t % 4 == 0, 2, 1)
    z = t >> a
    strips = np.zeros_like(z)
    while True:
        msk = (z & 7) == 5
        if not msk.any(): break
        z[msk] = (z[msk] - 1) >> 2; strips[msk] += 1
    tau = a + 2 * strips
    cv = c[idx]
    order = np.argsort(z, kind='stable')
    zs, ts, cs = z[order], tau[order], cv[order]
    brk = np.flatnonzero(np.diff(zs)) + 1
    starts = np.concatenate(([0], brk)); ends = np.concatenate((brk, [len(zs)]))
    acc = np.zeros(DMAX + 1)
    multi = np.flatnonzero(ends - starts >= 2)
    for g in multi:
        s0, e0 = starts[g], ends[g]
        tt = ts[s0:e0]; cc = cs[s0:e0]
        for i in range(len(tt)):
            for k in range(i + 1, len(tt)):
                D = abs(int(tt[i]) - int(tt[k]))
                if D <= DMAX:
                    acc[D] += (-1) ** (int(tt[i]) + int(tt[k])) * int(cc[i]) * int(cc[k])
    return acc / N

for N in [int(a) for a in sys.argv[1:]] or [2 ** 16, 2 ** 18, 2 ** 20]:
    acc = run(N)
    print('N = %8d  depth-one meet pairs, signed sum / N by lag D:' % N)
    print('   ' + '  '.join('D=%d:%+.4f' % (D, acc[D]) for D in range(1, 21) if abs(acc[D]) > 1e-9))
