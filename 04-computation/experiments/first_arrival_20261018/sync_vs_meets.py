#!/usr/bin/env python3
"""Synchronous cross terms (THM-4667) versus equal-depth meets (THM-4690 anatomy), mac-mini-2026-10-18-sixcolour.
Moving every source one transport step at a time, E_d - E_(d-1) = 2 * (sum over pairs whose paths first meet after
exactly d steps on both sides). This script computes E_d directly (signed synchronous transport, every source moved d
steps; the ROOT self-loop keeps 1 at 1) and the equal-depth meet sums split by literal clock lag, and compares."""
import sys
import numpy as np
from collections import defaultdict
sys.path.insert(0, '.')
src = open('pair_series_meets.py').read()
ns = {}; exec(src.split("N, KMAX, VF =")[0], ns)
c0_array, step = ns['c0_array'], ns['step']
N, DMAX = int(sys.argv[1]), int(sys.argv[2])
c = c0_array(N); idx = np.flatnonzero(c); x = 2 * idx + 1; cv = c[idx]
ys = [x]; ts = [np.zeros_like(x)]
for k in range(DMAX):
    z, tau = step(ys[-1]); z = np.where(ys[-1] == 1, 1, z); tau = np.where(ys[-1] == 1, 2, tau)   # 1 -> 1 with a = 2
    ys.append(z); ts.append(ts[-1] + tau)
def energy(d):
    sgn = np.where(ts[d] % 2 == 0, 1, -1)
    u, inv = np.unique(ys[d], return_inverse=True)
    acc = np.zeros(len(u), dtype=np.int64); np.add.at(acc, inv, sgn * cv)
    return float((acc.astype(np.float64) ** 2).sum())
E = [energy(d) for d in range(DMAX + 1)]
print('N = %d: synchronous E_d/N = %s' % (N, ', '.join('%.6f' % (e / N) for e in E)))
for d in range(1, DMAX + 1):
    # equal-depth meets at depth d: group by (y_d), branch y_(d-1); pairs from different branches
    key = ys[d]; br = ys[d - 1]; sg = np.where(ts[d] % 2 == 0, 1, -1); tt = ts[d]
    o = np.lexsort((br, key)); key, br, sg, tt, cc = key[o], br[o], sg[o], tt[o], cv[o]
    brk = np.flatnonzero(np.diff(key)) + 1; st = np.concatenate(([0], brk)); en = np.concatenate((brk, [len(key)]))
    lag = defaultdict(float); tot = 0.0
    for g in np.flatnonzero(en - st >= 2):
        s0, e0 = st[g], en[g]
        for i in range(s0, e0):
            for j in range(i + 1, e0):
                if br[i] == br[j]: continue
                v = float(sg[i] * sg[j] * cc[i] * cc[j]); tot += v; lag[abs(int(tt[i]) - int(tt[j]))] += v
    print('  d = %d: (E_d - E_(d-1))/N = %+.6f ; 2 * equal-depth meet sum / N = %+.6f ; by lag: %s' % (d, (E[d] - E[d - 1]) / N, 2 * tot / N,
          '  '.join('D%d %+.4f' % (D, 2 * v / N) for D, v in sorted(lag.items()) if abs(v) / N >= 0.001)))
