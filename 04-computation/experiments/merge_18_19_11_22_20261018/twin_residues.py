#!/usr/bin/env python3
"""Same-clock 'clock twins' of the signed ROOT transport (mac-mini-2026-10-18-sixcolour): pairs of sources x1, x2 whose
paths first meet after the same number k of steps with the same literal clock t. They sit at the same scale, so
x1 - x2 is a fixed integer Delta for each pair of backward words. Because c_0 has zero mean over residues mod 3 at
fixed chain level, a twin pair contributes about +E[c_0^2] if 3 | Delta and -E[c_0^2]/2 otherwise. This script
measures, for the full transport (V = {1}) at depth k <= KMAX, the number of twin pairs, the share with 3 | Delta, and
the sums of c_0(x1)c_0(x2) split by Delta mod 3. Usage: python3 twin_residues.py N KMAX"""
import sys
import numpy as np
from collections import defaultdict
sys.path.insert(0, '../first_arrival_20261018')
src = open('../first_arrival_20261018/pair_series_meets.py').read()
ns = {}; exec(src.split("N, KMAX, VF =")[0], ns)
c0_array, step = ns['c0_array'], ns['step']
N, KMAX = int(sys.argv[1]), int(sys.argv[2])
c = c0_array(N); idx = np.flatnonzero(c); x = 2 * idx + 1; cv = c[idx]
ys = [x]; ts = [np.zeros_like(x)]
for k in range(KMAX):
    z, tau = step(ys[-1]); z = np.where(ys[-1] == 1, 1, z); tau = np.where(ys[-1] == 1, 2, tau)
    ys.append(z); ts.append(ts[-1] + tau)
for k in range(1, KMAX + 1):
    key = ys[k]; br = ys[k - 1]; tt = ts[k]
    o = np.lexsort((tt, key)); K_, B_, T_, X_, C_ = key[o], br[o], tt[o], x[o], cv[o]
    # groups with equal (node, clock)
    brk = np.flatnonzero((np.diff(K_) != 0) | (np.diff(T_) != 0)) + 1
    st = np.concatenate(([0], brk)); en = np.concatenate((brk, [len(K_)]))
    npairs = [0, 0, 0]; sums = [0, 0, 0]
    for g in np.flatnonzero(en - st >= 2):
        s0, e0 = st[g], en[g]
        for i in range(s0, e0):
            for j in range(i + 1, e0):
                if B_[i] == B_[j]: continue                     # met earlier
                dm = int((X_[i] - X_[j]) % 3)
                npairs[dm] += 1; sums[dm] += int(C_[i]) * int(C_[j])
    tot = sum(sums)
    print('depth %d: twin pairs %d (3|Delta: %d = %.3f, else %d); sum c0c0/N: 3|Delta %+.4f, Delta=1 mod 3 %+.4f, Delta=2 mod 3 %+.4f, total %+.4f'
          % (k, sum(npairs), npairs[0], npairs[0] / max(1, sum(npairs)), npairs[1] + npairs[2], sums[0] / N, sums[1] / N, sums[2] / N, tot / N))
