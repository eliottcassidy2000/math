"""Tilted cross-energies N_t(c,c',D) = #{(u,u') in X_t^2, u != u': g(u)+c u - g(u') - c' u' = D}.
   Report sup_D N * (c + c' + 1) / 2^t over small tilt pairs."""
import numpy as np, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oalib import gtable

t = int(sys.argv[1]) if len(sys.argv) > 1 else 12
CMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 24
g = gtable(t)
u = np.arange(1 << (t - 1), 1 << t, dtype=np.int64)
gu = g[u]
best = []
worst = 0.0
for c in range(CMAX):
    A = gu + c * u
    for cp in range(CMAX):
        B = gu + cp * u
        # all differences A_i - B_j, excluding i == j
        lo = min(A.min() - B.max(), 0)
        diffs = (A[:, None] - B[None, :]).ravel()
        cnt = np.bincount(diffs - lo)
        # remove diagonal i == j
        dd = A - B - lo
        np.subtract.at(cnt, dd, 1)
        s = cnt.max()
        val = s * (c + cp + 1) / 2 ** t
        best.append((val, c, cp, int(np.argmax(cnt) + lo)))
        worst = max(worst, val)
best.sort(reverse=True)
print('t', t, 'max over tilts < %d of sup_D N_off*(c+c\'+1)/2^t = %.3f' % (CMAX, worst))
print('top 12 (value, c, c\', argmax D):', [(round(v, 3), c, cp, D) for v, c, cp, D in best[:12]])
print('c=c\'=0:', [b for b in best if b[1] == 0 and b[2] == 0])
