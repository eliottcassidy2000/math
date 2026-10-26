# Audit of HYP-9391 (c): size of the fitted low-order terms {1, q, q^2, a84} (and |b|^2) in the lane's projection
# ([307, 4999], weight q^-5), as normalised contributions c.col(q)/q^(5/2) at the primes q >= 2000: the HYP says they
# "contribute less than about 0.05 there".  Uses part_c_analysis.py's fit.  Writes lowblock_contrib.out.
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = ['x', sys.argv[1], sys.argv[2], sys.argv[3]]
src = open(os.path.join(HERE, 'part_c_analysis.py')).read().split('# ---------- 1. reproduce')[0]
exec(src)
qA = allq[(allq >= 307) & (allq <= 4999)]
cA = fit(qA)
B, Y = design(qA)
w = qA.astype(float) ** -2.5
sel = qA >= 2000
lines = [f"{sel.sum()} primes in [2000, 4999]; normalised fitted contributions (max |.| over these primes)"]
for i, nm in enumerate(['Y02', 'Y05', 'Y09', 'Y18']):
    parts = {n: np.max(np.abs(B[sel, j] * cA[j, i] * w[sel])) for j, n in enumerate(NAMES)}
    low = np.max(np.abs((B[sel][:, [0, 1, 2, 4]] @ cA[[0, 1, 2, 4], i]) * w[sel]))
    lowb = np.max(np.abs((B[sel][:, [0, 1, 2, 4, 8]] @ cA[[0, 1, 2, 4, 8], i]) * w[sel]))
    lines.append(f"  {nm}: " + ', '.join(f"{n}: {v:.3f}" for n, v in parts.items()) + f"; block {{1,q,q^2,a84}}: {low:.3f}; with |b|^2: {lowb:.3f}")
print('\n'.join(lines))
open(os.path.join(HERE, 'lowblock_contrib.out'), 'w').write('\n'.join(lines) + '\n')
