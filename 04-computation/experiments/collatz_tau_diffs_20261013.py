"""
collatz_tau_diffs_20261013.py -- HYP-9258 tests (ii) and (iii):
  (ii) the commonest nonzero differences tau(n+h) - tau(n) for h = 1, 2, 3, 4 (are they in the lattice <8, 19> of
       convergent denominators of log_2 3?), to N = 2^27;
  (iii) E[(-1)^tau] in dyadic blocks up to N = 2^29 (int16 storage), with the block-to-block spread.
Run: python3 collatz_tau_diffs_20261013.py
"""
import time
import numpy as np
from numba import njit

@njit(cache=True)
def tau_array16(N):
    tau = np.zeros(N + 8, dtype=np.int16)
    for n in range(2, N + 8):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
    return tau

t0 = time.time()
N = 1 << 29
tau = tau_array16(N)
print(f'tau to 2^29 in {time.time()-t0:.1f}s (max {int(tau.max())})')
M = 1 << 27
print('# (ii) commonest nonzero differences tau(n+h) - tau(n), n <= 2^27; share of nonzero differences in <8,19> (|d| <= 200)')
lattice = set()
for a in range(0, 30):
    for b in range(0, 15):
        v = 8 * a + 19 * b
        if 0 < v <= 200:
            lattice.add(v)
# also allow differences of the two generators (a or b negative)
for a in range(-30, 30):
    for b in range(-15, 15):
        v = 8 * a + 19 * b
        if 0 < v <= 200:
            lattice.add(v)
for h in (1, 2, 3, 4):
    d = tau[1 + h:M + 1 + h].astype(np.int32) - tau[1:M + 1].astype(np.int32)
    nz = d[d != 0]
    vals, cnt = np.unique(nz, return_counts=True)
    order = np.argsort(-cnt)[:10]
    top = [(int(vals[i]), int(cnt[i])) for i in order]
    absd = np.abs(nz)
    small = absd[absd <= 200]
    in_lat = np.isin(small, np.array(sorted(lattice)))
    print(f'h={h}: break density {len(nz)/M:.5f}; top {top}; share of |d|<=200 in <8,19>: {in_lat.mean():.4f} (|d|<=200 covers {len(small)/len(nz):.4f})')
print('# (iii) E[(-1)^tau] by dyadic block [2^(j-1), 2^j) and cumulative')
t = np.where(tau % 2 == 1, -1, 1).astype(np.int8); t[0] = 0
cum = 0.0; tot = 0
for j in range(16, 30):
    lo, hi = 1 << (j - 1), 1 << j
    blk = float(np.mean(t[lo:hi]))
    cum = float(np.mean(t[1:hi]))
    print(f' 2^{j:2d}: block mean {blk:+.6f}  cumulative mean {cum:+.6f}')
print('done')
