"""
collatz_tau_oscillation_20261013.py -- (a) the semigroup test for tau-differences (nonnegative combinations of 8 and 19
below the Frobenius number 125, where membership is non-trivial); (b) the dyadic-block oscillation of E[(-1)^tau]:
block means B_j over [2^(j-1), 2^j), split into even n (forced: tau(2n) = tau(n) + 1 gives -B_(j-1)) and odd n (O_j),
and quarter-octave means to see the shape in log_2 n.
Run: python3 collatz_tau_oscillation_20261013.py
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
print(f'tau to 2^29 in {time.time()-t0:.1f}s')
# (a) semigroup <8,19>_{>=0} below 125
semi = sorted({8 * a + 19 * b for a in range(0, 20) for b in range(0, 10) if 0 < 8 * a + 19 * b <= 125})
print('semigroup elements <= 125:', semi)
M = 1 << 27
for h in (1, 2, 3, 4):
    d = tau[1 + h:M + 1 + h].astype(np.int32) - tau[1:M + 1].astype(np.int32)
    nz = np.abs(d[d != 0])
    small = nz[nz <= 125]
    share_small = len(small) / len(nz)
    in_semi = np.isin(small, np.array(semi)).mean()
    vals, cnt = np.unique(small, return_counts=True)
    order = np.argsort(-cnt)[:14]
    top = [(int(vals[i]), int(cnt[i])) for i in order]
    non = [(int(v), int(c)) for v, c in zip(vals, cnt) if v not in semi]
    non.sort(key=lambda x: -x[1])
    print(f'h={h}: |d|<=125 covers {share_small:.4f} of nonzero differences; share of those in <8,19>_>=0: {in_semi:.5f}; top |d|: {top}; commonest non-semigroup values: {non[:6]}')
# (b) block oscillation
t = np.where(tau % 2 == 1, -1, 1).astype(np.int8); t[0] = 0
print('# dyadic blocks: B_j (all), E_j (even n), O_j (odd n), check E_j = -B_(j-1)')
prevB = None
for j in range(16, 30):
    lo, hi = 1 << (j - 1), 1 << j
    blk = t[lo:hi]
    B = float(blk.mean()); E = float(blk[0::2].mean()); O = float(blk[1::2].mean())
    chk = f'{E + prevB:+.6f}' if prevB is not None else '   n/a   '
    print(f' 2^{j:2d}: B={B:+.6f} E={E:+.6f} O={O:+.6f}  E_j + B_(j-1) = {chk}')
    prevB = B
print('# quarter-octave means of (-1)^tau, j = 22..29')
for j in range(22, 30):
    lo = 1 << (j - 1)
    qs = []
    for q in range(4):
        a = int(lo * 2 ** (q / 4)); b = int(lo * 2 ** ((q + 1) / 4))
        qs.append(float(t[a:b].mean()))
    print(f' 2^{j:2d}: ' + ' '.join(f'{v:+.5f}' for v in qs))
print('done')
