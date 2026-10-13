"""
collatz_merge_lags_20261013.py -- equal-time merging of the orbits of n and n+h, and what it does to the
two-point correlations of the total-stopping-time parity.

Part A (FINITE-EXACT, periodic): for each residue r mod 2^K, follow the T-orbits of r and r+h (as the residue
  classes they represent) for K steps and record the first k <= K at which T^k(n) = T^k(n+h) holds for the whole
  class: the parity words w, w' of n and n+h give T^k(n) = (3^a n + c)/2^k, T^k(n+h) = (3^a' (n+h) + c')/2^k, and
  the class-wide identity needs a = a' and c = 3^a h + c'.  The density of classes merged by depth K is the integer
  density of n with an equal-time merge by depth K (periodicity), and it tends to 1 (THM-4581 and its offsets).
Part B (NUMERICAL): for n <= N follow both integer orbits until they coincide at equal time or both reach 1;
  record merged fraction, merge-time distribution, and the correlation of (-1)^tau split into merged pairs
  (always +1) and unmerged pairs.
Run: python3 collatz_merge_lags_20261013.py <K> <N_log2>
"""
import sys, time
import numpy as np
from numba import njit, prange

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
NL = int(sys.argv[2]) if len(sys.argv) > 2 else 26
LAGS = (1, 2, 3, 4, 5, 6, 8, 16)

@njit(cache=True)
def merge_depth_periodic(K, h):
    """first class-wide equal-time merge depth of (r, r+h) for every residue r mod 2^K (0 = none by K)."""
    M = 1 << K
    out = np.zeros(M, dtype=np.int8)
    for r in range(M):
        n = r; m = r + h
        a = 0; ap = 0
        # c, c' tracked exactly via the actual integers: T^k(n) = (3^a n + c)/2^k, so c = 2^k T^k(n) - 3^a n
        for k in range(1, K + 1):
            if n & 1:
                n = (3 * n + 1) >> 1; a += 1
            else:
                n >>= 1
            if m & 1:
                m = (3 * m + 1) >> 1; ap += 1
            else:
                m >>= 1
            if a == ap and n == m:
                # equal values at equal time with equal 3-exponents: class-wide identity c = 3^a h + c'
                out[r] = k
                break
    return out

@njit(cache=True)
def tau_array(N):
    tau = np.zeros(N + 1, dtype=np.int32)
    for n in range(2, N + 1):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
    return tau

@njit(cache=True)
def merge_integers(N, h, tau):
    """for 1 <= n <= N: time of the first equal-time coincidence of the two orbits, or -1 if one of them reaches 1 first."""
    out = np.full(N + 1, -1, dtype=np.int32)
    for n in range(1, N + 1):
        x = n; y = n + h; k = 0
        while True:
            if x == y:
                out[n] = k; break          # equal-time coincidence before either orbit reached 1
            if x == 1 or y == 1:
                break                      # one orbit reached 1 first: no merge before the end
            x = (3 * x + 1) >> 1 if (x & 1) else (x >> 1)
            y = (3 * y + 1) >> 1 if (y & 1) else (y >> 1)
            k += 1
            if k > 100000:
                break
    return out

print(f'# Part A: class-wide equal-time merge of (n, n+h) by depth K={K} (periodic census over 2^K residues)')
for h in LAGS:
    t0 = time.time()
    d = merge_depth_periodic(K, h)
    merged = float(np.mean(d > 0))
    hist = np.bincount(d[d > 0], minlength=K + 1)
    firsts = ' '.join(f'{k}:{hist[k]/len(d):.4f}' for k in range(1, 13) if hist[k])
    print(f'h={h:2d}: merged by depth {K}: {merged:.5f}; first-merge depth masses {firsts} ... ({time.time()-t0:.1f}s)', flush=True)

print(f'\n# Part B: integer orbits to N=2^{NL}: equal-time coincidence before either orbit reaches 1 (parity of tau then agrees automatically)')
N = 1 << NL
tau = tau_array(N + max(LAGS))
t = np.where(tau % 2 == 1, -1, 1).astype(np.int8)
for h in LAGS:
    t0 = time.time()
    mt = merge_integers(N, h, tau)
    merged = mt >= 0
    pm = float(np.mean(merged[1:]))
    prod = t[1:N + 1] * t[1 + h:N + 1 + h]
    c_all = float(np.mean(prod))
    c_merged = float(np.mean(prod[merged[1:]])) if merged[1:].any() else float('nan')
    c_unm = float(np.mean(prod[~merged[1:]])) if (~merged[1:]).any() else float('nan')
    eq_tau = float(np.mean((tau[1:N + 1] == tau[1 + h:N + 1 + h])[merged[1:]])) if merged[1:].any() else float('nan')
    mean_mt = float(np.mean(mt[1:][merged[1:]])) if merged[1:].any() else float('nan')
    print(f'h={h:2d}: P(merge)={pm:.5f}  corr(all)={c_all:+.5f}  corr(merged)={c_merged:+.5f}  corr(unmerged)={c_unm:+.5f}  '
          f'tau equal on merged: {eq_tau:.4f}  mean merge time {mean_mt:.2f}  ({time.time()-t0:.1f}s)', flush=True)
    for NB in (1 << 20, 1 << 23):
        pmB = float(np.mean(merged[1:NB + 1])); cB = float(np.mean(prod[:NB]))
        print(f'       N=2^{int(np.log2(NB))}: P(merge)={pmB:.5f} corr={cB:+.5f}')
print('done')
