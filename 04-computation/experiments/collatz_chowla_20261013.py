"""
collatz_chowla_20261013.py -- two-point correlations of Collatz orbit functionals along (n, n+h).

Objects (T(n) = n/2 for even n, (3n+1)/2 for odd n):
  kappa(n) = coefficient stopping time (Terras): least k with 3^{a_k} < 2^k, a_k = number of odd steps among the
             first k; it depends only on n mod 2^k and equals the stopping time sigma(n) for all n outside a
             finite set in each residue class.
  s(n) = (-1)^kappa(n)        (stopping-time parity; sigma-measurable)
  tau(n) = number of T-steps from n to 1 (total stopping time; NOT sigma-measurable)
  omega(n) = number of odd steps on the way from n to 1 (the Collatz 'Liouville' count)
  t(n) = (-1)^tau(n), lam(n) = (-1)^omega(n)

Part 1 (FINITE-EXACT, depth k): over all residues r mod 2^k, the truncated s_k (0 if kappa > k) and the exact
  correlation E_k(h) = 2^-k sum_r s_k(r) s_k(r+h); the undetermined mass P(kappa > k) bounds the truncation error.
Part 2 (NUMERICAL): over n <= N the integer averages of s, t, lam and their correlations at lags h = 1..H.
Run: python3 collatz_chowla_20261013.py <kmax> <N_log2> <H>
"""
import sys, time
import numpy as np
from numba import njit, prange

kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
NL = int(sys.argv[2]) if len(sys.argv) > 2 else 28
H = int(sys.argv[3]) if len(sys.argv) > 3 else 40

LOG2_3 = np.log2(3.0)

@njit(cache=True)
def kappa_all(k):
    """coefficient stopping time of every residue class mod 2^k (0 means > k)."""
    M = 1 << k
    out = np.zeros(M, dtype=np.int8)
    for r in range(M):
        n = r
        a = 0
        kap = 0
        for j in range(1, k + 1):
            if n & 1:
                n = (3 * n + 1) >> 1
                a += 1
            else:
                n >>= 1
            # 3^a < 2^j  <=>  a*log2(3) < j
            if a * LOG2_3 < j:
                kap = j
                break
        out[r] = kap
    return out

@njit(parallel=True, cache=True)
def corr_periodic(s, hs):
    M = s.shape[0]
    out = np.zeros(hs.shape[0])
    for i in prange(hs.shape[0]):
        h = hs[i]
        acc = 0
        for r in range(M):
            acc += s[r] * s[(r + h) & (M - 1)]
        out[i] = acc / M
    return out

@njit(cache=True)
def total_stopping(N):
    """tau[n] = T-steps to reach 1, omega[n] = odd steps, for 1 <= n <= N (memoised below n)."""
    tau = np.zeros(N + 1, dtype=np.int32)
    om = np.zeros(N + 1, dtype=np.int32)
    for n in range(2, N + 1):
        m = n
        steps = 0
        odd = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
                odd += 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
        om[n] = odd + om[m]
    return tau, om

@njit(parallel=True, cache=True)
def corr_int(x, hs):
    N = x.shape[0]
    out = np.zeros(hs.shape[0])
    for i in prange(hs.shape[0]):
        h = hs[i]
        acc = 0
        for n in range(1, N - h):
            acc += x[n] * x[n + h]
        out[i] = acc / (N - 1 - h)
    return out

hs = np.arange(1, H + 1, dtype=np.int64)
print(f'# Part 1: periodic (Haar) correlations of stopping-time parity, truncated at depth k; lags 1..{H}')
print('k  P(kappa>k)   E_k[s]      ' + ' '.join(f'h={h}' for h in (1, 2, 3, 4, 5, 6, 8, 16, 32)))
rows = {}
for k in range(12, kmax + 1):
    t0 = time.time()
    kap = kappa_all(k)
    undet = float(np.mean(kap == 0))
    s = np.where(kap == 0, 0, np.where(kap % 2 == 1, -1, 1)).astype(np.int8)
    mean_s = float(np.mean(s))
    c = corr_periodic(s, hs)
    rows[k] = (undet, mean_s, c)
    sel = [1, 2, 3, 4, 5, 6, 8, 16, 32]
    print(f'{k:2d} {undet:.6f} {mean_s:+.6f} ' + ' '.join(f'{c[h-1]:+.5f}' for h in sel if h <= H) + f'   ({time.time()-t0:.1f}s)', flush=True)
np.save('collatz_chowla_20261013_periodic.npy', np.array([rows[k][2] for k in sorted(rows)]))

print(f'\n# Part 2: integer averages to N = 2^{NL} of s (stopping-time parity), t (total-stopping-time parity), lam (odd-step parity)')
N = 1 << NL
t0 = time.time()
tau, om = total_stopping(N)
print(f'total stopping times computed in {time.time()-t0:.1f}s; max tau {int(tau.max())} at n={int(np.argmax(tau))}; max omega {int(om.max())}', flush=True)
# sigma-parity on integers: kappa with full depth (use tau machinery: stopping time = first drop below n)
@njit(cache=True)
def stopping_parity(N):
    out = np.zeros(N + 1, dtype=np.int8)
    for n in range(2, N + 1):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        out[n] = -1 if (steps & 1) else 1
    return out
s_int = stopping_parity(N)
t_int = np.where(tau % 2 == 1, -1, 1).astype(np.int8); t_int[0] = 0
l_int = np.where(om % 2 == 1, -1, 1).astype(np.int8); l_int[0] = 0
for name, x in (('s (sigma parity)', s_int), ('t (tau parity)', t_int), ('lam (odd-step parity)', l_int)):
    m = float(np.mean(x[1:]))
    c = corr_int(x, hs)
    print(f'{name}: mean {m:+.6f}; correlations h=1..{H}:')
    print('   ' + ' '.join(f'{v:+.5f}' for v in c))
    # block means to see N-dependence
    for NB in (1 << 20, 1 << 24, 1 << 26):
        if NB < N:
            cb = corr_int(x[:NB + 1], hs[:8])
            print(f'   N=2^{int(np.log2(NB))}: mean {float(np.mean(x[1:NB+1])):+.6f} corr h=1..8 ' + ' '.join(f'{v:+.5f}' for v in cb))
np.save('collatz_chowla_20261013_int.npy', np.vstack([corr_int(s_int, hs), corr_int(t_int, hs), corr_int(l_int, hs)]))
print('done')
