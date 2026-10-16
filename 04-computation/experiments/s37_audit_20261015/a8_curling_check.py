#!/usr/bin/env python3
"""S37 audit, part 8: curling numbers recomputed by DIRECT block comparison (no incremental period table).
For every start over {2,3} of length n <= NMAX, c(S) = number of appended terms up to and including the first 1.
Compared with the session's output and with OEIS A094004 (a(n) = longest length before the 1 is produced, starting
terms included), fetched 2026-10-15: 1, 4, 5, 8, 9, 14, 15, 66, 68, 70, 123, 124, 125, 132, 133, 134, 135, 136,
138, 139, 140, 142, 143, 144, 145.  Gijswijt's sequence first occurrences checked to 1000 terms.
"""
import sys, time
import numpy as np
import numba as nb

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
CAP = 4000
A094004 = [1, 4, 5, 8, 9, 14, 15, 66, 68, 70, 123, 124, 125, 132, 133, 134, 135, 136, 138, 139, 140, 142, 143, 144, 145]
SESSION_MAX = [1, 3, 3, 5, 5, 9, 9, 59, 60, 61, 113, 113, 113, 119, 119, 119, 119, 119, 120, 120, 120, 121, 121, 121, 121]
SESSION_MEAN = [1.000, 1.750, 1.750, 2.125, 2.219, 2.500, 2.570, 3.273, 3.383, 3.529, 3.610, 3.666, 3.668, 3.696, 3.705,
                3.717, 3.721, 3.728, 3.730, 3.734, 3.735, 3.737, 3.738, 3.739, 3.739]


@nb.njit(cache=False)
def curl(s, L):
    best = 1
    for p in range(1, L // 2 + 1):
        k = 1
        while (k + 1) * p <= L:
            same = True
            a = L - (k + 1) * p
            b = L - p
            for i in range(p):
                if s[a + i] != s[b + i]:
                    same = False
                    break
            if not same:
                break
            k += 1
        if k > best:
            best = k
    return best


@nb.njit(parallel=True, cache=False)
def sweep(n, cap):
    total = 1 << n
    nch = 96
    best = np.zeros(nch, np.int64)
    ssum = np.zeros(nch, np.int64)
    fails = np.zeros(nch, np.int64)
    for c in nb.prange(nch):
        s = np.zeros(cap + n + 8, np.int64)
        lo = total * c // nch
        hi = total * (c + 1) // nch
        for st in range(lo, hi):
            for i in range(n):
                s[i] = 3 if (st >> i) & 1 else 2
            L = n
            steps = 0
            while True:
                k = curl(s, L)
                steps += 1
                if k == 1:
                    break
                if steps >= cap:
                    fails[c] += 1
                    break
                s[L] = k
                L += 1
            ssum[c] += steps
            if steps > best[c]:
                best[c] = steps
    return best.max(), ssum.sum(), fails.sum()


t0 = time.time()
g = np.zeros(1100, np.int64)
g[0] = 1
for L in range(1, 1000):
    g[L] = curl(g, L)
first = {}
for i in range(1000):
    first.setdefault(int(g[i]), i + 1)
print(f'Gijswijt first occurrences (1000 terms, direct curling): {first}', flush=True)
allok = True
for n in range(1, NMAX + 1):
    t1 = time.time()
    mx, tot, fails = sweep(n, CAP)
    mean = tot / (1 << n)
    ok = (mx == SESSION_MAX[n - 1]) and abs(mean - SESSION_MEAN[n - 1]) < 6e-4 and (n + mx - 1 == A094004[n - 1]) and fails == 0
    allok &= ok
    print(f'n = {n:2d}: max c = {mx:3d} (session {SESSION_MAX[n - 1]:3d}), mean {mean:.4f} (session {SESSION_MEAN[n - 1]:.3f}), '
          f'n + c - 1 = {n + mx - 1} (A094004 {A094004[n - 1]}), cap failures {fails}  {"ok" if ok else "MISMATCH"}  '
          f'[{time.time() - t1:.1f} s]', flush=True)
print(f'all rows agree: {allok}; total {time.time() - t0:.1f} s')
