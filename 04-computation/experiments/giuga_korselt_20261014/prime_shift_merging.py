#!/usr/bin/env python3
"""Collatz merging of n and n + h restricted to primes and to twin primes (THM-4627 (2) corollary).

Terras map T(n) = n/2 (even), (3n+1)/2 (odd).  (n, n+h) MERGE if T^k(n) = T^k(n+h) with equal numbers of odd steps for
some k (before either orbit reaches 1; pairs that reach 1 first count as non-merged, as in THM-4627's integer census).
Vectorized iteration (int64; trajectories below 10^7 stay far below 2^63).

Reports the merging fraction for: all odd n <= N, primes p <= N, lesser twin primes p <= N (p, p+2 prime), with h = 2,
and h = 1 for all n and primes.  Corollary (PROVED from THM-4627 (2)(a): the merge set by time K is periodic mod 2^K
with density -> 1; primes equidistribute over odd classes mod 2^K): relative density 1 among primes for every h.
Usage: python3 prime_shift_merging.py [N=10**7]
"""
import sys
import time
import numpy as np

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 7


def sieve(n):
    s = np.ones(n + 3, dtype=bool)
    s[:2] = False
    for p in range(2, int((n + 2) ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return s


def merge_fraction(starts, h, maxsteps=3000):
    a = starts.astype(np.int64).copy()
    b = a + h
    ca = np.zeros(len(a), dtype=np.int32)
    cb = np.zeros(len(a), dtype=np.int32)
    alive = np.ones(len(a), dtype=bool)
    merged = np.zeros(len(a), dtype=bool)
    for _ in range(maxsteps):
        idx = np.nonzero(alive)[0]
        if len(idx) == 0:
            break
        x, y = a[idx], b[idx]
        ox, oy = (x & 1) == 1, (y & 1) == 1
        x = np.where(ox, (3 * x + 1) >> 1, x >> 1)
        y = np.where(oy, (3 * y + 1) >> 1, y >> 1)
        ca[idx] += ox
        cb[idx] += oy
        a[idx], b[idx] = x, y
        m = (x == y) & (ca[idx] == cb[idx])
        merged[idx[m]] = True
        done = m | (x == 1) | (y == 1) | (x == 2) | (y == 2)
        alive[idx[done]] = False
    return merged.mean(), int(alive.sum())


def main():
    t0 = time.time()
    s = sieve(N + 2)
    odd = np.arange(3, N + 1, 2, dtype=np.int64)
    primes = np.nonzero(s[:N + 1])[0].astype(np.int64)
    primes = primes[primes > 2]
    twins = primes[s[primes + 2]]
    print(f"N = {N}: odd n {len(odd)}, odd primes {len(primes)}, lesser twin primes {len(twins)}  [{time.time() - t0:.1f}s]", flush=True)
    for name, arr, h in [('all odd n', odd, 2), ('primes p', primes, 2), ('twin primes (p, p+2)', twins, 2),
                         ('all odd n', odd, 1), ('primes p', primes, 1)]:
        fr, unfinished = merge_fraction(arr, h)
        print(f"   h = {h}: {name:22s} merging fraction {fr:.5f}  (unfinished {unfinished})  [{time.time() - t0:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
