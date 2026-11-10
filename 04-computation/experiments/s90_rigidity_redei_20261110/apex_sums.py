#!/usr/bin/env python3
"""Size of the exponential sums S(theta) = sum over configurations of exp(2 pi i theta H) for one and two apexes
over a transitive chain (README section 5, program item; HEURISTIC data, no theorem).

Family F(k, n): a transitive chain of n vertices plus k apexes with arbitrary cuts (for k = 2 the second apex also
chooses its arc to the first), all labelled configurations counted: |F(1, n)| = 2^n, |F(2, n)| = 2^n * 2^(n+1).
Every count is at most k! (k+1)^n (a Hamiltonian path orders the apexes and distributes the chain vertices among
the k+1 gaps, increasing inside each gap).

Method for k = 2: for each first-apex pattern the core C = chain + apex has the exact ear quadratic
E(x) = sum_{i != j} (1 - x_i) x_j P_ij + sum_j x_j s_j + sum_i (1 - x_i) e_i (checked in ear_weyl_identity.py);
the script evaluates it on all cuts and also re-counts a sample of configurations by direct subset DP.
Output: for random theta, |S| / |F| and the exponent log|S| / log|F| (1/2 = square-root cancellation).
Run:  python3 apex_sums.py > apex_sums.out
"""
import math
import sys

import numpy as np
from numba import njit


@njit(cache=True)
def core_tables(m, beats, end, start):
    size = 1 << m
    for S in range(size):
        for v in range(m):
            end[S, v] = 0
            start[S, v] = 0
    for v in range(m):
        end[1 << v, v] = 1
        start[1 << v, v] = 1
    for mask in range(1, size):
        for v in range(m):
            c = end[mask, v]
            if c:
                for w in range(m):
                    if not (mask >> w) & 1 and (beats[v] >> w) & 1:
                        end[mask | (1 << w), w] += c
            c2 = start[mask, v]
            if c2:
                # start[S][v]: paths of C[S] starting at v; extend at the front by a vertex w with w -> v
                for w in range(m):
                    if not (mask >> w) & 1 and (beats[w] >> v) & 1:
                        start[mask | (1 << w), w] += c2


@njit(cache=True)
def direct_count(m, beats):
    size = 1 << m
    cnt = np.zeros((size, m), dtype=np.int64)
    for v in range(m):
        cnt[1 << v, v] = 1
    for mask in range(1, size):
        for v in range(m):
            c = cnt[mask, v]
            if c:
                for w in range(m):
                    if not (mask >> w) & 1 and (beats[v] >> w) & 1:
                        cnt[mask | (1 << w), w] += c
    tot = 0
    for v in range(m):
        tot += cnt[size - 1, v]
    return tot


@njit(cache=True)
def two_apex_sums(n, thetas, check_every):
    m = n + 1                       # core = chain (vertices 0..n-1) + first apex (vertex n)
    size = 1 << m
    full = size - 1
    end = np.zeros((size, m), dtype=np.int64)
    start = np.zeros((size, m), dtype=np.int64)
    beats = np.zeros(m, dtype=np.int64)
    P = np.zeros((m, m), dtype=np.int64)
    s = np.zeros(m, dtype=np.int64)
    e = np.zeros(m, dtype=np.int64)
    nt = len(thetas)
    re = np.zeros(nt)
    im = np.zeros(nt)
    maxH = 0
    mism = 0
    checked = 0
    b2 = np.zeros(m + 1, dtype=np.int64)
    for u in range(1 << n):
        for i in range(m):
            beats[i] = 0
        for i in range(n):
            for j in range(i + 1, n):
                beats[i] |= 1 << j
            if (u >> i) & 1:
                beats[n] |= 1 << i
            else:
                beats[i] |= 1 << n
        core_tables(m, beats, end, start)
        for i in range(m):
            s[i] = start[full, i]
            e[i] = end[full, i]
            for j in range(m):
                P[i, j] = 0
        for S in range(1, full):
            T = full ^ S
            for i in range(m):
                a = end[S, i]
                if a:
                    for j in range(m):
                        b = start[T, j]
                        if b:
                            P[i, j] += a * b
        for x in range(size):
            val = 0
            for i in range(m):
                if (x >> i) & 1:
                    val += s[i]
                else:
                    val += e[i]
                    for j in range(m):
                        if (x >> j) & 1:
                            val += P[i, j]
            if val > maxH:
                maxH = val
            for t in range(nt):
                ph = 2.0 * math.pi * ((thetas[t] * val) % 1.0)
                re[t] += math.cos(ph)
                im[t] += math.sin(ph)
            if check_every > 0 and ((u * size + x) % check_every) == 0:
                for i in range(m):
                    b2[i] = beats[i]
                b2[m] = 0
                for i in range(m):
                    if (x >> i) & 1:
                        b2[m] |= 1 << i
                    else:
                        b2[i] |= 1 << m
                checked += 1
                if direct_count(m + 1, b2) != val:
                    mism += 1
    return re, im, maxH, mism, checked


def one_apex_counts(n):
    N = 1 << n
    idx = np.arange(N, dtype=np.int64)
    H = np.zeros(N, dtype=np.int64)
    prefix = np.zeros(N, dtype=np.int64)
    prev = np.zeros(N, dtype=np.int64)
    for t in range(n):
        bit = (idx >> (n - 1 - t)) & 1
        if t == 0:
            H += bit
        else:
            H += (prev == 0) & (bit == 1)
        H += np.where(bit == 0, prefix, 0)
        prefix = 2 * prefix + bit
        prev = bit
    H += (prev == 0)
    return H


def main():
    rng = np.random.default_rng(20261110)
    thetas = rng.random(16)
    print("one apex, |F| = 2^n, 16 random theta: median and max of |S|/|F|, exponent log|S|/log|F| of the median")
    for n in (12, 14, 16, 18, 20, 22):
        H = one_apex_counts(n).astype(np.float64)
        vals = np.array([abs(np.exp(2j * np.pi * ((th * H) % 1.0)).sum()) for th in thetas])
        F = 2.0 ** n
        print(f"  n = {n:2d}: median {np.median(vals)/F:.5f}  max {vals.max()/F:.5f}  exponent {math.log(np.median(vals))/math.log(F):.3f}"
              f"   (square root would be {F ** -0.5:.5f})")
    print("two apexes, |F| = 2^n * 2^(n+1): same statistics; max count against the bound 2 * 3^n")
    ok = True
    for n in (6, 8, 10, 11):
        re, im, maxH, mism, checked = two_apex_sums(n, thetas, 9973)
        vals = np.sqrt(re * re + im * im)
        F = 2.0 ** (2 * n + 1)
        ok &= (mism == 0) and (maxH <= 2 * 3 ** n)
        print(f"  n = {n:2d}: median {np.median(vals)/F:.6f}  max {vals.max()/F:.6f}  exponent {math.log(np.median(vals))/math.log(F):.3f}"
              f"   (square root would be {F ** -0.5:.6f});  max count {maxH} <= {2 * 3 ** n}; direct recounts {checked}, mismatches {mism}")
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
