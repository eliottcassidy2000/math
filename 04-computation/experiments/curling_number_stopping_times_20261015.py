#!/usr/bin/env python3
"""Curling numbers, Gijswijt's sequence and the Curling Number Conjecture as a 2-3 termination problem.

The curling number of a finite sequence S is the largest k such that S = X Y^k with Y nonempty.  Gijswijt's sequence
(OEIS A090822) is b(1) = 1, b(n+1) = curl(b(1..n)).  The Curling Number Conjecture (Chaffin, Linderman, Sloane, Wilks,
J. Integer Seq. 16 (2013)) states that every finite starting sequence, extended by repeatedly appending its curling
number, eventually produces a 1; it suffices to consider starts over {2, 3}.  For a start S over {2,3} define its curling
stopping time c(S) = number of appended terms up to and including the first appended 1.

Incremental curling: for every period p keep l_p = length of the longest suffix with period p; appending x keeps
l_p + 1 if x equals the symbol p places back, else resets to min(p, length).  curl = max_p floor(l_p / p).
Session opus-2026-10-15-S37.  Usage: python curling_number_stopping_times_20261015.py NMAX GIJ
"""
import sys, time
import numpy as np
import numba as nb

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
GIJ = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
CAP = 4000


@nb.njit(cache=True)
def append(buf, lp, L, x):
    # append x at position L; lp[p] = longest suffix with period p (O(1) update per p)
    buf[L] = x
    Ln = L + 1
    for p in range(1, Ln + 1):
        if p <= L and buf[L - p] == x:
            lp[p] = lp[p] + 1
        else:
            lp[p] = p if p < Ln else Ln
    return Ln


@nb.njit(cache=True)
def stopping_time(start_bits, n, cap, buf, lp):
    # start: symbols 2/3 from bits; returns appended count up to and including the first 1 (or -1 at cap)
    L = 0
    for i in range(lp.shape[0]):
        lp[i] = 0
    for i in range(n):
        x = 3 if (start_bits >> i) & 1 else 2
        L = append(buf, lp, L, x)
    steps = 0
    while steps < cap:
        c = 1
        for p in range(1, L // 2 + 1):
            k = lp[p] // p
            if k > c:
                c = k
        steps += 1
        if c == 1:
            return steps
        L = append(buf, lp, L, c)
    return -1


@nb.njit(parallel=True, cache=True)
def sweep(n, cap):
    total = 1 << n
    hist = np.zeros((64, cap + 2), dtype=np.int64)
    best = np.zeros(64, dtype=np.int64)
    arg = np.zeros(64, dtype=np.int64)
    nchunk = 64
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        buf = np.zeros(cap + n + 8, dtype=np.int64)
        lp = np.zeros(cap + n + 8, dtype=np.int64)
        for s in range(c * chunk, min(total, (c + 1) * chunk)):
            st = stopping_time(s, n, cap, buf, lp)
            if st < 0:
                hist[c, cap + 1] += 1
                continue
            hist[c, st] += 1
            if st > best[c]:
                best[c] = st
                arg[c] = s
    return hist.sum(axis=0), best, arg


def curl(seq):
    L = len(seq)
    best = 1
    for p in range(1, L // 2 + 1):
        k = 1
        while (k + 1) * p <= L and seq[L - (k + 1) * p:L - k * p] == seq[L - p:L]:
            k += 1
        best = max(best, k)
    return best


def gijswijt(m):
    b = [1]
    while len(b) < m:
        b.append(curl(b))
    return b


if __name__ == '__main__':
    t0 = time.time()
    g = gijswijt(GIJ)
    first = {}
    for i, v in enumerate(g, start=1):
        first.setdefault(v, i)
    print(f'Gijswijt sequence, first {GIJ} terms: first occurrences {dict(sorted(first.items()))}; '
          f'initial terms {g[:30]}  [{time.time() - t0:.1f} s]')
    # brute-force check of the incremental curling routine on random starts
    rng = np.random.default_rng(5)
    buf = np.zeros(CAP + 64, dtype=np.int64)
    lp = np.zeros(CAP + 64, dtype=np.int64)
    bad = 0
    for _ in range(300):
        n = int(rng.integers(1, 16))
        s = int(rng.integers(0, 1 << n))
        st = stopping_time(s, n, CAP, buf, lp)
        seq = [3 if (s >> i) & 1 else 2 for i in range(n)]
        steps = 0
        while True:
            c = curl(seq)
            steps += 1
            if c == 1 or steps > CAP:
                break
            seq.append(c)
        bad += (st != steps)
    print(f'incremental vs brute-force curling stopping times on 300 random starts: {bad} mismatches')
    for n in range(1, NMAX + 1):
        t1 = time.time()
        hist, best, arg = sweep(n, CAP)
        nz = np.nonzero(hist[:CAP + 1])[0]
        mx = int(nz.max())
        c = int(np.argmax(best))
        s = int(arg[c])
        start = ''.join('3' if (s >> i) & 1 else '2' for i in range(n))
        mean = float((np.arange(CAP + 1) * hist[:CAP + 1]).sum() / hist[:CAP + 1].sum())
        print(f'n = {n:2d}: all {1 << n} starts over {{2,3}} reach 1 (cap failures {int(hist[CAP + 1])}); '
              f'max curling stopping time {mx} (e.g. start {start}); mean {mean:.3f}  [{time.time() - t1:.1f} s]', flush=True)
