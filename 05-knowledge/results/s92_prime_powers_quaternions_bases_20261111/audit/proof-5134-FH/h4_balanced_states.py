#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (H)(iv): balanced states X_b(h) = y(h), i.e. 1 + b psi(h) = 2h + 1, i.e. b psi(h) = 2h.

Own code.  q(h) is computed in two independent ways:
  (R) the recursion q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2)   (table, h < 2^24);
  (S) the closed sum q(h) = sum over one-bits i of h of ceil(floor(h / 2^(i+1)) / 2)   (blocks, h < 2^30).
psi(h) = h - q(h) is also checked against the recursion psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1.
"""
import sys
import numpy as np
from numba import njit, prange


def q_table(nbits):
    q = np.zeros(1 << nbits, dtype=np.int64)
    for k in range(nbits):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return q


@njit(cache=True)
def q_sum(h):
    tot = 0
    i = 0
    while (h >> i) != 0:
        if (h >> i) & 1:
            H = h >> (i + 1)
            tot += (H + 1) >> 1
        i += 1
    return tot


@njit(parallel=True, cache=True)
def scan_block(lo, hi, b, out, maxout):
    """balanced states of base b in [lo, hi): b (h - q) = 2 h.  returns count; first maxout stored (unordered per thread, sorted later)."""
    nthreads = 64
    counts = np.zeros(nthreads, dtype=np.int64)
    step = (hi - lo + nthreads - 1) // nthreads
    buf = np.zeros((nthreads, maxout), dtype=np.int64)
    for t in prange(nthreads):
        a = lo + t * step
        e = min(hi, a + step)
        c = 0
        for h in range(a, e):
            q = q_sum(h)
            if b * (h - q) == 2 * h:
                if c < maxout:
                    buf[t, c] = h
                c += 1
        counts[t] = c
    n = 0
    for t in range(nthreads):
        for j in range(min(counts[t], maxout)):
            if n < maxout:
                out[n] = buf[t, j]
                n += 1
    return counts.sum()


@njit(parallel=True, cache=True)
def min_ratio_block(lo, hi):
    """min over lo <= h < hi of (2 psi(h) - h)  (must be >= 1 for h >= 1), and the h attaining the least psi/h."""
    nthreads = 64
    mins = np.full(nthreads, 1 << 60, dtype=np.int64)
    step = (hi - lo + nthreads - 1) // nthreads
    for t in prange(nthreads):
        a = lo + t * step
        e = min(hi, a + step)
        m = 1 << 60
        for h in range(a, e):
            v = 2 * (h - q_sum(h)) - h
            if v < m:
                m = v
        mins[t] = m
    return mins.min()


def main():
    P = print
    P("THM-5134 (H)(iv): balanced states (auditor proof-5134-FH)")
    NB = 24
    q = q_table(NB)
    h = np.arange(1 << NB, dtype=np.int64)
    psi = h - q
    # recursion of psi
    hh = h[: 1 << (NB - 1)]
    ok_psi = bool(np.all(psi[2 * hh] == psi[hh] + hh)) and bool(np.all(psi[2 * hh + 1] == psi[hh] + hh // 2 + 1))
    P(f"psi = h - q satisfies psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1 for h < 2^{NB - 1}: {ok_psi}")
    # closed sum agrees with the recursion on a sample
    samp = np.random.default_rng(92).integers(0, 1 << NB, 200000)
    ok_sum = all(q_sum(int(v)) == int(q[int(v)]) for v in samp[:20000])
    P(f"closed sum for q agrees with the recursion (20000 random h < 2^{NB}): {ok_sum}")
    P(f"q(h) < h/2, i.e. 2 psi(h) > h, for every 1 <= h < 2^{NB}: {bool(np.all(2 * psi[1:] > h[1:]))};  "
      f"min of 2 psi(h) - h over 1 <= h < 2^{NB}: {int((2 * psi[1:] - h[1:]).min())}")
    # where is q(h)/h largest
    r = q[1:] / h[1:]
    j = int(np.argmax(r)) + 1
    P(f"sup of q(h)/h over 1 <= h < 2^{NB}: {r.max():.9f} at h = {j} = {bin(j)}  (below 1/2; the bound h/2 is not attained)")
    P("")
    P(f"balanced states below 2^{NB} (table):")
    for b in range(2, 13):
        sol = h[b * psi == 2 * h]
        extra = ""
        if b == 2:
            extra = f"  all 0 or powers of two: {all(int(v) & (int(v) - 1) == 0 for v in sol.tolist())}"
        P(f"  b = {b:2d}: {len(sol):2d} states: {sol.tolist() if len(sol) <= 14 else sol[:14].tolist()}{extra}")
    sol3 = h[3 * psi == 2 * h]
    P("  b = 3 in detail: h, binary, psi, X = y:")
    for v in sol3.tolist():
        P(f"    h = {v:9d} = {bin(v)[2:]:>24s}  psi = {int(psi[v]):9d}  X = y = {2 * v + 1}  = {factor_str(2 * v + 1)}")
    P("")
    P("extension by the closed sum (independent of the table), blocks up to 2^30:")
    out = np.zeros(4096, dtype=np.int64)
    for b in (3, 4, 5):
        allsol = []
        tot = 0
        for k in range(0, 31):
            lo = 0 if k == 0 else 1 << (k - 1)
            hi = 1 if k == 0 else 1 << k
            c = scan_block(lo, hi, b, out, 4096)
            tot += int(c)
            allsol += sorted(out[: int(c)].tolist())
        P(f"  b = {b}: {tot} balanced states below 2^30: {allsol}")
        if b == 3:
            for kk in (24, 26, 28, 30):
                P(f"      below 2^{kk}: {sum(1 for v in allsol if v < (1 << kk))}")
    mn = min(min_ratio_block(1 << (k - 1), 1 << k) for k in range(1, 31))
    P(f"  min of 2 psi(h) - h over 1 <= h < 2^30: {int(mn)}   (psi(h) > h/2 for every 1 <= h < 2^30: {mn >= 1})")


def factor_str(n):
    f = []
    d = 2
    m = n
    while d * d <= m:
        while m % d == 0:
            f.append(d)
            m //= d
        d += 1
    if m > 1:
        f.append(m)
    return " * ".join(str(x) for x in f) if len(f) > 1 else (str(n) + " (prime)" if n > 1 else "1")


if __name__ == "__main__":
    main()
