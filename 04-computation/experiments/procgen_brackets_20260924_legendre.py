#!/usr/bin/env python3
"""Part 4: Legendre-type statements for odd-square brackets, and the (absent) link to Collatz orbits.

B_m = ((2m-1)^2, (2m+1)^2] splits at (2m)^2 into two Legendre intervals and at (2m-1)2m, (2m)^2, 2m(2m+1)
into four Oppermann intervals.  Legendre (a prime between consecutive squares) gives >= 2 primes per bracket,
Oppermann >= 4; both OPEN.  Brocard (>= 4 primes between squares of consecutive primes) is a statement about
unions of brackets, since for odd primes p < p' the interval (p^2, p'^2] is B_((p+1)/2) u ... u B_((p'-1)/2).

(a) counts per bracket / half / quarter for m <= MMAX (vectorised segmented sieve)      FINITE-EXACT
(b) Brocard unions in range                                                             FINITE-EXACT
(c) interval-length bookkeeping: bracket 8m ~ 4 sqrt(x) against x^0.525 (BHP) and sqrt(x) log x (RH)
(d) negative control: prime counts of brackets against Collatz stopping-time statistics  FINITE (no signal)
Usage: python3 <this> [MMAX]  (default 30000: odd squares to 3.6e9, ~20 s, < 200 MB)
"""
import math, sys, time
import numpy as np

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
LIM = (2 * MMAX + 1) ** 2
print("=" * 100)
print("PART 4. Legendre-type statements for odd-square brackets")
print("=" * 100)
t0 = time.time()
r = math.isqrt(LIM) + 1
base = np.ones(r + 1, dtype=bool); base[:2] = False
for i in range(2, math.isqrt(r) + 1):
    if base[i]: base[i * i::i] = False
small = np.flatnonzero(base)
cnt = np.zeros((MMAX + 2, 4), dtype=np.int64)       # quarters
seg = 1 << 24; start = 2
while start <= LIM:
    end = min(start + seg, LIM + 1)
    a = np.ones(end - start, dtype=bool)
    for p in small:
        p = int(p)
        if p * p >= end: break
        f = max(p * p, ((start + p - 1) // p) * p)
        a[f - start::p] = False
    q = np.flatnonzero(a).astype(np.int64) + start
    s = np.floor(np.sqrt(q.astype(np.float64))).astype(np.int64)
    s -= (s * s > q); s += ((s + 1) * (s + 1) <= q)                 # exact isqrt
    m = (s + 1) // 2
    m = np.where((2 * m + 1) ** 2 < q, m + 1, m)                    # (2m-1)^2 < q <= (2m+1)^2
    m = np.where((2 * m - 1) ** 2 >= q, m - 1, m)
    keep = (m >= 1) & (m <= MMAX); q = q[keep]; m = m[keep]
    quarter = (q > (2 * m - 1) * (2 * m)).astype(np.int64) + (q > (2 * m) ** 2) + (q > (2 * m) * (2 * m + 1))
    np.add.at(cnt, (m, quarter), 1)
    start = end
C = cnt[1:MMAX + 1]; tot = C.sum(1); half1 = C[:, 0] + C[:, 1]; half2 = C[:, 2] + C[:, 3]
M = np.arange(1, MMAX + 1)
print(f"\n(a) brackets m = 1..{MMAX} (odd squares up to {LIM:.3e}; {time.time()-t0:.1f}s)")
print(f"    min primes per bracket (m >= 2): {tot[1:].min()} at m = {int(M[1:][tot[1:].argmin()])};  per Legendre half: {min(half1.min(), half2.min())};"
      f"  per Oppermann quarter (m >= 2): {C[1:].min()}")
print(f"    first 12 brackets, quarter counts: {C[:12].tolist()}")
low = [(int(mm), int(tot[mm - 1])) for mm in M if tot[mm - 1] <= 6]
print(f"    brackets with <= 6 primes: {low[:20]}")
ratio = tot / (4 * M / np.log(2 * M + 1))
print(f"    count / (4m/log(2m+1)): mean {ratio[100:].mean():.4f}, min {ratio[100:].min():.4f} (m >= 101)")
qmin = C.min(axis=0)
print(f"    minimum over m of each quarter separately: {qmin.tolist()}; quarter minima for m >= 100: {C[99:].min(axis=0).tolist()}")

# ------------------------------------------------------------------------------------------ (b)
print("\n(b) Brocard's conjecture as a statement about unions of brackets")
odd_pr = [int(p) for p in small if 3 <= p <= 2 * MMAX + 1]
worst = None
for p, pp in zip(odd_pr, odd_pr[1:]):
    lo_m, hi_m = (p + 1) // 2, (pp - 1) // 2
    c = int(tot[lo_m - 1:hi_m].sum())
    if worst is None or c < worst[0]:
        worst = (c, p, pp, hi_m - lo_m + 1)
print(f"    min primes in (p^2, p'^2] over consecutive odd primes p' <= {2*MMAX+1}: {worst[0]} (p = {worst[1]}, p' = {worst[2]}, {worst[3]} bracket(s))")
print("    twin primes p' = p + 2 give a single bracket, so '>= 4 primes per bracket' (Oppermann) => Brocard, while")
print("    Legendre alone gives only >= 2 there.")

# ------------------------------------------------------------------------------------------ (c)
print("\n(c) Interval lengths at x ~ (2m)^2: bracket 8m = 4 sqrt(x); Legendre interval 2 sqrt(x) + 1")
print("    Baker-Harman-Pintz (2001): a prime in [x - x^0.525, x] for all large x (x_0 not explicit).  x^0.525 <= 4 sqrt(x)")
print(f"    iff x <= 4^40 = {4.0**40:.3e}; beyond that BHP intervals are longer than brackets, so BHP never implies")
print("    'a prime in every bracket'.  Under RH (Cramer) gaps are O(sqrt(x) log x): still a log too long.")
print("    Known: Legendre verified to n^2 ~ 2e19 (maximal prime gaps; Oliveira e Silva-Herzog-Pardi 2014, per the")
print("    Wikipedia summary, not read), hence >= 2 primes in every bracket below 2e19; exceptional-set bounds for")
print("    Legendre: Bazzanella (Arch. Math. 75 (2000) 29-34), and O(N^eps) exceptions under Lindeloef (Bazzanella 2013).")
print("    Cube analogue PROVED: Ingham (1937) primes between large consecutive cubes; Dudek (2016) explicit.")

# ------------------------------------------------------------------------------------------ (d)
print("\n(d) Negative control: Collatz statistics per bracket against prime counts (m <= 1000, n <= 4.0e6)")
MM = 1000; NN = (2 * MM + 1) ** 2
st = np.zeros(NN + 1, dtype=np.int32)                 # total stopping time of the shortcut map T
for n in range(2, NN + 1):
    v, k = n, 0
    while v >= n:
        v = v // 2 if v % 2 == 0 else (3 * v + 1) // 2
        k += 1
    st[n] = k + st[v]
mean_st = np.array([st[(2 * m - 1) ** 2 + 1:(2 * m + 1) ** 2 + 1].mean() for m in range(1, MM + 1)])
pc = tot[:MM].astype(float)
# detrend both against log m (quadratic fit) and correlate the residuals
x = np.log(np.arange(1, MM + 1))
def resid(y):
    c = np.polyfit(x[10:], y[10:], 2); return y[10:] - np.polyval(c, x[10:])
rr = np.corrcoef(resid(mean_st), resid(pc / (4 * np.arange(1, MM + 1) / np.log(2 * np.arange(1, MM + 1) + 1))))[0, 1]
print(f"    corr(detrended mean total stopping time of B_m, detrended normalised prime count of B_m), 11 <= m <= {MM}: {rr:+.4f}"
      f"  (|r| < {2/math.sqrt(MM-10):.3f} is noise)")
print("    No provable link is known: bracket statements see only archimedean size and primality, T sees 2-adic digits;")
print("    the pairing family of part 2 keeps every bracket and every prime fixed while changing tree-ness.")
