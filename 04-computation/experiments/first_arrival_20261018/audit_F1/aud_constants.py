#!/usr/bin/env python3
"""audit F1: independent exact evaluation of E_0(N) = sum_x c_0(x)^2 and SIB(N) = sum_{x=1 mod 8, 2x+1<=4N} c_0(x) c_0(2x+1)
by interval counting (no enumeration of x), valid for astronomically large N; checks THM-4690(4).

Method (own, different from source_ledger_constants.py):
  * For a seed x, c_0(x) = sum of K_N over its B-chain x, Bx, ... <= 4N (B(x) = 4x+1). The predicates
    "B^j(x) < 2N" and "B^j(x) <= 4N" are monotone in x, so (L(x), ends(x)) and, for the sibling z = 2x+1, (L(z), ends(z)) are
    constant on the integer intervals between the thresholds floor((2N-1-(4^j-1)/3)/4^j), floor((4N-(4^j-1)/3)/4^j) and the
    corresponding thresholds for 2x+1.
  * On each interval c_0 is evaluated from the chain definition at a representative x (exact), and the seeds of each residue
    class are counted exactly (seeds = x mod 8 in {1,3,7}; x mod 3 = r; sibling heads x = 1 mod 8).
  * Brute force for small N cross-checks the interval method exactly.
Outputs: E_0 - 34N/7 and SIB - 5N/7 for N up to 10^40 and for many random N; the sibling-case exception count.
"""
from fractions import Fraction as Fr
import random, math, sys

def chain_c0(x, N):
    """c_0(x) from the definition: sum over chain elements m <= 4N of w(m) + thr(m)[m < 2N]."""
    c, m = 0, x
    while m <= 4 * N:
        r = m % 3
        c += (2 if r == 2 else -1) + ((1 if r != 1 else -2) if m < 2 * N else 0)
        m = 4 * m + 1
    return c

def count_res(lo, hi, mod, residues):
    """number of integers x in [lo, hi] with x % mod in residues."""
    if hi < lo: return 0
    tot = 0
    for r in residues:
        # count x = r (mod mod) in [lo, hi]
        first = lo + ((r - lo) % mod)
        if first <= hi: tot += (hi - first) // mod + 1
    return tot

def thresholds_x(N, scale=1, shift=0):
    """thresholds t (in x) where a predicate on z = scale*x + shift changes between t and t+1."""
    ts = set()
    j = 0
    while True:
        p = 4 ** j; off = (p - 1) // 3
        a = (2 * N - 1 - off)                     # B^j(z) <= 2N-1  <=>  z <= a // p
        b = (4 * N - off)                         # B^j(z) <= 4N    <=>  z <= b // p
        if a < 0 and b < 0: break
        for zt in (a // p if a >= 0 else -1, b // p if b >= 0 else -1):
            # z <= zt  <=>  scale*x + shift <= zt  <=>  x <= (zt - shift) // scale
            if zt >= 0: ts.add((zt - shift) // scale)
        j += 1
    return ts

SEED_RES24 = {r: [x for x in range(24) if x % 2 == 1 and x % 8 in (1, 3, 7) and x % 3 == r] for r in range(3)}
HEAD_RES24 = {r: [x for x in range(24) if x % 8 == 1 and x % 3 == r] for r in range(3)}

def E0_interval(N):
    M = 4 * N
    ts = sorted(t for t in thresholds_x(N) if 0 <= t <= M) + [M]
    E, lo = 0, 1
    for t in ts:
        if t < lo: continue
        hi = t
        for r in range(3):
            # representative: smallest seed in [lo, hi] with x = r mod 3 (c_0 constant on the interval for fixed r)
            n = count_res(lo, hi, 24, SEED_RES24[r])
            if n == 0: continue
            first = min(lo + ((q - lo) % 24) for q in SEED_RES24[r])
            c = chain_c0(first, N)
            E += n * c * c
        lo = hi + 1
    return E

def SIB_interval(N):
    M = 4 * N
    xmax = (M - 1) // 2                                  # 2x+1 <= 4N
    ts = thresholds_x(N) | thresholds_x(N, 2, 1)
    ts = sorted(t for t in ts if 0 <= t <= xmax) + [xmax]
    S, lo = 0, 1
    exA = exB = 0
    for t in ts:
        if t < lo: continue
        hi = t
        for r in range(3):
            n = count_res(lo, hi, 24, HEAD_RES24[r])
            if n == 0: continue
            first = min(lo + ((q - lo) % 24) for q in HEAD_RES24[r])
            S += n * chain_c0(first, N) * chain_c0(2 * first + 1, N)
        lo = hi + 1
    return S

def brute(N):
    M = 4 * N
    c = {}
    for m in range(1, M + 1, 2):
        x = m
        while x % 8 == 5: x = (x - 1) // 4
        r = m % 3
        c[x] = c.get(x, 0) + (2 if r == 2 else -1) + ((1 if r != 1 else -2) if m < 2 * N else 0)
    E = sum(v * v for v in c.values())
    S = sum(c.get(x, 0) * c.get(2 * x + 1, 0) for x in range(1, M + 1, 8) if 2 * x + 1 <= M)
    return E, S, max(abs(v) for v in c.values())

def main():
    out = []
    # 1. interval method == brute force
    bad = 0
    Ns = list(range(1, 400)) + [1000, 4096, 7910, 10007, 65536, 66463]
    for N in Ns:
        Eb, Sb, cmax = brute(N)
        Ei, Si = E0_interval(N), SIB_interval(N)
        if (Eb, Sb) != (Ei, Si) or cmax > 3: bad += 1; out.append('MISMATCH N=%d brute (%d,%d) interval (%d,%d) max|c0| %d' % (N, Eb, Sb, Ei, Si, cmax))
    out.append('interval method vs brute force: %d values of N (1..399 and 6 larger), mismatches %d' % (len(Ns), bad))
    # 2. huge N
    out.append('%-44s %-28s %-24s %s' % ('N', 'E_0 - 34N/7', 'SIB - 5N/7', '(E_0-34N/7)/log2 N'))
    for k in list(range(3, 41, 3)) + [60, 100]:
        for N in (10 ** k, 2 ** int(k * 3.3219), 10 ** k + 12345):
            dE = Fr(E0_interval(N)) - Fr(34 * N, 7); dS = Fr(SIB_interval(N)) - Fr(5 * N, 7)
            out.append('%-44d %+-28.4f %+-24.4f %+.3f' % (N, float(dE), float(dS), float(dE) / math.log2(N)))
    # 3. random N: worst |E_0 - 34N/7| / log2 N
    rnd = random.Random(1018)
    wE = wS = 0.0; argE = argS = None; mE = mS = 0.0
    for _ in range(3000):
        N = rnd.randrange(2, 10 ** rnd.randrange(2, 30))
        dE = abs(Fr(7 * E0_interval(N) - 34 * N, 7)); dS = abs(Fr(7 * SIB_interval(N) - 5 * N, 7))   # exact (no float rounding)
        if float(dE) / math.log2(N) > wE: wE, argE = float(dE) / math.log2(N), N
        if float(dS) / math.log2(N) > wS: wS, argS = float(dS) / math.log2(N), N
        mE, mS = max(mE, float(dE)), max(mS, float(dS))
    out.append('3000 random N (2 .. 10^29), exact: max |E_0 - 34N/7|/log2 N = %.3f (N=%d), max |SIB - 5N/7|/log2 N = %.3f (N=%d); max |E_0 - 34N/7| = %.2f, max |SIB - 5N/7| = %.2f' % (wE, argE, wS, argS, mE, mS))
    print('\n'.join(out))
    return 1 if bad else 0

if __name__ == '__main__':
    sys.exit(main())
