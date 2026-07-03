#!/usr/bin/env python3
"""
Regime-C, done right (mac-mini-2026-07-03-S21): COVERING families with 7 near-equal far runners.
A family is COVERING (the hard LRC case) iff every q in {2..14} divides SOME speed -> t=1/q fails
for all q<=14. NON-covering families are trivially lonely at t=1/q. So regime C = COVERING +
7-near-equal-far.

KEY constraint: 7 consecutive far {w1..w1+6} cover q=2..7 (any 7 consecutive contain a multiple).
For q=8..14, the 6 near runners (<=22) must cover the ones the far miss. 6 near cover <=6 of the 7
moduli {8..14}; the 7th must be covered by the far -> CONSTRAINS w1 mod q. So regime-C covering
families are a THIN, structured set. Explore: do they exist? what are their lonely times / q?
"""
from fractions import Fraction as F
from math import gcd

R = F(1, 14)

def dmin(x):
    x = x % 1
    return min(x, 1 - x)

def is_covering(speeds):
    """covering iff every q in 2..14 divides some speed."""
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def lonely_at(speeds, a, q):
    return all(dmin(F(v * a, q)) >= R for v in speeds)

def smallest_lonely_q(speeds, qmax=4000):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and lonely_at(speeds, a, q):
                return q, a
    return None, None

def M_scan(speeds, N=60000):
    best, bt = F(0), F(0)
    for k in range(1, N):
        t = F(k, N); m = min(dmin(v * t) for v in speeds)
        if m > best: best, bt = m, t
    return best, bt

if __name__ == "__main__":
    print("=" * 82)
    print("COVERING regime-C families: near={8,9,10,11,12,13} (cover 8-13) + far={w1..w1+6} (14|w1 covers 14)")
    print("=" * 82)
    near = [8, 9, 10, 11, 12, 13]  # covers q=8,9,10,11,12,13 (and 2,3,4,6)
    # need q=14 covered: far must contain a multiple of 14 -> w1..w1+6 contains 14k
    print(f"{'w1':>6} {'covering?':>9} {'M(scan)':>10} {'>=1/14?':>8} {'smallest lonely q':>18} {'a':>5} {'q vs 183?'}")
    found = []
    for w1 in range(23, 600):
        far = list(range(w1, w1 + 7))
        speeds = near + far
        if not is_covering(speeds):
            continue
        q, a = smallest_lonely_q(speeds, qmax=3000)
        found.append((w1, q, a))
        if len(found) <= 20 or (q and q > 100):
            M, _ = M_scan(speeds, N=20000)
            qstr = f"{q}" if q else ">3000"
            hint = ""
            if q:
                hint = f"q={q} (183=3*61; q%61={q%61}, q%183={q%183}, q%7={q%7})"
            print(f"{w1:>6} {str(is_covering(speeds)):>9} {float(M):>10.5f} {str(M>=R):>8} {qstr:>18} {a if a else '-':>5} {hint}")
    print(f"\n# covering regime-C families found (w1 in 23..599): {len(found)}")
    if found:
        qs = [q for _, q, _ in found if q]
        print(f"# smallest-lonely-q distribution: min={min(qs)}, max={max(qs)}, all<=183? {all(q<=183 for q in qs)}")
        big = [(w,q) for w,q,a in found if q and q > 14]
        print(f"# families needing q>14 (genuinely past the sieve): {len(big)} -> {big[:15]}")
