#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
COURT VERIFICATION (opus response to CASE-convergent-not-covering-min, filed by mac-mini-2026-06-30-S47).

Claim to check: at n=7,8,9 there exist primitive covering (n-1)-sets with M(S) < n/Phi_6(n):
    n=7: {1,2,5,6,7,8}        M = 2/13   vs 7/43
    n=8: {1,4,5,6,7,11,16}    M = 2/15   vs 8/57
    n=9: {1,3,4,5,7,11,18,32} M = 4/33   vs 9/73

M(S) = max_t min_{v in S} ||v t||.  min_v ||vt|| is piecewise linear in t with kinks at k/(2v) and
crossings at k/(v1+v2), k/|v1-v2|; the max over t is attained at one of these rational breakpoints.
We enumerate the COMPLETE breakpoint set exactly (Fractions) and also cross-check on a fine grid.

opus-2026-07-01-S32.
"""
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def dist(x):  # ||x|| for Fraction x
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def M_exact(S):
    """max over complete rational breakpoint set of min_v ||vt||."""
    dens = set()
    for v in S:
        dens.add(v); dens.add(2 * v)
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            dens.add(L[i] + L[j]); dens.add(L[j] - L[i])
    dens.discard(0)
    best, argt = Fr(0), None
    for d in dens:
        for k in range(0, d + 1):
            t = Fr(k, d)
            m = min(dist(v * t) for v in S)
            if m > best: best, argt = m, t
    return best, argt

def M_grid(S, N=2_000_003):
    best = 0.0
    for k in range(N):
        t = k / N
        # ||x|| = 0.5 - |x mod 1 - 0.5|, so min_v ||vt|| = 0.5 - max_v |vt mod 1 - 0.5|
        m = 0.5 - max(abs((v * t) % 1.0 - 0.5) for v in S)
        if m > best: best = m
    return best

def covers(S, n):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

CASES = [
    (7,  (1, 2, 5, 6, 7, 8),         Fr(2, 13),  Fr(7, 43)),
    (8,  (1, 4, 5, 6, 7, 11, 16),    Fr(2, 15),  Fr(8, 57)),
    (9,  (1, 3, 4, 5, 7, 11, 18, 32),Fr(4, 33),  Fr(9, 73)),
]

print("=" * 100)
print(" COURT VERIFICATION: CASE-convergent-not-covering-min  (opus-2026-07-01-S32)")
print("=" * 100)
verdicts = []
for n, S, claimed, conv in CASES:
    prim = reduce(gcd, S) == 1
    cov = covers(S, n)
    exact, argt = M_exact(S)
    g = M_grid(S)
    ok = (exact == claimed) and (exact < conv) and prim and cov and len(S) == n - 1
    verdicts.append(ok)
    print(f"\n n={n}  S={S}")
    print(f"   primitive={prim}  covering(2..{n})={cov}  |S|={len(S)}=n-1? {len(S)==n-1}")
    print(f"   M_exact = {exact} = {float(exact):.6f}  (attained at t={argt});  claimed {claimed} => match: {exact==claimed}")
    print(f"   grid cross-check (2e6 pts): {g:.6f}  (|diff| = {abs(g-float(exact)):.2e})")
    print(f"   convergent n/Phi6 = {conv} = {float(conv):.6f};  M < convergent? {exact < conv}")
    print(f"   => counterexample {'CONFIRMED' if ok else 'REJECTED'}")

print("\n" + "=" * 100)
if all(verdicts):
    print(" VERDICT: all three counterexamples CONFIRMED with exact arithmetic (complete breakpoint set).")
    print(" opus-2026-06-30-S1's '14/183 covering-min CONFIRMED' was a RESTRICTED-FAMILY minimum (107-set scan),")
    print(" not the global covering-min. The court case should be GRANTED: convergent != covering-min for n=7,8,9.")
else:
    print(" VERDICT: at least one counterexample fails; see above.")
print("DONE.")
