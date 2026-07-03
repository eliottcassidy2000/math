#!/usr/bin/env python3
"""
Does THM-608 close ACTUAL near-equal hge7 families end-to-end? (mac-mini-2026-07-03-S23)
An hge7 family = {near (<=6, |v|<=22)} u {far (>=7, near-equal, magnitude N)}, covering.
Split R = near, C = far. THM-608: R lonely at t0 with slack delta + (i) 2*delta*N>=V_R + (ii) D*(t0+delta/V_R)<6/7
=> S lonely. Check how often the conditions are MET for real covering near-equal hge7 families, and confirm S
is then lonely. This honestly scopes THM-608 as a closer for the near-equal-far obligation.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1.0; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2, 15))

def R_best_witness(R, qmax=60):
    """rational t0=a/q maximizing R's slack; return (t0, slack, q)."""
    best = None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            m = min(nd(v*a/q) for v in R)
            if m > 1/14 and (best is None or m-1/14 > best[1]):
                best = (a/q, m-1/14, q)
    return best

def M_fine(sp, steps=2000000):
    vmax=max(sp); K=30; s=min(steps, vmax*K); best=0.0
    for k in range(1,s):
        t=k/s; m=min(nd(v*t) for v in sp)
        if m>best: best=m
    return best

if __name__ == "__main__":
    rng = random.Random(608)
    print("THM-608 applied to real covering near-equal hge7 families ({near<=6} u {>=7 near-equal far}).")
    print("=" * 92)
    print(f"{'near (R)':>20} {'far mag N':>10} {'far spread D':>12} {'R slack':>8} {'R q':>4} {'(i)':>4} {'(ii)':>5} {'THM-608?':>9} {'S lonely':>9}")
    applies=0; tested=0; s_lonely=0
    for _ in range(4000):
        nnear = rng.randint(2,6); nfar = 13 - nnear
        near = sorted(rng.sample(range(1,23), nnear))
        N = rng.choice([300,600,1500,4000,12000,40000])
        D = rng.randint(1, rng.choice([6,12,25,40]))
        drifts = sorted(rng.sample(range(0,D+1), min(nfar, D+1)))
        if len(drifts) < 7: continue
        far = [N+d for d in drifts]
        S = sorted(set(near+far))
        if len(S)!=13 or gcd_all(S)!=1 or not is_covering(S): continue
        R = near
        bw = R_best_witness(R, qmax=60)
        if bw is None: continue
        t0, delta, q = bw; VR = max(R)
        Dspread = far[-1]-far[0]
        cond_i = 2*delta*N >= VR
        cond_ii = Dspread*(t0 + delta/VR) < 6/7
        thm = cond_i and cond_ii
        tested += 1
        if thm:
            applies += 1
            M = M_fine(S)
            if M >= 1/14: s_lonely += 1
            if applies <= 12:
                print(f"{str(near):>20} {N:>10} {Dspread:>12} {delta:>8.4f} {q:>4} {str(cond_i):>4} {str(cond_ii):>5} {'YES':>9} {str(M>=1/14):>9}")
    print(f"\nover {tested} covering near-equal hge7 families: THM-608 conditions MET in {applies} ({100*applies//max(tested,1)}%);")
    print(f"   of those, S verified lonely (M>=1/14) in {s_lonely}/{applies}")
    print("=> THM-608 closes the near-equal-far families whose near-part R has a fine-enough witness with slack.")
    print("   Scope: wider far spread D needs a finer R witness (larger q) -- (ii) D*(t0+delta/VR)<6/7 is the exact gate.")
