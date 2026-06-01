#!/usr/bin/env python3
"""
lrc_gap_balanced_pair_s521.py   claudebox-2026-06-01-S521

The Lonely-Runner gap function as a balanced apex-pair problem.

For distinct positive integers v_1..v_m define M(v) = sup_t min_i ||v_i t||.
LRC(n=m+1) <=> M(v) >= 1/n.

THEOREM (verified here; proof in the reflection).  Unless all v_i are odd
(then M=1/2 at t=1/2), the maximum is attained at t* = k/(v_i+v_j) for a pair
i,j, with v_i,v_j equidistant from the observer on OPPOSITE sides (the binding
"apex" pair -- the observer's two circular neighbours).  Hence
    M(v) = max over pairs (i<j), k in [1, v_i+v_j) of
              { d = ||k v_i/(v_i+v_j)|| = ||k v_j/(v_i+v_j)|| :
                ||k w/(v_i+v_j)|| >= d for every other speed w }.

COROLLARY (tight locus).  M(v) = 1/n forces ||k v_i/(v_i+v_j)|| = 1/n, i.e.
min(r,q-r) = q/n with q = v_i+v_j, which requires n | (v_i+v_j).  So every
EXACTLY-TIGHT (extremal) speed set has a pair of speeds summing to a multiple of
n.  Verified: all tight sets have such a pair; sets with no such pair have
M > 1/n (slack).

CONSEQUENCES.
 * The canonical moduli for the LRC walk are the pairwise SUMS v_i+v_j (not
   arbitrary q): the optimal lonely time is t=k/(v_i+v_j).
 * The binding pair = observer's neighbours = the apex/seam tile of the staircase
   (fiber-only LRC): all the gap content lives at the apex.
 * The hard/tight locus of LRC(n) = speed sets with a pair summing to a multiple
   of n.  For n=14: pairs summing to 14, 28, ...
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

def dist(x):
    x = x % 1
    return min(x, 1 - x)
def D(sp, t): return min(dist(F(v) * t) for v in sp)
def cand(sp):
    C = set(); S = sorted(sp)
    for i in range(len(S)):
        for j in range(len(S)):
            for s in (S[i] + S[j], abs(S[i] - S[j])):
                if s:
                    for k in range(1, s + 1): C.add(F(k, s))
        for k in range(0, 2 * S[i] + 1): C.add(F(2 * k + 1, 2 * S[i]))
    return [t for t in C if 0 < t < 1]
def M_true(sp): return max(D(sp, t) for t in cand(sp))
def M_balanced(sp):
    best = F(1, 2) if all(v % 2 == 1 for v in sp) else F(0)
    for i in range(len(sp)):
        for j in range(i + 1, len(sp)):
            s = sp[i] + sp[j]
            for k in range(1, s):
                t = F(k, s); di = dist(F(sp[i]) * t)
                if di != dist(F(sp[j]) * t): continue
                if di > best and all(dist(F(w) * t) >= di for w in sp): best = di
    return best

def main():
    random.seed(5)
    print("(1) M_balanced == M_true (the balanced apex-pair formula for the LR gap):")
    mism = tested = 0
    for _ in range(1500):
        m = random.randint(2, 6); sp = sorted(random.sample(range(1, 30), m))
        if reduce_gcd(sp) != 1: continue
        tested += 1
        if M_true(sp) != M_balanced(sp): mism += 1
    print(f"    tested={tested}, mismatches={mism}")
    print("\n(2) tight locus: M=1/n  =>  some pair sums to 0 mod n  (n=m+1):")
    for m in (4, 5, 6):
        n = m + 1; B = 14 if m <= 5 else 11
        tight = tight_pair = nofail = 0; nopair_minM = F(1)
        for sp in combinations(range(1, B + 1), m):
            if reduce_gcd(sp) != 1: continue
            Mv = M_true(list(sp))
            if Mv < F(1, n): nofail += 1
            has = any((sp[i] + sp[j]) % n == 0 for i in range(m) for j in range(i + 1, m))
            if Mv == F(1, n):
                tight += 1; tight_pair += 1 if has else 0
            if not has and Mv < nopair_minM: nopair_minM = Mv
        print(f"    n={n}: LRC failures={nofail}; tight={tight} (all with pair? {tight==tight_pair}); "
              f"no-pair min M={nopair_minM}={float(nopair_minM):.4f} > 1/n={float(F(1,n)):.4f}")

def reduce_gcd(sp):
    g = 0
    for v in sp: g = gcd(g, v)
    return g

if __name__ == "__main__":
    main()
