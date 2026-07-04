#!/usr/bin/env python3
"""
The confinement = a covering-min piece (mac-mini-2026-07-03-S34).
Elementary facts pinning where q*>14 tight families live:
 (1) primitive tight (M=1/14) => covers {2..13} [THM-523 q-witness: else hides at 1/q>1/14].
 (2) if it also MISSES 14 => tight at t=1/14 => q*=14 (elementary).
 => contrapositive: a primitive tight family with q*>14 must COVER 14 => primitive COVERING with M=1/14,
    forbidden by the covering-min (M>=14/183, HYP-4060). So confinement <=> "no primitive tight covering
    family" = a covering-min piece, NOT an independent gap. Even block covers 14 but is IMPRIMITIVE (loophole).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_and_pts(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0); pts = []
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            t = F(a, q); m = min(nd(v*t) for v in sp)
            if m > best: best = m; pts = [t]
            elif m == best: pts.append(t)
    return best, sorted(set(pts))
def covers(sp, q): return any(v % q == 0 for v in sp)

if __name__ == "__main__":
    print(f"{'family':>11} {'M':>8} {'prim':>5} {'cov{2..13}':>10} {'has 14|v':>9} {'tight@1/14':>11} {'q*':>6}")
    for name, S in [('AP', list(range(1,14))), ('GW', [1,2,3,4,5,6,7,8,9,10,11,13,24]),
                    ('even block', [2*i for i in range(1,14)])]:
        M, pts = M_and_pts(S)
        cov213 = all(covers(S, q) for q in range(2,14))
        has14 = covers(S, 14)
        tat = (min(nd(v*F(1,14)) for v in S) == F(1,14))
        qs = sorted(set(t.denominator for t in pts))
        print(f"{name:>11} {float(M):>8.5f} {str(reduce(gcd,S)==1):>5} {str(cov213):>10} {str(has14):>9} {str(tat):>11} {str(qs):>6}")
    print("\n=> primitive tight => covers {2..13}; miss-14 => q*=14 (elementary); q*>14 => covers 14 =>")
    print("   primitive covering with M=1/14 => forbidden by covering-min. Confinement <=> covering-min piece.")
    print("   GAP-A (non-covering tight = {AP,GW}? three-gap) + GAP-B (covering-min) are the two real gaps.")
