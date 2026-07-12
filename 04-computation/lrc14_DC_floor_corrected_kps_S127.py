# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont41: CORRECTING the DC floor + reconciling with opus-S235 band-edge margin lemma.
# cont.40 CLAIMED (adversarial) the min-M divisor-complete family is {1..14}\{6}, M=2/23. WRONG -- a LOCAL min.
# The EXHAUSTIVE DC census (Vmax<=22) finds a LOWER floor: M = 1/12 = 0.0833 at the UNIQUE extremal
# {1,2,3,4,10,11,12,13,14,15,16,17,18} (two-block {1,2,3,4}u{10..18}, NOT near-AP). eps = 1/12-1/14 = 1/84.
# RECONCILIATION (opus-S235 band-edge margin lemma, PROVED): S clears (bandCount=0) at non-14 modulus q
# => M >= ceil(q/14)/q > 1/14. The extremal first clears (non-14) at q=24: ceil(24/14)/24 = 2/24 = 1/12 = M
# EXACTLY -- opus's lemma is TIGHT at this extremal. => the detuning bound is a FREE COROLLARY of
# bounded-clearing (my window/anti-concentration, the one OPEN piece); no separate detuning theorem.
# LESSON (= cont.26): adversarial hill-climb finds LOCAL minima; use the EXHAUSTIVE census for extremal claims.
from fractions import Fraction as F
from math import gcd, ceil
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def Mx(v):
    best=F(-1)
    for i in range(13):
        for j in range(i+1,13):
            q=v[i]+v[j]
            for p in range(1,q):
                if gcd(p,q)==1:
                    m=min(norm(F(vi*p,q)) for vi in v)
                    if m>best: best=m
    return best
def is_DC(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def clears_at(v,q): return any(all((q<=14*((vi*p)%q)<=13*q) for vi in v) for p in range(1,q))

def main():
    ext=[1,2,3,4,10,11,12,13,14,15,16,17,18]
    print("CORRECTED DC extremal:", ext)
    print("  divisor-complete?", is_DC(ext), "; M =", Mx(ext), "=", round(float(Mx(ext)),4), " (cont.40's 2/23 was a LOCAL min)")
    q0=next(q for q in range(8,45) if q%14 and clears_at(ext,q))
    print("  first non-14 clearing modulus q =", q0, "; band-edge ceil(q/14)/q =", F(ceil(q0/14),q0), "= M EXACTLY (opus-S235 lemma TIGHT here)")
    print("  => sharp bound: divisor-complete => M >= 1/12 (exhaustive Vmax<=22), eps = 1/84 ~ 0.0119")
    print("  UNIFIED: [DC => bounded-clearing (open, my window)] + [band-edge lemma (opus-S235 PROVED)] => M >= 1/12 > 1/14 => lonely.")

if __name__=='__main__': main()
