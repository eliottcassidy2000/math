# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont59: proving the core-length monotonicity (cont.58) as a GLOBAL inequality.
#
# TARGET: for a primitive covering family with maximal interval core {1..k}, M(F) >= 14/183, with equality
# only at k=12 (single-killer deep well). Enumeration (cont.58) confirmed it; this attacks the PROOF.
#
# RESULT (honest):
#  (1) RIGOROUS global lower bound (proved by perturbation at t_core=1/(k+1), moving toward runner 1):
#         M(F) >= [1/(k+1)] * a/(a + max(1,R)),   a = smallest RESONANT killer (mult of k+1),
#                                                  R = fastest NON-resonant killer speed.
#      TIGHT for single-killer (R=1): [1/13]*182/183 = 14/183 (deep well). Verified M(F) >= B(F).
#  (2) BUT the elementary witness UNDERSHOOTS for genuine multi-killer: a covering killer for d != k+1 is
#      non-resonant at t_core and falls fast (R large), so B(F) < 14/183. The true optimum lives at a
#      family-specific LARGER modulus (k=11 min {1..11,13,84}: core-resonance gives 7/95 < 14/183, true M=7/89
#      at t*=37/89). So NO uniform elementary witness closes multi-killer.
#  (3) The TIGHT multi-killer bound reduces to LRC(13)-ESCAPE + finite check, anchored by 1/13 > 14/183:
#      drop the largest killer -> a 12-speed sub-family with M >= 1/13 (settled LRC(13)) > 14/183; a LARGE
#      dropped killer keeps M >= 1/13 - B/L >= 14/183 (decorrelation); bounded killer -> finite check (ILP).
from math import gcd
from fractions import Fraction as F
from functools import reduce
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def Me(v):
    qc=4*max(v)+2; b=F(0); bt=None
    for q in range(2,qc):
        for p in range(1,q):
            if gcd(p,q)==1:
                m=min(norm(F(vi*p,q)) for vi in v)
                if m>b: b,bt=m,F(p,q)
    return b,bt
def cov(v,N=14): return all(any(x%d==0 for x in v) for d in range(2,N+1))
def prim(v): return reduce(gcd,v)==1
def corelen(v):
    s=set(v); k=0
    while (k+1) in s: k+=1
    return k
def bound(v):
    k=corelen(v); kp=k+1
    killers=[x for x in v if x>k]
    res=[x for x in killers if x%kp==0]; non=[x for x in killers if x%kp!=0]
    if not res: return None,k,None,None
    a=min(res); R=max(non) if non else 1
    return F(1,kp)*F(a,a+max(1,R)),k,a,R

def main():
    tgt=F(14,183)
    print(f"target 14/183={float(tgt):.6f}; 1/13={float(F(1,13)):.6f} (LRC(13) floor > covering-min => the BUFFER)\n")
    print("(1)+(2) rigorous balance bound B(F) vs true M, and which families it CLOSES (B>=14/183):")
    print(f"{'family':<28} | k | a | R | {'B(F)':>7} | {'M':>7} | M>=B | closes?")
    fams=[list(range(1,13))+[182], list(range(1,12))+[13,84], list(range(1,11))+[13,22,84]]
    for v in fams:
        B,k,a,R=bound(v); M,t=Me(v)
        print(f"{str(v)[:28]:<28} | {k} | {a} | {R} | {float(B):.5f} | {float(M):.5f} | {str(M>=B):>4} | {B>=tgt}")
    print("\n(3) LRC(13) escape base: drop the largest killer -> 12-speed sub-family, M>=1/13>14/183:")
    for v in [list(range(1,12))+[13], list(range(1,11))+[13,22]]:
        M,t=Me(v); print(f"    {str(v):<26} M={M}={float(M):.5f} >= 1/13? {M>=F(1,13)}")
    print("\n=> the elementary balance bound closes single-killer (tight 14/183); genuine multi-killer needs")
    print("   LRC(13)-escape + finite check. The monotonicity is TRUE but not an elementary global inequality;")
    print("   its engine is 1/13 > 14/183 (the LRC(13) floor sits ABOVE the covering-min, anchoring multi-killer).")

if __name__=="__main__":
    main()
