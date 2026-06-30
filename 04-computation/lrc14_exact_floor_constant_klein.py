#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The exact floor constant: inf R' over the covering family, vs 1/(2 zeta(2)) (klein-S14).

R' = m_S/(m_R m_Q), S=R u 14Q. HYP-3571 scan: inf R'=0.344 at R={1..13}\{7}, Q={1,2}. Here: the EXACT
rational R' at the binding case, a broader scan for the true inf, and comparison to the set-independent
1/(2 zeta(2))=3/pi^2 (mac-mini Gamma_0(14) claim) and the apex atom 4cos^2(3pi/7).
"""
import sys, os, math, itertools, random
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
M = __import__("lrc14_floor_CV_sheetcount_bound_macmini_20260629")
lonely_set, measure = M.lonely_set, M.measure

def Rprime(R,Q):
    R=tuple(sorted(set(x for x in R if x%14!=0))); Q=tuple(sorted(set(Q)))
    mR=measure(lonely_set(R)); mQ=measure(lonely_set(Q))
    if mR==0 or mQ==0: return None,mR,mQ
    S=tuple(sorted(set(R)|set(14*q for q in Q)))
    mS=measure(lonely_set(S))
    return mS/(mR*mQ), mR, mQ

if __name__=="__main__":
    Z2=F(1)  # placeholder; use float for zeta(2)
    zeta2=math.pi**2/6; bound=1/(2*zeta2)
    atom=4*math.cos(3*math.pi/7)**2
    print("="*78); print(" Exact floor constant inf R' (klein-S14)"); print("="*78)
    print(f" set-indep targets: 1/(2 zeta2) = 3/pi^2 = {bound:.6f};  apex atom 4cos^2(3pi/7) = {atom:.6f}")

    # binding case exact
    R=[x for x in range(1,14) if x!=7]; Q=[1,2]
    rp,mR,mQ = Rprime(R,Q)
    print(f"\n binding case R={{1..13}}\\{{7}}, Q={{1,2}}:")
    print(f"   m_R = {mR} = {float(mR):.6f};  m_Q = {mQ} = {float(mQ):.6f}")
    print(f"   R' = {rp} = {float(rp):.6f}")
    print(f"   R' vs 1/(2 zeta2): {float(rp):.4f} vs {bound:.4f}  (margin {float(rp)-bound:+.4f})")

    # broader scan for the true inf (size-valid coverings)
    rng=random.Random(7); fam=set()
    def add(Rs):
        Rs=tuple(sorted(set(x for x in Rs if 1<=x<=13)))
        if 2<=len(Rs)<=12: fam.add(Rs)
    base=list(range(1,14))
    for k in range(2,13): add(range(1,k+1))
    for x in base: add([y for y in base if y!=x])
    for x,y in itertools.combinations(base,2): add([z for z in base if z not in(x,y)])
    for _ in range(1200): add(rng.sample(range(1,14), rng.randint(2,12)))
    best=None
    for Rs in fam:
        Q=list(range(1,14-len(Rs)+1))
        rp2,_,_=Rprime(Rs,Q)
        if rp2 is None: continue
        if best is None or rp2<best[0]: best=(rp2,Rs,Q)
    print(f"\n scan ({len(fam)} coverings): inf R' = {best[0]} = {float(best[0]):.6f}")
    print(f"   at R={best[1]} Q={best[2]}")
    print(f"   inf R' >= 1/(2 zeta2)={bound:.4f}? {float(best[0])>=bound}  (margin {float(best[0])-bound:+.4f})")
    # is the exact value clean? show numerator/denominator
    print(f"   exact inf R' = {best[0]}  (num {best[0].numerator}, den {best[0].denominator})")
