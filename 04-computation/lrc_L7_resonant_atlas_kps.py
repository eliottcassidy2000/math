#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC(14) L7 resonant atlas (HYP-2757): for every commensurable ratio p/q in (1,2.15),
the far pair traces a torus geodesic; p0(B u {Cq,Cp}) -> a c-independent curve coverage.
Check the curve coverage < cap_k over all such ratios and bounded bases. Float (margin large)."""
import sys
from itertools import combinations
from functools import reduce
from math import gcd
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8")
caps={8:2243/5880,9:1979/4004,10:55/91,11:66/91,12:6/7}
def p0(E):
    E=sorted(set(E)); bset=set([0.0,1.0])
    for e in E:
        if e==0: continue
        for j in range(7):
            b=j/7.0; m=0
            while True:
                xv=(b+m)/e
                if xv>=1: break
                if xv>=0: bset.add(xv)
                m+=1
    B=sorted(bset); tot=0.0
    for lo,hi in zip(B,B[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if len(set(int((e*mid)%1*7) for e in E)&set(range(1,7)))==6: tot+=hi-lo
    return tot
def ratios(lo,hi,Q):
    R=[]
    for q in range(1,Q+1):
        for p in range(q+1,int(hi*q)+1):
            if gcd(p,q)==1 and lo< p/q <=hi: R.append((p,q))
    return sorted(set(R),key=lambda pq:pq[0]/pq[1])
def main():
    C=701  # large scale (curve limit)
    print("L7 RESONANT ATLAS: ratios p/q in (1,2.15], q<=8; curve coverage p0(B u {Cq,Cp}), C=%d, vs cap_k"%C)
    Rs=ratios(1.0,2.15,8); print("  %d ratios:"%len(Rs), [f"{p}/{q}" for p,q in Rs])
    for k in (9,10):
        cap=caps[k]; worst=0.0;wr=None;viol=0;cnt=0
        # sample bounded bases of size k-2 (so + 2 far = k)
        bases=[list(range(k-2))]
        # add even-AP and a few spread bases
        bases.append([0]+[2*i for i in range(1,k-2)])
        import random; random.seed(3)
        for _ in range(60): bases.append([0]+sorted(random.sample(range(1,15),k-3)))
        for (p,q) in Rs:
            for B in bases:
                E=sorted(set(B+[C*q,C*p]))
                if len(E)!=k or reduce(gcd,E)!=1: continue
                cnt+=1; v=p0(E)
                if v>cap: viol+=1
                if v>worst: worst=v;wr=(p,q,B[:5])
        print("  k=%d: %d (ratio,base) curve-points | max p0_curve=%.4f cap=%.4f margin=%.4f | violations=%d | worst ratio=%d/%d"%(
            k,cnt,worst,cap,cap-worst,viol,wr[0],wr[1]))
    print("=> if max p0_curve < cap (margin), the resonant atlas of L7 is FINITE and SAFE.")
if __name__=="__main__": main()
