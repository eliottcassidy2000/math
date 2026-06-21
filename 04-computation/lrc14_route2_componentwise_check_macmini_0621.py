#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- where does componentwise (sorted) T-vector dominance hold/fail?
The 12-vector is (T+(a),T-(a))_{a=1..6}. consec maximizes the SUM (WIN). Question:
does consec also componentwise-dominate (sorted T desc) every rival? i.e. is consec's
sorted T-vector >= every rival's sorted T-vector? If yes at all k, WIN-max is trivial.
The header of the old script claims k=10 (0..7,8,12) breaks componentwise. Verify.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)
def sector_of(e,a,y):
    pos=F(e*a)+F(7*e)*y; return (pos.numerator//pos.denominator)%7
def covered(E,a,y): return len({sector_of(e,a,y) for e in E})==7
def breakpoints(E,a):
    bps={F(0),-HALF,HALF}
    for e in E:
        if e==0: continue
        lo_val=F(7*e)*(-HALF)+F(e*a); hi_val=F(7*e)*(HALF)+F(e*a)
        lo_i=min(lo_val,hi_val);hi_i=max(lo_val,hi_val); m=lo_i.numerator//lo_i.denominator
        while m<=hi_i.numerator//hi_i.denominator+1:
            y=F(m-e*a,7*e)
            if -HALF<=y<=HALF: bps.add(y)
            m+=1
    return sorted(bps)
def window_TpTm(E,a):
    bps=breakpoints(E,a); ivals=list(zip(bps,bps[1:])); Tp=F(0)
    for lo,hi in ivals:
        if hi<=0: continue
        lo2=max(lo,F(0))
        if covered(E,a,(lo2+hi)/2): Tp=hi
        else:
            if lo2==F(0): Tp=F(0)
            break
    Tp=min(Tp,HALF); Tm=F(0)
    for lo,hi in reversed(ivals):
        if lo>=0: continue
        hi2=min(hi,F(0))
        if covered(E,a,(lo+hi2)/2): Tm=-lo
        else:
            if hi2==F(0): Tm=F(0)
            break
    Tm=min(Tm,HALF); return Tp,Tm
def Tvec(E): 
    v=[]
    for a in range(1,7):
        Tp,Tm=window_TpTm(E,a); v+=[Tp,Tm]
    return v
def WIN(E): return sum(Tvec(E))
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))
def dominates(u,v):  # sorted desc, u >= v componentwise
    su=sorted(u,reverse=True); sv=sorted(v,reverse=True)
    return all(x>=y for x,y in zip(su,sv))

if __name__=="__main__":
    for k,span in [(8,14),(9,16),(10,18)]:
        C=consec(k); Cv=Tvec(C); Cwin=WIN(C)
        bank=[(0,)+c for c in itertools.combinations(range(1,span+1),k-1)]
        bank=[E for E in bank if is_full_residue(E)]
        ndom=0; nbeat=0; nondom=[]
        for E in bank:
            if tuple(E)==tuple(C): continue
            E=list(E); w=WIN(E)
            if w>Cwin: nbeat+=1
            if dominates(Cv,Tvec(E)): ndom+=1
            else: nondom.append(E)
        print(f"k={k} span<={span}: shapes={len(bank)} WIN-beaters={nbeat} "
              f"consec-componentwise-dominates {ndom}/{len(bank)-1}")
        if nondom:
            print(f"   componentwise FAILS for {len(nondom)} shapes; examples:")
            for E in nondom[:3]:
                print(f"     {E}: WIN={float(WIN(E)):.6f} (<consec {float(Cwin):.6f}={Cwin})")
                print(f"        consec Tsorted={[str(x) for x in sorted(Cv,reverse=True)]}")
                print(f"        rival  Tsorted={[str(x) for x in sorted(Tvec(E),reverse=True)]}")
