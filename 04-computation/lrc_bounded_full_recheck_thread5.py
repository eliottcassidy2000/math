#!/usr/bin/env python3
"""THREAD 5: INDEPENDENT recheck of the BOUNDED span<=14 finite check.
The prompt says 'BOUNDED span<=14'. Confirm: over ALL primitive E with 0 in E, max(E)<=14,
|E|=k, k=8..12, p0(E)<=cap_k with positive margin, and consec is binding. Whole-circle
breakpoint p0 (the canonical method). Confirm the partition span<=14 / span>14 is exhaustive
(no E escapes both)."""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
ALL_INNER=0b1111110
def sector_of(p): return int((p%1)*7)
def p0(E):
    E=sorted(set(int(x) for x in E)); nz=[e for e in E if e]
    if not nz: return F(0)
    bps=set([F(0),F(1)])
    for e in nz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); num=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2; mask=0
        for e in nz: mask|=1<<(sector_of(e*mid))
        if (mask&ALL_INNER).bit_count()==6: num+=hi-lo
    return num
def primitive(E):
    nz=[e for e in E if e]
    return reduce(gcd,nz)==1 if nz else False
caps={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
print("THREAD 5 independent BOUNDED span<=14 recheck (whole-circle breakpoint p0):")
for k in (8,9,10,11,12):
    cap=caps[k]; viol=0; mx=F(-1); arg=None; nb=0; minmargin=F(10); minarg=None
    consec=tuple(range(k)); pc=p0(consec)
    for combo in itertools.combinations(range(1,15),k-1):
        E=(0,)+combo
        if not primitive(E): continue
        nb+=1
        v=p0(E)
        if v>cap: viol+=1
        if v>mx: mx=v; arg=E
        margin=cap-v
        if margin<minmargin: minmargin=margin; minarg=E
    print(f"  k={k}: primitive sets={nb} cap={float(cap):.5f} maxp0={float(mx):.5f}@{arg} "
          f"VIOL={viol} | binding margin={float(minmargin):.6f}@{minarg} "
          f"consec_k p0={float(pc):.5f} (consec is max? {arg==consec})")
