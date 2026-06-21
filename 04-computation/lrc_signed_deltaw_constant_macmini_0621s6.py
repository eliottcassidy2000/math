#!/usr/bin/env python3
"""
lrc_signed_deltaw_constant_macmini_0621s6.py  (mac-mini-2026-06-21-S6)
THE HONEST SIGNED CONSTANT for Delta_w (against HYP-2784's 125x-lossy absolute wall).
Delta_w = p0(B u {w}) - Phi(B), Phi(B)=p0(B)+p1(B)/7, B=consec_{k-1} (the binding base).
Compute EXACT signed Delta_w over w=15..150; report:
 - Delta_w * w  (the honest signed 'constant' -- is it bounded? by what? vs absolute 15.8)
 - max positive Delta_w / margin
 - the low-mode (n<=13) signed Fourier head vs the exact (head captures? tail size?)
 - sqrt7-scaling probe.
"""
import itertools, sys, math, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def p1(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        if 7-len(set(sector_of(e*((x0+x1)/2)) for e in E))==1: tot+=x1-x0
    return tot
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(5,7),12:F(6,7)}
print(f"{'k':>3}{'|B|':>4}{'maxDw*w':>10}{'argw':>6}{'maxDw/marg':>12}{'absConst15.8/w@w*':>10}")
for k in (8,9,10,11,12):
    B=tuple(range(k-1)); Phi=measS7(B)+p1(B)/7; cap=caps[k]; margin=cap-(measS7(B)+p1(B)/7)
    best=F(-1); bestw=0; bestrel=F(-1)
    for w in range(15,151):
        dw=measS7(B+(w,))-Phi
        if dw*w>best: best=dw*w; bestw=w
        if dw/margin>bestrel: bestrel=dw/margin
    print(f"{k:>3}{k-1:>4}{float(best):>10.4f}{bestw:>6}{float(bestrel):>12.4f}")
print("\n=> 'maxDw*w' is the HONEST signed constant. If it is bounded ~O(1) and SMALL (<< absolute 15.8),")
print("   the signed bound Delta_w <= C_signed/w closes the binding window. Compare to absolute V*6/49.")
print(f"   sqrt(7)={math.sqrt(7):.4f}  (the Weil/Gauss-sum scale; 6/49 was the absolute constant)")
