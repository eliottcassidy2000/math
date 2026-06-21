#!/usr/bin/env python3
"""gK8 concentration binding = SINGLE-FAR (r=1). Verify max_{far>=15} L_yK8(B u {far}) <= L_yK8(consec_k)
for bounded/dilated bases B (|B|=k-1), via the far-window sup (periodicity tail beyond). This + bounded
finite check + (r>=2 lower) => gK8 concentration. The single-far is the periodicity domain (THM-563)."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(2)
def sector_of(p): return int((p%1)*7)
def missdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if t<=6: q[t]+=x1-x0
    return q
def Lyk8(q): return 10*q[0]+q[3]+10*q[6]
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
for k in (8,9,10):
    consk=tuple(range(k)); Lbound_max=Lyk8(missdist(consk))
    # also include other bounded bases' single-far -- check sup over far for many bounded bases B (|B|=k-1)
    print(f"k={k}: L_yK8(consec_k)={float(Lbound_max):.4f}  10cap={float(10*caps[k]):.4f}")
    worst_sf=F(0); wcfg=None; N=0
    # bases: consec_{k-1}, even-AP, random bounded
    bases=[tuple(range(k-1)), tuple(2*i for i in range(k-1)), tuple([0]+list(range(2,k)))]
    for _ in range(25): bases.append(tuple(sorted([0]+random.sample(range(1,15),k-2))))
    for B in bases:
        for far in range(15, 200):
            L=Lyk8(missdist(B+(far,))); N+=1
            if L>worst_sf: worst_sf=L; wcfg=(B,far,float(L))
    print(f"  max single-far L_yK8 [far in 15..200, {len(bases)} bases] = {float(worst_sf):.4f} @ B={wcfg[0]} far={wcfg[1]}")
    print(f"  single-far <= bounded max? {worst_sf < Lbound_max}  margin to 10cap = {float(10*caps[k]-worst_sf):.4f}")
