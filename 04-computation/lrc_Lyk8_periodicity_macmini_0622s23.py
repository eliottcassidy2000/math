#!/usr/bin/env python3
"""L_yK8(consec_{k-1} u {far}) is PERIODIC in far (THM-563: base fixed) => sup_far = finite period-max.
Verify the swap inequality sup_far L_yK8(consec_{k-1} u {far}) < L_yK8(consec_k) as a finite periodic max
=> the gK8 concentration's single-far binding case is RIGOROUS by periodicity (mac-mini method)."""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b):return a*b//gcd(a,b)
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
for k in (8,9,10):
    consec=tuple(range(k)); Lc=Lyk8(missdist(consec))
    base=tuple(range(k-1))  # consec_{k-1}
    L=1
    for e in base:
        if e>0:L=lcm(L,e)
    P=7*L
    # L_yK8(base u {far}) periodic in far with period P; sup over one period
    mx=F(-10); mf=0; per_ok = (Lyk8(missdist(base+(15,))) == Lyk8(missdist(base+(15+P,))))
    step = 1 if P<=6000 else max(1,P//6000)
    for far in range(15,15+P,step):
        v=Lyk8(missdist(base+(far,)))
        if v>mx: mx=v; mf=far
    print(f"k={k}: L_yK8(consec_k)={float(Lc):.5f}={Lc}  sup_far L_yK8(consec_{k-1}u far)={float(mx):.5f} at far={mf}  periodic? {per_ok}  swap-decrease? {mx<Lc}  margin={float(Lc-mx):.4f}")
print("\n=> if sup_far < L_yK8(consec_k) for all k: the gK8 concentration single-far swap is a FINITE")
print("   periodic max (THM-563 method) => consec maximizes L_yK8 over single-far configs RIGOROUSLY.")
