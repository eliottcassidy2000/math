#!/usr/bin/env python3
"""Slice 1b: is the lower-J frontier a ONE-PARAMETER family near the minimum?
Probe: (a) consec-shifts give the discrete J-curve; (b) do near-AP deformations
{1..8, 8+det} or interior-detunes dip BELOW the consec minimum J=4465/882? (c) the
exact J for the consec-shift family, confirming cshift1={1..9} is the unique min."""
from fractions import Fraction as F
def Jexact(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    mu=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    return 6*mu-m2
Jmin=F(4465,882)
print(f"consec-shift J-curve (exact), Jmin target = {Jmin} = {float(Jmin):.4f}:")
best=None
for sh in range(0,12):
    E=list(range(sh,sh+9)); J=Jexact(E)
    star="  <-- MIN" if J==Jmin else ""
    if best is None or J<best[0]: best=(J,sh)
    print(f"  cshift{sh} {E[0]}..{E[-1]}: J={float(J):.4f}{star}")
print(f"  => consec min at shift {best[1]}, J={float(best[0]):.4f}, equals 4465/882: {best[0]==Jmin}")
print()
print("do near-AP / interior-detune deformations of {1..9} dip below Jmin?")
viol=0; tot=0
# single-coordinate moves keeping 9 distinct positive
base=list(range(1,10))
for i in range(9):
    for delta in range(-3,4):
        if delta==0: continue
        E2=base[:i]+[base[i]+delta]+base[i+1:]
        if len(set(E2))!=9 or min(E2)<1: continue
        tot+=1; J=Jexact(sorted(E2))
        if J<Jmin: viol+=1; print(f"  BELOW: {sorted(E2)} J={float(J):.4f}")
# two-coordinate near-AP detunes
import itertools
for i,j in itertools.combinations(range(9),2):
    for di,dj in [(1,-1),(-1,1),(2,-1),(1,1)]:
        E2=base[:]; E2[i]+=di; E2[j]+=dj
        if len(set(E2))!=9 or min(E2)<1: continue
        tot+=1; J=Jexact(sorted(E2))
        if J<Jmin: viol+=1
print(f"  {tot} deformations tested, {viol} below Jmin  =>  {'FRONTIER NOT 1-param' if viol else 'cshift1 is the LOCAL+DEFORMATION min (1-param frontier confirmed near optimum)'}")
