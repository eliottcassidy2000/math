#!/usr/bin/env python3
"""k=9,10 L_y(deg3) adversarial check: does consec maximize the THM-534 deg-3 functional?
   Does any shape breach L_y>cap (which would break THM-534's per-E certificate route)?"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)
def fm(E,R=6):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); S=[F(0)]*(R+1)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; L=x1-x0
        hit=set(int(7*e*xm)%7 for e in E)
        free=6-len([s for s in hit if s!=0])
        for r in range(R+1): S[r]+=L*comb(free,r)
    return S
Ly=lambda S: F(1)-F(13,18)*S[1]+F(4,9)*S[2]-F(1,6)*S[3]
cap={9:F(1979,4004),10:F(55,91)}
def gen(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g==1: out.append(E)
    return out
for k in [9,10]:
    ck=cap[k]; AP=tuple(range(k)); apLy=Ly(fm(AP))
    for maxE in [k+5,k+7]:
        shapes=gen(k,maxE)
        mx=(apLy,AP); beaters=0; overcap=[]
        for E in shapes:
            v=Ly(fm(E))
            if v>mx[0]: mx=(v,E)
            if v>apLy: beaters+=1
            if v>ck: overcap.append((float(v),E))
        print(f"k={k} maxE<={maxE} ({len(shapes)} sets): Ly(AP)={float(apLy):.5f} cap={float(ck):.4f} "
              f"consec_max={mx[1]==AP} beaters={beaters} OVERCAP={len(overcap)}", overcap[:3] if overcap else "")
print("DONE")
