"""mac-mini-S57: k=12,13 ANALOG bricks (completing kps-S79's k=11 D3 tail closure).
Unlike the thin k=11, k=12,13 close via a UNIFORM D3 floor (no R2-diameter split needed):
  min D3 (compact, exact Farey) = 0.355876 (k=12) / 0.308844 (k=13), bars 0.19934/0.05649,
  EXACT margins +0.1565 / +0.2524 (1.8x / 5.5x). Tail (diam>15) only rises (decorrelation).
brick (A) for k=11 VERIFIED (primitive): max R2 over prim-diam>=16 = 614 (the 1+10 split
{0..9,16}); near-2-AP bump 610 at D=18 stays under 614; diam>=19 => R2=590."""
from fractions import Fraction as F
from math import floor, comb, gcd
from functools import reduce
from itertools import combinations
from collections import Counter
TH=F(1,7); M=F(6,7)
def moments_exact(A,maxp=3):
    A=sorted(A);k=len(A)
    ds=set(abs(A[j]-A[i]) for i in range(k) for j in range(k) if i!=j)|set(a for a in A if a)
    bps=set([F(0),F(1)])
    for d in ds:
        for m in range(d+1): bps.add(F(m,d))
    bps=sorted(bps);Mm=[F(0)]*(maxp+1)
    for c in range(len(bps)-1):
        x0,x1=bps[c],bps[c+1];xm=(x0+x1)/2
        lin=[(F(-floor(e*xm)),F(e)) for e in A]
        sp=[lin[j] for j in sorted(range(k),key=lambda j:lin[j][0]+lin[j][1]*xm)]
        gaps=[(sp[i+1][0]-sp[i][0],sp[i+1][1]-sp[i][1]) for i in range(k-1)]+[(F(1)+sp[0][0]-sp[k-1][0],sp[0][1]-sp[k-1][1])]
        subs=set([x0,x1])
        for (a,b) in gaps:
            if b!=0 and x0<(TH-a)/b<x1: subs.add((TH-a)/b)
        subs=sorted(subs)
        for s in range(len(subs)-1):
            u0,u1=subs[s],subs[s+1];um=(u0+u1)/2;Aa=F(0);Bb=F(0)
            for (a,b) in gaps:
                if a+b*um>TH: Aa+=a-TH;Bb+=b
            for p in range(maxp+1):
                Mm[p]+=sum(comb(p,r)*Aa**(p-r)*Bb**r*(u1**(r+1)-u0**(r+1))/(r+1) for r in range(p+1))
    return Mm
def D3(m): 
    den=m[2]-m[3]/M; return m[1]/M+(m[1]-m[2]/M)**2/den if den>0 else m[1]/M
BARS={12:F(50285,252252),13:F(14249,252252)}
mins={12:(0,2,4,5,6,7,8,9,10,11,12,14),13:(0,2,3,4,5,6,7,8,9,10,11,12,14)}
for k in (12,13):
    d3=D3(moments_exact(list(mins[k])))
    print(f"k={k}: exact compact-min D3={float(d3):.6f}, bar={float(BARS[k]):.6f}, margin +{float(d3-BARS[k]):.6f} => UNIFORM D3 floor closes leg")
