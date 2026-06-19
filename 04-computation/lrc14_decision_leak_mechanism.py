#!/usr/bin/env python3
# Why does the nearest competitor {0..7,9} lose to consec at k=9? Show the signed moment/p-vector
# decomposition: which p_t mass moves, and how g weights it.
from fractions import Fraction as F
from math import gcd
from functools import reduce

def Nat(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi<=lo: continue
        t=Nat(E,(lo+hi)/2); p[t]+=hi-lo
    return p
g=[F(-(t-2)*(t-3)*(t-6),36) for t in range(7)]
C=list(range(9)); A=[0,1,2,3,4,5,6,7,9]
pC=dist(C); pA=dist(A)
print("g(t) =", [str(x) for x in g])
print(f"{'t':>2} {'p_consec':>12} {'p_comp':>12} {'dp':>14} {'g(t)':>8} {'g.dp':>14}")
tot=F(0)
for t in range(7):
    dp=pA[t]-pC[t]; contrib=g[t]*dp; tot+=contrib
    print(f"{t:>2} {float(pC[t]):>12.6f} {float(pA[t]):>12.6f} {float(dp):>14.6f} {float(g[t]):>8.4f} {float(contrib):>14.6f}")
print(f"\ndL_y = sum g.dp = {tot} = {float(tot):.6f}  (negative => consec wins)")
print(f"L_y(consec)={float(sum(g[t]*pC[t] for t in range(7))):.6f}  L_y(comp)={float(sum(g[t]*pA[t] for t in range(7))):.6f}")
# cap closure check
print(f"\ncap_9 ~ 0.49426 (THM-535 lower-bounded, paper-stated). L_y(consec)=0.492876 < cap: {0.492876<0.49426}")
