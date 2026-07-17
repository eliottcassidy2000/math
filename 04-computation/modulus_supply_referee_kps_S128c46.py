#!/usr/bin/env python3
"""modulus_supply_referee_kps_S128c46.py -- THE MODULUS SUPPLY: discrete (grid-q) B5 vs
continuous B5: |B5_q - B5| <= E(V)/q with E(V) = total arc-endpoint count summed over the
S_1..S_5 interval systems.  Referee the rate and the q0 formula on dissociated packets."""
import sys
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)
def depth_mu(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0: ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    mu=[F(0)]*(len(speeds)+1); depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    return mu
def B5_cont(V):
    mu=depth_mu(V); n=len(V)
    return mu[0]-sum(comb(d-1,5)*mu[d] for d in range(6,n+1))
def B5_grid(V,q):
    n=len(V); tot=0
    for p in range(q):
        t=F(p,q)
        d=0
        for v in V:
            x=(v*t)%1
            if min(x,1-x)<F(1,14): d+=1
        tot+= 1 if d==0 else 0
        if d>=6: tot-=comb(d-1,5)
    return F(tot,q)
cert=[307,425,541,671,800,944,1087,1413,1943,2147,2570,3056,3310]
small=[31,45,58,71,89,102,119,133,146,159,171,185,198]
for name,V in [("certified",cert),("small-13",small)]:
    b5=B5_cont(V)
    E=2*sum(V)  # endpoint count of the union system (each speed v contributes 2v endpoints)
    print("%s: continuous B5 = %.6f ; endpoint budget E = %d"%(name,float(b5),E))
    for q in (1009, 4001, 16001):
        bq=B5_grid(V,q)
        err=abs(float(bq-b5))
        print("   q=%d: B5_q = %+.6f  |err| = %.2e  vs E/q = %.2e  (rate OK: %s)"%(
            q,float(bq),err,E/q,err<=E/q))
    q0=int(E/float(b5))+1 if b5>0 else None
    if q0:
        bq=B5_grid(V,q0+1)
        print("   q0 = E/B5 = %d: B5_{q0+1} = %+.6f > 0: %s (THE MODULUS SUPPLY)"%(q0,float(bq),bq>0))
print("DONE")
