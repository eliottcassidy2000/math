#!/usr/bin/env python3
"""live_floor_bridge_referee_kps_S128c47.py -- the weight-1 sampling bridge:
liveCount(q) >= q*(mu0 - E/q) with E = 2*sum(v).  Referee on positive-measure packets
(incl. structured ones with B5 < 0 where the B5-route fails but the live-route works)
and the tight control (mu0 = 0: no false floor)."""
import sys
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)
def mu0_of(V):
    ev=[]
    for v in V:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0: ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    mu=F(0); depth=0; last=F(0)
    for x,d in ev:
        if depth==0: mu+=x-last
        depth+=d
        if depth==0: last=x
    if True and F(1)>last and depth==0: mu+=F(1)-last
    return mu
def live(V,q):
    c=0
    for p in range(q):
        t=F(p,q); ok=True
        for v in V:
            x=(v*t)%1
            if min(x,1-x)<F(1,14): ok=False; break
        if ok: c+=1
    return c
packs=[("geometric(B5<0,mu0>0)",[5*2**i for i in range(13)]),
       ("deepwell(B5<0,mu0>0)",list(range(1,13))+[182]),
       ("tight(mu0=0)",list(range(1,14)))]
for name,V in packs:
    m0=mu0_of(V); E=2*sum(V)
    print("%s: mu0 = %.6f, E = %d"%(name,float(m0),E))
    for q in (2003, 8009):
        lv=live(V,q)
        floor=float(q*(m0-F(E,q)))
        print("   q=%d: live = %d  >=  q(mu0-E/q) = %.1f : %s"%(q,lv,floor,lv>=floor))
    if m0>0:
        q0=int(2*E/float(m0))+1
        lv=live(V,q0)
        print("   q0 = 2E/mu0 = %d: live = %d > 0: %s (THE LIVE FLOOR FIRES)"%(q0,lv,lv>0))
print("DONE")
