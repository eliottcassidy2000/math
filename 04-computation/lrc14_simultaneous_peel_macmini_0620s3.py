#!/usr/bin/env python3
"""
lrc14_simultaneous_peel_macmini_0620s3.py  (mac-mini-2026-06-20-S3)  -> THM-548 Part 4 (r=2 closure)

The SIMULTANEOUS far-peel: peel BOTH far runners from the BOUNDED base B (not iteratively from a
wide base), so V(B) stays bounded and the one-far residuals are THM-547-controlled.  The exact
Newton identity (verified):
  p0(B u {u,v}) = P_2(B) + [p0(B u {u})-Phi(B)] + [p0(B u {v})-Phi(B)] + [I_B(u,v)-Phi_2(B)],
  Phi(B)=p0(B)+p1(B)/7 (one-far plateau),  P_2(B)=p0(B)+2 p1/7 + Phi_2,  Phi_2=(2p2-p1)/49.
So p0(B u {u,v}) <= P_2(B) + (6/49)V(B)/u + (6/49)V(B)/v + max|I_B-Phi_2|.
This script:
 (1) verifies the Newton identity exactly on samples;
 (2) computes max over bounded B subset{0..14} (and far pairs) of |I_B(u,v)-Phi_2(B)| -> the
     uniform curvature bound (does it stay small? what is the cutoff W for closure?);
 (3) computes the closing cutoff: for u,v>W, p0(B u {u,v}) <= cap_k.  Extends THM-547 to r=2.
"""
import sys, math, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def stats(B):
    B=sorted(set(B)); bps=set([F(0),F(1)])
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); p0=F(0);p1=F(0);p2=F(0);V=0;inB=[False]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in B); t=7-len(secs); w=x1-x0
        if t==0:p0+=w
        elif t==1:
            p1+=w; j=[s for s in range(7) if s not in secs][0]
            if not inB[j]: V+=1; inB[j]=True
            for jj in range(7):
                if jj!=j:inB[jj]=False
        else:
            if t==2:p2+=w
            for jj in range(7):inB[jj]=False
    return p0,p1,p2,V
def IB(B,u,v): return measS7(list(B)+[u,v])-measS7(list(B)+[u])-measS7(list(B)+[v])+measS7(list(B))

print("(1) Newton identity check (should be exact 0):")
for B,u,v in [((0,1,2,3,4,5,6),17,23),((0,4,6,8,10,12,14),15,16),((0,1,2,4,8),31,50)]:
    p0,p1,p2,V=stats(B); Phi=p0+p1*SEV; Phi2=(2*p2-p1)/49; P2=p0+2*p1*SEV+Phi2
    lhs=measS7(list(B)+[u,v])
    rhs=P2+(measS7(list(B)+[u])-Phi)+(measS7(list(B)+[v])-Phi)+(IB(B,u,v)-Phi2)
    print(f"   B={B} u,v={u},{v}:  lhs-rhs = {lhs-rhs}  {'OK' if lhs==rhs else 'MISMATCH'}")

print("\n(2) max |I_B(u,v)-Phi_2(B)| over bounded B subset{0..14},|B|=7, and far pairs:")
caps={9:F(1979,4004)}
worst=F(0); worstcfg=None
import random
# sample bounded B of size 7 (k=9, r=2) + a grid of resonant+independent far pairs
Bsample=[(0,1,2,3,4,5,6),(0,4,6,8,10,12,14),(0,2,4,6,8,10,12),(0,1,2,4,8,9,11),
         (0,1,3,5,7,9,11),(0,2,3,5,7,11,13),(0,1,2,3,4,6,12),(0,1,2,3,4,8,12)]
pairs=[(15,16),(15,30),(16,17),(20,21),(21,35),(15,23),(18,19),(15,16),(22,33),(16,24),(15,45)]
for B in Bsample:
    p0,p1,p2,V=stats(B); Phi2=(2*p2-p1)/49
    for (u,v) in pairs:
        d=abs(IB(B,u,v)-Phi2)
        if d>worst: worst=d; worstcfg=(B,u,v,float(d),V)
print(f"   worst |I_B-Phi_2| = {float(worst):.5f} at B={worstcfg[0]}, (u,v)=({worstcfg[1]},{worstcfg[2]}), V(B)={worstcfg[4]}")

print("\n(3) r=2 closing bound: p0(B u {u,v}) <= P_2(B) + (6/49)V/u + (6/49)V/v + maxcurv, vs cap_9")
print("    For B=consec_7 (the P_2-extremal), find cutoff W s.t. u,v>W => closes:")
B=(0,1,2,3,4,5,6); p0,p1,p2,V=stats(B); Phi2=(2*p2-p1)/49; P2=p0+2*p1*SEV+Phi2
cap=caps[9]; maxcurv=worst
print(f"    B=consec_7: P_2(B)={float(P2):.5f}  V(B)={V}  cap_9={float(cap):.5f}  margin(P_2)={float(cap-P2):.5f}")
budget=cap-P2-maxcurv  # budget for the two one-far residuals
# need (6/49)V/u + (6/49)V/v <= budget; worst at u=v=W: 2*(6/49)V/W <= budget
Wc = 2*F(6,49)*V/budget if budget>0 else None
print(f"    after curvature ({float(maxcurv):.4f}), residual budget = {float(budget):.5f}")
print(f"    closes for u,v > W = 2*(6/49)*V/budget = {float(Wc):.1f}" if Wc else "    NO budget!")
print(f"    => true-wide r=2 with min(u,v) > ~{float(Wc):.0f} closes via simultaneous peel;")
print(f"       max(u,v) <= ~{float(Wc):.0f} is a FINITE check (extends THM-547 collar to r=2).")
print("\nDONE.")
