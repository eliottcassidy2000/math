#!/usr/bin/env python3
"""
lrc14_boundary_value_macmini_0620s3.py  (mac-mini-2026-06-20-S3)  -> THM-548

Two checks for the true-wide region:

(1) SATURATION of the consecutive-pair curvature I_B(u,u+1) as u->infinity.
    The pair {u,u+1} has fixed internal relation (1,-1) => I_B(u,u+1) should converge to a
    FINITE limit (bounded), not diverge.  Verify u=20,40,...,2560.

(2) The FULLY-DECORRELATED BOUNDARY VALUE P_r(B) = lim p0(B u F) as r INDEPENDENT far runners ->infinity:
       P_r(B) = sum_M meas{B misses exactly M} * P(r iid uniform runners cover all of M)
              = sum_{t=0}^{6} ( sum_{|M|=t} meas(A_M) ) * c_t(r),
       c_t(r) = sum_{i=0}^{t} (-1)^i C(t,i) (1 - i/7)^r     [prob r runners hit all t sectors].
    This is the Fatou-type boundary limit (main term).  Verify p0(B u F) -> P_r(B) for a dissociated
    (Sidon-like) far set F, and check P_r(B) <= cap_k - margin for bounded B (the main-term safety).
"""
import sys, math, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

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
        xm=(x0+x1)/2
        if len(set(sector_of(e*xm) for e in E))==7: tot+=x1-x0
    return tot
def miss_profile(B):
    """meas{B misses exactly t sectors}, t=0..6, exact."""
    B=sorted(set(B)); bps=set([F(0),F(1)])
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); prof=[F(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; t=7-len(set(sector_of(e*xm) for e in B))
        if t<=6: prof[t]+=x1-x0
    return prof   # prof[t]=meas{exactly t missed}
def I_B(B,u,v):
    return measS7(list(B)+[u,v])-measS7(list(B)+[u])-measS7(list(B)+[v])+measS7(list(B))

print("="*80)
print("(1) SATURATION of consecutive-pair curvature I_B(u,u+1), B=consec_8")
print("="*80)
B=(0,1,2,3,4,5,6,7)
prev=None
for u in [20,40,80,160,320,640,1280,2560]:
    ib=I_B(B,u,u+1); d='' if prev is None else f'  (delta={float(ib-prev):+.6f})'
    print(f"   u={u:>5}: I_B(u,u+1) = {float(ib):.7f}{d}"); prev=ib
print("   => converges to a finite limit (bounded) if deltas -> 0.")

print("\n"+"="*80)
print("(2) FULLY-DECORRELATED BOUNDARY VALUE P_r(B) (Fatou limit), and p0(B u F) for dissociated F")
print("="*80)
def c_t(t,r):  # P(r iid uniform-on-7 runners hit all t given sectors)
    return sum((-1)**i*math.comb(t,i)*(1-i/7)**r for i in range(t+1))
def P_r(B,r):
    prof=miss_profile(B)
    return sum(float(prof[t])*c_t(t,r) for t in range(7))
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(5,7),12:F(6,7),13:1}
# dissociated far sets (primes-ish, no small relations) of size r
dissoc={1:[23],2:[23,41],3:[23,41,83],4:[23,41,83,167],5:[23,41,83,167,331]}
for B in [(0,1,2,3,4,5,6,7),(0,4,6,8,10,12,14),(0,1,2,4,8,9,11)]:
    prof=miss_profile(B)
    print(f"\nB={B}  miss-profile (t=0..6) = {[float(x) for x in prof]}")
    for r in range(1,6):
        k=len(B)+r
        Pr=P_r(B,r)
        F_=dissoc[r]; p0act=float(measS7(list(B)+F_))
        cap=float(caps.get(k,1))
        print(f"   r={r} (k={k}): P_r(B)={Pr:.5f}  p0(B u dissocF)={p0act:.5f}  cap_{k}={cap:.5f}  "
              f"margin(P_r)={cap-Pr:+.5f} {'OK' if Pr<cap else 'OVER!'}")
print("\nP_r(B) = the boundary value as r independent far runners ->infinity (main term).")
print("If P_r(B) <= cap_k with margin for all bounded B, true-wide reduces to bounding resonance corrections.")
print("\nDONE.")
