#!/usr/bin/env python3
"""beat_law_kps_S128c54.py -- kind-pasteur S128 cont.54.
TEST (data only, no pre-written conclusion): is the killer-block discrepancy q(K) governed
by the BEAT scale 1/d_min (d_min = smallest internal difference) rather than 1/max(K)?
Kernel fact (triangle inequality): ||v1 t||<L and ||v2 t||<L  =>  ||d t|| < 2L, d = v2-v1.
Scan PAIRS K={v, v+d} with v fixed and d swept, and report q(K), q*d, q*v."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
RHO=F(1,14)
def safe_intervals(B):
    pts=set([F(0),F(1)])
    for x in B:
        for k in range(0,x+1):
            for off in (-RHO,RHO):
                t=F(k,1)/x+off/x
                if 0<=t<=1: pts.add(t)
    pts=sorted(pts); out=[]
    for a,b in zip(pts,pts[1:]):
        if b<=a: continue
        mid=(a+b)/2; ok=True
        for x in B:
            y=x*mid; fy=y-int(y); fy=fy+1 if fy<0 else fy
            if min(fy,1-fy)<RHO: ok=False; break
        if ok:
            if out and out[-1][1]==a: out[-1]=(out[-1][0],b)
            else: out.append((a,b))
    return out
def measure(iv): return sum(b-a for a,b in iv)
def q_disc(B):
    iv=safe_intervals(B); dd=measure(iv)
    pts=[F(0)]
    for a,b in iv: pts.append(a); pts.append(b)
    pts.append(F(1)); pts=sorted(set(pts))
    H=F(0); vals=[F(0)]
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2
        inS=any(lo<=mid<=hi for lo,hi in iv)
        H+=(b-a)*((1 if inS else 0)-dd); vals.append(H)
    return max(vals)-min(vals), dd
# verify the kernel containment exactly on a sample
print("kernel check: {||v1 t||<1/14 and ||v2 t||<1/14} subset {||d t||<2/14}")
v1=200
bad=0
for d in (1,3,7,15,40):
    v2=v1+d
    # sample the overlap region on a fine rational grid
    ok=True
    Q=20000
    for p in range(Q):
        t=F(p,Q)
        def nd(x):
            fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
        if nd(v1*t)<RHO and nd(v2*t)<RHO:
            if not (nd(d*t) < 2*RHO): ok=False; break
    print("   d=%3d: containment holds on grid: %s"%(d,ok))
    if not ok: bad+=1
print()
print("q(K) for pairs K={v, v+d}, v=200 (data only):")
print("  %-5s %-12s %-10s %-10s %-8s"%("d","q(K)","q*d","q*v","delta(K)"))
for d in (1,2,3,5,8,13,21,34,55,89,144):
    K=[v1, v1+d]
    q,dd=q_disc(K)
    print("  %-5d %-12.4e %-10.4f %-10.4f %-8.4f"%(d,float(q),float(q)*d,float(q)*v1,float(dd)))
print()
print("same sweep at v=400 (does q track d or v?):")
print("  %-5s %-12s %-10s %-10s"%("d","q(K)","q*d","q*v"))
for d in (1,3,8,21,55,144):
    K=[400,400+d]
    q,dd=q_disc(K)
    print("  %-5d %-12.4e %-10.4f %-10.4f"%(d,float(q),float(q)*d,float(q)*400))
print("DONE")
