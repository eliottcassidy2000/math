#!/usr/bin/env python3
"""clustered_killer_gluing_kps_S128c53.py -- kind-pasteur S128 cont.53.
CLOSE THE CLUSTERED-KILLER CASE via THM-933 two-block gluing.
S = P (core, slow) u K (killer block, fast, internally CLUSTERED so balance-iteration fails).
THM-933 (BG) with m=2:   mu(S_rho(P) n S_rho(K)) >= d_P*d_K - q(K)*M(P),  M(P)=sum of core speeds.
Compute EXACTLY: delta(P), delta(K) (safe-set measures at rho=1/14), q(K) (discrepancy
sup H - inf H), and test positivity. Does the actual gap (min K > 13 max P) suffice?"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(933)
RHO=F(1,14)
def safe_intervals(B):
    """exact S_rho(B) = {t in [0,1): ||x t|| >= rho for all x in B} as sorted (lo,hi) list."""
    pts=set([F(0),F(1)])
    for x in B:
        for k in range(0,x+1):
            for off in (-RHO,RHO):
                t=F(k,1)/x + off/x
                if 0<=t<=1: pts.add(t)
    pts=sorted(pts)
    out=[]
    for a,b in zip(pts,pts[1:]):
        if b<=a: continue
        mid=(a+b)/2
        ok=True
        for x in B:
            y=x*mid; fy=y-int(y); fy=fy+1 if fy<0 else fy
            if min(fy,1-fy) < RHO: ok=False; break
        if ok:
            if out and out[-1][1]==a: out[-1]=(out[-1][0],b)
            else: out.append((a,b))
    return out
def measure(iv): return sum(b-a for a,b in iv)
def q_disc(B):
    """q(B) = sup H - inf H, H(u)=int_0^u (1_S - delta). Piecewise linear; extrema at endpoints."""
    iv=safe_intervals(B); d=measure(iv)
    # walk breakpoints accumulating H
    pts=[F(0)]
    for a,b in iv: pts.append(a); pts.append(b)
    pts.append(F(1)); pts=sorted(set(pts))
    H=F(0); vals=[F(0)]
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2
        inS=any(lo<=mid<=hi for lo,hi in iv)
        H += (b-a)*((1 if inS else 0) - d)
        vals.append(H)
    return max(vals)-min(vals), d
print("== THM-933 two-block gluing on CLUSTERED-killer families ==")
print("   condition for positivity: d_P*d_K > q(K)*M(P)   (M(P)=sum core speeds)")
rows=[]
for trial in range(14):
    j=random.choice([2,3])
    core=sorted(random.sample(range(1,12), 13-j))
    m=max(core)
    base=random.randint(13*m+1, 13*m+40)
    # CLUSTERED killers: all within a factor 2 of each other (balance-iteration fails here)
    K=sorted(set(base+random.randint(0,base//2) for _ in range(j)))
    if len(K)!=j: continue
    S=sorted(core+K)
    if len(set(S))!=13: continue
    dP=measure(safe_intervals(core))
    qK,dK=q_disc(K)
    MP=sum(core)
    lhs=dP*dK; rhs=qK*MP
    rows.append((float(lhs-rhs),j,tuple(core),tuple(K),float(dP),float(dK),float(qK),MP))
    print("   j=%d coreMax=%2d minK=%5d: d_P=%.4f d_K=%.4f q(K)=%.2e M(P)=%3d | d_P d_K=%.4f  q M(P)=%.4f  %s"%(
        j,m,min(K),float(dP),float(dK),float(qK),MP,float(lhs),float(rhs),
        "POSITIVE (closes)" if lhs>rhs else "negative"))
pos=[r for r in rows if r[0]>0]
print()
print("   %d/%d clustered configs close by (BG) at the actual gap (minK > 13 maxP)"%(len(pos),len(rows)))
if rows:
    worst=min(rows)
    print("   worst margin: %.4f (j=%d, core max %d, minK %d)"%(worst[0],worst[1],max(worst[2]),min(worst[3])))
    # what gap WOULD be needed? q(K) ~ c/minK, so needed minK ~ q(K)*minK*M(P)/(d_P d_K)
    need=[]
    for marg,j,core,K,dP,dK,qK,MP in rows:
        c = qK*min(K)          # q(K)*minK is roughly scale-free
        need.append(c*MP/(dP*dK))
    print("   implied gap requirement min(K) > c*M(P)/(d_P d_K): median %.1f x (vs actual min K)"%(sorted(need)[len(need)//2]))
print("DONE")
