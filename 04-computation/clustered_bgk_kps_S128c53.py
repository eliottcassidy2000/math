#!/usr/bin/env python3
"""clustered_bgk_kps_S128c53.py -- the SHARPER (BG-K) form of THM-933 on clustered killers.
(BG-K):  mu >= d_P*d_K - q(K)*K_P,  K_P = #circular components of S_rho(P)  (<< M(P)).
Crude (BG) used M(P)=sum of core speeds (~66) and FAILED. Supply exact component counts."""
import sys, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(933)
RHO=F(1,14)
def safe_intervals(B):
    pts=set([F(0),F(1)])
    for x in B:
        for k in range(0,x+1):
            for off in (-RHO,RHO):
                t=F(k,1)/x + off/x
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
def ncomp_circular(iv):
    """circular component count: merge wrap-around if 0 and 1 endpoints both in S."""
    if not iv: return 0
    n=len(iv)
    if iv[0][0]==F(0) and iv[-1][1]==F(1) and n>1: n-=1
    return n
def q_disc(B):
    iv=safe_intervals(B); d=measure(iv)
    pts=[F(0)]
    for a,b in iv: pts.append(a); pts.append(b)
    pts.append(F(1)); pts=sorted(set(pts))
    H=F(0); vals=[F(0)]
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2
        inS=any(lo<=mid<=hi for lo,hi in iv)
        H += (b-a)*((1 if inS else 0)-d); vals.append(H)
    return max(vals)-min(vals), d
print("== (BG-K): mu >= d_P*d_K - q(K)*K_P  with EXACT core component count K_P ==")
rows=[]
for trial in range(14):
    j=random.choice([2,3])
    core=sorted(random.sample(range(1,12),13-j)); m=max(core)
    base=random.randint(13*m+1,13*m+40)
    K=sorted(set(base+random.randint(0,base//2) for _ in range(j)))
    if len(K)!=j: continue
    if len(set(core+K))!=13: continue
    ivP=safe_intervals(core); dP=measure(ivP); KP=ncomp_circular(ivP)
    qK,dK=q_disc(K); MP=sum(core)
    lhs=dP*dK; rhs_K=qK*KP; rhs_M=qK*MP
    ok = lhs>rhs_K
    rows.append((float(lhs-rhs_K),ok,KP,MP))
    print("   j=%d minK=%5d: d_P=%.4f d_K=%.4f q(K)=%.2e | K_P=%2d (vs M(P)=%3d) | d_P d_K=%.4f  q*K_P=%.4f  %s"%(
        j,min(K),float(dP),float(dK),float(qK),KP,MP,float(lhs),float(rhs_K),
        "CLOSES" if ok else "negative"))
good=[r for r in rows if r[1]]
print()
print("   (BG-K) closes %d/%d clustered configs at the ACTUAL gap (minK > 13 maxP)"%(len(good),len(rows)))
if rows:
    print("   component count K_P is %.1fx smaller than the crude cap M(P) on average"%(
        sum(r[3] for r in rows)/max(1,sum(r[2] for r in rows))))
    print("   best margin %.4f | worst %.4f"%(max(r[0] for r in rows), min(r[0] for r in rows)))
print("DONE")
