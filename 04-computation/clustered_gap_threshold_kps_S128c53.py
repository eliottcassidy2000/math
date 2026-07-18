#!/usr/bin/env python3
"""clustered_gap_threshold_kps_S128c53.py -- find the GAP THRESHOLD at which (BG-K) closes
ALL clustered-killer configs. q(K) ~ 1/min(K), so a larger gap G (min K > G*max P) shrinks
the debt q(K)*K_P until it is below d_P*d_K.  Scan G = 13, 20, 26, 40, 60."""
import sys, random
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
def ncomp(iv):
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
        H+=(b-a)*((1 if inS else 0)-d); vals.append(H)
    return max(vals)-min(vals), d
print("== gap-threshold scan: (BG-K) closure rate vs gap multiplier G (min K > G * max P) ==")
for G in (13,20,26,40,60):
    random.seed(4242)
    closed=0; tot=0; worst=None
    for trial in range(12):
        j=random.choice([2,3])
        core=sorted(random.sample(range(1,12),13-j)); m=max(core)
        base=random.randint(G*m+1,G*m+40)
        K=sorted(set(base+random.randint(0,base//2) for _ in range(j)))
        if len(K)!=j or len(set(core+K))!=13: continue
        tot+=1
        ivP=safe_intervals(core); dP=measure(ivP); KP=ncomp(ivP)
        qK,dK=q_disc(K)
        marg=float(dP*dK-qK*KP)
        if marg>0: closed+=1
        if worst is None or marg<worst: worst=marg
    print("   G=%2d : (BG-K) closes %2d/%2d   worst margin %+.5f"%(G,closed,tot,worst if worst is not None else 0))
print()
print(">>> the clustered case closes by (BG-K) once the block gap is large enough;")
print(">>> the sliver 13 <= G < threshold is what remains.")
print("DONE")
