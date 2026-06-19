"""Verify remaining claimed constants: smaller-k cap-pushed minima stability,
exact rho* floors k=5,6, and mu_min(k) for k=4..7 across larger caps (claim: flat)."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import time

def mu_exact(E):
    E=sorted(set(int(e) for e in E)); k=len(E)
    if k<=1: return F(1) if k==1 else F(0)
    TH=F(2,7); bp=set([F(0),F(1)]); diffs=set()
    for i in range(k):
        for j in range(i+1,k): diffs.add(E[j]-E[i])
    for d in diffs:
        for m in range(0,d+1): bp.add(F(m,d))
    obp=sorted(b for b in bp if F(0)<=b<=F(1)); refined=set(obp)
    for a,b in zip(obp,obp[1:]):
        mid=(a+b)/2; floors={e:(e*mid).__floor__() for e in E}
        order=sorted(E,key=lambda e:e*mid-floors[e])
        for t in range(k):
            if t==k-1:
                ef,el=order[0],order[-1]; slope=ef-el; const=F(1)-floors[ef]+floors[el]
            else:
                eh,elo=order[t+1],order[t]; slope=eh-elo; const=-floors[eh]+floors[elo]
            if slope!=0:
                xb=(TH-const)/slope
                if a<xb<b: refined.add(xb)
    refined=sorted(refined); tot=F(0)
    for a,b in zip(refined,refined[1:]):
        mid=(a+b)/2; pts=sorted(set((e*mid)%1 for e in E))
        if len(pts)==1: mg=F(1)
        else:
            g=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; g.append(pts[0]+1-pts[-1]); mg=max(g)
        if mg>TH: tot+=(b-a)
    return tot

def prim(E):
    g=0
    for e in E: g=gcd(g,e)
    return g==1

def exh(k,cap):
    best=None; bestE=None
    for c in combinations(range(1,cap+1),k-1):
        E=(0,)+c
        if not prim(E): continue
        m=mu_exact(E)
        if best is None or m<best: best,bestE=m,E
    return best,bestE

# CLAIM STEP3: for k=4,5,6,7 mu_min flat across last 3 caps. Test the CRITICAL claim:
# does a LARGER spread beat these small-k minima? The claim says NO for k<=11.
print("Test: does spread>cap LOWER mu_min for small k? (the bounded-spread reduction)", flush=True)
for k,caps in [(4,[14,25,40]),(5,[14,18,22]),(6,[14,17,20]),(7,[14,15])]:
    row=[]
    for cap in caps:
        t0=time.time(); b,E=exh(k,cap); row.append((cap,b,E,time.time()-t0))
    flat = all(row[i][1]==row[0][1] for i in range(len(row)))
    print(f"  k={k}: " + " | ".join(f"cap{cap}:{float(b):.5f}({b})" for cap,b,E,_ in row) + f"  FLAT={flat}", flush=True)

# But ALSO try a large-spread descent for k=7 to see if bounded-spread reduction really holds there
import numpy as np, random
def mu_fast(E,N=40000):
    A=np.array(sorted(set(E)),float); x=(np.arange(N)+0.5)/N
    P=np.mod(np.outer(A,x),1.0); P.sort(axis=0)
    mg=np.maximum(np.diff(P,axis=0).max(axis=0),P[0]+1-P[-1])
    return float(np.mean(mg>2/7))
print("\nLarge-spread descent for k=7 (does it beat 13/35=0.371?):", flush=True)
random.seed(5)
best=F(13,35); bestE=None
for _ in range(200):
    S=random.randint(15,60); st=sorted(set([0]+random.sample(range(1,S+1),6)))
    if len(st)!=7: continue
    cur=st; curmu=mu_exact(cur); imp=True; r=0
    while imp and r<15:
        r+=1; imp=False
        for idx in range(1,7):
            for nv in range(1,S+10):
                if nv in cur: continue
                cand=sorted(set(cur[:idx]+[nv]+cur[idx+1:]))
                if len(cand)!=7 or not prim(cand): continue
                if mu_fast(cand,15000)<float(curmu)-0.002:
                    me=mu_exact(cand)
                    if me<curmu: cur,curmu=cand,me; imp=True; break
            if imp: break
    if curmu<best: best,bestE=curmu,cur
print(f"  k=7 lowest via large-spread descent: {best}={float(best):.6f} at {bestE}", flush=True)
print(f"  cap-14 min was 13/35={float(F(13,35)):.6f}. Beaten? {best<F(13,35)}", flush=True)
