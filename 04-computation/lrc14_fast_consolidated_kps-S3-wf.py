#!/usr/bin/env python3
"""Fast consolidated adversarial check (smaller trial counts for quick signal)."""
from fractions import Fraction as F
from math import gcd
import random
H=F(1,14)
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def safe_components(A,h=H):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append((sc[0][1])+(1-sc[-1][0]))
    return max(ws)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_primitive(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def case_of(S):
    S=sorted(set(S)); k=sum(1 for v in S if v>13)
    if k<=1: return "S1"
    if S[-1]<13*S[0]: return "S2"
    return "S3"
def best_margin(S):
    S=sorted(set(S)); bm=F(-1); bv=None
    for v in S:
        A=[u for u in S if u!=v]
        m=Wwidth(A)*7*v
        if m>bm: bm=m; bv=v
    return bm,bv

# ---- A: arc-width implication, 30000 trials ----
rng=random.Random(11); viol=0; tested=0
for _ in range(30000):
    nA=rng.randint(3,12)
    A=sorted(set(rng.randint(1,200) for _ in range(nA)))
    if len(A)<2: continue
    v=rng.randint(1,200)
    if v in A: continue
    W=Wwidth(A); thr=F(1,7*v)
    if W>thr:
        tested+=1
        if not safe_components(A+[v]): 
            viol+=1
            if viol<=5: print("IMPLICATION VIOLATION",A,v)
print(f"[A] arc-width implication: tested {tested}, violations {viol}")

# ---- B: S3 hunt, 15000 trials, ALWAYS compute M ----
def random_covering_13set(rng,max_speed=600):
    for _ in range(150):
        S=set()
        small=set(rng.sample(range(1,14),rng.randint(3,8)))
        S|=small
        V=rng.randint(20,max_speed); spread=rng.randint(2,50)
        need=13-len(S)
        if need<2: continue
        win=list(range(V,V+spread+1)); rng.shuffle(win)
        for w in win[:need]: S.add(w)
        while len(S)<13: S.add(rng.randint(1,13) if rng.random()<0.4 else V+rng.randint(0,spread+15))
        S=sorted(set(S))[:13]
        if len(S)!=13: continue
        S=sorted(S)
        if is_primitive(S) and is_covering(S): return S
    return None
rng=random.Random(999); n_s3=0; cfail=0; Mfail=0; worst=F(10**9); ws=None; wM=None; minM=F(1); minM_set=None
for _ in range(15000):
    S=random_covering_13set(rng)
    if S is None: continue
    if case_of(S)!="S3": continue
    n_s3+=1
    bm,bv=best_margin(S)
    M=Mval(S)
    if M<minM: minM=M; minM_set=S[:]
    if bm<worst: worst=bm; ws=S[:]; wM=M
    if bm<=1:
        cfail+=1
        if M<H: 
            Mfail+=1; print("!!! M<1/14:",S,M)
print(f"[B] S3 sets={n_s3}  C-fails(margin<=1)={cfail}  M<1/14={Mfail}")
print(f"    worst best-margin={worst}={float(worst):.5f} at {ws} (M={float(wM):.5f})")
print(f"    min M over S3={minM}={float(minM):.6f} at {minM_set}  (1/14={float(H):.6f})")

# ---- C: v=max determinism failures, 20000 trials ----
rng=random.Random(7); vmaxfail=0; ex=None
for _ in range(20000):
    S=random_covering_13set(rng)
    if S is None: continue
    if case_of(S)!="S3": continue
    S=sorted(S); vmax=S[-1]
    A=[u for u in S if u!=vmax]
    mv=Wwidth(A)*7*vmax
    if mv<=1:
        vmaxfail+=1
        if ex is None:
            bm,bv=best_margin(S)
            ex=(S[:],float(mv),bv,float(bm),Mval(S))
print(f"[C] v=max determinism fails (margin_vmax<=1): {vmaxfail}")
if ex: print(f"    example {ex[0]} vmax-margin={ex[1]:.5f} other-v={ex[2]} best-margin={ex[3]:.5f} M={float(ex[4]):.6f}")
