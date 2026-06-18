#!/usr/bin/env python3
"""
Probe: covering S3 sets where the CLUSTER (not the small part) carries the hard large
divisors (8,9,11,13,14...). Small part is low {1..7}-ish. This is a structurally different
S3 regime; test whether C(S) can fail or margin dips below the claimed floor.

Construction: small part = subset of {1,2,3,4,5,6,7}; cluster L of large speeds near V0 whose
union supplies multiples of every q in 2..14 not covered by the small part. Vmax>=13*Vmin.
"""
import sys
from fractions import Fraction as F
from math import gcd
import random
def out(*a): print(*a); sys.stdout.flush()
H=F(1,14)
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
def best_margin(S):
    S=sorted(set(S)); best=F(-1); arg=None
    for v in S:
        m=Wwidth([u for u in S if u!=v])*7*v
        if m>best: best=m; arg=v
    return best,arg
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
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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

rng=random.Random(31337)
out("Hard-cluster S3 hunt (small part low, cluster carries hard divisors):")
tested=0; cfail=0; worst=F(10**9); ws=None; cfx=[]; mfx=[]
for _ in range(300000):
    sp_size=rng.randint(3,7)
    small=sorted(rng.sample(range(1,8), sp_size))
    csize=13-len(small)
    if csize<2 or csize>10: continue
    V0=rng.choice([14,28,42,84,140,280,420,2520])*rng.randint(1,30)
    spread=rng.randint(csize, rng.choice([20,45,90,180,360]))
    # build cluster greedily to cover all q in 2..14 not covered by small
    need=[q for q in range(2,15) if not any(s%q==0 for s in small)]
    cl=set()
    # try to place a multiple of each needed q inside [V0,V0+spread]
    ok=True
    for q in need:
        ms=[x for x in range(V0,V0+spread+1) if x%q==0]
        if not ms: ok=False; break
        cl.add(rng.choice(ms))
    if not ok: continue
    while len(cl)<csize and len(cl)<spread+1:
        cl.add(V0+rng.randint(0,spread))
    if len(cl)!=csize: continue
    S=sorted(set(small)|cl)
    if len(S)!=13: continue
    if not is_primitive(S) or not is_covering(S): continue
    if case_of(S)!="S3": continue
    tested+=1
    bm,arg=best_margin(S)
    if bm<worst: worst=bm; ws=(S[:],arg)
    if bm<=1:
        cfail+=1; M=Mval(S)
        if len(cfx)<10: cfx.append((S[:],bm,arg,M))
        if M<H and len(mfx)<10: mfx.append((S[:],M))
out(f"tested={tested}  C-failures(margin<=1)={cfail}")
out(f"min best_margin={worst}={float(worst):.6f} at {ws}")
for s,b,a,M in cfx[:10]: out("  C-FAIL",s,float(b),"arg",a,"M",float(M))
for s,M in mfx[:10]: out("  M<1/14",s,M,float(M))
out("DONE3")
