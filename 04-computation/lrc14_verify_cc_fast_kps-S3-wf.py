#!/usr/bin/env python3
"""
LEAN fast decisive check. Focus on the criterion C(S) (best_margin>1) only, over a LARGE
number of covering+primitive S3 sets at modest V0 (fast). Report min best_margin and any
C-failure. Also a small batch with Mval to confirm M>=1/14 directly on the lowest-margin sets.
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

small_pool=[
 [7,8,9,10,11,12,13],[5,7,8,9,10,11,12,13],[6,7,8,9,10,11,12,13],
 [4,7,8,9,10,11,12,13],[1,7,8,9,10,11,12,13],[2,7,8,9,10,11,12,13],
 [9,10,11,12,13],[8,9,10,11,12,13],[11,12,13],[1,2,7,8,9,10,11,12,13],
 [1,2,3,8,9,10,11,12,13],[3,5,7,8,9,10,11,12,13],
]
rng=random.Random(20260618)
worst=F(10**9); ws=None; tested=0; cfail=0
low_sets=[]  # keep lowest-margin sets for Mval check
for _ in range(200000):
    small=rng.choice(small_pool); csize=13-len(small)
    if csize<2: continue
    V0=14*rng.randint(1,40)  # V0 up to 560 -> fast
    spread=rng.randint(csize, rng.choice([14,25,45]))
    c14=[x for x in range(0,spread+1) if (V0+x)%14==0]
    if not c14: continue
    offs=set([rng.choice(c14)])
    while len(offs)<csize: offs.add(rng.randint(0,spread))
    S=sorted(set(small+[V0+o for o in sorted(offs)]))
    if len(S)!=13 or not is_primitive(S) or not is_covering(S): continue
    if case_of(S)!="S3": continue
    tested+=1
    bm,arg=best_margin(S)
    if bm<=1: cfail+=1; out("  !!! C-FAIL",S,float(bm),"arg",arg,"M",float(Mval(S)))
    if bm<worst:
        worst=bm; ws=(S[:],arg)
        low_sets.append((bm,S[:]))
        low_sets=sorted(low_sets,key=lambda x:x[0])[:25]
out(f"FAST: tested covering primitive S3 = {tested}")
out(f"C(S) failures (best_margin<=1) = {cfail}")
out(f"min best_margin = {worst} = {float(worst):.6f} at {ws}")
out("Lowest-25 margin sets: verifying M(S)>=1/14 directly:")
mbad=0
for bm,S in low_sets[:25]:
    M=Mval(S)
    flag="" if M>=H else " <<< M<1/14!!!"
    if M<H: mbad+=1
    out(f"  bm={float(bm):.5f} M={M}={float(M):.6f} >=1/14?{M>=H}{flag}  S={S}")
out(f"M<1/14 among lowest-margin sets: {mbad}")
out("DONEFAST")
