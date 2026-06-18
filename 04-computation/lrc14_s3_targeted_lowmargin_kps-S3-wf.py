#!/usr/bin/env python3
"""
Targeted adversarial hunt for LOW best-margin covering+primitive S3 13-sets.
Goal: find best-margin <= 1 (C(S) fails) or M(S) < 1/14 (true LRC counterexample).
Strategy: systematically construct S3 sets that MINIMIZE the criterion margin:
  - small part forced to supply small divisors (q=2..7 via {2,3,4,5,6,7} etc.)
  - large cluster supplies the big divisors (8..14) and is tight (small internal spread)
  - vary V0 across many scales, vary cluster offsets, vary Vmax, ensure covering+primitive.
"""
from fractions import Fraction as F
from math import gcd
import itertools, random

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

# Targeted construction: small part hits 2..7, cluster hits 8..14 with multiples.
# Cluster near V0; choose cluster elements that are multiples of 8,9,10,11,12,13,14 close together.
# We'll search V0 over scales and pick cluster multiples in a window.
results=[]
worst=F(10**9); worst_set=None; worst_M=None
Mfail=[]

rng=random.Random(2024)
# Approach: pick a base B, cluster = multiples of various q's in [B, B+spread].
for trial in range(300000):
    # small part: random subset of {1..13} of size 4..8 guaranteeing 2..7 coverage where possible
    sp_size=rng.randint(3,8)
    small=set(rng.sample(range(1,14),sp_size))
    # large window
    V0=rng.choice([20,30,40,50,75,100,150,200,300,500,800,1300,2000,3500,6000])
    spread=rng.randint(10,60)
    # cluster: pick up to (13 - len(small) -1) larges in [V0, V0+spread], prefer covering-completing multiples
    need_size=13-len(small)
    if need_size<2: continue
    window=list(range(V0, V0+spread+1))
    rng.shuffle(window)
    cluster=set()
    # greedily add to complete missing divisors
    def missing(S):
        return [q for q in range(2,15) if not any(v%q==0 for v in S)]
    cur=set(small)
    for w in window:
        if len(cluster)>=need_size: break
        cluster.add(w); cur.add(w)
    S=sorted(small|cluster)
    if len(S)!=13:
        # pad/trim
        S=sorted(set(S))
        while len(S)<13: S.append(S[-1]+rng.randint(1,5)); S=sorted(set(S))
        S=sorted(set(S))[:13]
    if len(S)!=13: continue
    if not is_primitive(S): continue
    if not is_covering(S): continue
    if case_of(S)!="S3": continue
    bm,bv=best_margin(S)
    if bm<worst:
        worst=bm; worst_set=S[:]; 
        worst_M=Mval(S)
    if bm<=1:
        M=Mval(S)
        results.append((S[:],bm,M))
        if M<H: Mfail.append((S[:],M))

print(f"Targeted search done. S3 covering+primitive sets with best-margin<=1: {len(results)}")
print(f"Worst (min) best-margin found: {worst} = {float(worst):.6f}")
print(f"   set={worst_set}  M={worst_M}={float(worst_M):.6f}  M>=1/14? {worst_M>=H}")
if results:
    print("Sets with margin<=1 (C fails) -- check M:")
    for s,bm,M in results[:15]:
        print(f"   {s} margin={float(bm):.4f} M={M}={float(M):.6f} M>=1/14? {M>=H}")
if Mfail:
    print("!!!! TRUE LRC COUNTEREXAMPLES (M<1/14):")
    for s,M in Mfail: print("   ",s,M,float(M))
else:
    print("No M<1/14 found (no LRC counterexample in this search).")
