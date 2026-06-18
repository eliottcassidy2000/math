#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_boundary_kps-S3-wf

PIN DOWN the exact boundary: where does the K_j-ruler witness fail, and is the failure region
ALREADY covered by a PROVED lemma (S2 Lemma 1 global arc, or single-gap Collapse Lemma)?

The ruler-w0 witness failed on S=[4,7,8,9,11,12,28,32,36,40,44,48,52] (Vmax/Vmin=52/4=13, the S2/S3
boundary).  Hypothesis: the union of {S2 Lemma 1} u {single-gap Collapse Lemma} u {K_j-ruler}
covers ALL S3.  Test which lemma certifies each S3 set, and whether the UNION is 100%.

LEMMA inventory (all PROVED, elementary):
  L1 (S2 global arc): if 13*Vmin > Vmax, the arc (1/(14Vmin),13/(14Vmax)) is safe for ALL of S.
     Width = (13 Vmin - Vmax)/(14 Vmin Vmax) > 0.  => M(S)>=1/14 directly.  (No removal needed.)
  L_collapse (single-gap): exists j with J_j=[(14j+1)/(14 w0),(14j+13)/(14 wT)] P-safe, width>thresh.
  L_ruler (K_j multi-gap): exists cluster member r and period j of r with a >thresh window safe
     for all of S\{Vmax}.  (Generalizes L_collapse.)

We test, for a large S3 sample, the indicator of EACH lemma and report:
  - coverage of each lemma alone
  - coverage of the UNION
  - any set covered by NONE (the true residual of the angle) -- examine those by exact M(S).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def teeth_intervals(u,h=C):
    iv=[]
    for j in range(0,u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merged_intervals(iv):
    iv=sorted(iv); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return merged
def danger_merged(A,h=C):
    iv=[]
    for u in A: iv+=teeth_intervals(u,h)
    return merged_intervals(iv)
def safe_from_danger(merged):
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def safe_components(A,h=C): return safe_from_danger(danger_merged(A,h))
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)
def cand(S):
    S=sorted(set(S)); Cs=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cs.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cs.add(F(k,d)); k+=1
    Cs.add(F(1,2)); return Cs
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b

def L1_global(S):
    Vmin=min(S); Vmax=max(S)
    return 13*Vmin>Vmax  # global safe arc exists

def cluster_PL(S):
    Vmax=max(S)
    P=[v for v in S if v<=13]; L=sorted(v for v in S if v>13); Lp=[v for v in L if v!=Vmax]
    return P,Lp,Vmax

def L_collapse(S):
    P,Lp,Vmax=cluster_PL(S)
    if not Lp: return False
    w0=min(Lp); wT=max(Lp); thresh=F(1,7*Vmax)
    GP=safe_components(P) if P else [(F(0),F(1))]
    if 13*w0<=wT: return False
    jmax=(13*w0-wT-1)//(14*(wT-w0)) if wT>w0 else w0-1
    for j in range(0,max(0,jmax+1)):
        p=F(14*j+1,14*w0); q=F(14*j+13,14*wT)
        if q<=p or q>=1 or q-p<=thresh: continue
        for a,b in GP:
            if a<=p and q<=b: return True
    return False

def L_ruler(S):
    """K_j ruler over ALL cluster members (not just w0)."""
    P,Lp,Vmax=cluster_PL(S)
    if not Lp: return False
    thresh=F(1,7*Vmax)
    GA=safe_components(list(P)+list(Lp))
    for ruler in Lp:
        for j in range(0,ruler):
            lo=F(14*j+1,14*ruler); hi=F(14*j+13,14*ruler)
            if hi>1: hi=F(1)
            if lo>=hi: continue
            for a,b in GA:
                x=max(a,lo); y=min(b,hi)
                if y-x>thresh: return True
    return False

def make_S3(rng, spread_max, vlo, vhi):
    for _ in range(30000):
        small=set(rng.sample(range(1,14), rng.randint(3,9)))
        V0=rng.randint(vlo,vhi); s=rng.randint(1,spread_max)
        clsize=rng.randint(2, max(2,13-len(small)))
        cluster=set(); tries=0
        while len(cluster)<clsize and tries<300:
            cluster.add(V0+rng.randint(0,s)); tries+=1
        S=small|cluster
        while len(S)<13:
            S.add(rng.randint(V0,V0+s) if rng.random()<0.5 else rng.randint(1,13))
        if len(S)!=13: continue
        S=set(list(S))
        if len(S)!=13 or gcd_all(S)!=1 or not is_covering(S): continue
        Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
        if k>=2 and Vmax>=13*Vmin: return sorted(S)
    return None

def run(n=500, seed=401, spread_max=60, vlo=14, vhi=500, tag=""):
    rng=random.Random(seed)
    got=0; attempts=0
    cnt={'L1':0,'collapse':0,'ruler':0,'union':0}
    none_cases=[]
    while got<n and attempts<n*60:
        attempts+=1
        S=make_S3(rng,spread_max,vlo,vhi)
        if S is None: continue
        got+=1
        a=L1_global(S); b=L_collapse(S); c=L_ruler(S)
        if a: cnt['L1']+=1
        if b: cnt['collapse']+=1
        if c: cnt['ruler']+=1
        if a or b or c: cnt['union']+=1
        else: none_cases.append(S)
    flush(f"\n  [{tag}] S3 sets={got}")
    flush(f"    L1 (global arc, 13Vmin>Vmax): {cnt['L1']}/{got}")
    flush(f"    L_collapse (single-gap):      {cnt['collapse']}/{got}")
    flush(f"    L_ruler (K_j all members):    {cnt['ruler']}/{got}")
    flush(f"    UNION of the three:           {cnt['union']}/{got}")
    flush(f"    covered by NONE:              {len(none_cases)}")
    for S in none_cases[:6]:
        M=Mval(S)
        flush(f"      NONE: {S}  M={M}={float(M):.5f}  >=1/14? {M>=C}")
    return none_cases

if __name__=="__main__":
    flush("="*70)
    flush("BOUNDARY: union of proved lemmas vs full S3")
    flush("="*70)
    r1=run(500, 401, 60, 14, 500, "general")
    r2=run(400, 402, 14, 14, 500, "tight spread<=14")
    r3=run(400, 403, 120,14, 600, "loose spread<=120")
    r4=run(400, 404, 45, 14, 80,  "small V0 (S2/S3 boundary)")
    allnone=r1+r2+r3+r4
    flush(f"\n  TOTAL covered-by-NONE across regimes: {len(allnone)}")
