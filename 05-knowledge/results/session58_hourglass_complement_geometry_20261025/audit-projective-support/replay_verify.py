#!/usr/bin/env python3
"""Independent integer/support/projective checks; validation survives -O."""
from itertools import product, combinations
from math import prod, isqrt, gcd, comb
from pathlib import Path
import json

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gp(p):
    for b in range(1,isqrt(p)+1):
        a=isqrt(p-b*b)
        if a*a+b*b==p:return (a,b)
    raise RuntimeError('No split prime representation')
def gaps(ps):
    out={}
    for j in product((-1,0,1),repeat=len(ps)):
        if not any(j) or next(v for v in j if v)!=1:continue
        z=(1,0)
        for p,v in zip(ps,j):
            if v==0: z=(p*z[0],p*z[1]);continue
            a,b=gp(p); w=(a,v*b);z=mul(z,mul(w,w))
        d=abs(2*z[0]*z[1]);require(d>0 and d not in out,'positive unique Gaussian gap')
        out[d]=j
    return out

def brute(e):
    return {e*e-a*a for a in range(1,e) if isqrt(2*e*e-a*a)**2==2*e*e-a*a}
rows=[]
for ps in ((),(5,),(13,),(5,13),(5,17),(5,13,17),(5,13,17,29)):
    e=prod(ps);gs=gaps(ps);k=len(ps)
    require(set(gs)==brute(e),'independent endpoint census')
    hist={}
    for d,j in gs.items():
        S=tuple(i for i,v in enumerate(j) if v);hist[S]=hist.get(S,0)+1
        t=prod(p for i,p in enumerate(ps) if i not in S)
        require(gcd(e*e,d)==t*t,'exact truncated gcd')
        small=gaps(tuple(ps[i] for i in S))
        require(d//(t*t) in small and d%(t*t)==0,'full unused-prime division')
        # Every coordinate subspace deletion, not only maximal stripping.
        zero=[i for i,v in enumerate(j) if not v]
        for mask in range(1<<len(zero)):
            T={zero[a] for a in range(len(zero)) if mask>>a&1};u=prod(ps[i] for i in T)
            childps=tuple(p for i,p in enumerate(ps) if i not in T)
            childj=tuple(v for i,v in enumerate(j) if i not in T)
            child=gaps(childps)
            require(d%(u*u)==0 and child.get(d//(u*u))==childj,'exact index-preserving division')
    for r in range(1,k+1):
        for S in combinations(range(k),r): require(hist[S]==2**(r-1),'support stratum size')
    require(len(gs)==(3**k-1)//2,'projective size')
    rows.append({'center':e,'primes':ps,'gap_count':len(gs),'support_size_counts':{r:sum(c for S,c in hist.items() if len(S)==r) for r in range(1,k+1)}})
fieldrows=[]
for k in range(1,6):
    for q in (2,3,5,7):
        reps=set()
        for v in product(range(q),repeat=k):
            if not any(v):continue
            lead=next(a for a in v if a);inv=pow(lead,-1,q)
            reps.add(tuple((a*inv)%q for a in v))
        hist={}
        for v in reps:
            S=tuple(i for i,a in enumerate(v) if a);hist[S]=hist.get(S,0)+1
        require(len(reps)==sum(q**i for i in range(k)),'projective finite-field count')
        for r in range(1,k+1):
            for S in combinations(range(k),r):require(hist[S]==(q-1)**(r-1),'torus stratum count')
        require(sum(comb(k,r)*(q-1)**(r-1) for r in range(1,k+1))==len(reps),'polynomial identity')
        fieldrows.append({'k':k,'q':q,'points':len(reps)})
    supports=[S for r in range(1,k+1) for S in combinations(range(k),r)]
    require(len(supports)==2**k-1,'monoid topological supports')
    require(sum(len(S)==1 for S in supports)==k,'rank-zero strata')
gs=gaps((5,13)); byj={j:d for d,j in gs.items()}
a,b,c=(byj[j] for j in ((1,0),(0,1),(1,1)))
require((a,b,c)==(4056,3000,3696),'explicit addition hostile')
require(all(x+y!=z for x,y,z in ((a,b,c),(a,c,b),(b,c,a))),'no transported Schur relation')
require(gcd(65**2,3000)==25 and 3000%125==0,'balanced excess does not change truncated gcd')
result={'status':'PASS','integer_centers':rows,'projective_fields':fieldrows,'hostile':{'center':65,'linearly_dependent_index_gaps':[a,b,c],'schur_relation':False},'boundary':{'k0':'empty','k2_closed_rank_zero_points':2,'k2_monoid_topological_points':3,'hom_points_not_identified_with_rank_zero':True}}
Path(__file__).with_name('verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
