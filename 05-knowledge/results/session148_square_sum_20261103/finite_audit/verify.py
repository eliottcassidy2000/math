#!/usr/bin/env python3
"""Independent exact audit of square-sum graph growth, paths, and controls."""
from math import isqrt
from functools import lru_cache
import json

def require(c,msg):
    if not c: raise ValueError(msg)

def by_pairs(n):
    return {(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if isqrt(a+b)**2 == a+b}

def by_birth(n):
    return {(k*k-j,j) for j in range(1,n+1) for k in range(isqrt(j)+1,isqrt(2*j-1)+1)}

def count_by_square(n):
    return sum(max(0,min(n,(k*k-1)//2)-max(1,k*k-n)+1) for k in range(2,isqrt(2*n-1)+1))

def calc(n,edges):
    adj=[[] for i in range(n)]
    for a,b in edges: adj[a-1].append(b-1);adj[b-1].append(a-1)
    target=(1<<n)-1
    @lru_cache(None)
    def dp(mask,v):
        if mask == target: return 1
        return sum(dp(mask|1<<w,w) for w in adj[v] if not mask & 1<<w)
    endpoints=[dp(1<<v,v) for v in range(n)]
    return sum(endpoints),endpoints

def check_path(n,edges,path):
    return sorted(path)==list(range(1,n+1)) and all(tuple(sorted((a,b))) in edges for a,b in zip(path,path[1:]))

def run():
    for n in range(1,501):
        e=by_pairs(n); eb=by_birth(n)
        require(e==eb, f'edge disagreement n={n}')
        require(len(e)==count_by_square(n),f'square-sum count n={n}')
        birth={x for x in e if x[1]==n}
        require(len(birth)==isqrt(2*n-1)-isqrt(n),f'birth n={n}')
    expected=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,0,0,0,0,0,6,0,20]
    rows=[]
    for n in range(1,26):
        e=by_pairs(n); oriented,ends=calc(n,e)
        require(oriented==expected[n-1],f'path mismatch {n}')
        require(n==1 or oriented%2==0,f'reversal parity {n}')
        rows.append({'n':n,'birth':isqrt(2*n-1)-isqrt(n),'edges':len(e),'edges_minus_n_minus_1':len(e)-(n-1),'oriented_paths':oriented,'endpoint_counts':{str(i+1):c for i,c in enumerate(ends) if c}})
    n=15;e=by_pairs(n);p=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    require(check_path(n,e,p),'valid witness rejected')
    used={tuple(sorted(x)) for x in zip(p,p[1:])}
    require(e-used=={(1,3)},'wrong unused edge')
    require({v for a,b in e for u,v in ((a,b),(b,a)) if u==8}=={1},'wrong 8 neighborhood')
    require({v for a,b in e for u,v in ((a,b),(b,a)) if u==10}=={6,15},'wrong 10 neighborhood')
    require(not check_path(n,e,p[:-1]+[8]),'hostile repeated vertex accepted')
    q=p.copy();q[1],q[2]=q[2],q[1]
    require(not check_path(n,e,q),'hostile nonsquare edge accepted')
    require(not check_path(n,e-{(1,15)},p),'hostile missing edge accepted')
    # A proposed false birth formula overcounts the self-loop exactly at square n.
    require(isqrt(2*16-1)-isqrt(16-1)!=isqrt(2*16-1)-isqrt(16),'square-self-loop negative control failed')
    primes=[x for x in range(2,98) if all(x%d for d in range(2,isqrt(x)+1))]
    require(len(primes)==25 and primes[-1]==97,'first 25 primes')
    prime_rows=[{'p':p,'birth':isqrt(2*p-1)-isqrt(p),'edges':len(by_pairs(p))} for p in primes]
    return {'checks':'PASS','scope':'n=1..500 edge enumeration; n=1..25 exact oriented Hamilton paths','rows':rows,'p15':p,'prime_rows':prime_rows,'hostile_controls':'PASS'}
if __name__=='__main__': print(json.dumps(run(),indent=2))
