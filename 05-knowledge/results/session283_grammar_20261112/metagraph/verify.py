#!/usr/bin/env python3
"""Independent literal-metagraph checks; standard library, no imported lane code."""
import json, math
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache

def gate(x,msg):
    if not x: raise RuntimeError(msg)

def pairs(n): return list(combinations(range(n),2))
@lru_cache(None)
def hist(n,mask):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs(n)):
        a[i][j]=(mask>>k)&1; a[j][i]=1-a[i][j]
    h=[0]*(n-1)
    for i,j in pairs(n): h[sum(a[w][i]!=a[w][j] for w in range(n) if w not in (i,j))]+=1
    return tuple(h)
def step(n,p):
    N=n-2; M=n*(n-1)//2
    out=[F(0)]*(N+1)
    for d,x in enumerate(p):
        out[d]+=x*(1-F(2*N,M))
        if d: out[d-1]+=x*F(2*d,M)
        if d<N: out[d+1]+=x*F(2*(N-d),M)
    return out
def kraw(N,s,d):
    return sum((-1)**j*math.comb(d,j)*math.comb(N-d,s-j)
               for j in range(max(0,s-(N-d)), min(d,s)+1))
def rank_spectral(n,p,t):
    N=n-2; M=n*(n-1)//2
    return n-sum(sum(p[d]*kraw(N,s,d) for d in range(N+1))*F(M-4*s,M)**t for s in range(N+1))/2**N

def main():
    counts={'one_step_tournaments':0,'full_walk_distributions':0,'poissonized_edge_laws':0}
    for n in range(3,6):
        M=n*(n-1)//2
        for x in range(1<<M):
            p=hist(n,x)
            direct=[F(sum(hist(n,x^(1<<k))[d] for k in range(M)),M) for d in range(n-1)]
            gate(direct==step(n,p),f'one step {n},{x}')
            gate(rank_spectral(n,p,0)==n-p[0], 'spectral zero')
            counts['one_step_tournaments']+=1
    for n in range(3,5):
        M=n*(n-1)//2
        for x in range(1<<M):
            dist={x:1}; denom=1
            for t in range(5):
                actual=F(sum(c*(n-hist(n,y)[0]) for y,c in dist.items()),denom)
                gate(actual==rank_spectral(n,hist(n,x),t),f'walk {n},{x},{t}')
                nd={}
                for y,c in dist.items():
                    for k in range(M): nd[y^(1<<k)]=nd.get(y^(1<<k),0)+c
                dist=nd;denom*=M
                counts['full_walk_distributions']+=1
            # Each arc flipped independently with probability1/3: correlation1/3.
            # A separator XOR has correlation1/9.
            actual=F(sum(2**(M-y.bit_count())*(n-hist(n,x^y)[0]) for y in range(1<<M)),3**M)
            a=F(5,9);c=F(4,9);N=n-2;p=hist(n,x)
            predicted=n-sum(p[d]*c**d*a**(N-d) for d in range(N+1))
            gate(actual==predicted,f'poissonization {n},{x}')
            counts['poissonized_edge_laws']+=1
    # Infinite-size stationary mean is an exact rational identity for every n.
    stationary={str(n):str(F(n)-F(n*(n-1)//2,2**(n-2))) for n in (3,4,5,10,20,100)}
    out={'status':'PASS','controls':counts,'stationary_expected_rank':stationary,
         'scope':'labelled random arc reversals; no isoclass or path-count coverage claim'}
    print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__': main()
