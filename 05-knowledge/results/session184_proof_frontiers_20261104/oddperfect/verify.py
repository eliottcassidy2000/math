"""Exact sparse-domain, no-wrap modular selection tests. Standard library only."""
from itertools import combinations, product
from math import isqrt
from pathlib import Path
import json, random

def require(x,msg):
    if not x: raise RuntimeError(msg)
def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def order(p,r):
    x=p%r;d=1
    while x!=1: x=x*p%r;d+=1
    return d

def rows(S,q):
    ans=[]
    for i,p in enumerate(S):
        ds={order(p,r) for r in S if r!=p}
        ns=sorted(n for n in ds if n>1 and (n%4==2 if p==q else n%2==1) and all(d in ds for d in range(2,n+1) if n%d==0))
        rr=[]
        for n in ns:
            v=(p**n-1)//(p-1); v2=0
            while v%2==0:v//=2;v2+=1
            w=[]
            for r in S:
                k=0
                while v%r==0:v//=r;k+=1
                w.append(k)
            if v==1 and v2==int(p==q):
                w[i]-=n-1;rr.append(tuple(w))
        ans.append(rr)
    return ans

def moduli(R,d):
    if any(not rr for rr in R): return None
    low=[sum(min(c[j] for c in rr) for rr in R) for j in range(d)]
    high=[sum(max(c[j] for c in rr) for rr in R) for j in range(d)]
    return [1+max(abs(lo),abs(hi)) for lo,hi in zip(low,high)]

def count_mod(R,m):
    d=len(m);state={(0,)*d:1}
    for rr in R:
        new={}
        for old,k in state.items():
            for c in rr:
                x=tuple((old[j]+c[j])%m[j] for j in range(d))
                new[x]=new.get(x,0)+k
        state=new
    return state.get((0,)*d,0)

def count_exact(R,d):
    return sum(all(sum(c[j] for c in cs)==0 for j in range(d)) for cs in product(*R))

def main():
    rng=random.Random(184); trials=0; positives=0
    for _ in range(1000):
        d=rng.randrange(1,4);R=[[tuple(rng.randrange(-3,4) for _ in range(d)) for _ in range(rng.randrange(1,4))] for _ in range(rng.randrange(1,5))]
        exact=count_exact(R,d); m=moduli(R,d)
        require(count_mod(R,m)==exact,'no-wrap mismatch'); positives+=exact>0;trials+=1
    # Fractional feasibility is not integer selection: one block can choose -1 or +1.
    R=[[(-1,),(1,)]]
    require(count_mod(R,[2])==0 and count_exact(R,1)==0,'integrality control')
    # A too-small modulus can admit a spurious zero. Our bound forbids this.
    R=[[(2,)]]
    require(count_mod(R,[2])==1 and count_mod(R,moduli(R,1))==0,'wrap control')
    # Genuine perfect numbers 6 and 28, target e_2 subtracted once.
    for R in [[[(-1,1)],[(1,-1)]], [[(-2,1)],[(2,-1)]]]:
        require(count_mod(R,moduli(R,2))==count_exact(R,2)==1,'positive integer control')
    primes=[p for p in range(3,24,2) if prime(p)]; cases=0;closed=0;prefixes=0;maxprefix=0
    for k in (2,3,4):
        for S in combinations(primes,k):
            for q in S:
                if q%4!=1:continue
                R=rows(S,q);cases+=1;closed+=sum(map(len,R))
                # Any subset of nonempty rows is a nontrivial test of the exact modular compiler.
                Rnon=[rr for rr in R if rr]
                for t in range(1,len(Rnon)+1):
                    Rp=Rnon[:t];m=moduli(Rp,k)
                    require(count_mod(Rp,m)==count_exact(Rp,k),'arithmetic prefix mismatch')
                    prefixes+=1;maxprefix=max(maxprefix,t)
                if all(R):
                    require(count_mod(R,moduli(R,k))==count_exact(R,k),'arithmetic full mismatch')
                else: require(count_exact(R,k)==0,'empty domain control')
    out={'status':'PASS','random_signed_catalogues':trials,'positive_random_catalogues':positives,'support_euler_cases':cases,'exact_closed_rows':closed,'arithmetic_nonempty_prefix_checks':prefixes,'maximum_prefix_sources':maxprefix,'controls':['fractional feasible but modular impossible','small-modulus false positive detected','two positive integer perfect-number balances'],'global_odd_perfect_nonexistence':'NOT PROVED'}
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
