#!/usr/bin/env python3
"""Independent exact fixed-degree projective response recurrence checks."""
from fractions import Fraction as F
from math import gcd, factorial
from collections import defaultdict, Counter
from functools import lru_cache
import json, argparse

def pow2(e):
    return F(2**e) if e >= 0 else F(1, 2**(-e))

def fixed_degree(N, projective=True):
    # Store only still-observable coordinates; length 1 is implicit for gauge states.
    lengths=list(range(N if N%2 else N-1,0,-2))
    states={(N,(0,)*len(lengths)):F(1)}
    counts=[]
    for stage,p in enumerate(lengths):
        remaining=lengths[stage+1:]
        nxt=defaultdict(F)
        for (r,h),w in states.items():
            active=[j for j in lengths[stage:] if j<=r]
            H=dict(zip(active,h))
            for k in range(r//p+1):
                rr=r-k*p
                e=k*(p-1)//2+p*k*(k-1)//2+k*H.get(p,0)
                if projective:
                    e+=rr*k
                    hh=tuple(H[j]+k*(gcd(p,j)-j) for j in remaining if j<=rr)
                else:
                    hh=tuple(H[j]+k*gcd(p,j) for j in remaining if j<=rr)
                nxt[(rr,hh)]+=w*pow2(e)/ (p**k*factorial(k))
        states=nxt
        counts.append(len(states))
    result=states[(0,())]
    return result,counts

def partitions(n,lo=1):
    if n==0:
        yield ()
    for p in range(lo,n+1,2):
        for tail in partitions(n-p,p):
            yield (p,)+tail

def direct(n):
    total=F()
    for lam in partitions(n):
        e=sum((p-1)//2 for p in lam)+sum(gcd(p,q) for i,p in enumerate(lam) for q in lam[:i])
        z=1
        for p,k in Counter(lam).items(): z*=p**k*factorial(k)
        total+=F(2**e,z)
    return total

@lru_cache(None)
def hierarchy(n,s):
    if n==0:return F(1)
    total=F()
    odds=range(1,n+1,2)
    S=dict(zip(odds,s))
    for p in odds:
        rr=n-p
        ss=tuple(S[j]+gcd(p,j)-j for j in range(1,rr+1,2))
        total+=pow2((p-1)//2+S[p]+rr)*hierarchy(rr,ss)
    return total/n

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-order',type=int,default=30);ap.add_argument('--benchmark',type=int,nargs='*',default=[40,60,100]);args=ap.parse_args()
    values=[]
    for n in range(args.max_order+1):
        a,c=fixed_degree(n);b,d=fixed_degree(n,False)
        require(a==b==direct(n),f'direct mismatch {n}')
        if n<=20: require(a==hierarchy(n,(0,)*((n+1)//2)),f'hierarchy mismatch {n}')
        require(a.denominator==1,f'integrality {n}');values.append(int(a))
    require(values[:11]==[1,1,1,2,4,12,56,456,6880,191536,9733056],'OEIS prefix')
    benchmarks=[]
    for n in args.benchmark:
        a,c=fixed_degree(n);b,d=fixed_degree(n,False)
        require(a==b,f'benchmark mismatch {n}')
        benchmarks.append(dict(n=n,projective_peak=max(c),raw_peak=max(d),projective_total=sum(c),raw_total=sum(d),value=str(a),projective_by_stage=c,raw_by_stage=d))
    # Explicit reachable collision immediately before processing 3 and 1.
    mu=(9,)*5;nu=(15,)*3
    require(sum(mu)==sum(nu)==45,'collision mass')
    hm=[sum(gcd(p,j) for p in mu) for j in (1,3)]
    hn=[sum(gcd(p,j) for p in nu) for j in (1,3)]
    require(hm!=hn and [hm[1]-3*hm[0]]==[hn[1]-3*hn[0]],'strict quotient')
    print(json.dumps(dict(status='PASS',checked_through=args.max_order,hierarchy_through=min(20,args.max_order),values=values,benchmarks=benchmarks,collision=dict(mu=mu,nu=nu,raw=[hm,hn],projective=[0,0])),indent=2))
if __name__=='__main__':main()
