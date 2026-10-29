#!/usr/bin/env python3
"""Compressed realization of inherited deletion cokernel classification.
Never enumerate tournament classes or expand the rank-size diagonal for atlas.
"""
import argparse
from collections import Counter
from math import isqrt
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'enumeration'))
from baseline87 import degree_bucket_dp


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def primes_to(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]


def valuation(d,p):
    k=0
    while d%p==0:
        d//=p
        k+=1
    return k


def compressed_smith(n,counts):
    """Return [(invariant factor, multiplicity)], including unit factors."""
    if not 1<=n<len(counts):
        raise ValueError('need count prefix through n and n>=1')
    h=[counts[r]-(counts[r-1] if r else 0) for r in range(n)]
    rank=counts[n-1]
    require(all(x>=0 for x in h),'nonnegative harmonic counts')
    runs={}
    breaks={0,rank}
    for p in primes_to(n):
        buckets=Counter()
        for r,mult in enumerate(h):
            if mult:
                buckets[valuation(n-r,p)]+=mult
        position=0
        runs[p]=[]
        for v,mult in sorted(buckets.items()):
            position+=mult
            runs[p].append((position,v))
            breaks.add(position)
        require(position==rank,'valuation run mass')
    histogram=[]
    positions=sorted(breaks)
    for start,end in zip(positions,positions[1:]):
        factor=1
        for p,segments in runs.items():
            exponent=next(v for stop,v in segments if stop>start)
            factor*=p**exponent
        if histogram and histogram[-1][0]==factor:
            histogram[-1]=(factor,histogram[-1][1]+end-start)
        else:
            histogram.append((factor,end-start))
    require(sum(m for _,m in histogram)==rank,'Smith rank')
    require(all(b%a==0 for (a,_),(b,_) in zip(histogram,histogram[1:])),'invariant factor divisibility')
    bound=1+sum(max(valuation(d,p) for d in range(1,n+1)) for p in primes_to(n))
    require(len(histogram)<=bound,'compressed interval bound')
    # Determinant-order valuation check without constructing the group order.
    for p in primes_to(n):
        require(sum(m*valuation(d,p) for d,m in histogram)==sum(h[r]*valuation(n-r,p) for r in range(n)), 'torsion order valuation')
    return histogram,bound


def verify_small(counts):
    # Independent direct expansion and independent per-prime array sorting.
    for n in range(1,11):
        diagonal=[n-r for r in range(n) for _ in range(counts[r]-(counts[r-1] if r else 0))]
        actual=[1]*len(diagonal)
        for p in primes_to(n):
            powers=sorted(p**valuation(d,p) for d in diagonal)
            actual=[a*b for a,b in zip(actual,powers)]
        expected=sorted(Counter(actual).items())
        require(compressed_smith(n,counts)[0]==expected,'direct expanded Smith agreement')
    expected_small={1:[1],2:[2],3:[3],4:[1,4],5:[1,1,1,10]}
    for n,expected in expected_small.items():
        flat=[d for d,m in compressed_smith(n,counts)[0] for _ in range(m)]
        require(flat==expected,'independent tournament-matrix Smith agreement')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-n',type=int,default=200)
    args=parser.parse_args()
    require(args.max_n>=10,'max-n must be at least 10 for verification')
    if hasattr(sys,'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    counts,_=degree_bucket_dp(args.max_n)
    verify_small(counts)
    rows=[]
    for n in [5,10,25,50,100,200]:
        if n>args.max_n:
            continue
        histogram,bound=compressed_smith(n,counts)
        torsion={str(p):str(sum(m*valuation(d,p) for d,m in histogram)) for p in primes_to(n)}
        torsion={p:e for p,e in torsion.items() if e!='0'}
        rows.append({'n':n,'rank':str(counts[n-1]),'rank_decimal_digits':len(str(counts[n-1])), 'compressed_factors':len(histogram),'interval_bound':bound,'invariant_factors':[{'factor':str(d),'multiplicity':str(m)} for d,m in histogram], 'torsion_order_prime_factorization':torsion})
    print(json.dumps({'status':'PASS','classification':'Inherited PADDING_FROBENIUS_PROOF.md section 5; compressed algorithmic realization only','count_prefix_sha256':hashlib.sha256('\n'.join(map(str,counts)).encode()).hexdigest(),'max_n':args.max_n,'verification':'Direct expanded valuation-sort through n=10; independent tournament-matrix Smith controls through n=5; all atlas rank/divisibility/valuation/bound checks','atlas':rows},indent=2))


if __name__=='__main__':
    main()
