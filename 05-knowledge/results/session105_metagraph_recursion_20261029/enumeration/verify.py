#!/usr/bin/env python3
"""Deterministic checks, all active under -O. Timings in benchmark.py."""
import json, hashlib
import sys
sys.set_int_max_str_digits(0)
from pathlib import Path
from baseline98 import direct, fixed_degree, partitions
from baseline87 import degree_bucket_dp
from integer_response import count, count_packed, exact_shift
from math import gcd


def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def main():
    prefix,_=degree_bucket_dp(200)
    values=[]
    for n in range(31):
        expected=direct(n)
        values.append(int(expected))
        require(expected==prefix[n],f'direct {n}')
        for tail in (1,3):
            for cache in (False,True):
                require(count(n,tail,cache)[0]==expected,f'tuple {n},{tail},{cache}')
        for odd in (False,True):
            require(count_packed(n,odd)[0]==expected,f'packed {n},{odd}')
        for mu in partitions(n):
            for p in range(3,n+1,2):
                past=[a for a in mu if a>p]
                for j in range(3,p+1,2):
                    deficit=sum(j-gcd(a,j) for a in past)
                    require(0<=deficit<=n,'packed carry bound')
    require(values[:11]==[1,1,1,2,4,12,56,456,6880,191536,9733056],'known prefix')
    samples=[]
    for n in (40,60,100,150,200):
        reference,oldstates=fixed_degree(n)
        a,stats=count_packed(n)
        full,fullstats=count_packed(n,False)
        require(a==full==reference==prefix[n],f'large {n}')
        require(stats['stages']==oldstates[:-1],f'exact state equivalence {n}')
        samples.append(dict(n=n,value=str(a),stats=stats,full_factorial_stats=fullstats))
    require(exact_shift(24,-3)==3,'negative shift')
    try:exact_shift(25,-3)
    except ArithmeticError:pass
    else:raise RuntimeError('nonexact negative shift not rejected')
    for n in (-1,True,1.5):
        try:count_packed(n)
        except ValueError:pass
        else:raise RuntimeError('invalid n not rejected')
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).parent.glob('*.py'))}
    print(json.dumps(dict(status='PASS',direct_through=30,independent_prefix_through=200,values=values,samples=samples,sha256=files),indent=2))
if __name__=='__main__':main()
