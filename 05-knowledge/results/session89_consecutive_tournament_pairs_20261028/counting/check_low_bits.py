#!/usr/bin/env python3
"""Bounded exact checks of the all-order positive-channel low-bit theorem.

No committed JSON receipt or external baseline is needed. This is a checksum,
not a replacement for complete counting or a performance benchmark.
"""
import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys

from channels import fused_prefix
from verify import partitions, partition_data, partition_count

CHANNELS=((1,0),(1,3),(2,0),(2,1),(3,2),(4,1))


def check(condition,message):
    if not condition:
        raise AssertionError(message)


def valuation2(value):
    if type(value) is not int or value <= 0:
        raise ValueError('positive integer required')
    return (value & -value).bit_length()-1


def totient(n):
    result=n
    p=2
    remaining=n
    while p*p <= remaining:
        if remaining%p == 0:
            result-=result//p
            while remaining%p == 0:
                remaining//=p
        p+=1
    if remaining>1:
        result-=result//remaining
    return result


def check_count(n,r,s,value):
    if n == 0:
        check(value==1,'empty count')
        return 'empty'
    if n == 2:
        check(value==1<<(r+2*s-1),'order-two count')
        return 'order_two'
    if n%2:
        base=r*(n-1)//2+s
        check(value%(1<<(base+1)) == 1<<base,
              f'odd valuation failed: {(n,r,s)}')
        return 'odd_exact'
    base=r*n//2+2*s
    parity=(totient(n)//2)%2
    check(value%(1<<(base+1)) == parity<<base,
          f'even normalized low bit failed: {(n,r,s)}')
    check((valuation2(value)==base) == bool(parity),
          f'even equality condition failed: {(n,r,s)}')
    return 'even_equal' if parity else 'even_strict'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound',type=int,default=80)
    args=parser.parse_args()
    if args.bound<2:
        raise ValueError('bound must be at least 2')
    values,stats=fused_prefix(args.bound,CHANNELS)
    outcomes={}
    for (r,s),prefix in zip(CHANNELS,values):
        for n,value in enumerate(prefix):
            outcome=check_count(n,r,s,value)
            outcomes[outcome]=outcomes.get(outcome,0)+1
    partition_bound=min(26,args.bound)
    term_checks=0
    for n in range(partition_bound+1):
        for r,s in CHANNELS:
            check(partition_count(n,r,s)==values[CHANNELS.index((r,s))][n],
                  'independent rational count mismatch')
            if n in (0,2):
                continue
            base=r*(n-1)//2+s if n%2 else r*n//2+2*s
            for parts in partitions(n):
                if not all(k%2 for k in parts):
                    continue
                e,k,z=partition_data(parts)
                v=r*e+s*k-valuation2(z)
                if n%2:
                    check(v==base if k==1 else v>=base+1,
                          'odd summand valuation bound')
                else:
                    leading=(k==2 and parts[0]!=parts[1] and gcd(*parts)==1)
                    check(v==base if leading else v>=base+1,
                          'even summand valuation bound')
                term_checks+=1
    for n in range(4,args.bound+1,2):
        pairs=sum(1 for a in range(1,n//2+1,2)
                  if a<n-a and gcd(a,n-a)==1)
        check(pairs==totient(n)//2,'coprime odd-pair count')
    # The tempting bound v2(k!)<=k-2 is false at even k=4.
    check(valuation2(24)==3>4-2,'factorial-bound hostile control')
    here=Path(__file__).resolve().parent
    result={'status':'PASS','optimization_level':sys.flags.optimize,
            'channels':CHANNELS,'exact_counter_bound':args.bound,
            'independent_partition_bound':partition_bound,
            'termwise_valuation_checks':term_checks,'outcomes':outcomes,
            'source_sha256':{name:sha256((here/name).read_bytes()).hexdigest()
                             for name in ('channels.py','verify.py','check_low_bits.py')}}
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
