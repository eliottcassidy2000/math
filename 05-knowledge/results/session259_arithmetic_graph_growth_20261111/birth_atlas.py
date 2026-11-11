#!/usr/bin/env python3
"""Exact numeric-birth atlas, one node at a time. No prime heuristic."""
import csv
import sys
from math import isqrt
from pathlib import Path


def factor(n):
    result = {}
    p=2
    while p*p<=n:
        while n%p==0:
            result[p]=result.get(p,0)+1
            n//=p
        p+=1
    if n>1:
        result[n]=result.get(n,0)+1
    return result


def weak_birth(n):
    f=factor(n)
    return len(f)-1 if n>1 and n%2 and all(e==1 for e in f.values()) else None


def extra_birth(n):
    hits=[]
    for p,e in factor(n).items():
        if e>=2:
            m=n//(p*p)
            f=factor(m)
            if m%2 and all(e==1 for e in f.values()):
                hits.append(len(f))
    if len(hits)>1:
        raise RuntimeError('nonunique marked square')
    return hits[0] if hits else None


def rows(bound):
    for n in range(1,bound+1):
        f=factor(n)
        tau=1
        for e in f.values(): tau*=e+1
        square=isqrt(n)**2==n
        w=weak_birth(n)
        e=extra_birth(n)
        wd=weak_birth(n//2) if n%2==0 else None
        ed=extra_birth(n//2) if n%2==0 else None
        events=[x for x in [w,e,wd,ed] if x is not None]
        if len(events)>1:
            raise RuntimeError('simultaneous reduced barcode events')
        yield {'n':n,'sum_strict':(n-1)//2,'sum_weak':n//2,
               'product_strict':(tau-2-int(square))//2 if n>=2 else 0,
               'product_weak':(tau-2+int(square))//2 if n>=2 else 0,
               'double_diagonal':int(n%2==0),'square_diagonal':int(square and n>=4),
               'weak_birth_degree':w,'weak_death_degree':wd,
               'strict_extra_birth_degree':e,'strict_extra_death_degree':ed,
               'ordinary_essential_H0_birth':int(n==2)}


def main():
    allrows=list(rows(128))
    p=Path(__file__).with_name('BIRTH_ATLAS.csv')
    with p.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=allrows[0])
        writer.writeheader();writer.writerows(allrows)
    for n in range(1,20001):
        events=[weak_birth(n),extra_birth(n)]
        if n%2==0:
            events += [weak_birth(n//2),extra_birth(n//2)]
        if sum(x is not None for x in events)>1:
            raise RuntimeError(n)
        if n%16==0 and any(x is not None for x in events):
            raise RuntimeError('v2>=4 event')
    print('FINITE-EXACT PASS: atlas n1..128; unique-event and v2>=4 silence through20000.')

if __name__=='__main__': main()
