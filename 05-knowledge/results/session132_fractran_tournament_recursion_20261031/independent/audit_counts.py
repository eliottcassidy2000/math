#!/usr/bin/env python3
"""Independent conjugacy-class Burnside via explicit signed pair permutations.
Uses every integer partition, not the odd-partition/gcd formulas under audit.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import factorial
from pathlib import Path
import json,sys
from audit_fractions import require
ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent


def partitions(n,least=1):
    if not n:yield ()
    for k in range(least,n+1):
        for p in partitions(n-k,k):yield (k,)+p


@lru_cache(None)
def fixes(parts):
    n=sum(parts);p=list(range(n));offset=0
    for length in parts:
        for i in range(length):p[offset+i]=offset+(i+1)%length
        offset+=length
    edges=set(combinations(range(n),2));plain=twisted=True;orbits=0
    while edges:
        first=next(iter(edges));e=first;swaps=length=0
        while e in edges:
            edges.remove(e);u,v=e;a,b=p[u],p[v]
            swaps^=a>b;length+=1;e=tuple(sorted((a,b)))
        require(e==first,'Bad pair orbit')
        plain &= not swaps
        twisted &= not(swaps^(length%2))
        orbits+=1
    return (2**orbits if plain else 0,2**orbits if twisted else 0)


stored=json.loads((ROOT/'counting/normal.json').read_text())['counts']
results=[];classes=0
for n in range(26):
    a=s=0
    for parts in partitions(n):
        z=1
        for d,m in Counter(parts).items():z*=d**m*factorial(m)
        multiplicity=factorial(n)//z
        ordinary,twisted=fixes(parts)
        require(bool(ordinary)==all(d%2 for d in parts),'Ordinary support mismatch')
        require(bool(twisted)==(parts.count(1)<=1 and all(d==1 or d%4==2 for d in parts)),'Twisted support mismatch')
        a+=multiplicity*ordinary;s+=multiplicity*twisted;classes+=1
    require(a%factorial(n)==s%factorial(n)==0,'Nonintegral Burnside result')
    a//=factorial(n);s//=factorial(n)
    require((a+s)%2==0,'Converse quotient parity')
    expected=stored[n]
    require((a,s,(a+s)//2)==tuple(expected[k] for k in ['A','S','M']),'Published count mismatch '+str(n))
    # Pair of tournaments and one/two vertex colors, computed independently.
    m=n//2;paired=0
    for parts in partitions(m):
        z=1
        for d,k in Counter(parts).items():z*=d**k*factorial(k)
        fixed,_=fixes(parts)
        paired+=(factorial(m)//z)*fixed**2*(2**(n%2))**len(parts)
    require(paired%factorial(m)==0 and paired//factorial(m)==s,'Folding cardinality mismatch')
    results.append(dict(n=n,A=a,S=s,M=(a+s)//2))
print(json.dumps(dict(status='PASS',checked_partition_classes=classes,counts=results),indent=2))
