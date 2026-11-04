#!/usr/bin/env python3
"""Independent modular-power and BFS oracle; checks are active under -O."""
import importlib.util
from itertools import combinations
from math import isqrt
from pathlib import Path
from fractions import Fraction
import json
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('subject',ROOT/'check_odd_order_support.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def require(c,msg):
    if not c:raise RuntimeError(msg)
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def admissible(p,r,q):
    odd=r-1
    while odd%2==0:odd//=2
    return pow(p,odd*(2 if p==q else 1),r)==1
def geom(p,a,r):
    return (a+1)%r if p%r==1 else ((pow(p,a+1,r)-1)*pow(p-1,-1,r))%r
def reach(S,E,p):
    seen={p}; todo=[p]
    while todo:
        v=todo.pop()
        for a,b in E:
            if a==v and b not in seen:seen.add(b);todo.append(b)
    return seen
P=tuple(p for p in range(3,44,2) if prime(p));count=0;paircount=0
for size in range(2,6):
    for S in combinations(P,size):
        for q in S:
            if q%4!=1:continue
            E={(p,r) for p in S for r in S if p!=r and admissible(p,r,q)}
            require(E==m.graph(S,q),'modular-power graph mismatch')
            a=m.realizing_exponents(S,q,E)
            require(a[q]%4==1 and all(v>0 and v%2==0 for p,v in a.items() if p!=q),'Euler parity')
            for p in S:
                for r in S:
                    if p==r:continue
                    require((geom(p,a[p],r)==0)==((p,r) in E),'geometric oracle realization')
                    paircount+=1
            R={p:reach(S,E,p) for p in S}
            require(R==m.closure(S,E),'BFS closure mismatch')
            G={(p,r) for p in S for r in S if p!=r and admissible(p,r,None)}
            RG={p:reach(S,G,p) for p in S}
            components={frozenset(r for r in S if r in RG[p] and p in RG[r]) for p in S}
            sinks=[C for C in components if all(not(a in C and b not in C) for a,b in G)]
            require(all(q in RG[p] for p in S)==(len(sinks)==1 and q in sinks[0]),'unique-sink criterion')
            # H strong connectivity implies the first-pass path condition.
            require(not all(len(R[p])==len(S) for p in S) or all(q in RG[p] for p in S),'first-pass necessity')
            count+=1
S=(3,5,23,31,43);E=m.graph(S,5);a=m.realizing_exponents(S,5,E)
require(a=={3:10,5:461,23:20,31:1154,43:2},'sample exponents')
require(all(reach(S,E,p)=={3,23,43} for p in [3,23,43]),'sample sink')
require(Fraction(3,2)*Fraction(5,4)*Fraction(23,22)*Fraction(31,30)*Fraction(43,42)==Fraction(30659,14784),'sample abundancy')
require(all(geom(p,a[p],r)==sum(pow(p,i,r) for i in range(a[p]+1))%r for p in S for r in S if p!=r),'sample direct summation')
print(json.dumps({'status':'PASS','independent_support_cases':count,'independent_pair_checks':paircount,'prime_universe':P,'sample_verified':True},sort_keys=True))
