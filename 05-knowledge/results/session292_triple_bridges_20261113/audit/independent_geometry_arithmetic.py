#!/usr/bin/env python3
"""Independent reviewer controls: edge-set quotient and four hidden primes."""
from itertools import combinations
from math import isqrt,lcm,prod
import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('square',ROOT/'square_sum/verify.py')
sq=importlib.util.module_from_spec(spec);spec.loader.exec_module(sq)
def check(p,msg):
    if not p:raise RuntimeError(msg)
tests=0
for A,B,C in combinations(range(1,16),3):
    for N in range(1,22):
        edge_sets=set()
        for x in range(1,N+1):
            v=[x,A-x,x+B-A,A+C-B-x,x+B-C,C-x]
            if len(set(v))!=6 or min(v)<1 or max(v)>N:continue
            edge_sets.add(tuple(sorted(tuple(sorted((v[j],v[(j+1)%6]))) for j in range(6))))
        check(len(edge_sets)==sq.hex_count(A,B,C,N),str((A,B,C,N)))
        tests+=1
for f in (sq.direct,sq.compiled):
    check(f([4],2,True)==[],'loop boundary')
    check(f([4,4],3,True)==[],'backtrack boundary')
P=[149,347,2111,15583]
check(all(all(p%d for d in range(2,isqrt(p)+1)) for p in P),'prime factors')
check(prod(P)==99*2**34+23,'factorization')
def order(p):
    x=2%p;d=1
    while x!=1:x=x*2%p;d+=1
    return d
O=[order(p) for p in P];L=lcm(61209720,*O)
check(O==[148,346,1055,1113] and L==30670763342539320,'orders and period')
small=[p for p in range(2,100) if all(p%d for d in range(2,isqrt(p)+1))]
for t in range(10):
    k=34+L*t
    check(all((99*pow(2,k,p)+23)%p==0 for p in P),'forced factors')
    check(all((99*pow(2,k,p)+23)%p!=0 for p in small),'small-prime survival')
print(json.dumps({'status':'PASS','independent_edge_set_checks':tests,'hidden_primes':P,
                  'orders':O,'surviving_four_factor_cylinder':[34,L]},indent=2))
