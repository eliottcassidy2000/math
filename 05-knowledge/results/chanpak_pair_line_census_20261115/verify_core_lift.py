#!/usr/bin/env python3
"""Read a generated exact order-through-nine catalogue on stdin."""
import sys
from collections import defaultdict
from functools import lru_cache

rows = []
for line in sys.stdin:
    n,h,p,code,sc=line.split()
    n,h,p=map(int,(n,h,p))
    rows.append((n,h,p+2*h,sc!="-" and n>1))
assert {n for n,h,q,s in rows} == set(range(1,10))
atoms=defaultdict(set)
for n,h,q,s in rows:
    if s:
        atoms[h].add((n,q))
@lru_cache(None)
def core(h):
    if h == 1:
        return {(0,1)}
    return {(n0+n1,q0*q1)
            for h0,values in atoms.items() if h%h0==0
            for n0,q0 in values for n1,q1 in core(h//h0)
            if n0+n1<=9}
exact=defaultdict(set)
for n,h,q,_ in rows:
    exact[h].add((n,q))
for h,vals in exact.items():
    lifts={(n, q*2**(n-s)) for s,q in core(h) for n in range(max(1,s),10)}
    assert lifts==vals,(h,lifts^vals)
assert core(1)=={(0,1)}
assert core(3)=={(3,12)}
assert core(5)=={(4,32)}
print("Exact literal lift comparisons:",len(rows),"across",len(exact),"Hamiltonian counts")
assert core(9)=={(5,84),(6,144)}
print("Primitive records for H=9:", sorted(core(9)))
