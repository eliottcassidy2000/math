#!/usr/bin/env python3
"""REFUTES HYP-2864 (bounded-D witness for all covering 13-sets) and mac-mini HYP-2876
('every 13-set has witness D<=41'). The family {1..11,13, lcm(2..X)} is covering with
minimal witness denominator = nextprime(X) -> UNBOUNDED. Loneliness holds but only via
the large-D / EQUIDISTRIBUTION witness. kind-mendel-2026-06-22-S7."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def lcmall(xs):
    r=1
    for x in xs: r=r*x//gcd(r,x)
    return r
def min_witness_D(S, Dmax):
    for D in range(2,Dmax+1):
        for a in range(1,D//2+1):
            if gcd(a,D)==1 and all(nrm(F(s*a,D))>=F(1,14) for s in S): return D,a
    return None
core=[1,2,3,4,5,6,7,8,9,10,11,13]
print("S_X = {1..11,13, lcm(2..X)}: primitive covering 13-set, min witness D = nextprime(X):")
for X in [20,30,45,60,80,100]:
    v=lcmall(range(2,X+1))
    while v in core or gall(core+[v])!=1: v+=lcmall(range(2,X+1))
    S=sorted(core+[v]); w=min_witness_D(S, Dmax=X+80)
    print(f"  X={X:3d}: covering={is_covering(S)}, primitive={gall(S)==1}, min witness D={w[0]} (>X: {w[0]>X})")
print("=> witness D UNBOUNDED over covering 13-sets; NO finite certificate basis.")
print("=> bounded-D (HYP-2864) and 'D<=41 for all' (HYP-2876) are FALSE; these sets are lonely")
print("   ONLY via large-D witnesses = equidistribution (the irreducible analytic crux).")
