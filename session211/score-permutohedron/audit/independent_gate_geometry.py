#!/usr/bin/env python3
"""Independent rational-coordinate checks of gate interior and full affine slice."""
from fractions import Fraction
from itertools import permutations
from math import comb
import json

def check(x,msg):
    if not x:raise RuntimeError(msg)
def slacks(x):
    n=len(x)
    return [sum(x[i] for i in range(n) if mask>>i&1)-comb(mask.bit_count(),2) for mask in range(1,(1<<n)-1)]
def embed(y,t):
    n=len(y);return tuple(a+1 for a in y)+(t,n+1-t)

interior=0;slices=0
for n in range(1,6):
    for y in permutations(range(n)):
        for t in (Fraction(0),Fraction(1,3),Fraction(1,2),Fraction(1)):
            x=embed(y,t)
            check(sum(x)==comb(n+2,2),'total')
            check(min(slacks(x))==t,'exact minimum slack')
            interior+=1
for n in range(2,7):
    v=(-1,)+tuple(range(1,n-1))+(n,)
    for y in permutations(v):
        check(sum(y)==comb(n,2),'slice total')
        check(min(slacks(y))==-1,'core relaxed slack')
        check(min(slacks(embed(y,1)))==0,'ambient boundary')
        slices+=1
print(json.dumps({'status':'PASS','core_universe':'all permutohedron vertices n1..5, t=0,1/3,1/2,1','core_vertex_time_checks':interior,'slice_universe':'all relaxed permutohedron vertices n2..6','slice_vertex_checks':slices,'note':'unnormalized Landau facet slacks, not Euclidean distances'},indent=2))
