#!/usr/bin/env python3
"""Exact controls for the common-pair affine-lift protected-edge obstruction."""
from math import isqrt
import json
checks=0

def check(ok,message):
    global checks
    checks+=1
    if not ok: raise RuntimeError(message)

def centered(z):
    q=(z+24)//49
    return q,z-49*q

for n in range(2,1002):
    extras={c:49*(n+1)+(-1)**(n+1)*c for c in range(-24,25)}
    check(set(extras.values())==set(range(49*n+25,49*n+74)), 'extra partition')
    at24={c for c,z in extras.items() if z<=49*n+24}
    at72={c for c,z in extras.items() if z<=49*n+72}
    check(at24==set(), 'residue 24 must use no extra blocks')
    check(len(at72)==48 and all(c in at72 for c in range(-23,24)),
          'residue 72 interior blocks must use S1')
    excluded=set(extras)-at72
    check(len(excluded)==1 and abs(next(iter(excluded)))==24, 'only extreme block excluded')
    # Same two distinct endpoints would retain the unordered pair of endpoint parities.
    parity0=sorted([1,n%2])
    parity1=sorted([1,(n+1)%2])
    check(parity0!=parity1, 'consecutive compatible paths cannot share both endpoints')

check(centered(482)==(10,-8), 'first seed centered decomposition')
check(centered(3362)==(69,-19), 'middle seed centered decomposition')
check(centered(6242)==(127,19), 'last seed centered decomposition')
check(482+3362==62**2 and 62%7!=0, 'noninternal seed edge')
check(3362+6242==98**2 and 69+127==14**2, 'second edge descent')
check(centered(69)==(1,20) and centered(127)==(3,-20), 'second descent')
check(1+3==2**2, 'terminal descent')
for z in range(25,100001):
    r,c=centered(z)
    check(r>=1 and -24<=c<=24 and z==49*r+c, 'centered decomposition range')

print(json.dumps({'status':'passed','checks':checks,'n_range':[2,1001],
 'seed_edge':[482,3362],'centered_quotients':[10,69],'centered_residues':[-8,-19],
 'scope':'one common parity-compatible pair, all intact affine blocks, residues 24 and 72'},
 sort_keys=True,indent=2))
