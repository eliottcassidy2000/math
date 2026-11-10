#!/usr/bin/env python3
"""No imports from researched code: all8 cyclic difference sets, all9! orders."""
from itertools import product, permutations
from collections import Counter
import json

def require(x,label):
    if not x: raise RuntimeError(label)

def paths(vertices, arcs):
    return sum(all((a,b) in arcs for a,b in zip(order,order[1:])) for order in permutations(vertices))

cases=[]
for signs in product([False,True],repeat=3):
    diffs={k if s else 7-k for k,s in zip([1,2,3],signs)}
    arcs={(a,b) for a in range(7) for b in range(7) if (b-a)%7 in diffs}
    h=paths(range(7),arcs)
    inflated=arcs | {(i,7) for i in range(7)} | {(7,8)} | {(8,i) for i in range(7)}
    full=paths(range(9),inflated)
    two=0
    for mask in range(1,127):
        left=[i for i in range(7) if mask>>i&1]; right=[i for i in range(7) if not mask>>i&1]
        two+=paths(left,arcs)*paths(right,arcs)
    require(full==3*h+two,'three runs or split cover')
    require(two==1638,'two cover value')
    require(full in [2163,2205],'inflation range')
    cases.append({'differences':sorted(diffs),'internal_H':h,'two_cover':two,'inflated_H':full})
require(Counter(c['inflated_H'] for c in cases)==Counter({2163:6,2205:2}),'multiplicities')
require(175**2+14==30639>2163,'multi orbit branch')
print(json.dumps({'status':'PASS','method':'difference-set construction; exhaustive permutations of all9 vertices for each of8 cases','cases':cases},indent=2))
