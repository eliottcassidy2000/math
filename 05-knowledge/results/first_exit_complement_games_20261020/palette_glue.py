#!/usr/bin/env python3
"""Compose the incoming sixteen-point palette gate at its actual output.
The source geometry and port triples are attributed in PALETTE-GLUE.md.
"""
import sys
sys.dont_write_bytecode=True
from geometry_games import add,neg,mul,norm,colouring,need
from fractions import Fraction as F
from itertools import combinations
import json

zero=(F(0),)*4;one=(F(1),F(0),F(0),F(0))
z=(F(1,2),F(1,2),F(0),F(0));rho=(F(5,6),F(0),F(1,6),F(0))
points=[zero,one,z,add(one,z),rho,mul(rho,z),mul(rho,add(one,z))]
groups=[[(1,3,2),(2,4,0),(5,6,4)],[(1,4,0),(2,6,3),(3,5,6)],[(1,5,0),(2,3,1),(4,6,5)]]
for group in groups:
    for i,j,k in group:points.append(add(add(points[i],points[j]),neg(points[k])))
base_edges=[(i,j) for i,j in combinations(range(16),2) if norm(add(points[i],neg(points[j])))==one]
need(len(base_edges)==31,'incoming gate contact count')
selected=None
for t in range(1,4*15*15+2):
    rot=(F(1-3*t*t,1+3*t*t),F(2*t,1+3*t*t),F(0),F(0))
    need(norm(rot)==one,'rotation norm')
    copied=[mul(rot,p) for p in points]
    if any(norm(add(a,neg(b))) in (zero,one) for a in points[1:] for b in copied[1:]):continue
    selected=t;break
need(selected is not None,'finite bad-rotation bound')
joined=points+copied[1:]
edges=[(i,j) for i,j in combinations(range(31),2) if norm(add(joined[i],neg(joined[j])))==one]
need(len(set(joined))==31 and len(edges)==62,'exact induced one-vertex union')
adj=[set() for _ in joined]
for i,j in edges:adj[i].add(j);adj[j].add(i)

def output_counts(A,B):
    cs={}
    for row,c in enumerate(A):
        for v in range(7+3*row,10+3*row):cs[v]=c
    for row,c in enumerate(B):
        for v in range(22+3*row,25+3*row):cs[v]=c
    need(all(cs[i]!=cs[j] for i,j in edges if i in cs and j in cs),'proper full boundary')
    def count():
        if len(cs)==31:return 1
        options=[(set(range(6))-{cs[j] for j in adj[i] if j in cs},i) for i in range(31) if i not in cs]
        allowed,i=min(options,key=lambda p:(len(p[0]),-len(adj[p[1]])))
        total=0
        for c in allowed:
            cs[i]=c;total+=count()
        cs.pop(i,None)
        return total
    out=[]
    for c in range(6):
        cs[0]=c;out.append(count());cs.pop(0)
    return out

rows=[]
for B in combinations(range(6),3):
    A=(0,1,2);got=output_counts(A,B)
    expected=[576 if c in A and c in B else 0 for c in range(6)]
    need(got==expected,'full-graph output count, independent of product formula')
    rows.append({'left_palette':A,'right_palette':B,'output_counts':got})
palettes=list(combinations(range(6),3))
for i in range(6):
    for j in range(6):
        need(sum(i in p and j in p for p in palettes)==(10 if i==j else 4),
             'six-dimensional Gram spectrum')
c4=colouring(31,edges,4);need(c4 is not None,'four-colour hostile')
print(json.dumps({'status':'FINITE-EXACT; conditional boundary intersection, not unconditional equality forcing',
 'rotation_parameter':selected,'rotation_coefficients':list(map(str,rot)),
 'points':[[str(x) for x in p] for p in joined],'edges':edges,
 'four_colouring':c4,'palette_rows':rows,'pairs_checked':465},indent=2,sort_keys=True))
