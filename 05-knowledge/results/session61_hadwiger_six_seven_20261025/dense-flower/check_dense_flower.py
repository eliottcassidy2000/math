#!/usr/bin/env python3
"""Recheck inherited norm-seven 21-point/57-edge density hostile exactly."""
from itertools import combinations
from collections import Counter
import json
def need(c,s):
 if not c:raise RuntimeError(s)
def mul(z,w):
 a,b=z;c,d=w
 return a*c-b*d,a*d+b*c-b*d
def add(z,w):return z[0]+w[0],z[1]+w[1]
def sub(z,w):return z[0]-w[0],z[1]-w[1]
def norm(z):return z[0]**2-z[0]*z[1]+z[1]**2
A=(3,1);B=(2,-1)
flower=[(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
P=sorted({add(mul(A,p),mul(B,q)) for p in [(0,0),(1,0),(1,1)] for q in flower})
need(norm(A)==norm(B)==7,'normseven multipliers')
need(len(P)==21,'twentyone distinct points')
E=[(i,j) for i,j in combinations(range(len(P)),2) if norm(sub(P[i],P[j]))==7]
need(len(E)==57,'fiftyseven exact edges')
colors=[(a+b)%3 for a,b in P]
need(all(colors[i]!=colors[j] for i,j in E),'proper threecoloring')
S=set(E)
tri=[(a,b,c) for a,b,c in combinations(range(21),3) if {(a,b),(a,c),(b,c)}<=S]
need(bool(tri),'triangle lowerbound')
print(json.dumps({'status':'PASS','points':P,'edges':E,'edge_count':len(E),'triangles':len(tri),'colors':colors,'chromatic_number':3,'degree_histogram':dict(sorted(Counter(sum(i in e for e in E) for i in range(21)).items())),'scope':'57-edge construction verified; global extremality is inherited literature, not computed here'},indent=2))
