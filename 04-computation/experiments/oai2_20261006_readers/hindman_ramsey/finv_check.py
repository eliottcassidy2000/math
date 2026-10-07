#!/usr/bin/env python3
"""Independent brute-force check of a fully gap-determined witness E on [t]^3:
(1) triangle-free: no x<y<z in [t]^3 with all three gaps in E; (2) every binary subgrid has a gap in E."""
import sys, itertools
import numpy as np
fn=sys.argv[1]; t=int(sys.argv[2]); n=3
E=set(tuple(map(int,l.split())) for l in open(fn) if l.strip())
pts=sorted(itertools.product(range(t),repeat=n))
# (1) triangle check via sum-freeness over realizable triples (explicit points)
bad=0
for x,y,z in itertools.combinations(pts,3):
    if tuple(y[i]-x[i] for i in range(3)) in E and tuple(z[i]-y[i] for i in range(3)) in E and tuple(z[i]-x[i] for i in range(3)) in E:
        bad+=1; break
print("triangle found" if bad else "triangle-free: OK")
# (2) all binary subgrids, explicit: a0<a1, b-pairs, c-pairs
inE=lambda d: d in E
pairs=list(itertools.combinations(range(t),2))
missed=0; total=0
# precompute row-level (height 2) subgrids as 4 points
rows=[]
for (b0,b1) in pairs:
    for (c0,c1) in pairs:
        for (c2,c3) in pairs:
            rows.append(((b0,c0),(b0,c1),(b1,c2),(b1,c3)))
def rowhit(X):
    for i in range(4):
        for j in range(i+1,4):
            if (0,X[j][0]-X[i][0],X[j][1]-X[i][1]) in E: return True
    return False
rh=[rowhit(X) for X in rows]
for (a0,a1) in pairs:
    h=a1-a0
    for i,X in enumerate(rows):
        for j,Y in enumerate(rows):
            total+=1
            if rh[i] or rh[j]: continue
            hit=False
            for x in X:
                for y in Y:
                    if (h,y[0]-x[0],y[1]-x[1]) in E: hit=True;break
                if hit: break
            if not hit:
                missed+=1
print(f"binary subgrids: {total}, missed (independent): {missed}")
