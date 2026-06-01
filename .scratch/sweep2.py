from lrc_coloring_core import *
import collections, itertools, random

# (A) Hunt for chi != omega over many random sets (m=4,5,6) and ALL cells.
random.seed(1)
found = []
tested = 0
for trial in range(400):
    m = random.choice([4,5,6])
    speeds = tuple(sorted(random.sample(range(1, 25), m)))
    n = m+1
    crit = critical_times(speeds, n)
    mids = cell_midpoints(crit)
    for t in mids:
        tested += 1
        k,_ = chromatic_number(speeds, n, t)
        w = clique_number_arc(speeds, n, t)
        if k != w:
            found.append((speeds, t, k, w))
            break
print("(A) random search: tested", tested, "cells, chi!=omega found in", len(found), "sets")
for f in found[:10]:
    print("   ", f)

# (B) Is max_t chi exactly n, achieved away from t=0? Where does chi=n occur?
print()
print("(B) cells achieving chi=n, are they only near collisions of ALL points?")
for speeds in [(1,2,3,4),(1,3,4,7),(2,3,5,7),(1,2,3,4,5),(1,3,4,5,9)]:
    n=len(speeds)+1
    crit=critical_times(speeds,n); mids=cell_midpoints(crit)
    full=[]
    for t in mids:
        k,_=chromatic_number(speeds,n,t)
        if k==n:
            # are all n points within 1/n pairwise? = K_n
            N,edges,pts=danger_graph(speeds,n,t)
            iskn = len(edges)==n*(n-1)//2
            full.append((t,iskn))
    print("  ",speeds,"chi=n cells:",len(full),"all K_n?",all(x[1] for x in full), [str(t) for t,_ in full])
