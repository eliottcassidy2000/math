#!/usr/bin/env python3
"""Produce unlabeled 6-vertex representatives by lexicographic canonicalization."""
from itertools import permutations
from collections import Counter
from certificate import graph,strong,direct
n=6
pairs=[(i,j) for i in range(n) for j in range(i+1,n)]

def code(a,p):
    return sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs))

classes={}
for mask in range(1<<10):
    a=graph(n,mask)
    if not strong(a): continue
    c=min(code(a,p) for p in permutations(range(n)))
    h=direct(a)
    if c in classes:
        if classes[c][0]!=h: raise RuntimeError('non-invariant path count')
        classes[c][1]+=1
    else: classes[c]=[h,1]
print('canonical representatives: code, H, fixed-path labelings')
for c,(h,w) in sorted(classes.items()): print(c,h,w)
print('classes=',len(classes),'frequency=',dict(sorted(Counter(v[0] for v in classes.values()).items())))
