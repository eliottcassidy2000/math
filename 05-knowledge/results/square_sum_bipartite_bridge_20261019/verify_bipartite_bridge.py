#!/usr/bin/env python3
from itertools import combinations
from math import isqrt
from collections import Counter
import json
checks=0

def ck(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def sq(n):return n>=0 and isqrt(n)**2==n

# Direct exact verification of the root's bounded witness.
A=(1,193,673);B=(3,168,483)
roots=[[isqrt(a+b) for b in B] for a in A]
ck(roots==[[2,13,22],[14,19,26],[26,29,34]],'witness roots')
for a in A:
    for b in B:ck(sq(a+b),'witness cross edge')
for side in (A,B):
    for a,b in combinations(side,2):ck(not sq(a+b),'witness within-part nonedge')
ck(len(set(A+B))==6,'witness injection')

# Affine cleanup, independently testing cross edges iff and all internal pairs.
for a in range(1,301):
    for b in range(1,301):
        ck(sq((25*a+1)+(25*b-1))==sq(a+b),'cross-square equivalence')
        ck(25*a+1!=25*b-1,'disjoint affine parts')
for a,b in combinations(range(1,301),2):
    ck(not sq(25*a+1+25*b+1),'left independent')
    ck(not sq(25*a-1+25*b-1),'right independent')

# Pair-intersection scan, separate from a common-neighbor-triple scanner.
# For a pair {u,v}, set S=N(u) intersect N(v). A third vertex w completes
# a K_3,3 subgraph exactly when it has at least three neighbors in S.
def scan(N):
    neighbors=[set() for _ in range(N+1)]
    for root in range(2,isqrt(2*N)+1):
        total=root*root
        for a in range(max(1,total-N), min(N,(total-1)//2)+1):
            b=total-a
            neighbors[a].add(b);neighbors[b].add(a)
    examined=eligible=0
    for u in range(1,N+1):
        for v in range(u+1,N+1):
            examined+=1
            S=neighbors[u]&neighbors[v]
            if len(S)<3:continue
            eligible+=1
            hit=Counter(w for x in S for w in neighbors[x] if w not in (u,v))
            for w,count in sorted(hit.items()):
                if count>=3:
                    right=tuple(sorted(S&neighbors[w])[:3])
                    left=(u,v,w)
                    ck(len(set(left+right))==6,'scan disjoint parts')
                    ck(all(sq(a+b) for a in left for b in right),'scan witness actual')
                    return {'N':N,'found':True,'left':left,'right':right,
                            'pairs_examined':examined,'eligible_pairs':eligible}
    return {'N':N,'found':False,'pairs_examined':examined,'eligible_pairs':eligible}

lower=scan(672);upper=scan(673)
ck(not lower['found'],'no K33 subgraph through 672')
ck(upper['found'],'K33 subgraph at 673')
print(json.dumps({'status':'PASS','checks':checks,'minimum_K33_height':673,
                  'explicit_induced_witness':{'A':A,'B':B,'roots':roots},
                  'independent_pair_intersection_scans':[lower,upper],
                  'scope':'Finite exact threshold and affine controls; arbitrary m existence uses the cited sumset theorem plus the proved affine cleanup.'},indent=2,sort_keys=True))
