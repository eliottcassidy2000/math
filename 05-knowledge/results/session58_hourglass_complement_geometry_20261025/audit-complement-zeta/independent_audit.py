#!/usr/bin/env python3
"""Independent matrix/order and elementary-abelian checks; exceptions survive -O."""
from itertools import combinations, permutations, product
from collections import Counter
from math import factorial
from pathlib import Path
import json, sys

def req(condition, message):
    if not condition: raise RuntimeError(message)
if '--negative-exception' in sys.argv: req(False,'injected negative exception')
rows=[]
for n in range(1,6):
    pairs=list(combinations(range(n),2)); orders=list(permutations(range(n)))
    def matrix(bits,tournament):
        a=[[0]*n for _ in range(n)]
        for bit,(i,j) in zip(bits,pairs):
            a[i][j]=bit; a[j][i]=1-bit if tournament else bit
        return a
    def word(a,o): return tuple(a[o[i]][o[j]] for i,j in pairs)
    def canon(a): return min(word(a,o) for o in orders)
    def opposite(a): return [[a[j][i] for j in range(n)] for i in range(n)]
    def graph_pair(bits):
        return min(canon(matrix(bits,False)),canon(matrix(tuple(1-b for b in bits),False)))
    reps={}; graphs=set(); gpairs=set()
    for bits in product((0,1),repeat=len(pairs)):
        a=matrix(bits,True); reps.setdefault(canon(a),a)
        graphs.add(canon(matrix(bits,False))); gpairs.add(graph_pair(bits))
    cpairs={}; fixed=0
    for rep,a in reps.items():
        c=canon(opposite(a)); fixed+=rep==c; cpairs.setdefault(min(rep,c),a)
    profiles=[]
    for a in cpairs.values():
        h=Counter(graph_pair(word(a,o)) for o in orders)
        req(sum(h.values())==factorial(n),'histogram mass')
        req(n<2 or all(v%2==0 for v in h.values()),'histogram parity')
        req(h==Counter(graph_pair(word(opposite(a),o)) for o in orders),'converse invariance')
        profiles.append(tuple(sorted(h.items())))
    req(len(set(profiles))==len(cpairs),'finite profile separation')
    row=[n,len(graphs),len(gpairs),len(reps),fixed,len(cpairs),len(set(profiles))]
    rows.append(row)
req(rows==[[1,1,1,1,1,1,1],[2,2,1,1,1,1,1],[3,4,2,2,2,2,2],[4,11,6,4,2,3,3],[5,34,18,12,8,10,10]],'census')
# Exhaust every unordered p=3 factor pair directly in the additive plane.
p=3; V=list(product(range(p),repeat=2)); Z=(0,0)
def add(x,y): return tuple((a+b)%p for a,b in zip(x,y))
def neg(x): return tuple(-a%p for a in x)
def sub(x,y): return add(x,neg(y))
def shift(B,t): return frozenset(add(x,t) for x in B)
subsets=list(map(frozenset,combinations(V,3))); frames=[]; classes=set(); badchecks=0
for B,D in combinations(subsets,2):
    if len({add(b,d) for b in B for d in D})!=9: continue
    frames.append((B,D))
    orbit=[]
    for t in V:
        Bt,Dt=shift(B,t),shift(D,neg(t))
        orbit.append(tuple(sorted((tuple(sorted(Bt)),tuple(sorted(Dt))))))
        if Z in Bt:
            forbidden={tuple((2*z)%p for z in sub(d,b)) for b in Bt for d in Dt}
            req(forbidden==set(V),'normalized forbidden locus A/2'); badchecks+=1
    req(len(set(orbit))==9,'free unordered chart')
    classes.add(min(orbit))
    req(bool(B&D),'base dimension two disjoint count zero')
req(len(frames)==270 and len(classes)==30,'independent frame census')
# Affine offset in F_3^3, using normalized plane factors and all 27 translations.
ambient=list(product(range(3),repeat=3)); C=(0,0,1); A={add((x,y,0),C) for x,y in V}; halfA={tuple(2*z%3 for z in a) for a in A}
for orbit in classes:
    B0,D0=map(lambda L:frozenset(L),orbit)
    b=min(B0); B0=shift(B0,neg(b)); D0=shift(D0,b)
    B=frozenset((x,y,0) for x,y in B0); D=frozenset((x,y,1) for x,y in D0)
    req({add(b,d) for b in B for d in D}==A,'affine sum support')
    bad={t for t in ambient if shift(B,t)&shift(D,neg(t))}
    req(bad==halfA,'affine normalized forbidden A/2')
    req(len(bad)==9 and len(ambient)-len(bad)==18,'extension free chart count')
# Nontrivial involution must not be replaced by a fixed set or free action indiscriminately.
req([sum((x+m)%2==x for x in range(2)) for m in range(1,5)]==[0,2,0,2],'two-cycle trace')
req(3**3-3==24 and 3**3-9==18,'odd-degree deletion hostile')
req(3*(2**3-4)!=0,'p=2 disjoint-formula hostile with 0 in A')
result={'census_columns':['n','graphs','graph_pairs','tournaments','self_converse','converse_pairs','profiles'],'census':rows,'direct_p3_frame_count':len(frames),'gauge_classes':len(classes),'normalized_bad_locus_checks':badchecks,'affine_27_point_charts':len(classes),'negative_exception_flag':'--negative-exception','all_checks':'PASS'}
Path(__file__).with_name('independent_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
