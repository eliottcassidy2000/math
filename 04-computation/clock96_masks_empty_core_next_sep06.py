#!/usr/bin/env python3
"""Exact c96 affine padded-mask obstruction and an actual c90 boundary.

All labels are kept in Z/cZ. Standard-library integers only. The finite
universe contains all individual masks at orders8,16,24,32,96 and exactly
the pair/triple projections used by the written proof, not six-mask covers.
"""
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import gcd
import json

GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def mask(c,q,start,unit):
    k=(q+6)//7
    residues={(start+unit*j)%q for j in range(k)}
    return sum(1<<j for j in range(c) if j%q in residues)

def masks(c,q):
    atlas={}
    for unit in range(q):
        if gcd(unit,q)!=1:continue
        for start in range(q):
            m=mask(c,q,start,unit)
            atlas.setdefault(m,(start,unit))
    return atlas

def residues8(m):
    return {r for r in range(8) if any(m&(1<<j) for j in range(r,96,8))}

atlas={q:masks(96,q) for q in (8,16,24,32,96)}
expected={8:16,16:64,24:96,32:256,96:1536}
for q,rows in atlas.items():
    need(len(rows)==expected[q],'complete affine-mask count')
    for m in rows:
        need(m.bit_count()==(96//q)*((q+6)//7),'literal mask capacity')
        need(0<=m<1<<96,'shared label range')
# Only the pair inequalities used by the explicit trees.
pair_cases=0
for q,r,lower in [(96,8,2),(96,16,1),(32,24,1)]:
    minimum=97
    for a in atlas[q]:
        for b in atlas[r]:
            meet=(a&b).bit_count();minimum=min(minimum,meet)
            need(meet>=lower,'forced pair overlap')
            pair_cases+=1
    need(minimum==lower,'pair lower bound attained')
for q,multiple in [(16,6),(24,4)]:
    for e in atlas[8]:
        for a in atlas[q]:
            need((e&a).bit_count()%multiple==0,'intersection lattice quantum')
            pair_cases+=1
# All labelled triples E,D,A. Disjointness from E is kept literally.
shadow_cases=0
for e in atlas[8]:
    er=residues8(e);need(len(er)==2,'E shadow size')
    for d in atlas[16]:
        if e&d:continue
        dr=residues8(d)
        need(len(dr)==3 and not er&dr,'D shadow in six-class complement')
        for a in atlas[24]:
            if e&a:continue
            ar=residues8(a)
            need(len(ar)==4 and not er&ar,'A shadow in six-class complement')
            need(len(dr&ar)>=1,'seven slots in six shadow classes')
            need((d&a).bit_count()==2*len(dr&ar),'literal CRT multiplicity')
            need((d&a).bit_count()>=2,'coupled extra intersection')
            shadow_cases+=1
# Vertices preserve X,Y,A,B,D,E labels throughout.
# Each (i,j,w) is an edge with the proved lower intersection weight w.
A_trees={
 'ED':[(4,5,6),(0,5,2),(1,5,2),(2,5,0),(3,5,0)],
 'EA':[(2,5,4),(0,5,2),(1,5,2),(0,4,1),(2,3,0)],
 'EB':[(3,5,4),(0,5,2),(1,5,2),(0,4,1),(2,3,0)],
 'avoid':[(0,5,2),(1,5,2),(0,4,1),(4,2,2),(4,3,2)],
}
B_trees={
 'ED':[(4,5,6),(0,5,2),(1,2,1),(1,3,1),(1,5,0)],
 'EA':[(2,5,4),(0,5,2),(0,4,1),(1,2,1),(1,3,1)],
 'EB':[(3,5,4),(0,5,2),(0,4,1),(1,2,1),(1,3,1)],
 'avoid':[(0,5,2),(0,4,1),(4,2,2),(4,3,2),(1,2,1)],
}
case_bounds={}
for name,total,trees,bound in [('first',102,A_trees,93),('second',103,B_trees,95)]:
    row={}
    for case,tree in trees.items():
        need(len(tree)==5,'five tree edges')
        reached={0}
        while True:
            nxt=reached|{j for i,j,w in tree if i in reached}|{i for i,j,w in tree if j in reached}
            if nxt==reached:break
            reached=nxt
        need(len(reached)==6,'tree connectivity')
        # Labelwise Hunter inequality for every possible nonempty membership word.
        for membership in range(1,64):
            m=membership.bit_count()
            edges=sum(bool(membership&(1<<i)) and bool(membership&(1<<j)) for i,j,w in tree)
            need(edges<=m-1,'universal Hunter label inequality')
        row[case]=total-sum(w for i,j,w in tree)
        need(row[case]<=bound<96,'case excludes complete cover')
    case_bounds[name]=row

# Root's concurrently found c90 partition, with original labels and tail order.
c=90;gs=(2,2,2,3,3,6)
params=[(45,8,1),(45,23,1),(45,38,1),(30,3,1),(30,18,1),(15,0,1)]
cover=[mask(c,*p) for p in params]
need(all(m.bit_count()==g*((c//g+6)//7) for m,g in zip(cover,gs)),'c90 exact capacities')
need(sum(m.bit_count() for m in cover)==c,'c90 capacity saturation')
need(all(not(a&b) for a,b in combinations(cover,2)),'c90 disjoint masks')
need(reduce(int.__or__,cover,0)==(1<<c)-1,'c90 full labelled partition')
# Independent actual strict-danger construction, not padded-mask arithmetic.
num,den=126,1009
tails=(542,55082,25292,211773,30513,51126)
physical=[]
for w,g in zip(tails,gs):
    need(gcd(c,w)==g,'actual c90 tail gcd')
    need((w//g)%(c//g)==1,'actual affine unit')
    danger=0
    for j in range(c):
        phase=(w*(num+den*j))%(c*den)
        if 14*min(phase,c*den-phase)<c*den:danger|=1<<j
    physical.append(danger)
need(physical==cover,'actual common-phase masks equal padded partition')
need(all(14*min((v*num)%den,(-v*num)%den)>den for v in range(1,8)),
     'actual seven-speed body safe at y=126/1009')
row=tuple(c*v for v in range(1,8))+tails
need(len(set(row))==13 and gcd(*row)==1,'actual full row primitive and distinct')
need(all(14*min(v%(14*c),(-v)%(14*c))>=14*c for v in row),
     'full actual row safe at another time 1/(14c)')
# Hostile to forgetting shared labels: changing one mask address can destroy this cover.
shifted=cover[:-1]+[mask(c,15,1,1)]
need(reduce(int.__or__,shifted,0)!=(1<<c)-1,
     'independently readdressing one mask can destroy cover')

bank={'mask_counts':expected,'trees':{'first':A_trees,'second':B_trees},'bounds':case_bounds,
      'c90':{'gcds':gs,'params':params,'masks':cover,'phase':[num,den],'tails':tails}}
semantic=sha256(json.dumps(bank,sort_keys=True,separators=(',',':')).encode()).hexdigest()
atlas_digest=sha256(json.dumps({q:sorted(rows) for q,rows in atlas.items()},sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('STATUS: PROVED c96 padded-mask exclusion; actual c90 common-phase cover; LRC14 remains OPEN')
print('c96_affine_mask_counts',expected)
print('literal_pair_cases',pair_cases,'literal_disjoint_shadow_triples',shadow_cases)
print('c96_case_union_bounds',json.dumps(case_bounds,sort_keys=True))
print('c90_gcds',gs,'q_start_unit',params)
print('c90_actual_tail_speeds',tails,'body_phase',(num,den),'full_row_safe_time',(1,1260))
print('c90_mask_hex',[hex(m) for m in cover])
print('gates',GATES)
print('mask_atlas_sha256',atlas_digest)
print('semantic_sha256',semantic)
