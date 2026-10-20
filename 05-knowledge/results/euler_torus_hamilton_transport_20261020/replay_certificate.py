#!/usr/bin/env python3
"""Independent completeness count + direct replay of the stored face-flip forest."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
P=Path(__file__).resolve().parent

def check(condition,message):
    if not condition: raise RuntimeError(message)

def count_degree_coefficient(n=7,r=3):
    # For each unoriented edge ij, choose winner i or j: factor (x_i+x_j).
    coeff={(0,)*n:1}; maximum_states=1
    for i,j in combinations(range(n),2):
        nxt=Counter()
        for deg,c in coeff.items():
            for winner in (i,j):
                if deg[winner]==r: continue
                d=list(deg);d[winner]+=1;nxt[tuple(d)]+=c
        coeff=nxt;maximum_states=max(maximum_states,len(coeff))
    return coeff[(r,)*n],maximum_states

def omega(out):
    steps={1:(1,0),2:(-1,1),3:(0,1),4:(0,-1),5:(1,-1),6:(-1,0)}
    dx=dy=0
    for i,j in combinations(range(7),2):
        k=(j-i)%7 if out[i]>>j&1 else (i-j)%7
        x,y=steps[k];dx+=x;dy+=y
    check((dx+3*dy)%7==0,'closed total lattice chain')
    return ((dx+3*dy)//7,dy)

C=json.loads((P/'face_flip_certificate.json').read_text())
states=[tuple(t) for t in C['orientations']];faces=C['faces']
check(len(states)==len(set(states)),'state uniqueness')
for t in states:
    check(len(t)==7 and all(0<=x<128 and x.bit_count()==3 for x in t),'regular bitmask')
    check(all(not(t[i]>>i&1) for i in range(7)),'no loops')
    check(all(bool(t[i]>>j&1)!=bool(t[j]>>i&1) for i,j in combinations(range(7),2)),'exactly one orientation per edge')
total,maximum_states=count_degree_coefficient()
check(total==len(states)==2640,'independent polynomial completeness count')
expected_faces={tuple(sorted((i,(i+1)%7,(i+3)%7))) for i in range(7)}|{tuple(sorted(((i+1)%7,(i+3)%7,(i+4)%7))) for i in range(7)}
check({tuple(f) for f in faces}==expected_faces and len(faces)==14,'correct face set')
seen=set();root_of={};component_sizes=Counter();flux_roots={}
for child,parent,face_id in C['forest_child_parent_face']:
    check(0<=child<len(states) and child not in seen,'child occurs exactly once')
    if parent==-1:
        check(face_id==-1,'root marker')
        root=child;f=omega(states[child]);check(f not in flux_roots,'distinct root fluxes')
        flux_roots[f]=root
    else:
        check(parent in seen and 0<=face_id<len(faces),'previous parent and valid face')
        tri=faces[face_id];old=states[parent];new=list(old)
        check(all(sum(old[i]>>j&1 for j in tri)==1 for i in tri),'parent face is directed')
        for i,j in combinations(tri,2):new[i]^=1<<j;new[j]^=1<<i
        check(tuple(new)==states[child],'specified face actually gives child')
        root=root_of[parent]
        check(omega(states[child])==omega(states[root]),'tree flux conservation')
    root_of[child]=root;seen.add(child);component_sizes[root]+=1
check(len(seen)==len(states),'forest spans entire universe')
check(len(component_sizes)==25,'25 components')
check(Counter(component_sizes.values())=={1:6,14:6,42:6,273:6,660:1},'component size multiset')
result={'coefficient_count':total,'largest_polynomial_DP_layer':maximum_states,'forest_states':len(seen),
        'flux_components':len(component_sizes),'size_multiplicities':dict(sorted(Counter(component_sizes.values()).items())),
        'status':'PASS; universe completeness and every spanning-tree edge independently replayed'}
(P/'certificate_replay.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
