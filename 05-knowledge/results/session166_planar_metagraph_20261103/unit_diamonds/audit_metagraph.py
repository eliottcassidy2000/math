#!/usr/bin/env python3
from itertools import permutations,combinations
from functools import lru_cache
import json
from check import coords,norm2,ONE,demand
N=6;PAIRS=list(combinations(range(N),2));POSITIONS={p:i for i,p in enumerate(PAIRS)}
PERMS=list(permutations(range(N)))
def edge(mask,i,j):
    return bool(mask>>POSITIONS[(i,j)]&1) if i<j else not bool(mask>>POSITIONS[(j,i)]&1)
@lru_cache(None)
def canon(mask):
    return min(sum(1<<k for k,(i,j) in enumerate(PAIRS) if edge(mask,p[i],p[j])) for p in PERMS)
def strong(mask):
    for start in range(N):
        seen={start};todo=[start]
        while todo:
            v=todo.pop()
            for w in range(N):
                if w!=v and w not in seen and edge(mask,v,w):seen.add(w);todo.append(w)
        if len(seen)!=N:return False
    return True
names=list('OABCDEF');masks=[18,22,26,181,81,146,145]
cs=[canon(m) for m in masks]
demand(len(set(cs))==7,'classes collide');demand(all(strong(m) for m in masks),'not all strong')
expected={''.join(sorted(e)) for e in ['OA','OB','AB','AC','BC','OD','OE','DE','DF','EF','CF']}
actual=set();witnesses={}
for i,m in enumerate(masks):
    for j in range(i+1,len(masks)):
        flips=[PAIRS[b] for b in range(15) if canon(m^(1<<b))==cs[j]]
        if flips:
            label=''.join(sorted(names[i]+names[j]));actual.add(label);witnesses[label]=flips
demand(actual==expected,{'extra':sorted(actual-expected),'missing':sorted(expected-actual)})
p=coords(2);point_indices=[0,2,3,1,5,6,4];actual_geo=set()
for i,j in combinations(range(7),2):
    if norm2(p[point_indices[i]],p[point_indices[j]])==ONE:actual_geo.add(''.join(sorted(names[i]+names[j])))
demand(actual_geo==expected,'exact Moser geometry mismatch')
print(json.dumps({'vertex_names':names,'tournament_masks':masks,'canonical_masks':cs,'all_strong':True,'distinct_isomorphism_classes':True,'induced_single_arc_flip_edges':sorted(actual),'representative_arc_flip_witnesses':witnesses,'exact_planar_unit_edges':sorted(actual_geo),'field':'Q(sqrt(3),sqrt(11))','planar_coords_indices_in_diamond_chain':point_indices,'scope':'Simple single-arc-reversal metagraph on isomorphism classes, not the repository complement-line merged metagraph.'},indent=2,sort_keys=True))
