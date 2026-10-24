#!/usr/bin/env python3
"""Independent finite audit; no imports from, or generated data of, the producer.

The vector coordinates below use the first four bits of the even 6-bit
representative whose sixth bit is zero. They are not the producer's basis.
The alternating form is I+J rather than two standard symplectic blocks.
"""
from collections import Counter
from itertools import combinations, permutations, product
import json

gates = 0

def verify(condition, name):
    global gates
    gates += 1
    if not condition:
        raise ValueError(name)

def parity(x):
    return x.bit_count() % 2

def form(x,y):
    return parity(x & y) ^ (parity(x) & parity(y))

def rep(x):
    return x | (parity(x) << 4)

def coord(mask):
    if mask & 32:
        mask ^= 63
    return mask & 15

def mul(a,b):
    return tuple(a[t] for t in b)

def ctype(g):
    remaining = set(range(len(g)))
    out=[]
    while remaining:
        i = min(remaining)
        orbit = []
        while i not in orbit:
            orbit.append(i)
            i = g[i]
        remaining.difference_update(orbit)
        out.append(len(orbit))
    return tuple(sorted(out))

points=range(6)
edges=tuple(frozenset(e) for e in combinations(points,2))
edge_id={e:i for i,e in enumerate(edges)}
# Enumerate edge subsets, rather than recursively pair the first free point.
matchings=tuple(frozenset(t) for t in combinations(range(15),3)
                if all(sum(v in edges[i] for i in t)==1 for v in points))
match_id={m:i for i,m in enumerate(matchings)}
verify(len(matchings)==15,'15 perfect matchings')
for x,y,z in product(range(16),repeat=3):
    verify(form(x,y^z)==form(x,y)^form(x,z),'bilinear alternating form')
verify(all(form(x,x)==0 for x in range(16)),'alternating')
verify([x for x in range(16) if all(form(x,y)==0 for y in range(16))]==[0], 'nondegenerate')
lines=tuple(frozenset(t) for t in combinations(range(1,16),3)
            if t[0]^t[1]^t[2]==0 and form(t[0],t[1])==0)
verify(len(lines)==15,'15 Lagrangians')
verify(Counter(sum(v in L for L in lines) for v in range(1,16))=={3:15},'point degree 3')
edge_vector=tuple(coord(sum(1<<v for v in e)) for e in edges)
verify(len(set(edge_vector))==15 and 0 not in edge_vector,'all nonzero vectors represented by duads')
verify({frozenset(edge_vector[i] for i in M) for M in matchings}==set(lines),'matching-to-line bijection')
for i,j in combinations(range(15),2):
    verify((form(edge_vector[i],edge_vector[j])==0)==edges[i].isdisjoint(edges[j]),'orthogonality is duad disjointness')

neighbors=tuple(frozenset(j for j,N in enumerate(matchings) if i!=j and M&N)
                for i,M in enumerate(matchings))
verify(Counter(map(len,neighbors))=={6:15},'matching graph degree 6')
for i,j in combinations(range(15),2):
    verify(len(neighbors[i]&neighbors[j])==(1 if j in neighbors[i] else 3),'SRG common-neighbor counts')
pentads=tuple(frozenset(t) for t in combinations(range(15),5)
              if sorted(e for i in t for e in matchings[i])==list(range(15)))
verify(len(pentads)==6,'six one-factorizations')
pentad_id={p:i for i,p in enumerate(pentads)}
pairs=tuple(frozenset(i for i,P in enumerate(pentads) if m in P) for m in range(15))
verify(set(pairs)=={frozenset(t) for t in combinations(range(6),2)},'syntheme-to-pentad-duad bijection')
for i,j in combinations(range(15),2):
    verify((j in neighbors[i])==pairs[i].isdisjoint(pairs[j]),'KG(6,2) adjacency')

all_perms=tuple(permutations(points))
actions={}
outer={}
fingerprints={}
for g in all_perms:
    ea=tuple(edge_id[frozenset(g[v] for v in e)] for e in edges)
    ma=tuple(match_id[frozenset(ea[e] for e in M)] for M in matchings)
    pa=tuple(pentad_id[frozenset(ma[m] for m in P)] for P in pentads)
    actions[g]=ma
    outer[g]=pa
    for m in range(15):
        verify(pairs[ma[m]]==frozenset(pa[k] for k in pairs[m]),'outer-twisted equivariance')
        moved_vectors=frozenset(coord(sum(1<<g[v] for v in range(6) if rep(edge_vector[e])>>v&1)) for e in matchings[m])
        verify(moved_vectors==frozenset(edge_vector[e] for e in matchings[ma[m]]),'natural matching-Lagrangian equivariance')
    finger=(sum(i==ea[i] for i in range(15)),sum(i==ma[i] for i in range(15)),ctype(pa))
    typ=ctype(g)
    verify(typ not in fingerprints or fingerprints[typ]==finger,'class fingerprint well-defined')
    fingerprints[typ]=finger
verify(len(set(outer.values()))==720,'faithful onto six-point action')
# This checks every pair, not merely generator multiplication.
for g,h in product(all_perms,repeat=2):
    verify(outer[mul(g,h)]==mul(outer[g],outer[h]),'all 518400 homomorphism products')
trans=(1,0,2,3,4,5)
verify(fingerprints[ctype(trans)]==(7,3,(2,2,2)),'natural-action obstruction')
triple=(1,0,3,2,5,4)
verify(fingerprints[ctype(triple)][:2]==(3,7),'reversed obstruction')
dm=next(i for i,e in enumerate(edges) if e==frozenset([0,1]))
dstab=tuple(g for g in all_perms if frozenset(g[v] for v in edges[dm])==edges[dm])
mstab=tuple(g for g in all_perms if actions[g][0]==0)
verify(len(dstab)==len(mstab)==48,'two stabilizer sizes')
def orbit_sizes(group):
    return sorted(len(o) for o in {frozenset(g[v] for g in group) for v in range(6)})
verify(orbit_sizes(dstab)==[2,4] and orbit_sizes(mstab)==[6],'stabilizers not conjugate')

# Independently transcribed from Atlas section 4.14, not read from candidate JSON.
source_cycles=[[(8,9),(10,11),(12,13),(14,15)],
               [(1,9,8),(2,15,14),(4,12,13),(6,10,11)],
               [(1,4),(2,14,7,15),(3,10,6,11),(8,9,13,12)]]
gens=[]
for cycles in source_cycles:
    d=list(range(15))
    for cyc in cycles:
        for i in range(len(cyc)):
            d[cyc[i]-1]=cyc[(i+1)%len(cyc)]-1
    gens.append(tuple(d))
# Fixed-point closure over a whole set, not the producer's DFS.
group={tuple(range(15))}
while True:
    enlarged=group|{mul(g,h) for g in group for h in gens}
    if enlarged==group:
        break
    group=enlarged
verify(len(group)==24,'Bolza displayed action order 24')
orbits=sorted({tuple(sorted({g[v]+1 for g in group})) for v in range(15)})
verify(orbits==[(1,4,8,9,12,13),(2,7,14,15),(3,6,10,11),(5,)],'Bolza published orbit partition')
stabilizers=[sum(g[O[0]-1]+1==O[0] for g in group) for O in orbits]
verify(stabilizers==[4,6,6,24],'Bolza stabilizers')
verify(len({tuple(g[v-1]+1 for v in orbits[1]) for g in group})==24,'faithful four-point action gives S4')
verify(sum(map(len,orbits))==15,'weighted total 15, not 4')
verify(8000%5==0 and 8000%7==1728%7,'exceptional j collisions')
verify(8000==2**6*5**3 and 8000-1728==2**7*7**2,'only positive-prime j collisions')

result={
 'status':'PASS; finite actions checked independently, geometric quotients remain cited',
 'gates':gates,
 'universe':{'duads':15,'nonzero_vectors':15,'lagrangians':15,'flags':45,'synthemes':15,'pentads':6,'S6':720},
 'graph':{'SRG':[15,6,1,3],'KG6_2_isomorphism':'syntheme -> pair of containing pentads'},
 'fixed_point_characters':[{'cycle_type':key,'duads':value[0],'synthemes':value[1],'pentad_cycle_type':value[2]} for key,value in sorted(fingerprints.items())],
 'stabilizers':{'duad':{'order':48,'natural_orbits':[2,4]},'syntheme':{'order':48,'natural_orbits':[6]}},
 'bolza':{'order':24,'orbits':orbits,'stabilizers':stabilizers,'full_geometric_group':'CITED only in Atlas characteristic scope','weighted_degree':15,'edge_classes':4},
 'scope':{'Atlas_excluded_characteristics':[2,3,5],'p7':'three distinct loop classes of weights 4,4,1 plus product class of weight 6; two destination vertices, not four','p13_p29':'one four-kernel orbit target specializes from IV to V','j0_prime_divisors':[2,5],'j1728_prime_divisors':[2,7]},
 'hostile_claims_rejected':['natural S6-equivariant duad-syntheme bijection','same stabilizer order means conjugacy','syntheme graph is the Richelot graph','four edge classes always mean four destination vertices','a count of 15 identifies the Ogg prime set or selects modulus 15','finite permutation closure proves all characteristic specializations','Bolza CM bridge supplies an MRB value or defect formula']
}
print(json.dumps(result,indent=2,sort_keys=True))
