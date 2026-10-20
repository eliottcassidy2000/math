#!/usr/bin/env python3
"""Exact F5 determinant model and prior regular-tournament return graph."""
from itertools import product, permutations, combinations
from collections import Counter, deque
import json
p=5
checks=Counter()
def require(x,tag):
    checks[tag]+=1
    if not x:raise RuntimeError(tag)
def neg(v):return tuple(-a%5 for a in v)
def cls(v):return min(v,neg(v))
def det(v,w):return (v[0]*w[1]-v[1]*w[0])%5
def adj(v,w):return det(v,w)**2%5==1
X=tuple(product(range(5),repeat=2)); NZ=tuple(v for v in X if v!=(0,0)); V=tuple(sorted({cls(v) for v in NZ}))
N=cls((1,0));S=cls((2,0));U=tuple(cls((i,1)) for i in range(5));L=tuple(cls(((2*i-1)%5,2)) for i in range(5))
require(set((N,S)+U+L)==set(V),'twelve_distinct_vertices')
model={N:{*U},S:{*L}}
for i in range(5):
    model[U[i]]={N,U[(i-1)%5],U[(i+1)%5],L[i],L[(i+1)%5]}
    model[L[i]]={S,L[(i-1)%5],L[(i+1)%5],U[i],U[(i-1)%5]}
for v,w in product(V,repeat=2):
    require(adj(v,w)==(w in model[v]),'icosahedron_adjacency')
    for a,b in product((1,4),repeat=2):
        require(adj(tuple(a*x%5 for x in v),tuple(b*x%5 for x in w))==adj(v,w),'sign_quotient_well_defined')
    if det(v,w)==0:require(w in (v,cls(tuple(2*x%5 for x in v))),'zero_det_antipodal')
for v in V:
    counts=Counter(det(v,w)**2%5 for w in V)
    require(counts=={0:2,1:5,4:5},'per_vertex_determinant_counts')
    require(cls(tuple(4*x%5 for x in v))==v,'antipodal_involution')
# Literal matrix counts.
matrices=tuple(product(range(5),repeat=4))
nonzero_columns=[a for a in matrices if a[:2]!=(0,0) and a[2:]!=(0,0)]
dets=Counter(det(a[:2],a[2:]) for a in nonzero_columns)
require(len(matrices)==625 and len(nonzero_columns)==576,'625_to_576')
require(dets=={0:96,1:120,2:120,3:120,4:120},'matrix_determinant_counts')
paircounts=Counter(('same' if v==w else 'antipodal') if det(v,w)==0 else ('edge' if adj(v,w) else 'distance_two') for v,w in product(V,repeat=2))
require(paircounts=={'same':12,'antipodal':12,'edge':60,'distance_two':60},'144_pair_decomposition')
# Prior repository's explicit tournament-code certificate, independently evaluated.
codes={N:947,S:76}
codes.update(zip(U,(613,217,710,341,668)))
codes.update(zip(L,(313,682,355,410,806)))
E=tuple(combinations(range(5),2)); idx={e:i for i,e in enumerate(E)}
def beats(c,a,b):
    z=(c>>idx[tuple(sorted((a,b)))])&1
    return bool(z) if a<b else not bool(z)
def relabel(c,g):
    inverse={g[i]:i for i in range(5)}
    return sum(int(beats(c,inverse[a],inverse[b]))<<i for i,(a,b) in enumerate(E))
def even(g):return sum(g[i]>g[j] for i in range(5) for j in range(i+1,5))%2==0
for v,w in product(V,repeat=2):require(adj(v,w)==((codes[v]^codes[w]).bit_count()==6),'prior_return_graph_isomorphism')
for v in V:
    require(all(sum(beats(codes[v],a,b) for b in range(5) if b!=a)==2 for a in range(5)),'regular_tournament_codes')
    require(codes[cls(tuple(2*x%5 for x in v))]==1023^codes[v],'antipode_is_converse')
# Matrix is column-major a,b,c,d: [[a,c],[b,d]].
def matvec(m,v):return ((m[0]*v[0]+m[2]*v[1])%5,(m[1]*v[0]+m[3]*v[1])%5)
SL=[m for m in matrices if det(m[:2],m[2:])==1]
sl_actions={tuple(cls(matvec(m,v)) for v in V):m for m in SL}
require(len(SL)==120 and len(sl_actions)==60,'sl2_kernel_size_two')
revcode={c:v for v,c in codes.items()}
a5_actions={tuple(revcode[relabel(codes[v],g)] for v in V):g for g in permutations(range(5)) if even(g)}
require(len(a5_actions)==60,'a5_faithful')
require(set(sl_actions)==set(a5_actions),'sl2_mod_center_equals_a5_action')
gens=((1,0,1,1),(0,1,4,0))
generator_records=[]
for m in gens:
    a=tuple(cls(matvec(m,v)) for v in V)
    generator_records.append({'matrix_columns':m,'tournament_label_permutation':a5_actions[a]})
# Actual source anti-reversal transports to converse plus a fixed even relabelling.
reflection_perm=(1,0,2,4,3)
reflection_fixed=[]
for v in V:
    rv=cls((v[1],v[0]))
    require(codes[rv]==(1023^relabel(codes[v],reflection_perm)),'source_reflection_transport')
    if rv==v:reflection_fixed.append(v)
require(len(reflection_fixed)==4,'four_reflection_fixed_orbits')
require(sum(v[0]==v[1] for v in reflection_fixed)==2,'two_blue_fixed_orbits')
require(sum((v[0]+v[1])%5==0 for v in reflection_fixed)==2,'two_black_paired_fixed_orbits')
# Changing an allowed sign-equivariant field marking changes the edge predicate.
def cube(v):return tuple(x**3%5 for x in v)
require(adj((1,1),(2,1)) and not adj(cube((1,1)),cube((2,1))),'marking_dependence_hostile')
# Independently reconstruct and enumerate the actual balanced source paths.
source_states=range(5)
source_label=(1,2,3,4,0)
left={0,1,2};right={4,5,6}
def subset(i):return left if i==4 else {j+1 for j in range(2) if i>>j&1}
def source_path(i,k):
    I=subset(i);K=subset(k);J={6-x for x in left-K}
    return tuple(sorted(I)+sorted(J)+[3]+sorted(left-I)+sorted(right-J))
def source_beats(a,b):
    return (a<b) != (tuple(sorted((a,b))) in ((0,3),(3,6)))
source_rows=[]
source_paths=set()
for i,k in product(source_states,repeat=2):
    path=source_path(i,k);v=(source_label[i],source_label[k]);source_paths.add(path)
    require(len(path)==7 and set(path)==set(range(7)),'literal_source_vertex_coverage')
    require(all(source_beats(path[j],path[j+1]) for j in range(6)),'literal_source_path_legality')
    require(tuple(6-x for x in path[::-1])==source_path(k,i),'literal_source_reflection')
    source_rows.append({'states':[i,k],'field_vector':v,'path':path})
all_source_paths={q for q in permutations(range(7)) if all(source_beats(q[i],q[i+1]) for i in range(6))}
require(len(source_paths)==25 and source_paths==all_source_paths,'literal_source_full_path_census')
# The source path labels supplied by the balanced-fibre lane.
labels={'infinity':0,'empty':1,'one_two':4,'one':2,'two':3}
require(set(labels.values())==set(range(5)),'source_five_state_label_bijection')
print(json.dumps({'status':'PASS','total_checks':sum(checks.values()),'check_counts':dict(checks),'field_vertices':{'N':N,'S':S,'U':U,'L':L},'field_to_prior_tournament_code':[{'vector_mod_sign':v,'code':codes[v]} for v in V],'determinants_of_576_matrices':dict(dets),'quotient_144_pairs':dict(paircounts),'sl2_generators_to_a5':generator_records,'source_labels':labels,'literal_source_paths':source_rows,'source_reflection_tournament_permutation':reflection_perm,'source_reflection_fixed_vectors':reflection_fixed},sort_keys=True,indent=2))
