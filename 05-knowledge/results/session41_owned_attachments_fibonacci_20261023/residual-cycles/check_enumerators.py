#!/usr/bin/env python3
from itertools import product
from collections import Counter
from pathlib import Path
import hashlib,json
counts=Counter()
def need(v,label):
    if not v:raise RuntimeError(label)
    counts[label]+=1
pairs=list(product((0,1),repeat=2))
relations=[tuple(p for k,p in enumerate(pairs) if m>>k&1) for m in range(16)]
def mul(P,Q):
    R=Counter()
    for p,c in P.items():
        for q,d in Q.items():R[p+q]+=c*d
    return R
# A star is deliberately chosen: three children share the same root state.
# Observations have both cell and genuinely pair-dependent edge costs.
edges=((0,1),(0,2),(0,3))
def node(i,a):return (i+1)*a
def edge(j,a,b):return (-1)**j*(a+2*b+3*a*b)
for rels in product(relations,repeat=3):
    direct=Counter()
    for x in product((0,1),repeat=4):
        if all((x[i],x[j]) in R for (i,j),R in zip(edges,rels)):
            exponent=sum(node(i,a) for i,a in enumerate(x))+sum(edge(j,x[0],x[j]) for _,j in edges)
            direct[exponent]+=1
    recursive=Counter()
    for a in (0,1):
        P=Counter({node(0,a):1})
        for j,R in enumerate(rels,1):
            child=Counter()
            for aa,b in R:
                if aa==a:child[node(j,b)+edge(j,a,b)]+=1
            P=mul(P,child)
        recursive.update(P)
    need(direct==recursive,'all_star_coefficients_match')
    if direct:
        # Independent min/max DP, not derived by reading polynomial support.
        rootmin=[];rootmax=[]
        for a in (0,1):
            intervals=[]
            for j,R in enumerate(rels,1):
                vals=[node(j,b)+edge(j,a,b) for aa,b in R if aa==a]
                if not vals:break
                intervals.append((min(vals),max(vals)))
            if len(intervals)==3:
                rootmin.append(node(0,a)+sum(v[0] for v in intervals))
                rootmax.append(node(0,a)+sum(v[1] for v in intervals))
        need((min(rootmin),max(rootmax))==(min(direct),max(direct)),'minmax_certificate_exact')
        need((min(rootmin)==max(rootmax))==(len(direct)==1),'minmax_uniqueness_iff')
correct=Counter({0:1,2:1})
naive=mul(Counter({0:1,1:1}),Counter({0:1,1:1}))
need(naive==Counter({0:1,1:2,2:1}) and naive!=correct,'shared_owner_contraction_not_multiplication')
forced={(a,b):int(not(a and 0) and not(0 and b)) for a,b in pairs}
free={(a,b):sum(not(a and z) and not(z and b) for z in (0,1)) for a,b in pairs}
reconnected={(a,b):int(not(a and b)) for a,b in pairs}
need(tuple(forced[p] for p in pairs)==(1,1,1,1),'forced_zero_elimination_allows_all_pairs')
need(tuple(free[p] for p in pairs)==(2,1,1,1),'free_zero_elimination_retains_multiplicity')
need(forced!=reconnected and free!=reconnected,'zero_vertex_deletion_not_edge_contraction')
need({(1-t)+t for t in (0,1)}=={1},'coupling_rescues_sum_projection_uniqueness')
need(len({1-t for t in (0,1)})==2,'same_coupled_fiber_can_change_other_projection')
print(json.dumps({'status':'PASS','counts':dict(counts),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
