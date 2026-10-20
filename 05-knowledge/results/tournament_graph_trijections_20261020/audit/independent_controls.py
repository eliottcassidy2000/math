#!/usr/bin/env python3
"""Independent audit: no import from the producer verifier. Exact integer controls."""
import itertools as it
from collections import Counter, defaultdict, deque
import json, hashlib
from pathlib import Path
checks=0

def demand(ok, text):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(text)

def basis(vectors):
    b={}
    for x in vectors:
        while x:
            k=x.bit_length()-1
            if k in b: x ^= b[k]
            else:
                b[k]=x
                break
    return b

def reduce(x,b):
    for k in sorted(b,reverse=True):
        if x >> k & 1: x ^= b[k]
    return x

def setup(n):
    edges=list(it.combinations(range(1,n+1),2))
    bits={e:1<<i for i,e in enumerate(edges)}
    edge=lambda u,v: bits[tuple(sorted((u,v)))]
    cuts=[sum(w for (u,v),w in bits.items() if ((s>>(u-1))^(s>>(v-1)))&1) for s in range(1<<(n-1))]
    cb=basis(cuts)
    return edges,bits,edge,cuts,cb

out={'status':'PASS','generator_census':[], 'constraints':[], 'literal_controls':{}}
for n in range(2,9):
    edges,bits,edge,cuts,cb=setup(n)
    full=(1<<len(edges))-1
    gs=Counter()
    paths={}
    collision=None
    all_path_masks=[]
    for p in it.permutations(range(1,n+1)):
        if p[0]>p[-1]: continue
        pm=sum(edge(u,v) for u,v in zip(p,p[1:]))
        g=reduce(full^pm,cb)
        gs[g]+=2
        all_path_masks.append(pm)
        if collision is None and g in paths:
            prev=paths[g]
            prevmask=sum(edge(u,v) for u,v in zip(prev,prev[1:]))
            collision={'paths':[prev,p],'symmetric_difference_edges':[e for e,w in bits.items() if (pm^prevmask)&w]}
        paths[g]=p
    first=next(iter(gs))
    dim=len(edges)-n+1
    gr=len(basis(gs))
    dr=len(basis(g^first for g in gs))
    if n>=3:
        demand(0 not in gs,'zero generator')
        demand(gr==dim-(n%4==1),'generator rank')
        demand(dr==dim-(n%2==1),'difference rank')
        demand((gr==dr+1)==(n%4==3),'bipartite rank')
        demand(len(set(gs.values()))==1,'nonuniform weights')
        demand(len(basis(pm^all_path_masks[0] for pm in all_path_masks))==len(edges)-1,'ambient even span')
    else:
        demand(gs=={0:2},'n2 must be a loop aggregate')
    if n==8: demand(len(gs)==20160 and set(gs.values())=={2},'n8 injectivity')
    entry={'n':n,'dimension':dim,'rank':gr,'difference_rank':dr,'degree':len(gs),'weights':dict(Counter(gs.values()))}
    if n==7: entry['sharp_threshold_witness']=collision
    if n<=6:
        nodes={reduce(x,cb) for x in range(1<<len(edges))}
        comps=[]; colors={}; bipartite=True
        while nodes:
            start=nodes.pop(); found={start}; todo=deque([start]); colors[start]=0
            while todo:
                u=todo.popleft()
                for g in gs:
                    v=u^g
                    if v not in found:
                        found.add(v); nodes.discard(v); todo.append(v); colors[v]=colors[u]^1
                    elif colors[v]==colors[u]: bipartite=False
            comps.append(len(found))
        entry['component_sizes']=sorted(comps)
        entry['bipartite']=bipartite
        if n>=3:
            demand(len(comps)==(2 if n%4==1 else 1),'BFS components')
            demand(bipartite==(n%4==3),'BFS bipartite')
    out['generator_census'].append(entry)

# Exhaust every constraint/pattern for n<=4 and verify compatibility independently
# by solving edge equations through graph traversal, including every cycle edge.
for n in range(2,5):
    edges,bits,edge,cuts,cb=setup(n)
    full=(1<<len(edges))-1
    reps=sorted({reduce(x,cb) for x in range(full+1)})
    cases=0
    for F in range(full+1):
        fedges=[e for e,w in bits.items() if F&w]
        adjacency={u:[] for u in range(1,n+1)}
        for u,v in fedges: adjacency[u].append(v); adjacency[v].append(u)
        for aa in range(1<<len(fedges)):
            a=sum(bits[e] for i,e in enumerate(fedges) if aa>>i&1)
            image=0
            for x in reps:
                potentials={}; c=0; compatible=True
                for root in range(1,n+1):
                    if root in potentials: continue
                    c+=1; potentials[root]=0; todo=[root]
                    while todo:
                        u=todo.pop()
                        for v in adjacency[u]:
                            val=potentials[u]^bool((x^a)&edge(u,v))
                            if v not in potentials: potentials[v]=val; todo.append(v)
                            elif potentials[v]!=val: compatible=False
                actual=[x^cut for cut in cuts if ((x^cut)&F)==a]
                demand(bool(actual)==compatible,'cycle compatibility')
                demand(len(actual)==(1<<(c-1) if compatible else 0),'constraint multiplicity')
                image+=bool(actual)
            demand(image==1<<(len(edges)-len(fedges)-c+1),'constraint image count')
            cases+=1
    out['constraints'].append({'n':n,'constraint_pattern_cases':cases})

# Hamilton paths by subset dynamic programming, distinct from producer word tests.
def hp_count(n,x,edge):
    dp={(1<<(v-1),v):1 for v in range(1,n+1)}
    for mask in range(1,1<<n):
        for v in range(1,n+1):
            z=dp.get((mask,v),0)
            if not z: continue
            for w in range(1,n+1):
                if not mask>>(w-1)&1 and bool(v<w)^bool(x&edge(v,w)):
                    key=(mask|(1<<(w-1)),w); dp[key]=dp.get(key,0)+z
    return sum(dp.get(((1<<n)-1,v),0) for v in range(1,n+1))

edges,bits,edge,cuts,cb=setup(4)
x=edge(1,4); y=edge(2,4)|edge(3,4)
demand(reduce(x,cb)==reduce(y,cb),'literal hostile different cut class')
demand((hp_count(4,x,edge),hp_count(4,y,edge))==(5,1),'literal H hostile')
for a,b,c in it.product(range(2),repeat=3):
    x=a*edge(1,3)|b*edge(1,4)|c*edge(2,4)
    u,v,w=a,b^c,a^b
    y=u*edge(2,3)|v*edge(2,4)|w*edge(3,4)
    demand(reduce(x,cb)==reduce(y,cb),'literal chart class')
    demand((u,u^w,u^v^w)==(a,b,c),'inverse chart')
    shifted=(u^1)*edge(2,3)|v*edge(2,4)|w*edge(3,4)
    demand(reduce(x^edge(1,3)^edge(1,4)^edge(2,4),cb)==reduce(shifted,cb),'transported complement')
out['literal_controls']['order4_H_change']=[5,1]

# The smallest explicit aggregate-nonlifting hostile: original A3 is two K1,3s.
n=3
edges,bits,edge,cuts,cb=setup(n); full=(1<<len(edges))-1
adj={x:set() for x in range(full+1)}
for x in adj:
    for p in it.permutations(range(1,n+1)):
        if all(bool(u<v)^bool(x&edge(u,v)) for u,v in zip(p,p[1:])):
            adj[x].add(x^full^sum(edge(u,v) for u,v in zip(p,p[1:])))
remaining=set(adj); comps=[]
while remaining:
    root=remaining.pop(); todo=[root]; found={root}
    while todo:
        u=todo.pop()
        for v in adj[u]-found: found.add(v); remaining.discard(v); todo.append(v)
    comps.append(sorted(found))
demand(sorted(map(len,comps))==[4,4],'original A3 components')
demand(all(sorted(len(adj[x]) for x in c)==[1,1,1,3] for c in comps),'original A3 not two stars')
class_degrees=defaultdict(list)
for x in adj: class_degrees[reduce(x,cb)].append(len(adj[x]))
demand(all(sorted(v)==[1,1,1,3] for v in class_degrees.values()),'equitable hostile')
out['literal_controls']['order3_original_components']=comps
out['literal_controls']['order3_original_component_degrees']=[[1,1,1,3],[1,1,1,3]]
out['literal_controls']['order3_aggregate']='K2 with weight 6'
out['literal_controls']['order3_switching_class_degree_multisets']={str(k):sorted(v) for k,v in class_degrees.items()}
out['checks']=checks
out['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
print(json.dumps(out,indent=2,sort_keys=True))
