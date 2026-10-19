#!/usr/bin/env python3
"""Independent finite proof controls using directed edge sets and permutations."""
from itertools import combinations, permutations, product
from functools import lru_cache
import json
checks=0

def ck(value, label):
    global checks
    checks+=1
    if not value: raise RuntimeError(label)

@lru_cache(None)
def maps(n):return tuple(permutations(range(n)))

@lru_cache(None)
def embed(n,S,T):
    H=set(T)
    return sum(all((p[u],p[v]) in H for u,v in S) for p in maps(n))

def remove(n,S,v):
    vertices=[u for u in range(n) if u!=v]
    remap={u:i for i,u in enumerate(vertices)}
    return tuple(sorted((remap[u],remap[w]) for u,w in S if v not in (u,w)))

def source(n,T):return tuple(sorted(T+tuple((n,v) for v in range(n))))

def all_graphs(n, complete=False):
    pairs=list(combinations(range(n),2))
    for word in product((1,2) if complete else (0,1,2),repeat=len(pairs)):
        yield tuple(sorted((u,v) if x==1 else (v,u) for (u,v),x in zip(pairs,word) if x))

hosts=[tuple(all_graphs(n,True)) for n in range(5)]
classes=[]
for n,ts in enumerate(hosts):
    reps=[]
    for T in ts:
        if not any(embed(n,S,T) for S in reps):reps.append(T)
    classes.append(tuple(reps))
    for i,S in enumerate(reps):
        for j,T in enumerate(reps):ck(embed(n,S,T)%2==(i==j),'tournament indicator basis')
ck([len(c) for c in classes]==[1,1,1,2,4],'class counts')

# Both graph-level action formulas, including degree-zero padding.
graph_counts=[]
for n in range(5):
    graphs=tuple(all_graphs(n));graph_counts.append(len(graphs))
    for S in graphs:
        if n<4:
            for T in hosts[n+1]:
                left=sum(embed(n,S,remove(n+1,T,v)) for v in range(n+1))%2
                right=embed(n+1,S,T)%2  # label n is the newly isolated vertex
                ck(left==right,'padding action')
        if n:
            zero_indegree=[v for v in range(n) if all(w!=v for u,w in S)]
            for T in hosts[n-1]:
                left=embed(n,S,source(n-1,T))%2
                right=sum(embed(n-1,remove(n,S,v),T) for v in zero_indegree)%2
                ck(left==right,'source-deletion action')
        if n in (2,4) and any(all(v not in edge for edge in S) for v in range(n)):
            ck(any(embed(n,S,T)%2==0 for T in hosts[n]),'even-order isolate obstruction')

# Deletion rows and source evaluation are independently obtained by testing
# graph isomorphism through literal bijection enumeration, not canonical masks.
def class_id(n,T):
    for i,S in enumerate(classes[n]):
        if embed(n,S,T):return i
    raise RuntimeError('missing class')

D=[];H=[[]]
for n in range(4):
    D.append([sum(1<<i for i in range(len(classes[n]))
                  if sum(class_id(n,remove(n+1,T,v))==i for v in range(n+1))%2)
              for T in classes[n+1]])
for n in range(1,5):H.append([1<<class_id(n,source(n-1,T)) for T in classes[n-1]])
def apply(matrix,vector):return sum(((r&vector).bit_count()%2)<<i for i,r in enumerate(matrix))
ranks=[]
for n in range(4):
    images={apply(D[n],v) for v in range(1<<len(classes[n]))}
    ranks.append(len(images).bit_length()-1)
    for f in range(1<<len(classes[n])):
        dh=apply(D[n-1],apply(H[n],f)) if n else 0
        hd=apply(H[n+1],apply(D[n],f))
        ck(dh^hd==f,'contracting identity all functions')
        if n<3:ck(apply(D[n+1],apply(D[n],f))==0,'d squared all functions')
        if not apply(D[n],f):ck(dh==f,'explicit kernel primitive')
    ones=(1<<len(classes[n]))-1
    ck(apply(D[n],ones)==(((1<<len(classes[n+1]))-1) if n%2==0 else 0),'constant action')
    ck(ranks[n]==sum((-1)**(n-j)*len(classes[j]) for j in range(n+1)),'rank formula')
    if n%2:
        ck((1<<len(classes[n+1]))-1 not in images,'constant not in odd-degree image')
    else:
        for f in range(1<<len(classes[n])):
            df=apply(D[n],f)
            if df in (0,(1<<len(classes[n+1]))-1):
                c=bool(df)
                dh=apply(D[n-1],apply(H[n],f)) if n else 0
                ck(f==dh^(ones if c else 0),'constant decomposition')

# Direct positive-order sharpness controls, with every labelled tournament
# at orders 1,3,5; no use of Redei's theorem in these finite computations.
for m in (1,3,5):
    S=tuple((i,i+1) for i in range(m-2))
    for T in all_graphs(m,True):ck(embed(m,S,T)%2==1,'odd-order padded directed path')

print(json.dumps({'status':'PASS','checks':checks,'oriented_graph_counts':graph_counts,
                  'class_counts':[len(c) for c in classes],
                  'direct_ranks':ranks,'embedding_h_through_order':4,
                  'embedding_d_through_source_order':3,
                  'padded_path_all_labelled_host_orders':[1,3,5],
                  'scope':'Independent edge-set and bijection enumeration; no producer imports. Universal conclusions are supplied by the elementary proof.'},indent=2,sort_keys=True))
