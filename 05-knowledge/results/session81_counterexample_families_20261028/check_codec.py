#!/usr/bin/env python3
from itertools import combinations, permutations
from math import comb
import json

def encode(n,mask):
    N=comb(n,2)+8*n+1
    adj=[set() for _ in range(N)]
    def edge(a,b): adj[a].add(b);adj[b].add(a)
    for v in range(n):
        c,t,h=9*v,9*v+1,9*v+2
        edge(c,t);edge(c,h);edge(c,9*v+3)
        for j in (4,5): edge(t,9*v+j)
        for j in (6,7,8): edge(h,9*v+j)
    for k,(u,v) in enumerate(combinations(range(n),2)):
        if mask>>k&1:u,v=v,u
        edge(9*u+1,9*v+2)
    return adj

def decode(adj):
    core=[v for v in range(len(adj)) if len(adj[v])>1]
    roles={v:sum(len(adj[w])==1 for w in adj[v]) for v in core}
    centers=sorted(v for v in core if roles[v]==1)
    tails={c:next(v for v in adj[c] if roles.get(v)==2) for c in centers}
    heads={c:next(v for v in adj[c] if roles.get(v)==3) for c in centers}
    arcs={(i,j) for i,c in enumerate(centers) for j,d in enumerate(centers)
          if i!=j and heads[d] in adj[tails[c]]}
    if len(arcs)!=comb(len(centers),2):raise ValueError('decode arc count')
    return centers,arcs

def canon_graph(n,edges):
    pairs=list(combinations(range(n),2));index={p:i for i,p in enumerate(pairs)}
    return min(sum(1<<index[tuple(sorted((p[u],p[v])))] for u,v in edges) for p in permutations(range(n)))
checks=0
for n in range(1,6):
    for mask in range(1<<comb(n,2)):
        adj=encode(n,mask);N=len(adj)
        if sum(map(len,adj))//2!=N-1:raise ValueError('sparse count')
        if sum(not x for x in adj)!=comb(n,2)-n+1:raise ValueError('isolate count')
        # Every existing edge joins intrinsically different role categories,
        # giving an automorphism-invariant orientation (hence sign-evenness).
        role={v:('leaf' if len(adj[v])==1 else 'isolate' if not adj[v] else str(sum(len(adj[w])==1 for w in adj[v]))) for v in range(N)}
        allowed={frozenset(x) for x in [('leaf','1'),('leaf','2'),('leaf','3'),('1','2'),('1','3'),('2','3')]}
        if any(frozenset((role[u],role[v])) not in allowed for u in range(N) for v in adj[u]):raise ValueError('intrinsic orientation')
        original={(v,u) if mask>>k&1 else (u,v) for k,(u,v) in enumerate(combinations(range(n),2))}
        for p in (list(range(N)),list(reversed(range(N))),[(i+7)%N for i in range(N)]):
            shuffled=[set() for _ in range(N)]
            for u in range(N):shuffled[p[u]]={p[v] for v in adj[u]}
            centers,arcs=decode(shuffled)
            recovered={tuple(p[9*i] for i in arc) for arc in original}
            if {(centers[i],centers[j]) for i,j in arcs}!=recovered:raise ValueError('faithfulness roundtrip')
            checks+=1
counts=[]
for n in range(1,6):
    types={canon_graph(n,e) for e in combinations(list(combinations(range(n),2)),n-1)}
    counts.append(len(types))
if counts!=[1,1,1,3,6]:raise ValueError('independent sparse census')
print(json.dumps({'roundtrip_checks':checks,'all_labelled_tournaments_n_le':5,'three_relabellings_per_tournament':True,'intrinsic_aut_invariant_edge_orientation':True,'direct_sparse_orbits_n_1_to_5':counts},indent=2))
