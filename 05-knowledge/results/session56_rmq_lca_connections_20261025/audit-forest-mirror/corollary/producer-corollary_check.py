#!/usr/bin/env python3
"""Independent finite parity and weighted-complement checks; no dependencies."""
from collections import Counter
from itertools import combinations
import json

def require(ok, detail):
    if not ok: raise RuntimeError(detail)

def inspect(n, edges):
    m=len(edges); deg=[0]*n
    for u,v in edges: deg[u]+=1; deg[v]+=1
    counts=Counter(); forest_total=0
    for mask in range(1<<m):
        score=[0]*n; parent=list(range(n)); acyclic=True
        def root(u):
            while parent[u]!=u: u=parent[u]
            return u
        for j,(u,v) in enumerate(edges):
            score[v if mask>>j&1 else u]+=1
            if mask>>j&1:
                a,b=root(u),root(v)
                if a==b: acyclic=False
                else: parent[a]=b
        counts[tuple(score)]+=1
        forest_total+=acyclic
    all_even=all(d%2==0 for d in deg)
    require(forest_total==len(counts),('forest/score identity',n,edges))
    require((forest_total%2==1)==all_even,('forest parity',n,edges))
    fixed=[]; pairs=set(); weighted_edges=0
    for r,w in counts.items():
        mirror=tuple(d-a for d,a in zip(deg,r))
        require(counts[mirror]==w,('multiplicity symmetry',n,edges,r))
        if mirror==r:
            fixed.append(r)
            if m:
                require(w%2==0,('central parity',n,edges))
                weighted_edges+=w//2
        elif r<mirror:
            pairs.add((r,mirror)); weighted_edges+=w
    require(len(fixed)==int(all_even),('central score count',n,edges))
    if m: require(weighted_edges==1<<(m-1),('native edge total',n,edges))
    else: require(len(counts)==1 and not pairs,('edgeless case',n,edges))
    return {'vertices':n,'free_edges':m,'forests':forest_total,
            'fixed_score_vectors':len(fixed),'coarse_pairs':len(pairs),
            'native_distinct_state_pairs':weighted_edges}

def main():
    checked=0
    for n in range(1,6):
        possible=list(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            inspect(n,[e for j,e in enumerate(possible) if mask>>j&1]); checked+=1
    require(checked==1099,'graph census size')
    path=[]
    for n in range(3,8):
        row=inspect(n,[(u,v) for u,v in combinations(range(n),2) if v-u>1])
        require(row['fixed_score_vectors']==0,('path has fixed score',n))
        require(row['coarse_pairs']*2==row['forests'],('path matching',n))
        path.append(row)
    require([r['forests'] for r in path]==[2,8,52,478,5706],'path counts')
    print(json.dumps({'status':'all checks passed','simple_graphs_checked':checked,'path_cases':path},indent=2))

if __name__=='__main__':main()
