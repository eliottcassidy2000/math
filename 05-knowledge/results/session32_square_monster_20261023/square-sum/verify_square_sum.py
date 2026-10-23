#!/usr/bin/env python3
"""Exact, dependency-free square-sum certificates. No external or Git writes.
Reproduce: python3 verify_square_sum.py > certificates.json
"""
from math import isqrt
from itertools import combinations
from collections import Counter
import json

def check(b, message):
    if not b: raise RuntimeError(message)

def graph(n):
    return {v: {w for w in range(1,n+1) if w != v and isqrt(v+w)**2 == v+w} for v in range(1,n+1)}

def components(G, removed=()):
    unseen=set(G)-set(removed); ans=[]
    while unseen:
        root=min(unseen); unseen.remove(root); todo=[root]; block=[]
        while todo:
            v=todo.pop(); block.append(v)
            found=G[v]&unseen; unseen-=found; todo.extend(found)
        ans.append(sorted(block))
    return ans

def paths(G):
    n=len(G)
    if n==1:return [(1,)]
    leaves=[v for v in G if len(G[v])==1]
    starts=[min(leaves)] if leaves else sorted(G)
    found=[]
    # Only the prescribed starting leaf is fixed, never arbitrary degree-2 vertices.
    def dfs(path, unseen):
        if not unseen:
            if leaves or path[0]<path[-1]:found.append(tuple(path))
            return
        v=path[-1]
        for w in sorted(G[v]&unseen):
            dfs(path+[w], unseen-{w})
    for s in starts:dfs([s],set(G)-{s})
    return found

def kernel(G, endpoints):
    boundary={v for v in G if len(G[v])!=2}|set(endpoints)
    used=set(); chains=[]
    for u in sorted(boundary):
        for v in sorted(G[u]):
            if frozenset((u,v)) in used:continue
            path=[u,v];used.add(frozenset((u,v)))
            while path[-1] not in boundary:
                w=next(x for x in G[path[-1]] if x!=path[-2])
                used.add(frozenset((path[-1],w)));path.append(w)
            chains.append(path)
    check(sum(len(p)-1 for p in chains)==sum(map(len,G.values()))//2,'kernel edges lost')
    return boundary, chains

def kernel_paths(G):
    """Independently enumerate edge subsets of the endpoint-marked kernel.
    Degree-2 nonendpoints force long chains. Direct kernel edges remain optional.
    """
    n=len(G)
    if n==1:return [(1,)]
    if len(components(G))!=1:return []
    leaves={v for v in G if len(G[v])==1}
    result=[]
    for pair in combinations(sorted(G),2):
        if not leaves <= set(pair):continue
        B,chains=kernel(G,pair)
        mandatory=[p for p in chains if len(p)>2]
        optional=[p for p in chains if len(p)==2]
        target={v:(1 if v in pair else 2) for v in B}
        degree=Counter()
        for p in mandatory:degree[p[0]]+=1;degree[p[-1]]+=1
        if any(degree[v]>target[v] for v in B):continue
        # Need exactly n-1 original edges. This also limits the subset search.
        need=n-1-sum(len(p)-1 for p in mandatory)
        if not 0<=need<=len(optional):continue
        for choice in combinations(optional,need):
            selected=mandatory+list(choice);d=degree.copy()
            for p in choice:d[p[0]]+=1;d[p[-1]]+=1
            if any(d[v]!=target[v] for v in B):continue
            H={v:set() for v in G}
            for p in selected:
                for a,b in zip(p,p[1:]):H[a].add(b);H[b].add(a)
            if len(components(H))!=1:continue
            path=[pair[0]]
            while len(path)<n:
                nxt=H[path[-1]]-set(path)
                check(len(nxt)==1,'non-path selected');path.append(next(iter(nxt)))
            result.append(tuple(path))
    return result

def canonical(p):return min(p,p[::-1])

arrival=[]
for n in range(1,37):
    G=graph(n);old=graph(n-1)
    oldC=components(old);touch=sum(bool(set(C)&G[n]) for C in oldC)
    d=isqrt(2*n-1)-isqrt(n)
    check(len(G[n])==d,'arrival formula')
    check(len(components(G))==len(oldC)+1-touch,'component recurrence')
    beta=sum(map(len,G.values()))//2-n+len(components(G))
    oldbeta=sum(map(len,old.values()))//2-(n-1)+len(oldC)
    check(beta==oldbeta+d-touch,'cycle-rank recurrence')
    arrival.append(dict(n=n,neighbors=sorted(G[n]),degree=d,touched_components=touch,components=len(components(G)),cycle_rank=beta))

cutsets={18:[7],19:[7],20:[5,7],21:[5,7],22:[5,7,14]}
cuts=[]
for n,S in cutsets.items():
    C=components(graph(n),S);check(len(C)>len(S)+1,'invalid cut witness')
    cuts.append(dict(n=n,cut=S,remaining_components=C))

census=[]
for n in range(1,26):
    G=graph(n);P=paths(G);K=kernel_paths(G)
    check({canonical(p) for p in P}=={canonical(p) for p in K},f'independent kernel mismatch {n}')
    check(len(P)==len(K),f'kernel multiplicity mismatch {n}')
    census.append(dict(n=n,count=len(P),endpoint_counts={str(k):v for k,v in sorted(Counter(tuple(sorted((p[0],p[-1]))) for p in P).items())},paths=[list(p) for p in P]))
check([x['n'] for x in census if x['count']]==[1,15,16,17,23,25],'path support')
check(census[22]['count']==3 and census[24]['count']==10,'known nontrivial counts')
for n,a,b in [(15,1,10),(23,2,13),(25,11,24)]:
    for p in paths(graph(n)):
        i=p.index(n);check(0<i<n-1 and {p[i-1],p[i+1]}=={a,b},'compulsory arrival ear')

# Explicit zero-to-two cut transitions and user component certificates.
check(components(graph(12))==[[1,3,6,8,10],[2,7,9],[4,5,11,12]],'G12')
check(components(graph(13))==[[1,3,4,5,6,8,10,11,12,13],[2,7,9]],'G13')
# The G24 endpoint gate's adjacency facts and final proper cycle.
G=graph(24)
expected={18:{7},9:{7,16},11:{5,14},22:{3,14},2:{7,14,23},20:{5,16},4:{5,12,21},8:{1,17},17:{8,19},19:{6,17},10:{6,15},21:{4,15},24:{1,12}}
for v,N in expected.items():check(G[v]==N,f'G24 adjacency {v}')
cycle24=[1,8,17,19,6,10,15,21,4,12,24]
check(all(b in G[a] for a,b in zip(cycle24,cycle24[1:]+cycle24[:1])),'proper cycle')
# Finite direct checks support but do not replace the all-n proof in README.
for n in range(1,100001):
    d=isqrt(2*n-1)-isqrt(n)
    check((d==2)==(13<=n<=15 or 19<=n<=32 or 36<=n<=40 or 49<=n<=50),'degree-two classification')
    if n>=5:check(d>=1,'new vertex attachment')

print(json.dumps(dict(status='FINITE-EXACT, independent DFS/kernel agreement through25',arrival=arrival,degree_two_births_all=[*range(13,16),*range(19,33),*range(36,41),49,50],cut_certificates=cuts,proper_cycle24=cycle24,census=census),indent=2))
