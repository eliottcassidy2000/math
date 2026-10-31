#!/usr/bin/env python3
"""Independent brute adjacency, deletion invariants, and graph certificates.
No imports of the graph constructor/search under audit.
"""
from collections import Counter,deque
from functools import lru_cache
from itertools import combinations
from math import isqrt
from pathlib import Path
import json,sys
from audit_fractions import require
ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent
D=json.loads((ROOT/'prime_graphs/graphs.json').read_text())
W=json.loads((ROOT/'prime_graphs/witnesses.json').read_text())


def graph(n):
    edges={(a,b) for a,b in combinations(range(1,n+1),2) if isqrt(a+b)**2==a+b}
    adj={v:set() for v in range(1,n+1)}
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    return edges,adj


def components(adj,removed=(),deleted_edge=None):
    todo=set(adj)-set(removed);result=[]
    while todo:
        v=todo.pop();component={v};q=[v]
        while q:
            a=q.pop()
            for b in adj[a]&todo:
                if deleted_edge and set((a,b))==set(deleted_edge):continue
                todo.remove(b);q.append(b);component.add(b)
        result.append(component)
    return result


def girths(adj):
    shortest=odd=None
    # The shortest odd closed walk is a shortest odd cycle.
    for start in adj:
        dist={(start,0):0};q=deque([(start,0)])
        while q:
            v,p=q.popleft()
            for w in adj[v]:
                key=(w,1-p)
                if key not in dist:
                    dist[key]=dist[v,p]+1;q.append(key)
        if (start,1) in dist:odd=min(odd or 10**9,dist[start,1])
    for a in adj:
        for b in adj[a]:
            if a>b:continue
            dist={a:0};q=deque([a])
            while q and b not in dist:
                v=q.popleft()
                for w in adj[v]:
                    if {v,w}=={a,b}:continue
                    if w not in dist:dist[w]=dist[v]+1;q.append(w)
            if b in dist:shortest=min(shortest or 10**9,dist[b]+1)
    return shortest,odd


def validate_walk(seq,n,cycle=False,spanning=False):
    require(all(1<=v<=n for v in seq),'Out-of-range witness')
    require(len(set(seq))==len(seq),'Repeated vertex')
    if spanning:require(set(seq)==set(range(1,n+1)),'Nonspanning witness')
    pairs=list(zip(seq,seq[1:]))+([(seq[-1],seq[0])] if cycle else [])
    require(all(a!=b and isqrt(a+b)**2==a+b for a,b in pairs),'Nonsquare edge witness')


for row in D['all_graphs_1_through_97']:
    n=row['n'];edges,adj=graph(n);cc=components(adj)
    require(edges==set(map(tuple,row['edge_list'])),'Incorrect edges '+str(n))
    require(len(edges)==row['edges'] and len(cc)==row['components'],'Graph count mismatch')
    require(sorted(map(len,cc),reverse=True)==row['component_sizes'],'Component sizes')
    deg=Counter(map(len,adj.values()))
    require(dict(deg)=={int(k):v for k,v in row['degree_histogram'].items()},'Degree histogram')
    cuts=[v for v in adj if len(components(adj,[v]))>len(cc)]
    bridges=[list(e) for e in sorted(edges) if len(components(adj,deleted_edge=e))>len(cc)]
    require(cuts==row['articulation_points'] and bridges==row['bridges'],'Deletion invariant mismatch')
    require(girths(adj)==(row['girth'],row['odd_girth']),'Girth mismatch')
    oddadj={v:{w for w in neighbors if isqrt(v+w)%2} for v,neighbors in adj.items()}
    require(len(components(oddadj))==row['odd_root_components'],'Odd-root graph mismatch')
    require(sorted({isqrt(n+v) for v in adj[n]})==row['arrival_roots'],'Arrival roots mismatch')
for row in W['finite_hamilton']:
    for name in ['path','cycle']:
        witness=row[name]['witness']
        if witness:validate_walk(witness,row['n'],name=='cycle',True)
for row in W['prime_cycles']:validate_walk(row['cycle'],row['n'],True,True)
validate_walk(W['first_triangle'],30,True)
validate_walk(W['first_odd_cycle'],15,True)
# Ear decomposition certifies 2-connectivity of a spanning G31 subgraph.
base,*ears=W['ears31']
require(base[0]==base[-1],'Ear base not closed')
validate_walk(base[:-1],31,True);used=set(base)
for ear in ears:
    validate_walk(ear,31)
    require(ear[0] in used and ear[-1] in used and ear[0]!=ear[-1],'Ear endpoints')
    require(not(set(ear[1:-1])&used),'Ear interior not fresh')
    used.update(ear)
require(used==set(range(1,32)),'Ears not spanning')
# Supply-list consistency, with no guessed repair.
ps=[p for p in range(2,98) if all(p%d for d in range(2,isqrt(p)+1))]
supplied=[p for g in D['groups_as_supplied'].values() for p in g]
require(len(ps)==25 and ps==D['first_25_primes'],'Prime labels')
require(len(supplied)==25 and len(set(supplied))==24 and Counter(supplied)[71]==2 and set(ps)-set(supplied)=={67},'Supply duplicate/missing')
print(json.dumps(dict(status='PASS',graphs=97,prime_cycle_certificates=len(W['prime_cycles']),
                     ear_decomposition_vertices=len(used),duplicate=71,missing=67,
                     scope='Hamiltonian witness existence independently checked; exact Hamiltonian counts audited separately by source review.'),indent=2))
