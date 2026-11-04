#!/usr/bin/env python3
"""Independent audit: edge subsets, all-start vertex search, and frozen inputs.
Only transfer functions under test are imported from the producer. Square-sum
adjacency, path counting, forest extraction, witnesses, and checks are local.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
from math import isqrt
import importlib.util
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
P = ROOT / 'session190-squareports'
sys.path.insert(0, str(P))
spec = importlib.util.spec_from_file_location('tested_transfer', P / 'verify.py')
tested = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tested)

def check(condition, label):
    if not condition:
        raise RuntimeError(label)

def graph(V, E):
    G = {v: set() for v in V}
    for a,b in E:
        G[a].add(b)
        G[b].add(a)
    return G

def square(n, deleted=()):
    V = set(range(1,n+1)) - set(deleted)
    E = set()
    for q in range(2, isqrt(2*n)+1):
        for a in V:
            b = q*q-a
            if b in V and a < b:
                E.add((a,b))
    return graph(V,E)

def subsets(G):
    """No union-find and no transfer: traverse each selected edge graph."""
    V = sorted(G)
    E = [(a,b) for a in V for b in sorted(G[a]) if a < b]
    result = Counter()
    for mask in range(1 << len(E)):
        F = {v: [] for v in V}
        for i,(a,b) in enumerate(E):
            if mask & (1 << i):
                F[a].append(b)
                F[b].append(a)
        if any(len(F[v]) > 2 for v in V):
            continue
        unseen = set(V)
        pairs = []
        while unseen:
            first = min(unseen)
            stack = [first]
            comp = set()
            while stack:
                v=stack.pop()
                if v in comp:
                    continue
                comp.add(v)
                stack.extend(F[v])
            unseen -= comp
            ends = sorted(v for v in comp if len(F[v]) < 2)
            if len(comp) == 1:
                pairs.append((first,first))
            elif len(ends) == 2:
                pairs.append(tuple(ends))
            else:
                break
        else:
            result[tuple(sorted(pairs))] += 1
    return result

def paths(G):
    """Enumerate all ordered vertex paths; choose reversal by numeric endpoints.
    No forced-leaf start, subset memoization, or imported counting/pruning code.
    """
    V=sorted(G)
    C=Counter()
    if not V:
        return C
    if len(V)==1:
        return Counter({(V[0],V[0]):1})
    bit = {v:1<<i for i,v in enumerate(V)}
    allmask = (1<<len(V))-1
    A = {v:tuple((w,bit[w]) for w in sorted(G[v])) for v in V}
    for start in V:
        def visit(v, used):
            if used == allmask:
                if start < v:
                    C[start,v] += 1
                return
            for w,b in A[v]:
                if not used & b:
                    visit(w,used | b)
        visit(start,bit[start])
    return C

def decode(entries):
    return Counter({tuple(x['endpoints']):x['count'] for x in entries})

def check_path(G,path):
    check(len(path)==len(G) and set(path)==set(G), 'witness vertex set')
    check(all(b in G[a] for a,b in zip(path,path[1:])), 'witness edge')

def main():
    manifest = json.loads((P/'FROZEN_MANIFEST.json').read_text())
    for name,h in manifest['files'].items():
        check(hashlib.sha256((P/name).read_bytes()).hexdigest()==h, 'frozen hash '+name)
    R={'status':'PASS','frozen_manifest_sha256':hashlib.sha256((P/'FROZEN_MANIFEST.json').read_bytes()).hexdigest()}
    total=ports=additions=deletions=0
    profile_cache={}
    for n in range(6):
        E=list(combinations(range(n),2))
        for mask in range(1<<len(E)):
            G=graph(range(n),[e for i,e in enumerate(E) if mask & (1<<i)])
            F=subsets(G)
            check(F==tested.full_profile(G),('all-component profile',n,mask))
            check(paths(G)==tested.solve(G)[0],('pruned endpoint profile',n,mask))
            total+=1
            if n <= 4:
                profile_cache[n,mask] = (G,F)
            if n == 5:
                pairings=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
                X=[F[p] for p in pairings]
                Y=[]
                for chosen,ends in pairings:
                    H={v:set(nb) for v,nb in G.items()}
                    H[5]=set(chosen)
                    for v in chosen:H[v].add(5)
                    Y.append(paths(H)[ends])
                check(Y==[X[1]+X[2],X[0]+X[2],X[0]+X[1]],('internal-vertex probes',mask))
                ports+=3
    for (n,mask),(G,F) in profile_cache.items():
        for sm in range(1<<n):
            S={v for v in G if sm & (1<<v)}
            H={v:set(nb) for v,nb in G.items()}
            H[-1]=S.copy()
            for v in S:H[v].add(-1)
            check(tested.full_transfer(F,-1,S)==subsets(H),('addition',n,mask,sm))
            additions+=1
            for deleted in G:
                old={v:nb-{deleted} for v,nb in G.items() if v!=deleted}
                new={v:nb-{deleted} for v,nb in H.items() if v!=deleted}
                check(tested.full_transfer(subsets(old),-1,S-{deleted})==subsets(new),('deletion',n,mask,sm,deleted))
                deletions+=1
    R['exhaustive_graphs_including_empty']=total
    R['arbitrary_label_arrivals']=additions
    R['deletion_commutation_cases']=deletions
    R['four_port_queries_with_fifth_internal_vertex']=ports
    single=Counter({((0,0),):1})
    check(tested.full_transfer(single,1,{0})==Counter({((0,0),(1,1)):1,((0,1),):1}),'singleton multiplicity')
    two=Counter({((0,0),(1,1)):1})
    T=tested.full_transfer(two,2,{0,1})
    check(sum(T.values())==4 and T[((0,1),)]==1,'two-singleton multiplicity')
    edge=Counter({((0,1),):1})
    T=tested.full_transfer(edge,2,{0,1})
    check(sum(T.values())==3 and ((0,1),) not in T,'same-component cycle rejection')
    R['singleton_and_cycle_controls']=3
    for r in range(1,11):
        V=range(2*r+2)
        old=graph(V,[(2*i,2*i+1) for i in range(r+1)])
        C=subsets(old)
        check(min(map(len,C))==r+1,'horizon lower bound')
        C=Counter({s:c for s,c in C.items() if len(s)<=r+1})
        for j in range(1,r+1):
            C=tested.full_transfer(C,2*r+1+j,{2*j-1,2*j})
            C=Counter({s:c for s,c in C.items() if len(s)<=r-j+1})
        check(C==Counter({((0,2*r+1),):1}),('strict horizon witness',r))
    R['strict_horizon_witnesses']=10
    data=json.loads((P/'RESULTS.json').read_text())
    for row in data['square_endpoint_census']:
        G=square(row['n'])
        check(paths(G)==decode(row['endpoint_pairs']),('square census',row['n']))
        check(row['path_count']==sum(decode(row['endpoint_pairs']).values()),'path count')
        for w in row['witnesses']:check_path(G,w['path'])
        if row['n']>1:
            n=row['n']
            check(set(row['arrival_neighbors'])=={q*q-n for q in range(2,isqrt(2*n-1)+1) if n<q*q},'arrival formula')
    R['independent_square_endpoint_matrices']=30
    cards=0
    for entry in data['square_deletion_census']:
        n=entry['n']
        for card in entry['cards']:
            G=square(n,[card['deleted']])
            check(paths(G)==decode(card['endpoint_pairs']),('square deletion',n,card['deleted']))
            if card['one_path'] is not None:check_path(G,card['one_path'])
            cards+=1
    R['independent_deletion_endpoint_matrices']=cards
    S=json.loads((P/'STRUCTURAL.json').read_text())
    check_path(square(25,[3]),S['Q25_minus3_path'])
    G=square(26,[3])
    for v,neighbors in S['Q26_minus3_adjacency'].items():
        check(G[int(v)]==set(neighbors),('structural local adjacency',v))
    # The human argument's exact degree obstruction, independently represented.
    check(G[18]=={7} and G[22]=={14},'forced endpoints')
    check(G[9]=={7,16} and G[13]=={12,23} and G[26]=={10,23},'forced degree-two vertices')
    forced={frozenset(e) for e in [(18,7),(22,14),(9,7),(9,16),(13,12),(13,23),(26,10),(26,23)]}
    forced_degree=Counter(v for e in forced for v in e)
    available_for_two={v for v in G[2] if forced_degree[v]<2}
    check(available_for_two=={14},'internal vertex 2 has at most one usable neighbor')
    R['human_certificate']='Q25-3 path and eight Q26-3 adjacency lists; endpoints 18/22 force saturation at 7/23, leaving internal 2 only neighbor 14'
    R['finite_rows_25_to_30']=[{k:r[k] for k in ('n','path_count','endpoint_pair_count','endpoint_extensions','two_component_joins')} for r in data['square_endpoint_census'] if r['n']>=25]
    print(json.dumps(R,sort_keys=True,indent=2))

if __name__=='__main__':main()
