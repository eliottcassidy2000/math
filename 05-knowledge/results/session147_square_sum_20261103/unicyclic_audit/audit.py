#!/usr/bin/env python3
"""Independent audit: Prüfer generation, cut tests, Hamilton subset DP."""
from itertools import combinations, product
from math import isqrt
from collections import Counter
from pathlib import Path
import importlib.util

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def tree_from_prufer(code, n):
    deg=[1]*n
    for v in code: deg[v]+=1
    edges=[]
    for v in code:
        leaf=next(i for i,d in enumerate(deg) if d==1)
        edges.append((min(leaf,v),max(leaf,v)))
        deg[leaf]-=1; deg[v]-=1
    u,v=(i for i,d in enumerate(deg) if d==1)
    edges.append((u,v))
    return frozenset(edges)

def components(n, E):
    sets=[{i} for i in range(n)]
    for u,v in E:
        a=next(s for s in sets if u in s); b=next(s for s in sets if v in s)
        if a is not b:
            a.update(b); sets.remove(b)
    return sorted((tuple(sorted(s)) for s in sets), key=lambda s:s[0])

def pathsets(n,E):
    # DP by visited bitset and last vertex, storing edge sets; no permutations.
    states={(1<<v,v):{frozenset()} for v in range(n)}
    for size in range(1,n):
        nxt={}
        for (mask,last), paths in states.items():
            for u,v in E:
                if u==last: w=v
                elif v==last: w=u
                else: continue
                if mask>>w&1: continue
                bucket=nxt.setdefault((mask|1<<w,w),set())
                edge=(u,v)
                bucket.update(p|{edge} for p in paths)
        states=nxt
    return set().union(*states.values()) if states else set()

def verify(n,E):
    # Independently identify cycle edges by connectivity after single-edge removal.
    require(len(E)==n and len(components(n,E))==1, 'not connected unicyclic')
    C={e for e in E if len(components(n,E-{e}))==1}
    Vc=set().union(*(set(e) for e in C))
    deg=[sum(v in e for e in E) for v in range(n)]
    S={v for v in Vc if deg[v]==3}
    bad=any(deg[v]>=4 if v in Vc else deg[v]>=3 for v in range(n))
    D=set() if bad else {e for e in C if S.issubset(e)}
    closed=0 if bad else len(C) if not S else 2 if len(S)==1 else 1 if len(S)==2 and tuple(sorted(S)) in C else 0
    paths=pathsets(n,E)
    require(len(D)==closed and paths=={frozenset(E-{e}) for e in D},f'path mismatch n={n}, E={sorted(E)}')
    for p in paths:
        require(E-C<=p and len(E-p)==1 and next(iter(E-p)) in C, f'bridge or omission mismatch {E}')
    for r in range(n+1):
        for F0 in combinations(E,r):
            F=set(F0)
            require(len(components(n,E-F))==1+r-int(bool(F&C)),f'component count mismatch {E=} {F=}')
    return len(paths),len(C),len(S),bad

def expected_rejection(call, label):
    try: call()
    except (ValueError,RuntimeError,AssertionError): return
    raise RuntimeError('expected rejection did not occur: '+label)

def main():
    for n in range(3,7):
        allE=set(combinations(range(n),2)); graphs=set()
        for code in product(range(n),repeat=n-2):
            T=tree_from_prufer(code,n)
            for e in allE-T: graphs.add(T|{e})
        dist=Counter(); cycles=Counter(); bads=0
        for E in graphs:
            h,c,s,bad=verify(n,E); dist[h]+=1; cycles[c]+=1; bads+=bad
        expected={3:1,4:15,5:222,6:3660}[n]
        require(len(graphs)==expected, f'count mismatch n={n}')
        print(f'n={n}: graphs={len(graphs)} path_counts={dict(sorted(dist.items()))} cycle_lengths={dict(sorted(cycles.items()))} degree_obstructions={bads}')
    # Ensure intentional errors and upstream input rejection still fire under -O.
    expected_rejection(lambda: require(False,'intentional failure'), 'explicit require')
    spec=importlib.util.spec_from_file_location('original',Path(__file__).resolve().parent.parent / 'unicyclic/check.py')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    invalid=[(2,[(0,1)]),(3,[(0,1),(1,2)]),(3,[(0,0),(0,1),(1,2)]),(3,[(0,1),(1,0),(1,2)]),(3,[(0,1),(1,2),(0,3)]),(6,[(0,1),(1,2),(0,2),(3,4),(4,5),(3,5)])]
    for n,edges in invalid: expected_rejection(lambda:mod.theorem(n,edges),f'original validator {n} {edges}')
    expected_rejection(lambda: verify(4,frozenset({(0,1),(1,2),(2,3)})), 'non-unicyclic graph')
    print('explicit failure controls: intentional false condition, six malformed/disconnected upstream cases, and a non-unicyclic graph all rejected')
    for n in range(1,20):
        E=frozenset((i,j) for i,j in combinations(range(1,n+1),2) if isqrt(i+j)**2==i+j)
        E0=frozenset((i-1,j-1) for i,j in E)
        comps=components(n,E0)
        neighbors=sorted(i for i in range(1,n) if isqrt(i+n)**2==i+n)
        require(len(neighbors)==isqrt(2*n-1)-isqrt(n),f'd(n) formula failed {n}')
        rank=len(E)-n+len(comps)
        print(f'G{n:02}: edges={len(E):2} comps={len(comps)} rank={rank} older={neighbors}')
        if n<=13 and n>=2: require(len(comps)>1,f'G{n} unexpectedly connected')
        if n==14: require(len(comps)==1 and rank==0,'G14 not tree')
        if 15<=n<=18: require(len(comps)==1 and rank==1,f'G{n} not unicyclic')
        if n==19: require(len(comps)==1 and rank==2,'G19 not rank 2')
    print('all independent checks passed')
if __name__=='__main__': main()
