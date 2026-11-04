#!/usr/bin/env python3
"""Independent endpoint-kernel replay of every square-sum deletion card.
Method inherited from S32, reimplemented independently of verify.py's DFS.
"""
from collections import Counter
from itertools import combinations
from math import isqrt
from pathlib import Path
import json

def graph(n,delete):
    V=set(range(1,n+1))-{delete}
    return {v:{w for w in V if v!=w and isqrt(v+w)**2==v+w} for v in V}

def comps(G):
    todo=set(G);parts=[]
    while todo:
        seen=set();work=[min(todo)]
        while work:
            v=work.pop()
            if v in seen:continue
            seen.add(v);work.extend(G[v]-seen)
        todo-=seen;parts.append(seen)
    return parts

def kernel_decision(G):
    if len(comps(G))!=1:return False
    leaves={v for v in G if len(G[v])==1}
    for ends in combinations(sorted(G),2):
        if not leaves<=set(ends):continue
        B={v for v in G if len(G[v])!=2}|set(ends)
        used=set();mandatory=[];optional=[]
        for a in sorted(B):
            for b in sorted(G[a]):
                edge=tuple(sorted((a,b)))
                if edge in used:continue
                chain=[a,b];used.add(edge)
                while chain[-1] not in B:
                    opts=G[chain[-1]]-{chain[-2]}
                    if len(opts)!=1:raise RuntimeError('degree2 chain')
                    c=next(iter(opts));used.add(tuple(sorted((chain[-1],c))));chain.append(c)
                (mandatory if len(chain)>2 else optional).append(chain)
        if len(used)!=sum(map(len,G.values()))//2:raise RuntimeError('lost edge')
        target={v:(1 if v in ends else 2) for v in B};deg=Counter(v for p in mandatory for v in (p[0],p[-1]))
        if any(deg[v]>target[v] for v in B):continue
        needed=len(G)-1-sum(len(p)-1 for p in mandatory)
        if not 0<=needed<=len(optional):continue
        for chosen in combinations(optional,needed):
            ds=deg.copy()
            for p in chosen:ds[p[0]]+=1;ds[p[-1]]+=1
            if any(ds[v]!=target[v] for v in B):continue
            H={v:set() for v in G}
            for p in mandatory+list(chosen):
                for x,y in zip(p,p[1:]):H[x].add(y);H[y].add(x)
            if len(comps(H))==1:return True
    return False

r=json.loads(Path(__file__).with_name('RESULTS.json').read_text());out=[]
for row in r['square_deletion_profiles']:
    n=row['n'];expected=set(map(int,row['restoring_deletions']));actual=set()
    if kernel_decision(graph(n,None)):raise RuntimeError('base not negative')
    for v in range(1,n+1):
        if kernel_decision(graph(n,v)):actual.add(v)
    if actual!=expected:raise RuntimeError(('deletion mismatch',n,actual,expected))
    out.append({'n':n,'restoring_deletions':sorted(actual),'all_cards_checked':n})
print(json.dumps({'status':'PASS','method':'independent endpoint-aware chain-kernel edge-subset search','base_graphs':6,'deletion_cards':sum(x['all_cards_checked'] for x in out),'results':out},indent=2))
