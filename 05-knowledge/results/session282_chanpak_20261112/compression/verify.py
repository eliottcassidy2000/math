#!/usr/bin/env python3
"""Literal SCC neutral-suffix compression checks; exact, assertion independent."""
from pathlib import Path
from itertools import combinations
import importlib.util,json
ROOT=Path(__file__).resolve().parents[4]
s=importlib.util.spec_from_file_location('c',ROOT/'05-knowledge/results/session275_exact_selection_20261111/switches/carrier.py');c=importlib.util.module_from_spec(s);s.loader.exec_module(c)
def check(p,m):
    if not p:raise RuntimeError(m)
def induce(T,ids):return tuple(sum(1<<j for j,v in enumerate(ids) if T[u]>>v&1) for u in ids)
def sccs(D):
    n=len(D);reach=[D[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    todo=set(range(n));parts=[]
    while todo:
        u=min(todo);part={v for v in todo if reach[u]>>v&1 and reach[v]>>u&1};parts.append(sorted(part));todo-=part
    parts.sort(key=lambda P:-sum(1 for v in range(n) if v not in P and D[P[0]]>>v&1))
    for i,P in enumerate(parts):
        for R in parts[i+1:]:check(all(D[u]>>v&1 for u in P for v in R),'SCC ordinal orientation')
    return parts

def compress(T):
    D=induce(T,list(range(1,len(T))));parts=sccs(D)
    incoming={v-1 for v in range(1,len(T)) if T[v]&1}
    last=max((j for j,P in enumerate(parts) if incoming.intersection(P)),default=-1)
    removed={v+1 for j,P in enumerate(parts) if j>last and len(P)==1 for v in P}
    ids=[i for i in range(len(T)) if i not in removed]
    return induce(T,ids),len(removed)

def control(T):
    p=c.profile(T);R,nremoved=compress(T);q=c.profile(R)
    check(all(p[k]==q[k] for k in ('H','h','V','e')),'preserved exact receiver')
    h,V=p['h'],p['V'];b=0 if h==1 else h
    if V==h:B=b+1
    else:
        check(V>=2*h,'positive receiver gap')
        B=b+2+((V-h)//h).bit_length()-1
    check(len(R)<=B,'global witness bound')
    c.hp.cache_clear()
    return nremoved,p['Q']!=q['Q']
count=removed=changed=0
for n in range(1,6):
    for T in c.tournaments(n):
        a,b=control(T);count+=1;removed+=bool(a);changed+=b
# Expose neutral suffix deletion in nontrivial rooted examples at larger order.
for T in list(c.tournaments(4)):
    for tails in (1,2,3):
        n=len(T);m=n+tails
        # Append whole transitive singleton SCC suffix, all root-outgoing.
        R=tuple(T[i]|(((1<<m)-1)^((1<<n)-1)) for i in range(n))+tuple(((1<<m)-1)^((1<<(i+1))-1) for i in range(n,m))
        a,b=control(R);check(a>=tails,'explicit neutral suffix removed');count+=1;removed+=bool(a);changed+=b
check(changed>0,'full carrier not preserved hostile')
print(json.dumps(dict(status='PASS',literal_tournaments=count,nontrivial_compressions=removed,Q_changed_controls=changed,max_literal_order=7),indent=2,sort_keys=True))
