#!/usr/bin/env python3
"""Exact arithmetic enumeration, graph BFS, and signed-permutation metric checks.
Standard library only; every load-bearing check uses explicit exceptions.
"""
from collections import Counter, deque
from itertools import product
from math import factorial
from pathlib import Path
import hashlib,json

CHECKS=0

def check(condition, message):
    global CHECKS
    CHECKS+=1
    if not condition:
        raise RuntimeError(message)

def partition(n, edges):
    labels=list(range(n))
    for u,v in edges:
        a,b=labels[u],labels[v]
        if a==b:return None
        labels=[a if z==b else z for z in labels]
    return labels

def data(k,m):
    b=[i*(k+1) for i in range(1,k+1)]
    inv={x:i for i,x in enumerate(b,1)}
    fixed=[(inv[d],d+inv[d]) if d in inv else (0,d) for d in range(k+1,m+1)]
    p=partition(m+1,fixed)
    check(p is not None,'fixed forest acyclic')
    check(len(set(p))==k+1,'fixed component count')
    check(all((p[v]==p[0])==(v not in b) for v in range(m+1)),'hub exact')
    opts=[]
    for d in range(1,k+1):
        aa=[a for a in range(m+1-d) if p[a]!=p[a+d]]
        check(len(aa)==2*k,'actual nonloop choices')
        check(set(aa)=={a for x in b for a in (x-d,x)},'both signed contacts exact')
        opts.append(aa)
    return fixed,p,opts,b

def enumerate_arithmetic(k,m):
    fixed,p,opts,b=data(k,m)
    states=[]
    def visit(d,cur,code):
        if d==k+1:
            states.append(tuple(code));return
        for a in opts[d-1]:
            x,y=cur[a],cur[a+d]
            if x!=y:
                visit(d+1,[x if z==y else z for z in cur],code+[a])
    visit(1,p,[])
    check(len(states)==factorial(k)*2**k,'exact arithmetic state count')
    desc={}
    for code in states:
        # Independent terminal validation uses connectivity by adjacency BFS,
        # not the partition recursion.
        edges=fixed+[(a,a+d) for d,a in enumerate(code,1)]
        adj=[[] for _ in range(m+1)]
        for a,c in edges:adj[a].append(c);adj[c].append(a)
        seen={0};todo=[0]
        for v in todo:
            for w in adj[v]:
                if w not in seen:seen.add(w);todo.append(w)
        check(len(seen)==m+1 and len(edges)==m,'terminal tree')
        check(sorted(abs(a-c) for a,c in edges)==list(range(1,m+1)),'terminal colors')
        perm=[];sign=[]
        for d,a in enumerate(code,1):
            hit=[i for i,x in enumerate(b) if x in (a,a+d)]
            check(len(hit)==1,'one marked endpoint')
            i=hit[0];perm.append(i);sign.append(a==b[i])
        check(len(set(perm))==k,'permutation sidecar')
        desc[code]=(tuple(perm),tuple(sign))
    check(len(set(desc.values()))==len(states),'description injective')
    return states,desc

def graph(states,width):
    neighbors=[set() for _ in states]
    # Complete pair scan deliberately does not use signed-permutation moves.
    for i,x in enumerate(states):
        for j in range(i):
            if sum(a!=b for a,b in zip(x,states[j]))<=width:
                neighbors[i].add(j);neighbors[j].add(i)
    return neighbors

def component_sizes(neighbors):
    unseen=set(range(len(neighbors)));out=[]
    while unseen:
        seed=unseen.pop();todo=[seed]
        for a in todo:
            for c in neighbors[a]:
                if c in unseen:unseen.remove(c);todo.append(c)
        out.append(len(todo))
    return out

def predicted(x,y):
    p,e=x;t,h=y;k=len(p)
    inv={v:i for i,v in enumerate(t)}
    rel=[inv[v] for v in p];unseen=set(range(k));cycles=0
    while unseen:
        v=next(iter(unseen));cycles+=1
        while v in unseen:
            unseen.remove(v);v=rel[v]
    f=sum(p[i]==t[i] and e[i]!=h[i] for i in range(k))
    return k-cycles+(f+1)//2

def main():
    rows=[]
    for k in range(1,6):
        m=k*k+2*k;states,desc=enumerate_arithmetic(k,m)
        # A second arithmetic interval confirms extra hub labels leave the fiber.
        extra,extra_desc=enumerate_arithmetic(k,m+3)
        check(set(states)==set(extra),'padding preserves actual states')
        check(desc==extra_desc,'padding preserves sidecars')
        if k<=4:
            one=graph(states,1)
            sizes=component_sizes(one)
            check(sizes==[2**k]*factorial(k),'cube component count and size')
            check(all(len(ns)==k for ns in one),'one-color degree k')
            for i,ns in enumerate(one):
                for j in ns:
                    p,e=desc[states[i]];t,h=desc[states[j]]
                    check(p==t and sum(a!=b for a,b in zip(e,h))==1,'exact cube adjacency')
            two=graph(states,2)
            check(all(len(ns)==k+5*k*(k-1)//2 for ns in two),'width-two degree')
            diameter=0
            for i,code in enumerate(states):
                ds=[-1]*len(states);ds[i]=0;todo=deque([i])
                while todo:
                    a=todo.popleft()
                    for c in two[a]:
                        if ds[c]<0:ds[c]=ds[a]+1;todo.append(c)
                check(min(ds)>=0,'width-two connectivity')
                for j,got in enumerate(ds):
                    check(got==predicted(desc[code],desc[states[j]]),'exact all-pairs metric')
                diameter=max(diameter,max(ds))
            check(diameter==(1 if k==1 else k-1),'sharp diameter')
        else:
            sizes=None;diameter=None
        rows.append({'k':k,'m':m,'arithmetic_states':len(states),'one_color_components':factorial(k) if sizes is not None else 'not enumerated','two_color_diameter':diameter,'distance_pairs_checked':len(states)**2 if k<=4 else 0})
    # Large-parameter construction/rank checks, without exhaustive states.
    for k in range(1,31):
        for padding in (0,1,17):
            data(k,k*k+2*k+padding)
    files=['verify_factorial_fibers.py','GRACEFUL_FACTORIAL_FIBERS.md']
    hashes={f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files}
    receipt={'status':'PASS','exact_arithmetic_and_graph_checks':rows,'large_construction_k_max':30,'checks':CHECKS,'sha256':hashes,'limitations':'All-pairs graph distances independently enumerated only for k<=4; all-k theorem is proved in the note. No Graceful Tree Conjecture or global one-color closure is claimed.'}
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__':main()
