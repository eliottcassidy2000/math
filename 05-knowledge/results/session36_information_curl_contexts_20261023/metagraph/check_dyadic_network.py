#!/usr/bin/env python3
"""Actual fixed-path line graph: exact blue bicycle/critical-group experiments.
No repository imports; colors use literal reversal symmetry, not node self-converseness.
"""
import itertools as it, json, time, hashlib
from collections import Counter
from pathlib import Path
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE=Path(__file__).resolve().parent

def need(ok,msg='check failed'):
    if not ok:raise RuntimeError(msg)

def rank2(M):
    rows=[sum((int(v)&1)<<j for j,v in enumerate(row)) for row in M]
    basis={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p in basis:row^=basis[p]
            else:basis[p]=row;break
    return len(basis)

def components(n,edges):
    nei=[set() for _ in range(n)]
    for i,j,w in edges:
        if w and i!=j:nei[i].add(j);nei[j].add(i)
    out=[];seen=set()
    for a in range(n):
        if a in seen:continue
        group=[];todo=[a];seen.add(a)
        while todo:
            i=todo.pop();group.append(i)
            for j in nei[i]-seen:seen.add(j);todo.append(j)
        out.append(sorted(group))
    return out

def lap(n,edges):
    L=sp.zeros(n)
    for i,j,w in edges:
        if i!=j:L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
    return L

def tree(L):return 1 if L.rows==1 else int(L[:-1,:-1].det(method='domain-ge'))
def val2(z):
    need(z>0)
    return (z&-z).bit_length()-1

def census(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    pairs=list(it.combinations(range(n),2));full=(1<<len(free))-1
    reflection=[free.index((n-1-j,n-1-i)) for i,j in free]
    cache={}
    def canon(a):
        raw=sum(a[i][j]<<k for k,(i,j) in enumerate(pairs))
        if raw in cache:return cache[raw]
        d=[sum(row) for row in a]
        blocks=[[v for v in range(n) if d[v]==s] for s in sorted(set(d))]
        best=None
        for chunks in it.product(*(it.permutations(b) for b in blocks)):
            p=tuple(v for chunk in chunks for v in chunk)
            code=sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs))
            if best is None or code<best:best=code
        key=tuple(sorted(d)),best;cache[raw]=key;return key
    nodes=[];colors=[];triangles=[]
    for x in range(full+1):
        a=[[0]*n for _ in range(n)]
        for i in range(n-1):a[i][i+1]=1
        for k,(i,j) in enumerate(free):
            u,v=(i,j) if x>>k&1 else (j,i);a[u][v]=1
        nodes.append(min(canon(a),canon(list(map(list,zip(*a))))))
        colors.append(all((x>>k&1)==(x>>rk&1) for k,rk in enumerate(reflection)))
        triangles.append(sum((a[i][j] and a[j][k] and a[k][i]) or (a[j][i] and a[k][j] and a[i][k]) for i,j,k in it.combinations(range(n),3)))
    keys=sorted(set(nodes));ix={key:i for i,key in enumerate(keys)};node=[ix[c] for c in nodes]
    blue=Counter();black=Counter();r=len(keys);ct={}
    for x in range(full+1):
        if node[x] in ct:need(ct[node[x]]==triangles[x])
        ct[node[x]]=triangles[x]
        if x<(full^x):
            edge=tuple(sorted((node[x],node[full^x])))
            (blue if colors[x] else black)[edge]+=1
            need(colors[x]==colors[full^x])
    need(all(v%2==0 for v in black.values()),'black twin')
    edgeb=[(i,j,w) for (i,j),w in blue.items()]
    edgek=[(i,j,w//2) for (i,j),w in black.items()]
    Lb=lap(r,edgeb);Lk=lap(r,edgek);L=Lb+2*Lk
    cb=components(r,edgeb);ck=components(r,edgek);cg=components(r,edgeb+edgek)
    need(len(cg)==1,'finite connectedness')
    bcycle=r-rank2(Lb.tolist())-len(cb)
    rankcritical=r-1-rank2(L[:-1,:-1].tolist())
    need(rankcritical==len(cb)-1+bcycle)
    # Independent F2 incidence rank formula.
    actualblue=[(i,j) for (i,j),w in blue.items() for _ in range(w)]
    incidence=[[int(v==i)^int(v==j) for i,j in actualblue] for v in range(r)]
    need(rank2(incidence)==r-len(cb))
    # Intersection dim Cut cap Cycle via rank(C)-rank(C C^T).
    need(bcycle==rank2(incidence)-rank2(Lb.tolist()))
    nonSC=sum(not any(colors[x] for x in range(full+1) if node[x]==v) for v in range(r))
    # Actual harmonic current obeys every kernel congruence.
    current=L*sp.Matrix([ct[v] for v in range(r)])
    need(all(sum(int(current[v]) for v in component)%2==0 for component in cb))
    if n%2:
        need(all(i!=j and (ct[i]-ct[j])%2==1 for i,j,w in edgeb),'odd-order bipartite')
    out=dict(n=n,tilings=full+1,nodes=r,blue_lines=sum(blue.values()),black_lines=sum(black.values()),blue_components=len(cb),black_components=len(ck),blue_bicycle_dimension=bcycle,critical_two_rank=rankcritical,nonselfconverse_merged_nodes=nonSC,blue_component_sizes=sorted(map(len,cb)),blue_edges=edgeb,paired_black_edges=edgek,node_cyclic_triangles=ct)
    if n<=6:
        tau=tree(L);snf=smith_normal_form(L[:-1,:-1],domain=sp.ZZ)
        diag=[abs(int(snf[i,i])) for i in range(r-1)]
        need(all(z>0 for z in diag))
        need(sp.prod(diag)==tau)
        need(sum(z%2==0 for z in diag)==rankcritical)
        need(val2(tau)>=rankcritical)
        out.update(spanning_trees=tau,tree_valuation2=val2(tau),smith_invariants=diag)
    return out

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--max-n',type=int,default=6);args=p.parse_args()
    out=[]
    for n in range(3,args.max_n+1):
        t=time.time();a=census(n);out.append(a)
        print({k:v for k,v in a.items() if k not in ('blue_edges','paired_black_edges','node_cyclic_triangles','blue_component_sizes','smith_invariants')},'seconds',round(time.time()-t,2),flush=True)
        (HERE/'dyadic_network_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print('ALL CHECKS PASSED',flush=True)
