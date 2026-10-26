#!/usr/bin/env python3
import json
from pathlib import Path

def rank(rows):
    basis={}
    for r in rows:
        while r:
            k=r.bit_length()-1
            if k in basis:r^=basis[k]
            else:basis[k]=r;break
    return len(basis)

def run(n):
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)];N=1<<len(pairs)
    E=set()
    for x in range(N):
        a=[[False]*n for _ in range(n)]
        for i in range(n-1):a[i][i+1]=True
        for k,(i,j) in enumerate(pairs):a[i][j]=bool(x>>k&1)
        for i in range(n):
            for j in range(i+1,n):a[j][i]=not a[i][j]
        for i in range(n):
            for j in range(i+2,n):
                p=list(range(i))+list(range(i+1,j+1))+[i]+list(range(j+1,n))
                if not all(a[p[k]][p[k+1]] for k in range(n-1)):continue
                y=sum(int(a[p[u]][p[v]])<<k for k,(u,v) in enumerate(pairs))
                if x!=y:E.add(tuple(sorted((x,y))))
    L=[0]*N
    for x,y in E:
        L[x]^=(1<<x)^(1<<y);L[y]^=(1<<x)^(1<<y)
    incidence_rank=rank([(1<<x)^(1<<y) for x,y in E]); lr=rank(L)
    d={'n':n,'N':N,'simple_edges':len(E),'incidence_rank_F2':incidence_rank,'class_count':N-incidence_rank,'laplacian_rank_F2':lr,'laplacian_nullity_F2':N-lr}
    if lr<incidence_rank and N<=64:
        # Find smallest connected component with an additional null vector.
        neigh=[set() for _ in range(N)]
        for x,y in E:neigh[x].add(y);neigh[y].add(x)
        seen=set()
        for x in range(N):
            if x in seen:continue
            stack=[x];comp={x};seen.add(x)
            while stack:
                for y in neigh[stack.pop()]:
                    if y not in seen:seen.add(y);comp.add(y);stack.append(y)
            comp=sorted(comp);ix={v:i for i,v in enumerate(comp)};M=[]
            for v in comp:
                M.append(sum(((L[v]>>u)&1)<<i for i,u in enumerate(comp)))
            if rank(M)<len(comp)-1:
                witness=None
                for w in range(1,(1<<len(comp))-1):
                    if all((r&w).bit_count()%2==0 for r in M):witness=[comp[i] for i in range(len(comp)) if w>>i&1];break
                d['hostile']={'component':comp,'edges':sorted((a,b) for a,b in E if a in ix),'laplacian_kernel_nonconstant_support':witness,'component_laplacian_rows':M}
                break
    return d
rows=[run(n) for n in range(1,7)]
Path('binary_laplacian_hostile.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
