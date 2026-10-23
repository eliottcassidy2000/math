#!/usr/bin/env python3
"""Same-agent independent construction audit; imports no sibling checker."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations,product
from pathlib import Path
import json

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def positive_semidefinite(B):
    B=[list(row) for row in B]
    for k in range(len(B)):
        pivot=B[k][k]
        if pivot<0:return False
        if not pivot:
            if any(B[j][k] for j in range(k+1,len(B))):return False
            continue
        for i in range(k+1,len(B)):
            for j in range(k+1,len(B)):
                B[i][j]-=B[i][k]*B[k][j]/pivot
    return True

def iso(B,C,a,b):
    degreesB=list(map(sum,B));degreesC=list(map(sum,C));mapping={a:a,b:b}
    if any(degreesB[u]!=degreesC[v] for u,v in mapping.items()):return False
    if any(B[u][v]!=C[mapping[u]][mapping[v]] for u in mapping for v in mapping):return False
    def extend():
        if len(mapping)==8:return True
        candidates={u:[v for v in range(8) if v not in mapping.values()
                        and degreesB[u]==degreesC[v]
                        and all(B[u][j]==C[v][mapping[j]] for j in mapping)]
                    for u in range(8) if u not in mapping}
        u=min(candidates,key=lambda x:len(candidates[x]))
        for v in candidates[u]:
            mapping[u]=v
            if extend():return True
            del mapping[u]
        return False
    return extend()

def main():
    # Decode the inherited canonical integer, instead of copying the sibling matrix.
    A=[[0]*8 for _ in range(8)]
    for k,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(29558785>>(27-k))&1;A[j][i]=1-A[i][j]
    rows=[];automorphisms=0
    for p in permutations(range(8)):
        automorphisms+=all(A[u][v]==A[p[u]][p[v]] for u in range(8) for v in range(8))
        if not all(A[u][v] for u,v in zip(p,p[1:])):continue
        pos={u:i for i,u in enumerate(p)}
        K=[[(0 if u==v else A[u][v]^(abs(pos[u]-pos[v])!=1))
            for v in range(8)] for u in range(8)]
        rows.append((p,K))
    need(len(rows)==409 and automorphisms==1,'independent canonical integer source census')
    groups=defaultdict(list)
    for p,K in rows:
        key=(p[0],p[-1],tuple(sorted(map(sum,K))))
        for group in groups[key]:
            if iso(K,group[0][1],p[0],p[-1]):group.append((p,K));break
        else:groups[key].append([(p,K)])
    actual={frozenset(p for p,K in group) for bucket in groups.values() for group in bucket}
    saved=json.loads((Path(__file__).parent/'atlas_receipt.json').read_text())
    expected={frozenset(tuple(p) for p in chart['paths']) for chart in saved['charts_data']}
    need(actual==expected,'direct degree-pruned pointed isomorphisms reproduce entire atlas')
    need(Counter(map(len,actual))=={1:320,2:38,3:3,4:1},'atlas exact histogram')
    weighted=[]
    for r in range(1,5):
        N=1<<r;w=[1+(x%5) for x in range(N)];Z=sum(w);pi=[F(x,Z) for x in w]
        K=[[F(0) for _ in range(N)] for _ in range(N)]
        for x in range(N):
            for j in range(r):
                y=x^(1<<j);K[x][y]=F(1,2*r)*min(F(1),F(w[y],w[x]))
            K[x][x]=1-sum(K[x])
        need(all(pi[x]*K[x][y]==pi[y]*K[y][x] for x in range(N) for y in range(N)), 'actual weighted detailed balance')
        gap=F(min(w),r*max(w))
        B=[[(pi[x] if x==y else F(0))-pi[x]*K[x][y]
            -gap*((pi[x] if x==y else F(0))-pi[x]*pi[y])
            for y in range(N)] for x in range(N)]
        need(positive_semidefinite(B),'exact rational weighted Poincare matrix')
        weighted.append(dict(r=r,kappa=str(F(max(w),min(w))),lower_gap=str(gap)))
    # Opposite corners with heavy weights show why an unrestricted weight bound fails.
    eps=F(1,1000);w=[F(1),eps,eps,F(1)];f=[1,0,0,-1];Z=sum(w)
    energy=sum(min(w[x],w[x^(1<<j)])*F((f[x]-f[x^(1<<j)])**2,8*Z)
               for x in range(4) for j in range(2))
    variance=sum(w[x]*f[x]*f[x]/Z for x in range(4))
    need(energy/variance==eps/2,'small intermediate mass gives vanishing gap upper bound')
    mixed=[]
    core_paths=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2));g=(4,0,2,3,1,5,6,7)
    for r,M in product((1,2),(3,9,15)):
        offset=8*r;n=offset+2*M+2
        def arc(x,y):
            if x==y:return 0
            if x<offset or y<offset:
                return A[x%8][y%8] if x//8==y//8 else int(x<y)
            x-=offset;y-=offset
            if x==0 or y==2*M+1:return 1
            if y==0 or x==2*M+1:return 0
            cx,ux=divmod(x-1,M);cy,uy=divmod(y-1,M)
            if cx==cy:return int(0<(uy-ux)%M<=(M-1)//2)
            return int((ux-uy)%3==0) if cx==0 else int((ux-uy)%3!=0)
        B=[[arc(x,y) for y in range(n)] for x in range(n)]
        def make_path(bits,u,v):
            a=(u+3*v)%M;b=u
            return tuple(8*i+j for i in range(r) for j in core_paths[bits>>i&1])+(offset,)+tuple(offset+1+(a+j)%M for j in range(M))+tuple(offset+1+M+(b+2+j)%M for j in range(M))+(n-1,)
        base=make_path(0,0,0);bp={x:i for i,x in enumerate(base)}
        K=[[0 if x==y else B[x][y]^(abs(bp[x]-bp[y])!=1) for y in range(n)] for x in range(n)]
        states=0
        for bits,u,v in product(range(1<<r),range(M),range(M//3)):
            p=make_path(bits,u,v);positions={x:i for i,x in enumerate(p)}
            need(len(positions)==n and all(B[x][y] for x,y in zip(p,p[1:])), 'mixed actual Hamiltonian word')
            h=list(range(n))
            for i in range(r):
                if bits>>i&1:
                    for x in range(8):h[8*i+x]=8*i+g[x]
            for x in range(M):
                h[offset+1+x]=offset+1+(x+u+3*v)%M
                h[offset+1+M+x]=offset+1+M+(x+u)%M
            need(all(K[x][y]==(B[h[x]][h[y]]^(abs(positions[h[x]]-positions[h[y]])!=1))
                     for x in range(n) for y in range(n) if x!=y), 'mixed cross-orbit and retained-kernel partner transport')
            need(h[base[0]]==p[0] and h[base[-1]]==p[-1],'mixed actual boundary endpoints')
            states+=1
        mixed.append(dict(r=r,M=M,owners=n,states=states,gap_bound=str(F(2,2*r+M*M))))
    return dict(status='PASS: same-agent independent integer/isomorphism/rational-matrix paths',
                pointed_charts=len(actual),weighted=weighted,mixed_kernel_cases=mixed,
                weight_hostile_rayleigh=str(energy/variance),checks=CHECKS)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
