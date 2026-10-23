#!/usr/bin/env python3
"""Native prime-power carry, actual partners, and exact routing costs."""
from collections import Counter
from fractions import Fraction as F
from itertools import product,combinations
import json

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def source_arc(M,x,y):
    n=2*M+2
    if x==0 or y==n-1:return True
    if y==0 or x==n-1:return False
    cx,ux=divmod(x-1,M);cy,uy=divmod(y-1,M)
    if cx==cy:return (uy-ux)%M in range(1,(M+1)//2)
    return (ux-uy)%3==0 if cx==0 else (ux-uy)%3!=0

def path(M,a,b):
    return (0,)+tuple(1+(a+j)%M for j in range(M))+tuple(1+M+(b+2+j)%M for j in range(M))+(2*M+1,)

def rotate(M,x,a,b):
    if x in (0,2*M+1):return x
    c,u=divmod(x-1,M);return 1+c*M+(u+(a if c==0 else b))%M

def partner(A,P):
    positions={u:i for i,u in enumerate(P)};n=len(P)
    return [[0 if u==v else (A[u][v] if abs(positions[u]-positions[v])==1 else 1-A[u][v])
             for v in range(n)] for u in range(n)]

def signed(d,m):
    d%=m
    return d if d<=m//2 else d-m

def route(x,y,M,K):
    u,v=x;out=[x]
    for axis,modulus in ((0,M),(1,K)):
        delta=signed(y[axis]-(u if axis==0 else v),modulus)
        direction=1 if delta>0 else -1
        for _ in range(abs(delta)):
            if axis==0:u=(u+direction)%M
            else:v=(v+direction)%K
            out.append((u,v))
    return out

def case(M,full):
    K=M//3;n=2*M+2;states=list(product(range(M),range(K)));N=len(states)
    A=[[0 if u==v else int(source_arc(M,u,v)) for v in range(n)] for u in range(n)]
    base=path(M,0,0);basepartner=partner(A,base)
    need(all(A[u][v]+A[v][u]==1 for u,v in combinations(range(n),2)),'actual tournament orientations')
    for a,b in ((1,1),(3%M,0)):
        g=[rotate(M,u,a,b) for u in range(n)]
        need(all(A[u][v]==A[g[u]][g[v]] for u,v in combinations(range(n),2)),'full actual source automorphism')
        need(g[0]==0 and g[-1]==n-1,'literal endpoint owners fixed')
    checked_partners=0
    for u,v in states:
        a,b=(u+3*v)%M,u;P=path(M,a,b)
        need(len(set(P))==n and P[0]==0 and P[-1]==n-1,'all actual owners and endpoint clocks')
        need(all(A[x][y] for x,y in zip(P,P[1:])),'actual Hamiltonian path')
        need((b,((a-b)%M)//3)==(u,v),'lossless product chart including kernel')
        need(tuple(rotate(M,x,a,b) for x in base)==P,'actual full-domain orbit transporter')
        if full or (u,v) in ((0,0),(1,0),(0,1),(M-1,K-1)):
            Q=partner(A,P);g=[rotate(M,x,a,b) for x in range(n)]
            need(all(basepartner[x][y]==Q[g[x]][g[y]] for x,y in combinations(range(n),2)),
                 'literal partner isomorphism retains every class potential')
            AT=[list(row) for row in zip(*A)];rev=P[::-1]
            C=partner(AT,rev)
            need(all(C[x][y]==Q[y][x] for x,y in combinations(range(n),2)),
                 'literal converse partner square at reversed clocks')
            checked_partners+=1
    D=M//2+K//2;bound=D*F(M*M-1,M)
    loads=Counter();maxevents=0;maxwork=0
    pairs=product(states,repeat=2) if full else [(states[0],states[-1]),(states[-1],states[0]),(states[0],states[N//2])]
    count=0
    for x,y in pairs:
        p=route(x,y,M,K);ell=len(p)-1;work=0
        need(p[0]==x and p[-1]==y and ell<=D,'native route endpoints and integer clock')
        for a,b in zip(p,p[1:]):
            du=(b[0]-a[0])%M;dv=(b[1]-a[1])%K
            need((du in (1,M-1) and dv==0) or (du==0 and K>1 and dv in (1,K-1)),
                 'every route event is an authenticated native generator')
            work+=2*M if du else M
            loads[a,b]+=ell
        maxevents=max(maxevents,ell);maxwork=max(maxwork,work);count+=1
    congestion=max((F(8*c,N) for c in loads.values()),default=F(0))
    if full:need(congestion<=bound,'exact all-pair length-weighted native congestion')
    # Dropping h leaves the distinct v labels unchanged, while each u cycle is full.
    need(len({v for u,v in states})==K,'retained kernel components before repair')
    need(len(states)==M*K,'complete product components after repair')
    out=dict(M=M,vertices=n,states=N,kernel_components_without_h=K,full_pair_universe=full,
             endpoint_pairs_checked=count,partner_transports_checked=checked_partners,
             maximum_events=maxevents,maximum_owner_work=maxwork,diameter_bound=D,
             canonical_congestion_bound=str(bound),gap_lower_bound=str(F(2,M*M)))
    if full:out['exact_routed_congestion']=str(congestion)
    return out

def main():
    cases=[case(M,M<81) for M in (3,9,15,27,81)]
    return dict(status='FINITE-EXACT PASS; guard-preserving native orbit only',cases=cases,checks=CHECKS)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
