#!/usr/bin/env python3
"""arb_flip_law_referee_kps_S128c33.py -- kind-pasteur S128 cont.33.
Referee the arborescence flip laws:
  (T) transfer: Dtau_u = -F(u,v), Dtau_v = +F(u,v), F = det(L del {u,v})
  (B) bystander: Dtau_r = (e_u-e_v)^T adj(Ltil_r) (e_u+e_v), r not in {u,v}
on ALL tournaments n=4,5 and 300 random n=6..9. Exact integers throughout."""
import sys, random
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)
random.seed(12345)

def det_int(M):
    M=[r[:] for r in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            for i in range(k+1,n):
                if M[i][k]: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]

def lap(B,n):
    L=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]: L[v][v]+=1; L[v][u]-=1
    return L

def tau(B,n,r):
    L=lap(B,n); idx=[i for i in range(n) if i!=r]
    return det_int([[L[i][j] for j in idx] for i in idx])

def minor(L,rows,cols):
    ri=[i for i in range(len(L)) if i not in rows]
    ci=[j for j in range(len(L)) if j not in cols]
    return det_int([[L[i][j] for j in ci] for i in ri])

def adj_entry(L,idx,i,j):
    # adj(M)[i][j] = (-1)^{i+j} * minor_{j,i}(M) for M = L restricted to idx
    M=[[L[a][b] for b in idx] for a in idx]
    ii=idx.index(i); jj=idx.index(j)
    ri=[a for a in range(len(M)) if a!=jj]; ci=[b for b in range(len(M)) if b!=ii]
    if not ri: return 1  # 1x1 adjugate
    s=(-1)**(ii+jj)
    return s*det_int([[M[a][b] for b in ci] for a in ri])

def check(B,n):
    L=lap(B,n)
    t0=[tau(B,n,r) for r in range(n)]
    for u in range(n):
        for v in range(n):
            if u==v or not B[u][v]: continue
            B[u][v],B[v][u]=False,True
            t1=[tau(B,n,r) for r in range(n)]
            B[u][v],B[v][u]=True,False
            F=minor(L,{u,v},{u,v})
            assert t1[u]-t0[u]==-F, ("T-u",n,u,v,t1[u]-t0[u],-F)
            assert t1[v]-t0[v]==+F, ("T-v",n,u,v,t1[v]-t0[v],F)
            for r in range(n):
                if r in (u,v): continue
                idx=[i for i in range(n) if i!=r]
                # (e_u-e_v)^T adj (e_u+e_v) = adj[u][u]+adj[u][v]-adj[v][u]-adj[v][v]
                d=(adj_entry(L,idx,u,u)+adj_entry(L,idx,u,v)
                   -adj_entry(L,idx,v,u)-adj_entry(L,idx,v,v))
                assert t1[r]-t0[r]==d, ("B",n,r,u,v,t1[r]-t0[r],d)
    return True

cnt=0
for n in (4,5):
    pairs=list(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        B=[[False]*n for _ in range(n)]
        for k,(a,b) in enumerate(pairs):
            if (mask>>k)&1: B[a][b]=True
            else: B[b][a]=True
        check(B,n); cnt+=1
    print("n=%d: ALL %d tournaments pass (T)+(B)"%(n,1<<len(pairs)))
for n in (6,7,8,9):
    for _ in range(75):
        B=[[False]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1,n):
                if random.random()<0.5: B[a][b]=True
                else: B[b][a]=True
        check(B,n); cnt+=1
    print("n=%d: 75 random tournaments pass"%n)
print("TOTAL %d tournaments, zero failures. Flip laws EXACT."%cnt)
