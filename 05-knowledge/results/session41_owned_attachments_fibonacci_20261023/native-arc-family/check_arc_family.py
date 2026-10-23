#!/usr/bin/env python3
import itertools,json

def need(p,m):
    if not p:raise ValueError(m)
def double(A):
    n=len(A);M=[[0]*(2*n+1) for _ in range(2*n+1)]
    for i in range(n):
        M[n][i]=1;M[n+1+i][n]=1
        for j in range(n):
            M[i][j]=A[i][j];M[n+1+i][n+1+j]=A[j][i]
            M[i][n+1+j]=int(i==j or A[i][j]);M[n+1+j][i]=1-M[i][n+1+j]
    return M
reports=[]
for q in (3,7,11,19,23,31):
    squares=sorted({x*x%q for x in range(1,q)});L=q+1;n=4*q+3
    A=double(double([[int((j-i)%q in squares) for j in range(q)] for i in range(q)]))
    centers=(q,2*q+1,3*q+2)
    P=[3*q+2]+[i*L+x for x in range(q) for i in (0,1,3,2)]+[2*q+1,q]
    need(sorted(P)==list(range(n)) and all(A[u][v] for u,v in zip(P,P[1:])),'actual family path')
    G=[tuple(v if v in centers else L*(v//L)+(a*(v%L)+b)%q for v in range(n)) for a in squares for b in range(q)]
    need(all(all(A[u][v]==A[g[u]][g[v]] for u in range(n) for v in range(n)) for g in G),'actual family automorphisms')
    paths=[tuple(g[v] for v in P) for g in G]
    need(len(set(paths))==q*(q-1)//2,'distinct affine family')
    edges=[set(zip(Q,Q[1:])) for Q in paths]
    base=list(zip(P,P[1:]));hits=[sum(e in E for E in edges) for e in base]
    expected=[(q-1)//2]+[q*(q-1)//2 if k%4!=3 else q-1 for k in range(4*q-1)]+[(q-1)//2,q*(q-1)//2]
    need(hits==expected,'exact arc incidence multiplicities')
    need(sum(base[0] in E and base[-2] in E for E in edges)==1,'first last pins isolate')
    need(sum((Q[4],Q[5])==(P[4],P[5]) for Q in paths)==1,'one indexed flag isolates')
    if q>=7:need(min(hits)>1,'no unindexed arc isolates')
    reports.append({'q':q,'order':n,'orbit_size':len(G),'minimum_arc_survivors':min(hits),'two_arc_survivors':1,'indexed_flag_survivors':1})
print(json.dumps({'status':'PASS','family':reports},indent=2))
