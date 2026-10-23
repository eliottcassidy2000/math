#!/usr/bin/env python3
"""Standalone native order31 counterexample: arc membership is not a flag."""
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
A=double(double([[int((j-i)%7 in (1,2,4)) for j in range(7)] for i in range(7)]))
P=(23,0,8,24,16,1,9,25,17,2,10,26,18,3,11,27,19,4,12,28,20,5,13,29,21,6,14,30,22,15,7)
need(sorted(P)==list(range(31)) and all(A[u][v] for u,v in zip(P,P[1:])),'actual Hamiltonian path')
G=[tuple(v if v in (7,15,23) else 8*(v//8)+(a*(v%8)+b)%7 for v in range(31)) for a in (1,2,4) for b in range(7)]
need(all(all(A[u][v]==A[g[u]][g[v]] for u in range(31) for v in range(31)) for g in G),'actual source automorphisms')
paths=[tuple(g[v] for v in P) for g in G]
need(len(set(paths))==21,'distinct source-automorphic paths')
edges=[set(zip(Q,Q[1:])) for Q in paths];base=list(zip(P,P[1:]))
hits=[sum(e in E for E in edges) for e in base]
need(min(hits)>1,'no single arc occurrence identifies P')
pairs=[(e,f) for e,f in itertools.combinations(base,2) if sum(e in E and f in E for E in edges)==1]
need(len(pairs)==13 and pairs[0]==((23,0),(16,1)),'exact optimal two arc certificate')
i=4; e=base[i]
need(e==(16,1) and sum((Q[i],Q[i+1])==e for Q in paths)==1,'one indexed flag identifies P')
print(json.dumps({'status':'PASS','path':P,'source_automorphism_count':21,'arc_membership_counts':hits,'minimum_unindexed_arcs':2,'isolating_pair_count':len(pairs),'first_isolating_pair':pairs[0],'indexed_flag':{'position_zero_based':i,'arc':e,'survivors':1}},indent=2))
