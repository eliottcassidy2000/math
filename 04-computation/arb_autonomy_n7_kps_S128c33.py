#!/usr/bin/env python3
"""arb_autonomy_n7_kps_S128c33.py -- is the tau-drift autonomous (function of tau_tot alone) at n=7?
Exhaustive 2^15 tilings; saves the g_7 table."""
import sys, json
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
n=7
def det_int(M):
    M=[r[:] for r in M]; k_n=len(M); sign=1; prev=1
    for k in range(k_n-1):
        if M[k][k]==0:
            for i in range(k+1,k_n):
                if M[i][k]: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,k_n):
            for j in range(k+1,k_n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]
def tau_tot(B):
    L=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]: L[v][v]+=1; L[v][u]-=1
    s=0
    for r in range(n):
        idx=[i for i in range(n) if i!=r]
        s+=det_int([[L[i][j] for j in idx] for i in idx])
    return s
tiles=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
m=len(tiles); m2=comb(n,2)
g={}; splits=0
for t in range(1<<m):
    B=[[False]*n for _ in range(n)]
    for k in range(2,n+1): B[k-1][k-2]=True
    for i,(x,y) in enumerate(tiles):
        if (t>>i)&1: B[x-1][y-1]=True
        else: B[y-1][x-1]=True
    tt=tau_tot(B)
    if tt in g and g[tt]=='SPLIT': continue
    s=0
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]:
                B[u][v],B[v][u]=False,True
                s+=tau_tot(B)-tt
                B[u][v],B[v][u]=True,False
    d=F(s,m2)
    if tt in g:
        if g[tt]!=d: g[tt]='SPLIT'; splits+=1
    else: g[tt]=d
    if t%4096==0: print("  ...%d/32768 tau-values %d splits %d"%(t,len(g),splits),flush=True)
print("n=7: distinct tau_tot values %d ; AUTONOMY %s (splits %d)"%(len(g),"HOLDS" if splits==0 else "FAILS",splits))
json.dump({str(k):(str(v) if v!='SPLIT' else 'SPLIT') for k,v in g.items()},
          open(r"C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-math\f631d0eb-9f23-494b-bb86-e0501bc456e9\scratchpad\n7_tau_g7.json","w"))
print("DONE")
