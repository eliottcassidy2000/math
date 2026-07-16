#!/usr/bin/env python3
"""arb_n7_tuples_kps_S128c33.py -- per-tiling (H,c5,c7,tau_tot,tau-drift) at n=7, saved;
then sub-coordinate minimality for the tau-drift at n=7."""
import sys, json
from math import comb
from fractions import Fraction as F
from itertools import permutations, combinations
sys.stdout.reconfigure(line_buffering=True)
n=7
def det_int(M):
    M=[r[:] for r in M]; kn=len(M); sign=1; prev=1
    for k in range(kn-1):
        if M[k][k]==0:
            for i in range(k+1,kn):
                if M[i][k]: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,kn):
            for j in range(k+1,kn):
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
def ham(B):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1<<n):
        for v in range(n):
            c=dp[S][v]
            if not c or not (S>>v)&1: continue
            for u in range(n):
                if not (S>>u)&1 and B[v][u]: dp[S|1<<u][u]+=c
    return sum(dp[(1<<n)-1][v] for v in range(n))
def codd(B,Lc):
    tot=0
    for S in combinations(range(n),Lc):
        u=S[0]
        for perm in permutations(S[1:]):
            prev=u; ok=True
            for w in perm:
                if not B[prev][w]: ok=False; break
                prev=w
            if ok and B[prev][u]: tot+=1
    return tot
tiles=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
m=len(tiles); m2=comb(n,2)
fib={}
for t in range(1<<m):
    B=[[False]*n for _ in range(n)]
    for k in range(2,n+1): B[k-1][k-2]=True
    for i,(x,y) in enumerate(tiles):
        if (t>>i)&1: B[x-1][y-1]=True
        else: B[y-1][x-1]=True
    tt=tau_tot(B)
    H=ham(B); c5=codd(B,5); c7=codd(B,7)
    key=(H,c5,c7,tt)
    if key in fib: continue   # (tau,census) determines drift (engine Q4-pre: 0/427 splits) -> dedup safe
    s=0
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]:
                B[u][v],B[v][u]=False,True
                s+=tau_tot(B)-tt
                B[u][v],B[v][u]=True,False
    fib[key]=F(s,m2)
    if t%4096==0: print("  ...%d/32768 states %d"%(t,len(fib)),flush=True)
print("n=7 (H,c5,c7,tau) states:",len(fib))
json.dump({"%d,%d,%d,%d"%k:[v.numerator,v.denominator] for k,v in fib.items()},
          open(r"C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-math\f631d0eb-9f23-494b-bb86-e0501bc456e9\scratchpad\n7_tau_census_tuples.json","w"))
def test(name,pr):
    f={}; bad=0
    for k,v in fib.items():
        p=pr(k)
        if p in f and f[p]!=v: bad+=1
        f.setdefault(p,v)
    print("  %-14s fibers=%4d %s"%(name,len(f),"DETERMINES" if bad==0 else "FAILS(%d)"%bad))
print("sub-coordinate minimality at n=7:")
test("tau",         lambda k:(k[3],))
test("(tau,H)",     lambda k:(k[3],k[0]))
test("(tau,c5)",    lambda k:(k[3],k[1]))
test("(tau,c7)",    lambda k:(k[3],k[2]))
test("(tau,H,c5)",  lambda k:(k[3],k[0],k[1]))
test("(tau,H,c7)",  lambda k:(k[3],k[0],k[2]))
test("(tau,c5,c7)", lambda k:(k[3],k[1],k[2]))
test("(H,c5,c7)",   lambda k:(k[0],k[1],k[2]))
test("full",        lambda k:k)
print("DONE")
