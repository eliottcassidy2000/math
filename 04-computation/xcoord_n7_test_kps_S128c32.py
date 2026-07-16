#!/usr/bin/env python3
"""xcoord_n7_test_kps_S128c32.py -- does X = H - 6 c5 - w7 c7 fiber the n=7 drift? Save tuples."""
import sys, json
from math import comb
from fractions import Fraction as F
from itertools import permutations, combinations
sys.stdout.reconfigure(line_buffering=True)
n = 7
tiles = [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
tidx = {t:i for i,t in enumerate(tiles)}; m = len(tiles); m2 = comb(n,2)
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
def codd(B,L):
    tot=0
    for S in combinations(range(n),L):
        u=S[0]
        for perm in permutations(S[1:]):
            prev=u; ok=True
            for w in perm:
                if not B[prev][w]: ok=False; break
                prev=w
            if ok and B[prev][u]: tot+=1
    return tot
tup={}
for t in range(1<<m):
    tv=[(t>>i)&1 for i in range(m)]
    B=[[False]*n for _ in range(n)]
    for k in range(2,n+1): B[k-1][k-2]=True
    for (x,y),i in tidx.items():
        if tv[i]: B[x-1][y-1]=True
        else: B[y-1][x-1]=True
    key0=(ham(B), codd(B,5), codd(B,7))
    if key0 in tup: continue
    H0=key0[0]
    s=0
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]:
                B[u][v],B[v][u]=False,True; s+=ham(B)-H0; B[u][v],B[v][u]=True,False
    tup[key0]=F(s,m2)
    if t % 4096 == 0: print("  ...%d, tuples %d" % (t, len(tup)), flush=True)
print("tuples:", len(tup))
json.dump({"%d,%d,%d"%k: [v.numerator, v.denominator] for k,v in tup.items()}, open(r"C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-math\f631d0eb-9f23-494b-bb86-e0501bc456e9\scratchpad\n7_drift_tuples.json","w"))
# X-fibering test for w7 in candidates
for w7num, w7den in [(0,1),(6,1),(12,1),(18,1),(24,1),(30,1),(36,1),(2,1),(3,1),(4,1),(8,1),(10,1),(20,1)]:
    fib={}
    ok=True
    for (H,c5,c7),v in tup.items():
        X = F(H) - 6*c5 - F(w7num,w7den)*c7
        if X in fib and fib[X]!=v: ok=False; break
        fib[X]=v
    if ok:
        print("X = H - 6 c5 - %s/%s c7 FIBERS the n=7 drift (distinct X: %d)" % (w7num,w7den,len(fib)))
        break
else:
    print("no single-w7 linear X fibers n=7; (H,c5,c7) remains the minimal coordinate triple")
print("DONE")
