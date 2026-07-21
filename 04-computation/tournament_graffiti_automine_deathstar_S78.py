#!/usr/bin/env python3
"""
tournament_graffiti_automine_deathstar_S78.py  (HYP-8636)

Capstone: the TournamentGraffiti engine at full capacity. Computes the combined
classical+spectral invariant set over all tournaments n<=6 and SYSTEMATICALLY mines:
 (1) comparability DAG: which ordered pairs (X,Y) have X <= Y for ALL tournaments;
 (2) tight additive bounds X <= Y + c (smallest integer c), flagging TIGHT ones
     (equality achieved) = WOWII-shaped candidates;
 (3) forbidden values: which integer invariants skip values (ndev=2 is the known one).
Survivors are candidate theorems; tight+equality ones are the WOWII targets.
"""
from math import comb
from itertools import combinations
import numpy as np

def all_t(n):
    P = list(combinations(range(n), 2))
    for b in range(1 << len(P)):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(P):
            if b>>k&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def inv(A,n):
    s=[sum(A[i]) for i in range(n)]
    sumC2=sum(comb(x,2) for x in s); c3=comb(n,3)-sumC2
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        for v in range(n):
            c=dp[m][v]
            if c:
                for w in range(n):
                    if not(m>>w&1) and out[v]>>w&1: dp[m|1<<w][w]+=c
    H=sum(dp[full][v] for v in range(n))
    beta=0
    for m in range(1<<n):
        S=[v for v in range(n) if m>>v&1]
        if S and sorted(sum(A[v][w] for w in S) for v in S)==list(range(len(S))): beta=max(beta,len(S))
    kings=0
    for v in range(n):
        r=out[v]|(1<<v)
        for w in range(n):
            if out[v]>>w&1: r|=out[w]
        if r==full: kings+=1
    M=np.array(A,float); ev=np.linalg.eigvals(M)
    ndev=len({(round(z.real,4),round(z.imag,4)) for z in ev})
    disc=round(abs(np.linalg.det(np.eye(n)+M-M.T)))//2**(n-1)
    return dict(n=n,c3=c3,H=H,beta=beta,kings=kings,smax=max(s),smin=min(s),
                srange=max(s)-min(s),sumC2=sumC2,ndev=ndev,disc=disc)

print("computing invariants n=3..6 ...",flush=True)
DATA=[inv(A,n) for n in (3,4,5,6) for A in all_t(n)]
print(f"{len(DATA)} tournaments\n",flush=True)
KEYS=["c3","H","beta","kings","smax","smin","srange","sumC2","ndev","disc","n"]

print("="*70,"\n(1) COMPARABILITY DAG: X <= Y for ALL tournaments (nontrivial only)\n","="*70)
edges=[]
for X in KEYS:
    for Y in KEYS:
        if X==Y: continue
        if all(d[X]<=d[Y] for d in DATA):
            # skip if trivially dominated via a third key chain? just report direct
            tight=any(d[X]==d[Y] for d in DATA)
            edges.append((X,Y,tight))
for X,Y,t in edges:
    print(f"   {X:7s} <= {Y:7s}   {'(TIGHT: equality achieved)' if t else ''}")

print("\n"+"="*70,"\n(2) TIGHT ADDITIVE BOUNDS  X <= Y + c  (smallest c; flag WOWII-shaped)\n","="*70)
import itertools
found=[]
for X,Y in itertools.permutations(KEYS,2):
    if X=="n" or Y=="n": continue
    c=max(d[X]-d[Y] for d in DATA)          # smallest c with X<=Y+c always
    if 1<=c<=3:
        tight=any(d[X]==d[Y]+c for d in DATA)
        if tight: found.append((X,Y,c))
for X,Y,c in sorted(found,key=lambda z:z[2])[:25]:
    print(f"   {X:7s} <= {Y:7s} + {c}   (tight)")

print("\n"+"="*70,"\n(3) FORBIDDEN VALUES: integer invariants that SKIP values (n<=6 union)\n","="*70)
for K in ["ndev","disc","beta","kings","srange","c3"]:
    vals=sorted({d[K] for d in DATA}); full=list(range(min(vals),max(vals)+1))
    miss=[v for v in full if v not in vals]
    print(f"   {K:7s}: values {vals}"+(f"   FORBIDDEN: {miss}" if miss else "   (no gap)"))
print("\n   -> ndev forbids {2} (THM-1858, PROVED). Others: reported above.")
