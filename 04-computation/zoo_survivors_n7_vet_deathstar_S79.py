#!/usr/bin/env python3
"""
zoo_survivors_n7_vet_deathstar_S79.py

Vet the full-zoo survivor inequalities at n=7 (holds-small-breaks-larger discipline).
Random sample of n=7 tournaments; report which survivors HOLD vs BREAK, with the
first counterexample. The spectral-dominance cluster gamma,fas,dichr <= ndev is the
prize; H<=arb, disc<=H, fas<=c3 the others.
"""
import random
from math import comb
from itertools import combinations, permutations
import numpy as np
random.seed(20260721)

def rand_t(n):
    A=[[0]*n for _ in range(n)]
    for i,j in combinations(range(n),2):
        if random.random()<0.5: A[i][j]=1
        else: A[j][i]=1
    return A

def ndev(A,n):
    ev=np.linalg.eigvals(np.array(A,float))
    return len({(round(z.real,4),round(z.imag,4)) for z in ev})
def disc(A,n):
    M=np.array(A,float); return round(abs(np.linalg.det(np.eye(n)+M-M.T)))//2**(n-1)
def c3(A,n):
    s=[sum(A[i]) for i in range(n)]; return comb(n,3)-sum(comb(x,2) for x in s)
def H(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        for v in range(n):
            c=dp[m][v]
            if c:
                for w in range(n):
                    if not(m>>w&1) and out[v]>>w&1: dp[m|1<<w][w]+=c
    return sum(dp[full][v] for v in range(n))
def gamma(A,n):
    for k in range(1,n+1):
        for S in combinations(range(n),k):
            dom=set(S)
            for s in S:
                for w in range(n):
                    if A[s][w]: dom.add(w)
            if len(dom)==n: return k
    return n
def dichr(A,n):
    def acyc(S): return sorted(sum(A[v][w] for w in S) for v in S)==list(range(len(S)))
    for k in range(1,n+1):
        def bt(v,parts):
            if v==n: return True
            for p in parts:
                p.append(v)
                if acyc(p) and bt(v+1,parts): return True
                p.pop()
            return False
        if bt(0,[[] for _ in range(k)]): return k
    return n
def fas(A,n):
    best=comb(n,2)
    for perm in permutations(range(n)):
        pos={v:i for i,v in enumerate(perm)}
        b=sum(1 for i in range(n) for j in range(n) if A[i][j] and pos[i]>pos[j])
        if b<best: best=b
        if best==0: break
    return best
def arb(A,n):
    indeg=[sum(A[i][j] for i in range(n)) for j in range(n)]
    L=[[(indeg[j] if i==j else 0)-A[i][j] for j in range(n)] for i in range(n)]
    tot=0
    for r in range(n):
        M=[[L[i][j] for j in range(n) if j!=r] for i in range(n) if i!=r]
        tot+=round(abs(np.linalg.det(np.array(M,float))))
    return tot
def kings(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    k=0
    for v in range(n):
        r=out[v]|(1<<v)
        for w in range(n):
            if out[v]>>w&1: r|=out[w]
        if r==full: k+=1
    return k

SURV=[
 ("gamma <= ndev", lambda d: d['gamma']<=d['ndev']),
 ("fas <= ndev",   lambda d: d['fas']<=d['ndev']),
 ("dichr <= ndev", lambda d: d['dichr']<=d['ndev']),
 ("H <= arb",      lambda d: d['H']<=d['arb']),
 ("disc <= H",     lambda d: d['disc']<=d['H']),
 ("fas <= c3",     lambda d: d['fas']<=d['c3']),
]
n=7; R=12000
viol={name:None for name,_ in SURV}
fvals=set(); ndevvals=set(); kingsvals=set(); hamc=set()
for _ in range(R):
    A=rand_t(n)
    d=dict(ndev=ndev(A,n),disc=disc(A,n),c3=c3(A,n),H=H(A,n),gamma=gamma(A,n),
           dichr=dichr(A,n),fas=fas(A,n),arb=arb(A,n))
    for name,f in SURV:
        if viol[name] is None and not f(d): viol[name]=dict(d)
    ndevvals.add(d['ndev']); kingsvals.add(kings(A,n)); fvals.add(d['fas'])
print(f"n=7 survivor vetting ({R} random samples):")
for name,_ in SURV:
    v=viol[name]
    print(f"  {'HOLDS ' if v is None else 'BREAKS'} {name}"+("" if v is None else f"   counterexample {v}"))
print(f"\n  forbidden-value check n=7: ndev values {sorted(ndevvals)} (2 absent? {2 not in ndevvals})")
print(f"                            kings values {sorted(kingsvals)} (2 absent? {2 not in kingsvals})")
print(f"                            fas values {sorted(fvals)}")
