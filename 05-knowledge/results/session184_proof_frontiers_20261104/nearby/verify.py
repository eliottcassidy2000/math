#!/usr/bin/env python3
"""Dependency-free exact controls; this is not the proof of infinite statements."""
from itertools import combinations, product
from math import isqrt
import json

def require(ok, message):
    if not ok: raise RuntimeError(message)

def factors(n):
    out=[]; d=2
    while d*d<=n:
        while n%d==0: out.append(d); n//=d
        d+=1
    if n>1: out.append(n)
    return out

def direct(n, hole):
    weak=[]; strict=[]
    for d in range(2,isqrt(n)+1):
        if n%d: continue
        pair=(d,n//d); weak.append(pair)
        if pair[0]<pair[1]: strict.append(pair)
    if hole(n): return False,False
    return (bool(weak) and all(hole(a) or hole(b) for a,b in weak),
            bool(strict) and all(hole(a) or hole(b) for a,b in strict))

def prime_classifier(n,H):
    f=factors(n); k=len(f)
    w=(k==2 and any(p in H for p in f)) or (k==3 and all(p in H for p in f))
    t=(w and not(k==2 and f[0]==f[1])) or(k==4 and len(set(f))==1 and f[0] in H)
    return w,t

N=2500
facts={n:factors(n) for n in range(2,N+1)}
primes=[2,3,5,7,11,13]
prime_cases=0
for mask in range(1<<len(primes)):
    H={p for j,p in enumerate(primes) if mask>>j&1}
    for n in range(2,N+1):
        require(direct(n,H.__contains__)==prime_classifier(n,H),("prime",H,n))
        prime_cases+=1

rank_ideal_cases=0
for K in range(1,6):
    hole=lambda n: len(factors(n))<=K
    for n,f in facts.items():
        k=len(f)
        w=K+1<=k<=2*K+1
        t=(w and not(k==2 and len(set(f))==1)) or (k==2*K+2 and len(set(f))==1)
        require(direct(n,hole)==(w,t),("rank ideal",K,n))
        rank_ideal_cases+=1
    for p in (2,3,5):
        H={p**j for j in range(1,K+1)}
        require(direct(p**(2*K+1),H.__contains__)==(True,True),("weak sharp",K,p))
        require(direct(p**(2*K+2),H.__contains__)==(False,True),("strict sharp",K,p))

bounded_cases=0
# All subsets of a mixed-prime universe, including composites.
universe=[2,3,4,5,6,8,9,10,12]
for mask in range(1<<len(universe)):
    H={p for j,p in enumerate(universe) if mask>>j&1}
    K=max([len(factors(h)) for h in H]+[1])
    for n in range(2,401):
        w,t=direct(n,H.__contains__); f=facts[n]; k=len(f)
        require(not w or k<=2*K+1,("weak bound",H,n))
        require(not t or k<=2*K+2,("strict bound",H,n))
        if k==2*K+1:
            forced=all(d in H for d in range(2,n) if n%d==0 and len(factors(d))<=K)
            require(w==forced,("weak equality",H,n))
        if k==2*K+2:
            predicted=len(set(f))==1 and all(f[0]**j in H for j in range(1,K+1))
            require(t==predicted,("strict equality",H,n))
        bounded_cases+=1

# Dirichlet series as formal integer coefficient arrays, compared with direct fibres.
# Multiply by 6 to avoid all rational arithmetic.
M=1800
P=[0]*(M+1)
for n in range(2,M+1): P[n]=int(len(facts[n])==1)
def conv(a,b):
    c=[0]*(M+1)
    for i in range(2,M+1):
        if a[i]:
            for j in range(2,M//i+1):
                if b[j]: c[i*j]+=a[i]*b[j]
    return c
def dilate(a,r):
    b=[0]*(M+1)
    for n in range(2,M+1):
        if n**r>M: break
        b[n**r]=a[n]
    return b
profile_cases=0
for mode in ("all","finite","mod4","sparse"):
    H={n for n in range(2,M+1) if P[n] and (mode=="all" or mode=="finite" and n in {2,5,11} or mode=="mod4" and n%4==1 or mode=="sparse" and n in {2,17,257})}
    Q=[int(n in H) for n in range(M+1)]
    Q2=dilate(Q,2); Q3=dilate(Q,3); Q4=dilate(Q,4)
    QP=conv(Q,P); QQ=conv(Q,Q); QQQ=conv(QQ,Q); QQ2=conv(Q,Q2)
    for n in range(2,M+1):
        sw=6*QP[n]-3*QQ[n]+3*Q2[n]+QQQ[n]+3*QQ2[n]+2*Q3[n]
        st=sw-6*Q2[n]+6*Q4[n]
        w,t=direct(n,H.__contains__)
        require((sw,st)==(6*int(w),6*int(t)),("series",mode,n,sw,st,w,t))
        profile_cases+=1

# Explicit collisions/diagonals and infinite-hole-vs-finite-cage controls.
require(direct(4,{2}.__contains__)==(True,False),"square diagonal")
require(direct(8,{2}.__contains__)==(True,True),"cube repetition")
require(direct(16,{2}.__contains__)==(False,True),"fourth power strict")
require(direct(30,{2,3}.__contains__)==(False,False),"one retained prime saves triple")
require(direct(30,{2,3,5}.__contains__)==(True,True),"all deleted primes kill triple")
require(direct(36,{2,3}.__contains__)==(False,False),"unequal rank-two factors 4*9 survive")
print(json.dumps({"status":"PASS","prime_hole_target_cases":prime_cases,"rank_ideal_target_cases":rank_ideal_cases,"bounded_rank_target_cases":bounded_cases,"exact_series_coefficient_cases":profile_cases,"sharp_endpoint_checks":30,"hostiles":6},sort_keys=True,indent=2))
