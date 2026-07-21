#!/usr/bin/env python3
"""modular_prime_census_n7_boxeph_S197.py -- boxeph-2026-07-21-S197

n>=7 pattern hunt: the REDUCTION HIERARCHY stratifies tournaments three ways of increasing fineness:
  (1) order-join/SCC atoms = STRONG tournaments        (boxeph THM-1862/1926): 1,1,6,35,353 (n=3..7)
  (2) modular/substitution SEEDS = MODULAR-PRIME       (opus THM-1960): 1,1,1,0,3,15 (n=1..6), n=7 OPEN
  (3) circulant character-generated                    (boxeph THM-1955): 1,0,1,0,2 (n=3..7)
This completes opus's open 'seed census to n=7' and reports the strong-fraction and prime-fraction
trends into the hard-to-enumerate regime.

A HOMOGENEOUS SET / module M (2<=|M|<=n-1): every v outside M relates UNIFORMLY to all of M
(v beats all of M, or v loses to all of M). Modular-prime = has NO nontrivial module.
"""
from itertools import permutations, combinations

def canon(A,n,perms):
    best=None
    for p in perms:
        c=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or c<best: best=c
    return best
def iso_reps(nmax):
    reps={1:[[[0]]]}
    for n in range(2,nmax+1):
        perms=list(permutations(range(n))); seen=set(); out=[]
        for B in reps[n-1]:
            for pat in range(1<<(n-1)):
                A=[row[:]+[0] for row in B]+[[0]*n]
                for k in range(n-1):
                    if pat>>k&1: A[n-1][k]=1
                    else: A[k][n-1]=1
                c=canon(A,n,perms)
                if c not in seen:
                    seen.add(c); M=[[0]*n for _ in range(n)]; idx=0
                    for i in range(n):
                        for j in range(n):
                            if i!=j: M[i][j]=c[idx]; idx+=1
                    out.append(M)
        reps[n]=out
    return reps
def strong(A,n):
    R=[[1 if(i==j or A[i][j])else 0 for j in range(n)]for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j]=1
    return all(R[i][j] and R[j][i] for i in range(n) for j in range(n))
def is_module(A,n,M):
    Ms=set(M)
    for v in range(n):
        if v in Ms: continue
        vals=set(A[v][u] for u in M)
        if len(vals)>1: return False
    return True
def modular_prime(A,n):
    for k in range(2,n):
        for M in combinations(range(n),k):
            if is_module(A,n,list(M)): return False
    return True

print("computing iso reps up to n=7 ...", flush=True)
reps=iso_reps(7)
print("iso classes:", {n:len(reps[n]) for n in range(1,8)}, flush=True)
print()
print(" n : total  strong  modular-prime(seed)  reducible  strong-frac  prime-frac")
for n in range(1,8):
    tot=len(reps[n])
    st=sum(1 for A in reps[n] if strong(A,n)) if n>=1 else 0
    pr=sum(1 for A in reps[n] if modular_prime(A,n))
    red=tot-st
    print(" %d : %5d  %5d   %6d               %5d      %s        %s"
          % (n, tot, st, pr, red,
             ("%.3f"%(st/tot)) if tot else "-", ("%.3f"%(pr/tot)) if tot else "-"))
print()
print("=> opus seed census 1,1,1,0,3,15 (n=1..6) EXTENDED to n=7 (the entry above).")
print("   strong atoms (SCC) >= modular-prime seeds: substitution carves inside strong components.")
print("   Both fractions' trend into n=7 quantifies why the hard regime forces enumeration.")
