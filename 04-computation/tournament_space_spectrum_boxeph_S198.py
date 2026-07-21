#!/usr/bin/env python3
"""tournament_space_spectrum_boxeph_S198.py -- boxeph-2026-07-21-S198

Tournament space on n vertices as a SPECTRUM from a single point to a continuum.

Tournament space FIBERS over the score sequence (Landau's polytope of valid score sequences).
The spectral coordinate is the score SPREAD sigma^2 = Var(scores), ranging from
  sigma^2 = (n^2-1)/12   (the TRANSITIVE vertex, scores 0..n-1)   -- the single point
  sigma^2 = 0            (the REGULAR center, all scores (n-1)/2)  -- the continuum.
The FIBER (# iso classes with a given score sequence) is a SINGLETON at the transitive vertex and
swells to its maximum at the regular center, where the strongly-connected / modular-prime /
structurally-diverse tournaments live. This computes the fibration and shows structure runs
OPPOSITE to score spread.
"""
from itertools import permutations, combinations
from fractions import Fraction as Fr
from collections import defaultdict
from math import comb

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
def scores(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def strong(A,n):
    R=[[1 if(i==j or A[i][j])else 0 for j in range(n)]for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j]=1
    return all(R[i][j] and R[j][i] for i in range(n) for j in range(n))
def modular_prime(A,n):
    for k in range(2,n):
        for M in combinations(range(n),k):
            Ms=set(M)
            if all(len(set(A[v][u] for u in M))==1 for v in range(n) if v not in Ms):
                return False
    return True
def num_c3(A,n):
    sc=[sum(A[i]) for i in range(n)]; return comb(n,3)-sum(comb(s,2) for s in sc)
def variance(sq,n):
    m=Fr(n-1,2); return sum((Fr(s)-m)**2 for s in sq)/n

print("enumerating iso classes n<=7 ...", flush=True)
reps=iso_reps(7)
print("iso classes:", {n:len(reps[n]) for n in range(3,8)}, flush=True)

for n in range(3,8):
    fibers=defaultdict(list)
    for A in reps[n]:
        fibers[scores(A,n)].append(A)
    rows=[]
    for sq,cls in fibers.items():
        v=variance(sq,n)
        st=sum(1 for A in cls if strong(A,n)); mp=sum(1 for A in cls if modular_prime(A,n))
        c3s=[num_c3(A,n) for A in cls]
        rows.append((v, len(cls), st, mp, sq, min(c3s), max(c3s)))
    rows.sort(key=lambda r: float(r[0]))
    vmax=Fr(n*n-1,12)
    print("\n"+"="*92)
    print("n=%d : %d score sequences (Landau); %d iso classes; sigma^2 in [0, (n^2-1)/12=%s]" %
          (n, len(fibers), len(reps[n]), vmax))
    print("="*92)
    print("  sigma^2   fiber  strong  modprime  c3-range   score sequence")
    for v,sz,st,mp,sq,c3lo,c3hi in rows:
        tag=""
        if v==vmax: tag=" <- TRANSITIVE vertex (single point)"
        if v==0:    tag=" <- REGULAR center (continuum)"
        print("  %-7s  %4d   %4d    %5d     [%d..%d]   %s%s" %
              (str(v), sz, st, mp, c3lo, c3hi, sq, tag))
    # spectrum summary: correlation of fiber size with spread
    maxf=max(r[1] for r in rows); trans=[r for r in rows if r[0]==vmax][0]
    reg=[r for r in rows if r[0]==0]
    print("  --> max fiber = %d (near regular); transitive fiber = %d; %s" %
          (maxf, trans[1], ("regular fiber = %d (all strong=%d, all modprime=%d)"%(reg[0][1],reg[0][2],reg[0][3])) if reg else "no exactly-regular seq (n even)"))
