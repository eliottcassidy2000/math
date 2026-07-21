#!/usr/bin/env python3
"""tournament_space_spectrum_boxeph_S198.py -- boxeph-2026-07-21-S198

Tournament space on n vertices as a SPECTRUM from a single point to a continuum.

Tournament space FIBERS over the score sequence (Landau's polytope of valid score sequences).
The score SPREAD sigma^2 = Var(scores) ranges from the transitive maximum
(n^2-1)/12 down to epsilon_n=0 for odd n and 1/4 for even n.  It determines
c3 exactly, but it does NOT determine the score fiber's size or structure.
This audit reports the parity seam and same-variance structural collisions.

Tournament Analysis is deliberately not forced here: taking score fibers as
vertices and variance difference as the pairwise observable gives only a
transitive preorder (with genuine ties), while an arbitrary tie gauge would
destroy precisely the within-level information under audit.
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
    eps=Fr(0) if n%2 else Fr(1,4)
    vmax=Fr(n*n-1,12)
    c3max=n*(n*n-(1 if n%2 else 4))//24
    assert min(r[0] for r in rows)==eps
    assert max(r[6] for r in rows)==c3max
    assert all(lo==hi for _,_,_,_,_,lo,hi in rows)
    print("\n"+"="*92)
    print("n=%d : %d score sequences (Landau); %d iso classes; sigma^2 in [%s, (n^2-1)/12=%s]" %
          (n, len(fibers), len(reps[n]), eps, vmax))
    print("="*92)
    print("  sigma^2   fiber  strong  modprime  c3-range   score sequence")
    for v,sz,st,mp,sq,c3lo,c3hi in rows:
        tag=""
        if v==vmax: tag=" <- TRANSITIVE vertex (single point)"
        if v==eps:  tag=" <- BALANCED edge (%s)" % ("regular" if n%2 else "near-regular")
        print("  %-7s  %4d   %4d    %5d     [%d..%d]   %s%s" %
              (str(v), sz, st, mp, c3lo, c3hi, sq, tag))
    # Exact summary: the scalar coordinate has nontrivial structural fibers.
    maxf=max(r[1] for r in rows); trans=[r for r in rows if r[0]==vmax][0]
    balanced=[r for r in rows if r[0]==eps]
    maxlocs=sorted(set(r[0] for r in rows if r[1]==maxf))
    byvar=defaultdict(list)
    for row in rows: byvar[row[0]].append(row)
    varying=[]
    for v,group in sorted(byvar.items(), key=lambda item: float(item[0])):
        sigs=sorted(set((r[1],r[2],r[3]) for r in group))
        if len(sigs)>1: varying.append((v,sigs))
    print("  --> max fiber = %d at sigma^2=%s; transitive fiber = %d; balanced fiber = %d (strong=%d, modprime=%d)" %
          (maxf, ",".join(map(str,maxlocs)), trans[1], balanced[0][1], balanced[0][2], balanced[0][3]))
    print("  --> same-variance levels with differing (fiber,strong,modprime): %d%s" %
          (len(varying), ("; first=%s -> %s"%(varying[0][0],varying[0][1])) if varying else ""))

print("\nTOURNAMENT-ANALYSIS NOTE: variance comparison gives a transitive preorder with ties;")
print("forcing a tie gauge would erase the score-fiber information being measured, so no tournament fingerprint is reported.")
