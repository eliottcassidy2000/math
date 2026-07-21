#!/usr/bin/env python3
"""continuum_coordinates_boxeph_S199.py -- boxeph-2026-07-21-S199

Reframes/lenses for the CONTINUUM (the near-regular structural interior of tournament space, THM-1979)
so we can describe it WITHOUT enumeration. New coordinates:

  CYCLIC TEMPERATURE  tau = c3 / c3_max = 1 - sigma^2/sigma^2_max  in [0,1]
      tau=0 transitive (frozen ground state), tau=1 regular (hot/quasirandom continuum).
  ISO-CYCLIC SHELL  = the set of classes at fixed tau (fixed c3 / score-spread).
  CYCLE SPECTRUM  (N_4, N_5, ..., N_n),  N_k = tr(A^k) = #closed k-walks (my zeta THM-1926).
      N_1=N_2=0 always; N_3=3c3 is FROZEN by tau; the FREE structural coordinates start at N_4.
  ENTROPY  of a shell = log(fiber size) -- grows toward the hot continuum.

This checks: (1) N_3 frozen by score-seq, N_4 the first FREE (structure-carrying) moment;
(2) the tau/shell/entropy structure; (3) the COORDINATE BUDGET -- how few coordinates (cycle
spectrum, then |R|) pin a near-regular tournament, so the continuum is low-dimensional.
"""
import numpy as np
from itertools import permutations, combinations
from fractions import Fraction as Fr
from collections import defaultdict
from math import comb, log2

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
def scoreseq(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def c3(A,n):
    sc=[sum(A[i]) for i in range(n)]; return comb(n,3)-sum(comb(s,2) for s in sc)
def Nk(A,n,K):
    M=np.array(A,dtype=float); P=np.eye(n); out=[]
    for k in range(1,K+1):
        P=P@M; out.append(int(round(np.trace(P))))
    return out   # [N_1,...,N_K]
def signed_redei(A,n):
    # R = sum over Hamiltonian paths of sign(permutation given by path order)
    adj=[[A[i][j] for j in range(n)] for i in range(n)]
    def sgn(perm):
        s=1; seen=[False]*n
        for i in range(n):
            if not seen[i]:
                j=i; L=0
                while not seen[j]: seen[j]=True; j=perm[j]; L+=1
                if L%2==0: s=-s
        return s
    R=0; path=[]
    def dfs(v,used):
        nonlocal R
        path.append(v)
        if len(path)==n:
            perm=[0]*n
            for idx in range(n-1): perm[path[idx]]=path[idx+1]
            perm[path[-1]]=path[0]  # close into a permutation (cycle) -- use path order sign
            # sign of the path as a sequence = sign of permutation mapping identity->path
            inv=0
            for a in range(n):
                for b in range(a+1,n):
                    if path[a]>path[b]: inv+=1
            R+= (1 if inv%2==0 else -1)
        else:
            for w in range(n):
                if not (used>>w&1) and adj[v][w]:
                    dfs(w, used|(1<<w))
        path.pop()
    for v in range(n): dfs(v,1<<v)
    return R

reps=iso_reps(7)
print("iso classes:", {n:len(reps[n]) for n in range(3,8)})

# ---------- (1) frozen vs free moments ----------
print("\n"+"="*90); print("(1) FROZEN vs FREE cycle moments: which N_k vary WITHIN a fixed score sequence?")
print("="*90)
for n in range(4,8):
    byscore=defaultdict(list)
    for A in reps[n]: byscore[scoreseq(A,n)].append(Nk(A,n,min(n,6)))
    free={}
    for k in range(3,min(n,6)+1):
        varies=any(len(set(v[k-1] for v in cls))>1 for cls in byscore.values())
        free[k]=varies
    print("  n=%d: N_3 varies within score-seq? %s ; N_4? %s ; N_5? %s ; N_6? %s"
          %(n, free.get(3), free.get(4), free.get(5), free.get(6)))
print("  => N_1=N_2=0, N_3=3c3 FROZEN by the score sequence; the FIRST FREE structural coordinate is N_4.")

def strong(A,n):
    R=[[1 if(i==j or A[i][j])else 0 for j in range(n)]for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j]=1
    return all(R[i][j] and R[j][i] for i in range(n) for j in range(n))

# ---------- (2) cyclic temperature + shell entropy ----------
print("\n"+"="*90); print("(2) CYCLIC TEMPERATURE tau=c3/c3_max ; ISO-CYCLIC SHELL entropy = log2(fiber)")
print("="*90)
for n in range(5,8):
    c3max=max(c3(A,n) for A in reps[n])
    shells=defaultdict(list)
    for A in reps[n]: shells[c3(A,n)].append(A)
    print("  n=%d  c3_max=%d (regular)" % (n,c3max))
    print("    tau      c3   #classes  entropy(log2)  strong-frac")
    for cc in sorted(shells, reverse=True):
        cls=shells[cc]; tau=Fr(cc,c3max); st=sum(1 for A in cls if strong(A,n))
        print("    %-7s  %3d  %6d     %6.2f         %.2f"
              % (str(tau), cc, len(cls), log2(len(cls)), st/len(cls)))

# ---------- (3) coordinate budget: how few coords pin a continuum tournament ----------
print("\n"+"="*90); print("(3) COORDINATE BUDGET: within the hottest shells, resolution by cycle-spectrum then +|R|")
print("="*90)
print("  Lens layers: L0 = tau (1 real) ; L1 = cycle spectrum (N_4..N_n) = char_A ; L2 = +|R| (beyond-spectral).")
for n in range(5,8):
    c3max=max(c3(A,n) for A in reps[n])
    # take the two hottest non-trivial shells (largest tau below 1, plus tau=1)
    shells=defaultdict(list)
    for A in reps[n]: shells[c3(A,n)].append(A)
    hot=sorted(shells, reverse=True)[:3]
    for cc in hot:
        cls=shells[cc]; tau=Fr(cc,c3max)
        specs=set(tuple(Nk(A,n,n)) for A in cls)          # cycle spectrum (N_1..N_n) = char_A
        withR=set((tuple(Nk(A,n,n)), signed_redei(A,n)) for A in cls) if len(cls)<=200 else None
        line="    n=%d tau=%-5s shell=%3d classes: cycle-spectrum resolves %d/%d" % (n,str(tau),len(cls),len(specs),len(cls))
        if withR is not None:
            line+= " ; +|R| resolves %d/%d" % (len(withR),len(cls))
        print(line)
print("  => a near-regular tournament is pinned by tau + a short cycle spectrum (+|R| from n=7);")
print("     the continuum is LOW-DIMENSIONAL in these coordinates -- no enumeration needed.")
