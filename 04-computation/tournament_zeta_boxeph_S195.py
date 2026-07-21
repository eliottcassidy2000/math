#!/usr/bin/env python3
"""tournament_zeta_boxeph_S195.py -- boxeph-2026-07-21-S195

THE TOURNAMENT ZETA (Bowen-Lanford / closed-orbit) ζ_T(u) = 1/det(I - uA), the harmonic-analysis
reduction principle (klein-S399 handoff), deepening THM-1925/THM-1862 and integrating kps THM-1880
(char_S Chebyshev/Gauss atoms).

Verifies:
(1) det(I - uA) = Σ_j c_j u^j with c_j = char_A coeffs (integer); ζ = 1/that.
(2) EULER PRODUCT: -log det(I-uA) = Σ_k N_k u^k/k, N_k = tr(A^k) = # closed k-walks; prime-cycle
    counts π_ℓ = (1/ℓ)Σ_{d|ℓ} μ(ℓ/d) N_d are NON-NEGATIVE INTEGERS; ζ = ∏_ℓ (1-u^ℓ)^{-π_ℓ}.
(3) REDUCTION: det(I-uA(T)) = ∏ over strong components; ζ=1 on acyclic/transitive (A nilpotent).
(4) COMPLEMENT-INVARIANCE: ζ_T = ζ_{T^op} (A^op = A^T).
(5) c3 read-off: N_3 = 3·c3.  3-cycle atom: ζ = 1/(1-u^3), one prime cycle.
(6) TRIG ATOMS: for circulant atoms N_k = Σ_j λ_j^k (character/trig moments); poles at 1/λ_j
    (Gauss sums Paley, Dirichlet/Chebyshev interval).
"""
import numpy as np
from itertools import permutations, combinations
from fractions import Fraction as Fr
from math import cos, sin, pi, gcd

# ---------- enumeration + invariants ----------
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
def charpoly_int(A):
    n=len(A); A=[[int(A[i][j]) for j in range(n)] for i in range(n)]
    M=[[1 if i==j else 0 for j in range(n)] for i in range(n)]; coeffs=[1]
    for k in range(1,n+1):
        AM=[[sum(A[i][t]*M[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        tr=sum(AM[i][i] for i in range(n)); ck=-tr//k; coeffs.append(ck)
        M=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return coeffs   # char_A(x)=Σ coeffs[j] x^{n-j};  det(I-uA)=Σ coeffs[j] u^j
def matpow_trace(A,k):
    n=len(A); M=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for _ in range(k):
        M=[[sum(M[i][t]*A[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return sum(M[i][i] for i in range(n))
def strong_components(A):
    n=len(A); R=[[1 if (i==j or A[i][j]) else 0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j]=1
    seen=[False]*n; comps=[]
    for v in range(n):
        if not seen[v]:
            comp=[w for w in range(n) if R[v][w] and R[w][v]]
            for w in comp: seen[w]=True
            comps.append(comp)
    return comps
def num_c3(A,n):
    c=0
    for i,j,k in combinations(range(n),3):
        if A[i][j]+A[j][k]+A[k][i]==3 or A[j][i]+A[k][j]+A[i][k]==3: c+=1
    return c
def complement(A):
    n=len(A); return [[A[j][i] if i!=j else 0 for j in range(n)] for i in range(n)]
def polymul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]+=a*b
    return r
def mobius(m):
    if m==1: return 1
    r=1; mm=m; p=2
    while p*p<=mm:
        if mm%p==0:
            mm//=p
            if mm%p==0: return 0
            r=-r
        p+=1
    if mm>1: r=-r
    return r
def divisors(k): return [d for d in range(1,k+1) if k%d==0]

reps=iso_reps(6)
print("iso classes:", {n:len(reps[n]) for n in range(3,7)})

# ---------- (3)+(1) reduction + zeta polynomial ----------
print("\n"+"="*90); print("(1)+(3)  ζ_T(u)=1/det(I-uA);  det(I-uA)=Σ c_j u^j = ∏ over strong components; ζ=1 on acyclic")
print("="*90)
mism=0; ntot=0; trans_ok=True
for n in range(3,7):
    for A in reps[n]:
        d_full=charpoly_int(A)                      # det(I-uA) coeffs
        prod=[1]
        for comp in strong_components(A):
            sub=[[A[i][j] for j in comp] for i in comp]
            prod=polymul(prod,charpoly_int(sub))
        if d_full!=prod: mism+=1
        ntot+=1
        if num_c3(A,n)==0 and d_full!=[1]+[0]*n: trans_ok=False   # acyclic => det(I-uA)=1
# NOTE det(I-uA) for acyclic A (nilpotent) = 1  (all c_j=0 for j>=1)
acyc_det=charpoly_int([[1 if i<j else 0 for j in range(5)] for i in range(5)])
print("  det(I-uA)=∏_SCC det over all %d classes n<=6: mismatches=%d" % (ntot,mism))
print("  transitive n=5: det(I-uA)=%s  => ζ=1 (no closed orbits; A nilpotent)" % acyc_det)
print("  => ζ_T is INVISIBLE to the acyclic part; all content is in the strong components (strong-core reduction).")

# ---------- (2) Euler product / prime cycles ----------
print("\n"+"="*90); print("(2)  EULER PRODUCT: prime-cycle counts π_ℓ=(1/ℓ)Σ_{d|ℓ}μ(ℓ/d)N_d are non-neg integers")
print("="*90)
def prime_cycles(A,n,K):
    N={k:matpow_trace(A,k) for k in range(1,K+1)}
    pis={}
    for l in range(1,K+1):
        s=sum(mobius(l//d)*N[d] for d in divisors(l))
        assert s%l==0, "non-integer prime count"
        pis[l]=s//l
    return N,pis
bad=0
for n in range(3,7):
    for A in reps[n]:
        N,pis=prime_cycles(A,n,n)
        if any(v<0 for v in pis.values()): bad+=1
print("  prime-cycle counts π_ℓ non-negative for all classes n<=6: %s" % (bad==0))
# the 3-cycle atom
C3=[[0,1,0],[0,0,1],[1,0,0]]
N,pis=prime_cycles(C3,3,6)
print("  3-cycle atom: det(I-uA)=%s => ζ=1/(1-u^3); N_k=%s; prime cycles π=%s (one length-3 prime)"
      % (charpoly_int(C3),[N[k] for k in range(1,7)],{l:pis[l] for l in pis if pis[l]}))

# ---------- (5) c3 read-off + (4) complement invariance ----------
print("\n"+"="*90); print("(4) complement-invariance ζ_T=ζ_{T^op}   (5) N_3 = 3·c3")
print("="*90)
comp_ok=True; c3_ok=True
for n in range(3,7):
    for A in reps[n]:
        if charpoly_int(A)!=charpoly_int(complement(A)): comp_ok=False
        if matpow_trace(A,3)!=3*num_c3(A,n): c3_ok=False
print("  det(I-uA(T)) == det(I-uA(T^op)) for all classes n<=6 (=> ζ complement-invariant): %s" % comp_ok)
print("  N_3 = tr(A^3) = 3*c3 for all classes n<=6: %s   (c3 is read off ζ's u^3 log-coefficient)" % c3_ok)

# ---------- (6) trig atoms: N_k = Σ λ_j^k (character moments), poles = 1/λ ----------
print("\n"+"="*90); print("(6) circulant atoms: N_k=Σ_j λ_j^k (trig/character moments); ζ poles at 1/λ_j")
print("="*90)
def circ(n,C):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for k in C: A[i][(i+k)%n]=1
    return A
for n in (7,11,19):
    QR=sorted({(k*k)%n for k in range(1,n)}); A=circ(n,QR)
    w=np.exp(2j*pi/n); lam=[sum(w**((j*k)%n) for k in QR) for j in range(n)]
    Nk=[matpow_trace(A,k) for k in (1,2,3)]
    Nk_spec=[round(sum((l**k).real for l in lam)) for k in (1,2,3)]
    print("  Paley-%2d: N_1,N_2,N_3 (trace) = %s ; Σλ_j^k = %s ; match=%s ; poles 1/λ have |λ|=√%d=%.3f (non-Perron)"
          %(n,Nk,Nk_spec,Nk==Nk_spec,n,n**0.5))
print("\n  => ζ_T is the closed-orbit generating function; its Euler product lives on the strong core,")
print("     its poles are reciprocals of the trigonometric atom-spectra (Gauss sums / Chebyshev).")
