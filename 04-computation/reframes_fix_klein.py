#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Corrected reframes: (B) symmetric=Galilean-invariant; (D) Paley GRAPH (odd n,1mod4) vs TOURNAMENT (even n,3mod4)."""
import math
from fractions import Fraction as F
def dist0(x,D): x%=D; return min(x,D-x)
def gap_set(diffs,Dmax):  # gap of a difference set = symmetric (all-pairs) lonely measure
    best=F(0)
    for D in range(2,Dmax+1):
        for a in range(1,D):
            m=min(dist0(d*a,D) for d in diffs)
            if F(m,D)>best: best=F(m,D)
    return best

print("="*84); print(" (B-fix) SYMMETRIC (all-observers) M IS Galilean-invariant; one-observer M is NOT"); print("="*84)
S=[1,2,5,6,7,8]; V=[0]+S                       # observer 0 is just runner number n+1
def diffs(V): return sorted({abs(a-b) for i,a in enumerate(V) for b in V[i+1:]})
M_sym=gap_set(diffs(V),60)
print(f"   V={V} (observer = a runner). pairwise-difference set {diffs(V)}")
print(f"   symmetric M (every runner lonely from ALL others) = {M_sym}={float(M_sym):.5f}")
for w in [1,3,10]:
    Vw=[v-w for v in V]
    print(f"     shift V-> V-{w} = {Vw}: diffs {diffs(Vw)} -> M_sym={gap_set(diffs(Vw),60)}  (invariant: differences unchanged)")
print("   => the SYMMETRIC version is Galilean-invariant (lives on the difference set / the translation-")
print("      invariant danger-circulant). The one-observer M=2/13 is the asymmetric (distinguished-0) value;")
print("      the symmetric all-pairs value is SMALLER/different -- copying the observer to all n points is the")
print("      STRONGER (group-invariant) condition. Same machine, two analogues.")

print("\n"+"="*84); print(" (D-fix) odd n -> Paley GRAPH (Ramanujan); even n -> Paley TOURNAMENT (Redei/OCF)"); print("="*84)
def is_pr(q): return q>1 and all(q%d for d in range(2,int(q**.5)+1))
def ham_paths_tournament(q):
    QR={(x*x)%q for x in range(1,q)}
    T=[[1 if ((j-i)%q) in QR else 0 for j in range(q)] for i in range(q)]
    # sanity: tournament iff exactly one of i->j, j->i
    ok=all((T[i][j]+T[j][i]==1) for i in range(q) for j in range(q) if i!=j)
    dp={}
    for v in range(q): dp[(1<<v,v)]=1
    for mask in range(1<<q):
        for last in range(q):
            c=dp.get((mask,last),0)
            if not c: continue
            for nx in range(q):
                if not(mask&(1<<nx)) and T[last][nx]:
                    k=(mask|(1<<nx),nx); dp[k]=dp.get(k,0)+c
    full=(1<<q)-1
    return ok, sum(dp.get((full,v),0) for v in range(q))
for q in [7,11]:        # 7=GF(7) n=4, 11=GF(11) n=6 ; both 3 mod 4 -> tournaments
    n=(q+1)//2; ok,H=ham_paths_tournament(q)
    print(f"   2n-1={q} (n={n}, {q%4} mod 4): valid Paley TOURNAMENT={ok}; #Hamiltonian paths H={H}, Redei ODD? {H%2==1}")
for q in [13,17]:       # 1 mod 4 -> Paley GRAPH (undirected), Ramanujan
    n=(q+1)//2
    QR={(x*x)%q for x in range(1,q)}; k=len(QR)
    import numpy as np
    A=np.array([[1 if ((j-i)%q) in QR else 0 for j in range(q)] for i in range(q)],float)
    ev=sorted(np.linalg.eigvalsh(A)); nt=max(abs(l) for l in ev if abs(l-k)>1e-6)
    print(f"   2n-1={q} (n={n}, {q%4} mod 4): Paley GRAPH, deg {k}; max|nontriv eig|={nt:.3f} <= 2sqrt({k-1})={2*math.sqrt(k-1):.3f}? Ramanujan={nt<=2*math.sqrt(k-1)+1e-9}")
print("   => EVEN n: 2n-1=3 mod 4 = Paley TOURNAMENT -> Redei/OCF H-count (odd); ODD n: 2n-1=1 mod 4 =")
print("      Paley GRAPH -> Ramanujan/Ihara. The even/odd split AGAIN: tournament(OCF) vs graph(Ihara).")
