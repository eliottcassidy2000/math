#!/usr/bin/env python3
"""
klein-2026-07-21-S400.  The exact SCC-composition law for disc(T) (why disc is NOT multiplicative).

disc(T) = |det(I+K)|/2^{n-1},  K=A-A^T.  For the ordered sum T = T1 => T2 (all T1 beats all T2),
Schur complement gives:
        disc(T1 => T2) = disc(T1) * disc(T2) * (1 + s1*s2)/2,
        s_i := 1^T (I+K_i)^{-1} 1   (total inverse-response of block i).
Identity: (I+K)x=1, K skew => s = x^T(I+K)x = x^T x = ||x||^2 >= 0, and s <= n (Cauchy-Schwarz).

This script:
 (1) VERIFIES the composition law on all ordered sums T1=>T2 (n<=6).
 (2) computes s over all strong tournaments n<=6: range, and does s<1 hold? (=> sub-mult with singletons)
 (3) settles the direction: is disc(T) <= prod disc(SCC_i)  (sub-mult => H>=disc REDUCES to strong)
     or can it exceed (super-mult => reduction fails, need direct strong bound with slack)?
"""
import itertools, math
import numpy as np

def all_tournaments(n):
    pairs=list(itertools.combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def Kmat(A,n): return np.array([[A[i][j]-A[j][i] for j in range(n)] for i in range(n)],dtype=float)
def disc_from_K(K,n): return abs(np.linalg.det(np.eye(n)+K))/(2**(n-1))
def disc(A,n): return disc_from_K(Kmat(A,n),n)
def s_val(A,n):
    K=Kmat(A,n); x=np.linalg.solve(np.eye(n)+K, np.ones(n)); return float(np.ones(n)@x)

def strong(A,n):
    reach=[set([i]) for i in range(n)]
    for _ in range(n):
        for i in range(n):
            for j in list(reach[i]):
                for k in range(n):
                    if A[j][k]: reach[i].add(k)
    return all(len(reach[i])==n for i in range(n))

def ordered_sum(A1,n1,A2,n2):
    n=n1+n2; A=[[0]*n for _ in range(n)]
    for i in range(n1):
        for j in range(n1): A[i][j]=A1[i][j]
    for i in range(n2):
        for j in range(n2): A[n1+i][n1+j]=A2[i][j]
    for i in range(n1):
        for j in range(n2): A[i][n1+j]=1   # all T1 -> T2
    return A,n

# (1) verify composition law
print("="*76); print("(1) COMPOSITION LAW  disc(T1=>T2) = disc(T1)disc(T2)(1+s1 s2)/2"); print("="*76)
maxerr=0.0; checks=0
small={1:[[[0]]], 2:None, 3:None}
reps={}
for n in range(1,5):
    reps[n]=list(all_tournaments(n)) if n>=1 else []
for n1 in range(1,4):
    for n2 in range(1,4):
        for A1 in reps[n1]:
            for A2 in reps[n2]:
                A,n=ordered_sum(A1,n1,A2,n2)
                lhs=disc(A,n)
                s1=s_val(A1,n1) if n1>0 else 0; s2=s_val(A2,n2)
                rhs=disc(A1,n1)*disc(A2,n2)*(1+s1*s2)/2
                err=abs(lhs-rhs); maxerr=max(maxerr,err); checks+=1
print(f"  checked {checks} ordered sums (n<=6); max |lhs-rhs| = {maxerr:.2e}  "
      f"=> {'FORMULA VERIFIED' if maxerr<1e-8 else 'FORMULA WRONG'}")

# (2) s over strong tournaments
print(); print("="*76); print("(2) s = 1^T(I+K)^{-1}1 over STRONG tournaments (and singleton=1)"); print("="*76)
print(f"  singleton s = {s_val([[0]],1):.4f}")
for n in range(3,7):
    ss=[]
    seen=set()
    for A in all_tournaments(n):
        if not strong(A,n): continue
        key=tuple(A[i][j] for i in range(n) for j in range(n))
        # dedup by canon (small n ok)
        best=None
        for p in itertools.permutations(range(n)):
            k=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or k<best: best=k
        if best in seen: continue
        seen.add(best)
        ss.append(s_val(A,n))
    print(f"  n={n}: {len(ss)} strong classes; s in [{min(ss):.4f}, {max(ss):.4f}], mean {sum(ss)/len(ss):.4f}; "
          f"#with s>=1: {sum(1 for v in ss if v>=1-1e-9)}/{len(ss)}")

# (3) direction: disc(T) vs prod disc(SCC)
print(); print("="*76); print("(3) disc(T) vs prod disc(SCC_i): sub-multiplicative (<=) => reduction to strong holds"); print("="*76)
def sccs(A,n):
    reach=[set([i]) for i in range(n)]
    for _ in range(n):
        for i in range(n):
            for j in list(reach[i]):
                for k in range(n):
                    if A[j][k]: reach[i].add(k)
    comp=[];seen=set()
    for i in range(n):
        if i in seen: continue
        c=[k for k in range(n) if k in reach[i] and i in reach[k]]
        for v in c: seen.add(v)
        comp.append(c)
    return comp
def subAmat(A,S):
    S=sorted(S);m=len(S)
    return [[A[S[a]][S[b]] for b in range(m)] for a in range(m)],m
for n in range(3,7):
    ratios=[]
    for A in all_tournaments(n):
        comp=sccs(A,n)
        if len(comp)==1: continue
        dT=disc(A,n); dp=1.0
        for c in comp:
            sa,m=subAmat(A,c); dp*=disc(sa,m)
        ratios.append(dT/dp)
    print(f"  n={n}: reducible; disc(T)/prod disc(SCC) in [{min(ratios):.4f}, {max(ratios):.4f}]  "
          f"({'SUB-mult (<=1)' if max(ratios)<=1+1e-6 else 'can EXCEED 1 => SUPER-mult, reduction FAILS'})")
print("\nDONE.")
