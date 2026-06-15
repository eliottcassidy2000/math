#!/usr/bin/env python3
"""
skew_spectral_mechanism_monad.py -- WHY is the skew-spectrum of a tournament
SPECTRAL (a function of the ordinary spectrum)?  [monad-S7 finding]

S5/S6 background: per(xI+A) [unsigned all-length face] is NON-spectral (THM-506).
S7 finding (skew_even_face_monad.py): det(xI - S), S=A-A^T, is SPECTRAL --
it carries NO info beyond det(xI-A), verified exhaustively n<=6 incl. the
cospectral-but-different-H pairs at n=6.

Mechanism to test:  det(xI - S) = det((x-1)I - 2A) - bordered_det, where
  bordered_det(x) = det([[ (x-1)I-2A, 1],[1^T, 0]]) = -1^T adj((x-1)I-2A) 1.
det((x-1)I-2A) = 2^n charA((x-1)/2) is spectral.  So det(xI-S) spectral
<=> bordered_det spectral <=> walk-generating fn 1^T(zI-2A)^{-1}1 spectral
<=> the WALK COUNTS w_k = 1^T A^k 1 are spectral for tournaments.

Tests:
  (a) Is w_k = 1^T A^k 1 spectral (constant on cospectral classes)?  For each k.
  (b) Verify det(xI-S) = det((x-1)I-2A) - bordered_det identity exactly.
  (c) Reconcile: confirm the UNSIGNED even-real-cycle face I(Omega_even,2) IS
      non-spectral (splits the n=6 cospectral classes) -- so the SIGN is the
      whole difference (det side spectral, per/unsigned side not).
  (d) Express tr(S^{2t}) as a polynomial in tr(A^j) (find the formula t=1,2,3).

Author: monad-explorer-2026-06-15-S7
Pure-python exact integer arithmetic.
"""
from itertools import combinations
from fractions import Fraction

# ---- reuse exact int charpoly + det ----
def charpoly_int(M):
    n=len(M)
    def trace(X): return sum(X[i][i] for i in range(n))
    def matmul(X,Y): return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    def adddiag(X,s): return [[X[i][j]+(s if i==j else 0) for j in range(n)] for i in range(n)]
    Mprev=[[M[i][j] for j in range(n)] for i in range(n)]
    cc=[1,-trace(Mprev)]; cprev=cc[1]
    for k in range(2,n+1):
        Mk=matmul(M, adddiag(Mprev,cprev))
        num=-trace(Mk); assert num%k==0
        ck=num//k; cc.append(ck); Mprev=Mk; cprev=ck
    poly=[0]*(n+1)
    for k in range(n+1): poly[n-k]=cc[k]
    return poly

def det_int(M):
    n=len(M);A=[row[:] for row in M];sign=1;prev=1
    for k in range(n-1):
        if A[k][k]==0:
            sw=None
            for r in range(k+1,n):
                if A[r][k]!=0: sw=r;break
            if sw is None: return 0
            A[k],A[sw]=A[sw],A[k];sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                A[i][j]=(A[i][j]*A[k][k]-A[i][k]*A[k][j])//prev
        prev=A[k][k]
    return sign*A[n-1][n-1]

def tour_from_bits(bits,n):
    A=[[0]*n for _ in range(n)];pos=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>pos)&1:A[i][j]=1
            else:A[j][i]=1
            pos+=1
    return A
def skew(A):
    n=len(A);return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
def matvec(M,v):
    return [sum(M[i][j]*v[j] for j in range(len(v))) for i in range(len(M))]
def matmul(X,Y):
    n=len(X);m=len(Y[0]);p=len(Y)
    return [[sum(X[i][k]*Y[k][j] for k in range(p)) for j in range(m)] for i in range(n)]

def walk_counts(A, kmax):
    n=len(A); ones=[1]*n; w=[]; v=ones[:]
    w.append(sum(v))            # w_0 = 1^T I 1 = n
    for k in range(1,kmax+1):
        v=matvec(A,v); w.append(sum(v))
    return w   # w[k]=1^T A^k 1

def ham_path_count(A):
    n=len(A)
    if n==1:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    full=(1<<n)-1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if d==0:continue
            for w in range(n):
                if mask&(1<<w):continue
                if A[v][w]==1:dp[mask|(1<<w)][w]+=d
    return sum(dp[full][v] for v in range(n))

def even_face_at2(A):
    """I(Omega_even, 2) = sum over independent sets of even directed cycles of 2^{size}."""
    n=len(A)
    cyc=[]
    for start in range(n):
        stack=[(start,1<<start,1)]
        while stack:
            v,mask,length=stack.pop()
            for w in range(n):
                if A[v][w]!=1:continue
                if w==start and length>=3:
                    if length%2==0: cyc.append(mask)
                elif w>start and not(mask&(1<<w)):
                    stack.append((w,mask|(1<<w),length+1))
    m=len(cyc)
    acc=[1]
    def rec(idx,used,cnt):
        for j in range(idx,m):
            if cyc[j]&used:continue
            acc[0]+=(1<<(cnt+1))
            rec(j+1,used|cyc[j],cnt+1)
    rec(0,0,0)
    return acc[0]

def run(n, sample_limit=None):
    print(f"\n===== n={n} =====")
    total=1<<(n*(n-1)//2)
    if sample_limit is None or total<=sample_limit:
        bits_list=range(total); exh=True
    else:
        import random;random.seed(7);bits_list=[random.randrange(total) for _ in range(sample_limit)];exh=False
    kmax=n
    acp_to_wk={}     # acp -> set of tuple(w_1..w_kmax)
    acp_to_even={}   # acp -> set of I(Omega_even,2)
    acp_to_H={}
    for bits in bits_list:
        A=tour_from_bits(bits,n)
        acp=tuple(charpoly_int(A))
        wk=tuple(walk_counts(A,kmax)[1:])  # w_1..w_kmax
        acp_to_wk.setdefault(acp,set()).add(wk)
        acp_to_even.setdefault(acp,set()).add(even_face_at2(A))
        acp_to_H.setdefault(acp,set()).add(ham_path_count(A))
    # (a) walk counts spectral?
    wk_split=[a for a,s in acp_to_wk.items() if len(s)>1]
    print(f"  exhaustive={exh}; #cospectral(charA) classes total={len(acp_to_wk)}")
    print(f"  (a) walk counts w_1..w_{kmax} SPECTRAL? classes where w_k differs within charA: {len(wk_split)} "
          f"-> {'ALL w_k SPECTRAL' if not wk_split else 'NON-SPECTRAL w_k EXISTS'}")
    # per-k: find smallest k where w_k is non-spectral, if any
    for k in range(1,kmax+1):
        bad=0
        for a,s in acp_to_wk.items():
            vals=set(t[k-1] for t in s)
            if len(vals)>1:bad+=1
        if bad>0:
            print(f"      w_{k}: NON-spectral in {bad} classes")
    # (c) even face non-spectral?
    even_split=[a for a,s in acp_to_even.items() if len(s)>1]
    H_split=[a for a,s in acp_to_H.items() if len(s)>1]
    print(f"  (c) I(Omega_even,2) NON-spectral in {len(even_split)} classes; H non-spectral in {len(H_split)} classes")

if __name__=="__main__":
    for n in [3,4,5,6]:
        run(n)
    run(7, sample_limit=60000)
