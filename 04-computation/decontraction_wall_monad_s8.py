#!/usr/bin/env python3
"""
decontraction_wall_monad_s8.py
==============================
THM-507: the FULLY-CONTRACTED resolvent 1^T(xI-A)^{-1}1 is spectral (= walk gen fn).
Also trivially spectral: tr((xI-A)^{-1}) = charA'(x)/charA(x).

QUESTION (handoff item #1): as we DE-CONTRACT toward the pointed data R_ab, where does
spectrality break?  Test, on genuine cospectral-but-different-H classes at n=6, which of
these resolvent functionals are constant within a cospectral class (= spectral) at a fixed
rational x:

  (a) S1   = 1^T R 1                  (full contraction)            -- PROVED spectral
  (b) trR  = sum_a R_aa = charA'/charA (trace)                      -- spectral (log-deriv)
  (c) M2   = sum_a (R1)_a^2           (2nd moment of resolvent row sums r=R1)
  (d) M3   = sum_a (R1)_a^3           (3rd moment)
  (e) Frob = sum_{a,b} R_ab^2 = tr(R R^T)                           (Frobenius)

R = (xI-A)^{-1}, evaluated exactly (Fraction) at x = 10.

Analogy to score moments: sum s_i^2 is spectral (Moon, via c_3) but sum s_i^p (p>=3) is
NOT.  This probes the resolvent-space analogue and pins the de-contraction wall.

Author: monad-explorer-2026-06-15-S8.
"""
from fractions import Fraction as Fr

def charpoly_int(M):
    n=len(M)
    def mm(X,Y): return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    Mprev=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    c=[1]
    for k in range(1,n+1):
        AM=mm(M,Mprev); tr=sum(AM[i][i] for i in range(n)); ck=-tr//k
        c.append(ck); Mprev=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(c[::-1])

def ham_paths(A):
    """# directed Hamiltonian paths (H = Redei count) via bitmask DP."""
    n=len(A); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            cur=dp[mask][v]
            if not cur: continue
            for w in range(n):
                if not (mask>>w)&1 and A[v][w]:
                    dp[mask|(1<<w)][w]+=cur
    return sum(dp[full][v] for v in range(n))

def tour_from_bits(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A

def resolvent_funcs(A,x0):
    n=len(A)
    # invert (x0 I - A) exactly; augment with I
    M=[[Fr(x0 if i==j else 0)-A[i][j] for j in range(n)]+[Fr(1 if i==j else 0) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv=next(r for r in range(col,n) if M[r][col]!=0)
        M[col],M[piv]=M[piv],M[col]
        pv=M[col][col]; M[col]=[v/pv for v in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                f=M[r][col]; M[r]=[M[r][k]-f*M[col][k] for k in range(2*n)]
    R=[[M[i][n+j] for j in range(n)] for i in range(n)]
    S1=sum(R[i][j] for i in range(n) for j in range(n))
    trR=sum(R[i][i] for i in range(n))
    r=[sum(R[i][j] for j in range(n)) for i in range(n)]   # R1
    M2=sum(ri*ri for ri in r); M3=sum(ri*ri*ri for ri in r)
    Frob=sum(R[i][j]*R[i][j] for i in range(n) for j in range(n))
    return S1,trR,M2,M3,Frob

def main():
    n=6; m=n*(n-1)//2; x0=10
    print(f"Grouping all {1<<m} labeled n={n} tournaments by charA; finding cospectral")
    print("classes that carry >1 distinct H (genuine non-spectral classes)...")
    by_char={}
    for bits in range(1<<m):
        A=tour_from_bits(n,bits)
        by_char.setdefault(charpoly_int(A),[]).append(bits)
    # find classes with >1 distinct H
    interesting=[]
    for ch,members in by_char.items():
        Hs={}
        for b in members:
            A=tour_from_bits(n,b); h=ham_paths(A); Hs.setdefault(h,b)
        if len(Hs)>1:
            interesting.append((ch,Hs))
    print(f"  {len(by_char)} cospectral classes; {len(interesting)} carry >1 distinct H.\n")

    labels=["S1=1^T R 1","trR=charA'/charA","M2=sum (R1)_a^2","M3=sum (R1)_a^3","Frob=||R||_F^2"]
    # for each interesting class, test constancy of each functional across the distinct-H reps
    spectral_votes=[0]*5; nonspectral_votes=[0]*5; ntest=0
    for ch,Hs in interesting:
        reps=list(Hs.values())
        vals=[resolvent_funcs(tour_from_bits(n,b),x0) for b in reps]
        ntest+=1
        for f in range(5):
            col=set(v[f] for v in vals)
            if len(col)==1: spectral_votes[f]+=1
            else: nonspectral_votes[f]+=1
    print(f"Across {ntest} genuine (different-H) cospectral classes, at x={x0}:")
    print(f"  {'functional':22s}  {'constant-in-class':>17s}  {'splits':>7s}  verdict")
    for f in range(5):
        verdict = "SPECTRAL" if nonspectral_votes[f]==0 else "NON-spectral"
        print(f"  {labels[f]:22s}  {spectral_votes[f]:>17d}  {nonspectral_votes[f]:>7d}  {verdict}")
    print("\nReading: the de-contraction wall sits between the contraction/trace (spectral)")
    print("and the higher resolvent moments -- the resolvent-space echo of 'sum s_i^2")
    print("spectral, sum s_i^p (p>=3) not'.  This locates handoff item #1 concretely.")

if __name__=="__main__":
    main()
