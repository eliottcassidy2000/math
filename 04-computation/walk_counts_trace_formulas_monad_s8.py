#!/usr/bin/env python3
"""
walk_counts_trace_formulas_monad_s8.py   (pure python, exact -- no sympy)
========================================================================
Consequences of THM-507 (walk counts spectral, F(x)=prod(x+1+lam)/prod(x-lam)-1).

 (A) Explicit walk counts as  w_k = C(n,k+1) + (spectral cycle-count corrections).
     Hand-derived from log(F+1)=sum_m q_m/(m x^m), q_m = p_m - (-1)^m sum_j C(m,j)p_j,
     using the tournament facts p_1=tr A=0 and p_2=tr A^2=0 (no 2-cycles), p_3=3c_3:
        w_0 = C(n,1) = n
        w_1 = C(n,2)                       (# arcs)
        w_2 = C(n,3) + 2 c_3
        w_3 = C(n,4) + (2n-3) c_3
     Leading term C(n,k+1) = the TRANSITIVE-tournament walk count; the corrections are
     the spectral cycle counts (c_3=tr A^3/3, ...). This script VERIFIES w_2,w_3 exactly.

 (B) Recentred reciprocity: the THM-507 involution x -> -x-1 is centred at x = -1/2.
     With G(x)=1+F(x) = (-1)^n charA(-x-1)/charA(x), G(x)G(-x-1)=1 identically;
     equivalently g(t):=G(t-1/2) satisfies g(t)g(-t)=1.  The -1/2 is the SAME centre as
     the central-factorial / W(r)=tr M(r) structure (THM-055/059/080).  Exact check at
     several rational points.

 (C) Complement fixed-point: for A'=J-I-A, F_{A'}(x)=h/(1-h), h(x)=-F_A(-x-1)
     (Sherman-Morrison).  Tournament: A'=A^T, F_{A^T}=F_A, so this self-reference IS the
     reciprocity -- the tournament walk function is the complement-operation fixed point.
     Exact check at rational points.

Author: monad-explorer-2026-06-15-S8.
"""
from fractions import Fraction as Fr
from math import comb
import random

def matmul(X, Y):
    n=len(X); m=len(Y[0]); p=len(Y)
    return [[sum(X[i][k]*Y[k][j] for k in range(p)) for j in range(m)] for i in range(n)]

def charpoly_int(M):
    """det(xI-M) low->high coeffs (Faddeev-LeVerrier)."""
    n=len(M)
    Mprev=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    c=[1]
    for k in range(1,n+1):
        AM=matmul(M,Mprev)
        tr=sum(AM[i][i] for i in range(n))
        ck=-tr//k
        assert -tr%k==0
        c.append(ck)
        Mprev=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return c[::-1]   # [x^0..x^n]

def poly_eval(coeffs, x):
    s=Fr(0)
    for c in reversed(coeffs):
        s=s*x+c
    return s

def rand_tour(nn,rng):
    A=[[0]*nn for _ in range(nn)]
    for i in range(nn):
        for j in range(i+1,nn):
            if rng.random()<0.5: A[i][j]=1
            else: A[j][i]=1
    return A

def walk_counts(A,K):
    n=len(A); vk=[1]*n; w=[]
    for k in range(K+1):
        w.append(sum(vk))
        vk=[sum(A[i][j]*vk[j] for j in range(n)) for i in range(n)]
    return w

def F_at(A, x0):
    """F(x0)=1^T (x0 I - A)^{-1} 1 exactly via Fraction Gaussian elimination."""
    n=len(A)
    # augmented [x0 I - A | 1]
    Mat=[[Fr(x0 if i==j else 0)-A[i][j] for j in range(n)]+[Fr(1)] for i in range(n)]
    # solve (x0 I - A) y = 1 ; then F = 1^T y
    for col in range(n):
        piv=None
        for r in range(col,n):
            if Mat[r][col]!=0: piv=r; break
        Mat[col],Mat[piv]=Mat[piv],Mat[col]
        pv=Mat[col][col]
        Mat[col]=[v/pv for v in Mat[col]]
        for r in range(n):
            if r!=col and Mat[r][col]!=0:
                f=Mat[r][col]
                Mat[r]=[Mat[r][k]-f*Mat[col][k] for k in range(n+1)]
    y=[Mat[i][n] for i in range(n)]
    return sum(y)

def banner(s): print("\n"+"="*70+"\n "+s+"\n"+"="*70)

def main():
    rng=random.Random(20260615)

    banner("(A) verify w_2 = C(n,3)+2c_3 and w_3 = C(n,4)+(2n-3)c_3 (exact)")
    bad=0; tot=0
    for n in range(3,11):
        for _ in range(150):
            A=rand_tour(n,rng)
            w=walk_counts(A,3)
            # c_3 = tr(A^3)/3
            A2=matmul(A,A); A3=matmul(A2,A); trA3=sum(A3[i][i] for i in range(n))
            assert trA3%3==0; c3=trA3//3
            pred2=comb(n,3)+2*c3
            pred3=comb(n,4)+(2*n-3)*c3
            tot+=1
            if not (w[0]==n and w[1]==comb(n,2) and w[2]==pred2 and w[3]==pred3):
                bad+=1
    print(f"  checked {tot} random tournaments n=3..10 -- failures: {bad}")
    print("  => w_k = (transitive count C(n,k+1)) + spectral cycle-count corrections.")

    banner("(B) recentred reciprocity (1+F(x))(1+F(-x-1))=1, centre x=-1/2 (exact)")
    A=rand_tour(6,rng); n=len(A)
    cA=charpoly_int(A)
    pts=[Fr(2),Fr(3),Fr(-5,2),Fr(7,3),Fr(-1,2)+Fr(11,7)]  # avoid eigenvalues
    okB=True
    for x0 in pts:
        try:
            Fx=F_at(A,x0); Fy=F_at(A,-x0-1)
        except ZeroDivisionError:
            continue
        prod=(1+Fx)*(1+Fy)
        # also closed-form G(x)=(-1)^n charA(-x-1)/charA(x)
        G=((-1)**n)*poly_eval(cA,-x0-1)/poly_eval(cA,x0)
        if prod!=1 or (1+Fx)!=G: okB=False
        print(f"    x={x0}:  (1+F(x))(1+F(-x-1)) = {prod};  (1+F(x))-G(x) = {(1+Fx)-G}")
    print(f"  reciprocity + closed form hold at all tested points: {okB}")
    print("  fixed point of x->-x-1 is x=-1/2 (same -1/2 as central-factorial structure).")

    banner("(C) complement fixed-point  F_{A'}=h/(1-h), h=-F_A(-x-1), A'=A^T (exact)")
    At=[[A[j][i] for j in range(n)] for i in range(n)]
    okC=True
    for x0 in [Fr(2),Fr(5,2),Fr(-3),Fr(9,4)]:
        FA=F_at(A,x0); FAt=F_at(At,x0)
        h=-F_at(A,-x0-1)
        rhs=h/(1-h)
        if FAt!=rhs or FAt!=FA: okC=False
        print(f"    x={x0}:  F_A={FA}  F_A^T={FAt}  h/(1-h)={rhs}  (all equal: {FA==FAt==rhs})")
    print(f"  complement relation + tournament self-reference hold: {okC}")
    print("\n  Tournament walk function = fixed point of the graph-complement walk map;")
    print("  the reciprocity is its fixed-point equation.")

if __name__=="__main__":
    main()
