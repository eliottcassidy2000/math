#!/usr/bin/env python3
"""
S713 — Even/odd and the Pfaffian: the Pfaffian as the algebraic incarnation of the 2-adic seam.

Setup (canon THM-174/120/213): tournament T on n vtx, adjacency A (A_ij=1 iff i->j), skew-adjacency
S = A - A^T (S_ij = +1 iff i->j, -1 iff j->i, 0 diag). The Pfaffian Pf(S) is defined only in EVEN
dimension; det(S)=Pf(S)^2 (even n), det(S)=0 (odd n). THM-174: det(I+2A)=Pf(S)^2 (even), =(1^T w)^2
(odd), always a perfect square.

This session pins the EVEN/ODD structure the Pfaffian encodes, and adds what canon left implicit:

(A) PARITY LAW (the even-n sibling of Redei). Pf(S_T) is ALWAYS ODD for every tournament on even n.
    Reason: S = J - I (mod 2), Pf over GF(2) = #perfect-matchings(K_n) mod 2 = (n-1)!! mod 2 = 1.
    Verify across ALL tournaments n=4 (64) and n=6 (32768); cross-check |Pf| values vs THM-120 {1,3,5,7,9}.

(B) ODD-n KERNEL LADDER. For odd n, S is singular (rank n-1); the cofactor vector w_i=(-1)^i Pf(S_hat_i)
    spans ker S. Its entries are Pfaffians of the EVEN (n-1)-vertex deleted subtournaments => ALL ODD.
    So 1^T w = sum of n (odd-many) odd numbers = ODD, and det(I+2A)=(1^T w)^2 = an ODD perfect square.
    The even/odd seam is a LADDER: delete a vertex, odd<->even, Pf<->kernel. Verify n=3 (8), n=5 (1024).

(C) UNIFIED: det(I+2A) is always an ODD PERFECT SQUARE (=> = 1 mod 8). Even n: Pf^2, Pf odd.
    Odd n: (1^T w)^2, 1^T w odd. Verify mod 8.

(D) COMPLEMENT MOD-4 PHASE. T->T^op sends S->-S, Pf(-S)=(-1)^(n/2) Pf(S): n=0 mod4 invariant,
    n=2 mod4 sign-flip. The reciprocity/XNOR mod-4 phase, now on the Pfaffian. Verify + self-complement.

(E) Pf vs H (#Hamiltonian paths, Redei-odd). Both odd; explore the joint law / congruences (n=4, n=6 sample).

No numpy/sympy. Exact integer arithmetic.
"""
import random, itertools
from math import gcd

# ---------- Pfaffian (exact integer, recursive) ----------
def pf(M):
    n=len(M)
    if n==0: return 1
    if n%2==1: return 0
    total=0
    row0=M[0]
    for j in range(1,n):
        a=row0[j]
        if a==0: continue
        rest=[k for k in range(1,n) if k!=j]
        sub=[[M[r][c] for c in rest] for r in rest]
        total += ((-1)**(j-1))*a*pf(sub)
    return total

# ---------- integer determinant (Bareiss, fraction-free) ----------
def det_int(Min):
    M=[row[:] for row in Min]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            piv=None
            for i in range(k+1,n):
                if M[i][k]!=0: piv=i; break
            if piv is None: return 0
            M[k],M[piv]=M[piv],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]

# ---------- tournament enumeration ----------
def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def skew(A):
    n=len(A); return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]

def ident_plus_2A(A):
    n=len(A); return [[(1 if i==j else 0)+2*A[i][j] for j in range(n)] for i in range(n)]

def ham_paths(A):
    n=len(A); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if not d or not (mask>>v)&1: continue
            for w in range(n):
                if (mask>>w)&1: continue
                if A[v][w]: dp[mask|(1<<w)][w]+=d
    full=(1<<n)-1
    return sum(dp[full][v] for v in range(n))

def matvec(M,v):
    return [sum(M[i][k]*v[k] for k in range(len(v))) for i in range(len(M))]

def kernel_cofactor(S):
    n=len(S); w=[]
    for i in range(n):
        rest=[k for k in range(n) if k!=i]
        sub=[[S[r][c] for c in rest] for r in rest]
        w.append(((-1)**i)*pf(sub))
    return w

# ============================ RUN ============================
if __name__=="__main__":
    print("="*80)
    print("S713 — Even/Odd and the Pfaffian: the algebraic incarnation of the 2-adic seam")
    print("="*80)

    # (A) parity law, even n
    print("\n(A) PARITY LAW: Pf(S_T) is ALWAYS ODD (even n) — the Pfaffian sibling of Redei's H-odd")
    for n in (4,6):
        cnt=0; odd=0; absvals={}; detmatch=0
        for A in all_tournaments(n):
            S=skew(A); p=pf(S); cnt+=1
            if p%2!=0: odd+=1
            absvals[abs(p)]=absvals.get(abs(p),0)+1
            if det_int(S)==p*p: detmatch+=1
        print(f"  n={n}: {cnt} tournaments | Pf odd: {odd}/{cnt} | det(S)=Pf^2: {detmatch}/{cnt} | "
              f"|Pf| values: {dict(sorted(absvals.items()))}")
    print("  (n=6 |Pf| in {1,3,5,7,9} reproduces THM-120; all odd.)")

    # (B) odd-n kernel ladder
    print("\n(B) ODD-n KERNEL LADDER: w_i=(-1)^i Pf(S_hat_i) spans ker S; entries are even-subtournament")
    print("    Pfaffians => ALL ODD; 1^T w odd; det(I+2A)=(1^T w)^2.")
    for n in (3,5):
        cnt=0; ker_ok=0; all_odd=0; sum_odd=0; det_ok=0
        for A in all_tournaments(n):
            S=skew(A); w=kernel_cofactor(S); cnt+=1
            if all(x==0 for x in matvec(S,w)): ker_ok+=1
            if all(x%2!=0 for x in w): all_odd+=1
            s=sum(w)
            if s%2!=0: sum_odd+=1
            if det_int(ident_plus_2A(A))==s*s: det_ok+=1
        print(f"  n={n}: {cnt} | S.w=0: {ker_ok}/{cnt} | all w_i odd: {all_odd}/{cnt} | "
              f"1^Tw odd: {sum_odd}/{cnt} | det(I+2A)=(1^Tw)^2: {det_ok}/{cnt}")

    # (C) unified: det(I+2A) is an odd perfect square == 1 mod 8
    print("\n(C) UNIFIED: det(I+2A) is an ODD PERFECT SQUARE (=> = 1 mod 8) for all n")
    for n in (3,4,5,6):
        cnt=0; mod8=0; sq=0
        for A in all_tournaments(n):
            d=det_int(ident_plus_2A(A)); cnt+=1
            r=int(round(d**0.5))
            if r*r==d: sq+=1
            if d%8==1: mod8+=1
        print(f"  n={n}: det(I+2A) perfect square {sq}/{cnt} | =1 mod 8: {mod8}/{cnt}")

    # (D) complement mod-4 phase
    print("\n(D) COMPLEMENT MOD-4 PHASE: Pf(-S)=(-1)^(n/2) Pf(S)  (T->T^op)")
    for n in (4,6):
        cnt=0; ok=0; selfc=0
        for A in all_tournaments(n):
            S=skew(A); p=pf(S); pneg=pf([[-x for x in row] for row in S]); cnt+=1
            if pneg==((-1)**(n//2))*p: ok+=1
        phase="invariant (n=0 mod4)" if (n//2)%2==0 else "SIGN-FLIP (n=2 mod4)"
        print(f"  n={n}: Pf(-S)=(-1)^(n/2)Pf(S) holds {ok}/{cnt}; phase=(-1)^{n//2} => {phase}")

    # (E) Pf vs H
    print("\n(E) Pf vs H (#Hamiltonian paths, Redei-odd): joint law / congruences")
    # n=4 full
    for n in (4,):
        joint={}; both_odd=0; cnt=0; cong4={}
        for A in all_tournaments(n):
            S=skew(A); p=abs(pf(S)); h=ham_paths(A); cnt+=1
            if p%2 and h%2: both_odd+=1
            joint[(p,h)]=joint.get((p,h),0)+1
            cong4[(p%4,h%4)]=cong4.get((p%4,h%4),0)+1
        print(f"  n=4 ({cnt}): both odd {both_odd}/{cnt}; (|Pf|,H) joint {dict(sorted(joint.items()))}")
        print(f"        (|Pf| mod4, H mod4) distribution: {dict(sorted(cong4.items()))}")
    # n=6 sample
    rng=random.Random(0); N=4000; joint={}; both_odd=0
    pairs=[(i,j) for i in range(6) for j in range(i+1,6)]
    for _ in range(N):
        bits=rng.randrange(1<<15); A=[[0]*6 for _ in range(6)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        S=skew(A); p=abs(pf(S)); h=ham_paths(A)
        if p%2 and h%2: both_odd+=1
        joint[(p,h%8)]=joint.get((p,h%8),0)+1
    print(f"  n=6 (sample {N}): both odd {both_odd}/{N}; (|Pf|, H mod 8) joint {dict(sorted(joint.items()))}")
    print("="*80)
