#!/usr/bin/env python3
"""trig_reduction_tournament_boxeph_S194.py -- boxeph-2026-07-21-S194

THE REDUCTION IS A PRODUCT OF TRIGONOMETRIC FUNCTIONS (THM-1875). Three legs:

(1) SPECTRAL REDUCTION: char_A multiplicative under order-join T1|>T2 (block-triangular), so a
    tournament's adjacency spectrum is the disjoint union of its strong components' spectra. The
    spectral companion of THM-1862 / THM-1830.  [verified over iso classes n<=6]
(2) SINE-PRODUCT BRIDGE: Sum_T (-1)^{back(T)} prod x_k^{score_k} = prod_{i<j}(x_j - x_i)
    (mac-mini's signed involution engine = the Vandermonde). On the unit circle x_k=e^{i th_k}
    this is prod_{i<j} 2i sin((th_j-th_i)/2) e^{i(th_i+th_j)/2} -- a PRODUCT OF SINES; the
    transitive-core reduction IS a trigonometric factorization. At roots of unity = cyclotomic
    discriminant.  [verified numerically]
(3) TRIGONOMETRIC ATOMS: circulant (rotational) strong tournaments have eigenvalues that are
    explicit trig sums -- Gauss sums (Re=-1/2) for Paley (quadratic-residue connection set), and
    Dirichlet/Chebyshev-kernel values for interval connection sets.  [verified numerically]
"""
import numpy as np
from itertools import combinations, permutations, product
from math import cos, sin, pi, sqrt

np.set_printoptions(precision=4, suppress=True)

# ---------------- helpers ----------------
def from_bits(bits, n):
    A = np.zeros((n, n), dtype=int); idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits >> idx & 1: A[i, j] = 1
            else:               A[j, i] = 1
            idx += 1
    return A

def canon(A, n, perms):
    best = None
    for p in perms:
        c = tuple(A[p[i], p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or c < best: best = c
    return best

def iso_reps(nmax):
    reps = {1: [np.zeros((1,1), dtype=int)]}
    for n in range(2, nmax+1):
        perms = list(permutations(range(n))); seen=set(); out=[]
        for B in reps[n-1]:
            for pat in range(1 << (n-1)):
                A = np.zeros((n,n), dtype=int); A[:n-1,:n-1]=B
                for k in range(n-1):
                    if pat>>k&1: A[n-1,k]=1
                    else:        A[k,n-1]=1
                c = canon(A,n,perms)
                if c not in seen:
                    seen.add(c)
                    M=np.zeros((n,n),dtype=int); idx=0
                    for i in range(n):
                        for j in range(n):
                            if i!=j: M[i,j]=c[idx]; idx+=1
                    out.append(M)
        reps[n]=out
    return reps

def order_join(A1, A2):
    n1,n2=len(A1),len(A2); n=n1+n2; B=np.zeros((n,n),dtype=int)
    B[:n1,:n1]=A1; B[n1:,n1:]=A2; B[:n1,n1:]=1
    return B

def strong_components(A):
    n=len(A)
    R=[[1 if (i==j or A[i][j]) else 0 for j in range(n)] for i in range(n)]  # Warshall closure
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

# ============================================================ LEG 1
print("="*88); print("LEG 1  SPECTRAL REDUCTION: char_A(T1|>T2)=char_A(T1)*char_A(T2); spec(T)=U spec(SCC)")
print("="*88)
reps = iso_reps(6)
print("iso classes:", {n: len(reps[n]) for n in range(3,7)})
# (a) join multiplicativity of char poly
maxerr=0.0
for A1 in reps[3]:
    for A2 in reps[3]:
        B=order_join(A1,A2)
        e_join=np.sort_complex(np.linalg.eigvals(B))
        e_union=np.sort_complex(np.concatenate([np.linalg.eigvals(A1),np.linalg.eigvals(A2)]))
        maxerr=max(maxerr, np.max(np.abs(e_join-e_union)))
print("  (a) eig(T1|>T2) == eig(T1) U eig(T2) for all n=3 x n=3 joins:  max err = %.2e" % maxerr)
# (b) char_A(T) == product of char_A(strong components), EXACT integer arithmetic (Faddeev-LeVerrier)
def charpoly_int(A):
    """exact integer coeffs [1, c1, ..., cn] of det(xI - A) via Faddeev-LeVerrier."""
    n=len(A); A=[[int(A[i][j]) for j in range(n)] for i in range(n)]
    M=[[1 if i==j else 0 for j in range(n)] for i in range(n)]  # M0 = I
    coeffs=[1]
    for k in range(1,n+1):
        AM=[[sum(A[i][t]*M[t][j] for t in range(n)) for j in range(n)] for i in range(n)]  # A*M_{k-1}
        tr=sum(AM[i][i] for i in range(n))
        ck=-tr//k                        # exact
        coeffs.append(ck)
        M=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]  # M_k = A*M_{k-1}+ck I
    return coeffs
def polymul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]+=a*b
    return r
nchk=0; fails=0
for n in range(3,7):
    for A in reps[n]:
        comps=strong_components(A)
        full=charpoly_int(A)
        prod=[1]
        for comp in comps:
            sub=[[int(A[i][j]) for j in comp] for i in comp]
            prod=polymul(prod, charpoly_int(sub))
        if full!=prod: fails+=1
        nchk+=1
print("  (b) char_A(T) == prod char_A(strong components), EXACT over all %d classes n<=6:  mismatches = %d" % (nchk,fails))
print("  => the SPECTRAL atoms are the strongly-connected tournaments (block-triangular reducibility).")
print("     [numpy float eigvals err ~sqrt(3) on nilpotent transitive blocks is ill-conditioning, not a real gap]")

# ============================================================ LEG 2
print("\n"+"="*88); print("LEG 2  SINE-PRODUCT: Sum_T (-1)^{back} x^{score} = prod(x_j-x_i) = product of sines on circle")
print("="*88)
def signed_partition(x):
    n=len(x); m=n*(n-1)//2; total=0j
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<m):
        score=[0]*n; back=0
        for t,(i,j) in enumerate(pairs):
            if bits>>t&1:   score[j]+=1               # j beats i : winner j (larger index) = forward
            else:           score[i]+=1; back+=1      # i beats j : back-arc (smaller index wins)
        term=(-1)**back
        val=term
        for k in range(n): val*= x[k]**score[k]
        total+=val
    return total
def vandermonde(x):
    n=len(x); p=1+0j
    for i in range(n):
        for j in range(i+1,n): p*=(x[j]-x[i])
    return p
# generic complex point
for n in (3,4,5):
    x=np.array([complex(cos(1.3*k+0.2), sin(0.7*k-0.4))*(1.0+0.15*k) for k in range(n)])
    lhs=signed_partition(x); rhs=vandermonde(x)
    print("  n=%d  |Sum_T(-1)^back x^score - prod(x_j-x_i)| = %.2e" % (n, abs(lhs-rhs)))
# on the unit circle => product of sines
n=5; th=np.array([0.3,1.1,2.0,3.4,5.0]); x=np.exp(1j*th)
lhs=signed_partition(x)
sineprod=1+0j
for i in range(n):
    for j in range(i+1,n):
        sineprod*= 2j*sin((th[j]-th[i])/2)*np.exp(1j*(th[i]+th[j])/2)
print("  circle n=5: |signed_partition(e^{i th}) - prod 2i sin((th_j-th_i)/2) e^{i(th_i+th_j)/2}| = %.2e"%abs(lhs-sineprod))
# at n-th roots of unity => cyclotomic discriminant modulus n^{n/2}
for n in (5,7):
    x=np.array([np.exp(2j*pi*k/n) for k in range(n)])
    v=abs(vandermonde(x))
    print("  roots of unity n=%d: |prod(w^j - w^i)| = %.4f   (n^{n/2} = %.4f)  [=sqrt|disc(x^n-1)|]"%(n,v,n**(n/2)))

# ============================================================ LEG 3
print("\n"+"="*88); print("LEG 3  TRIGONOMETRIC ATOMS: circulant tournament eigenvalues as explicit trig sums")
print("="*88)
def circulant_tournament(n, C):
    A=np.zeros((n,n),dtype=int)
    for i in range(n):
        for k in C: A[i,(i+k)%n]=1
    return A
def is_tournament(A):
    n=len(A)
    return all(A[i,j]+A[j,i]==1 for i in range(n) for j in range(n) if i!=j)
# Paley (quadratic residues), n=7,11,19  (need n==3 mod 4 for a tournament)
for n in (7,11,19):
    QR=sorted({(k*k)%n for k in range(1,n)})
    A=circulant_tournament(n,QR)
    assert is_tournament(A), "Paley %d not a tournament"%n
    ev=np.linalg.eigvals(A.astype(float))
    reparts=sorted(set(round(e.real,4) for e in ev))
    # Gauss-sum eigenvalue: lam_j = sum_{k in QR} w^{jk}; for j!=0 Re = -1/2
    w=np.exp(2j*pi/n); lam=[sum(w**((j*k)%n) for k in QR) for j in range(n)]
    reps_nonperron=sorted(set(round(l.real,4) for j,l in enumerate(lam) if j!=0))
    print("  Paley-%d: Re(eigenvalues) = %s ; Gauss-sum Re for j!=0 = %s  (theory: -0.5)"
          %(n, reparts, reps_nonperron))
# interval connection set {1..m}, m=(n-1)/2 => Dirichlet-kernel eigenvalues
for n in (7,9,11):
    m=(n-1)//2; C=list(range(1,m+1)); A=circulant_tournament(n,C)
    ok=is_tournament(A)
    w=np.exp(2j*pi/n)
    lam=[sum(w**((j*k)%n) for k in C) for j in range(n)]
    # closed form: lam_j = sum_{k=1}^m w^{jk} = (Dirichlet_m at 2pi j/n minus 1)/... check Re
    reals=[round(l.real,4) for l in lam]
    # Dirichlet: sum_{k=-m}^{m} e^{ik t} = sin((2m+1)t/2)/sin(t/2); our sum_{k=1}^m e^{ikt}=(D-1)/2 real part
    dirich=[]
    for j in range(n):
        t=2*pi*j/n
        if j==0: dval=m
        else:    dval=( sin((2*m+1)*t/2)/sin(t/2) - 1)/2
        dirich.append(round(dval,4))
    match=all(abs(reals[j]-dirich[j])<1e-6 for j in range(n))
    print("  interval C={1..%d} n=%d (tournament=%s): Re(lam_j) == (Dirichlet_%d(2pi j/n)-1)/2 ? %s"
          %(m,n,ok,m,match))
print("\n  => circulant atoms' eigenvalues are Gauss sums (Paley) / Dirichlet-kernel values (interval);")
print("     via LEG 1 the spectrum of ANY reducible tournament is built from these trig atoms.")
