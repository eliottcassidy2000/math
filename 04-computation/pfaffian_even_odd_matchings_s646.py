#!/usr/bin/env python3
"""
S646 / HYP-2324 — Even/odd and the Pfaffian: WHY the tournament Pfaffian is odd.

Develops S645 (tournament discriminant det M = Pf^2 even n / 0 odd n).

THE PFAFFIAN as a signed sum over PERFECT MATCHINGS:
  Pf(M) = sum_{matchings pi of {1..2n}} sgn(pi) * prod_{(i,j) in pi} M_ij
  #matchings = (2n-1)!! = (2n-1)(2n-3)...3.1  (double factorial = product of ODDS).

KEY CHAIN (explains S645's 'Pfaffian always odd'):
  (2n-1)!! is ODD (product of odds)
   => Pf of a +-1 matrix = sum of (2n-1)!! many +-1 terms = sum of ODD-many odd terms = ODD
   => det = Pf^2 = ODD SQUARE.
And (2n-1)!! = 2^n * (1/2)^(n) [rising factorial of 1/2, S644] = the apex/sqrt(pi) constant.

PARITY MASTER SWITCH: even n <=> perfect matching exists <=> Pfaffian defined <=> det = square;
                      odd n  <=> NO matching <=> Pf = 0 <=> det = 0 <=> rank deficient.
Plus: a skew-symmetric matrix has EVEN RANK (always) -- the structural even/odd.
No external libs.
"""
from itertools import permutations, product
from math import prod

# ---------- perfect matchings & Pfaffian ----------
def perfect_matchings(elts):
    if not elts:
        yield []
        return
    a = elts[0]
    for i in range(1, len(elts)):
        b = elts[i]
        rest = elts[1:i] + elts[i+1:]
        for m in perfect_matchings(rest):
            yield [(a, b)] + m

def matching_sign(matching, n):
    # sign of the permutation (a1 b1 a2 b2 ...) reading the matching as a permutation of 0..2n-1
    perm = []
    for (a, b) in matching:
        perm += [a, b]
    # sign via inversions
    s = 1
    for i in range(len(perm)):
        for j in range(i+1, len(perm)):
            if perm[i] > perm[j]: s = -s
    return s

def pfaffian_via_matchings(M):
    n = len(M)
    if n % 2 == 1: return 0
    total = 0
    for mt in perfect_matchings(list(range(n))):
        total += matching_sign(mt, n) * prod(M[a][b] for (a, b) in mt)
    return total

def det_int(M):
    import copy
    M = [row[:] for row in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            piv=next((i for i in range(k+1,n) if M[i][k]!=0), None)
            if piv is None: return 0
            M[k],M[piv]=M[piv],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]

def double_factorial(m):
    r=1
    while m>1: r*=m; m-=2
    return r

print("="*68)
print("(A) #perfect matchings of 2n points = (2n-1)!! = product of ODDS (so ODD)")
print("="*68)
print("  2n | #matchings | (2n-1)!! | odd? | 2^n*(1/2)^(n) [rising 1/2, S644]")
for n in range(1,7):
    twoN=2*n
    cnt=sum(1 for _ in perfect_matchings(list(range(twoN))))
    df=double_factorial(twoN-1)
    # rising factorial of 1/2: (1/2)(3/2)...((2n-1)/2) = (2n-1)!!/2^n
    from fractions import Fraction
    rf=Fraction(1)
    for i in range(n): rf*= Fraction(2*i+1,2)
    check = (Fraction(2)**n * rf == df)
    print(f"  {twoN:2d} | {cnt:10d} | {df:8d} | {str(df%2==1):>4} | 2^{n}*(1/2)^({n})={float(2**n*rf):.4g} -> ={df}? {check}")
print("  -> the Pfaffian term-count (2n-1)!! is ODD; = 2^n * rising-factorial-of-1/2 (S644).")

print("\n" + "="*68)
print("(B) WHY the tournament Pfaffian is ODD: odd-many +-1 terms")
print("="*68)
# random-ish +-1 skew matrices (tournaments) and their Pfaffian parity
def skew_from_bits(bits, n):
    M=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            v=1 if bits[idx] else -1; idx+=1
            M[i][j]=v; M[j][i]=-v
    return M
for n in [2,4,6]:
    npairs=n*(n-1)//2
    allodd=True; pfset=set(); detok=True
    import itertools
    sample = list(itertools.product([0,1],repeat=npairs)) if npairs<=15 else None
    if sample is None:
        sample=[tuple(__import__('random').randint(0,1) for _ in range(npairs)) for _ in range(300)]
    for bits in sample:
        M=skew_from_bits(bits,n)
        pf=pfaffian_via_matchings(M); pfset.add(pf)
        if pf%2==0: allodd=False
        if det_int(M)!=pf*pf: detok=False
    print(f"  n={n}: Pfaffian values {sorted(set(abs(x) for x in pfset))}; all ODD? {allodd}; "
          f"det==Pf^2 ? {detok}")
print("  -> every +-1 skew Pfaffian is ODD (odd # of +-1 terms); det = Pf^2 exactly.")

print("\n" + "="*68)
print("(C) the PARITY MASTER SWITCH and even rank of skew matrices")
print("="*68)
print("  n | matching exists? | Pfaffian | det | rank | rank even?")
for n in range(1,8):
    M=skew_from_bits([1]*(n*(n-1)//2), n)   # the transitive tournament
    pf=pfaffian_via_matchings(M); d=det_int(M)
    # rank over rationals via the determinant pattern: skew rank is even; here n or n-1
    # compute rank by gaussian elim over Q
    from fractions import Fraction
    A=[[Fraction(x) for x in row] for row in M]; r=0; m=len(A)
    for col in range(m):
        piv=next((i for i in range(r,m) if A[i][col]!=0),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        for i in range(m):
            if i!=r and A[i][col]!=0:
                f=A[i][col]/A[r][col]; A[i]=[A[i][j]-f*A[r][j] for j in range(m)]
        r+=1
    me = "yes" if n%2==0 else "no"
    print(f"  {n} | {me:>15} | {pf:8d} | {d:4d} | {r:4d} | {r%2==0}")
print("  -> even n: matching exists, Pf!=0, det=Pf^2 square, full rank n (even).")
print("     odd n: NO perfect matching, Pf=0, det=0, rank n-1 (even). Skew rank ALWAYS even.")
print("     The even/odd of n is the master switch; the Pfaffian lives only on the even side,")
print("     and the n->n-2 Pfaffian recursion stays within one parity class (Mode B, S645).")
