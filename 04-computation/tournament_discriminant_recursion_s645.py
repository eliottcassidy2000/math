#!/usr/bin/env python3
"""
S645 / HYP-2323 — The tournament DISCRIMINANT (skew-adjacency det) and how tournament
properties change recursively as n grows.

TOURNAMENT DISCRIMINANT := det(M), M = skew-adjacency (M_ij=+1 if i->j, -1 if j->i).
Skew-symmetric ⟹ det = 0 (n odd) or det = Pf(M)^2 (n even, a perfect SQUARE). The
Pfaffian Pf = the Vandermonde/√-discriminant analogue (sign under odd relabel = An, S643).
Always a square ⟺ Aut(T)⊆An (S643) ⟺ the falling-factorial disc model (S644).

Computes, per n: the discriminant spectrum; the Pfaffian spectrum; and several other
tournament invariants and their RECURSION (n->n-1 parity flip; n->n-2 Pfaffian / Mode B).
No external libs (exact integer determinant + Pfaffian by recursion).
"""
from itertools import permutations, product
from fractions import Fraction

# ---------- exact integer determinant (Bareiss) ----------
def det_int(M):
    n=len(M); M=[row[:] for row in M]; sign=1; prev=1
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

# ---------- Pfaffian (recursive, exact) ----------
def pfaffian(M):
    n=len(M)
    if n==0: return 1
    if n%2==1: return 0
    # expand along first row
    idx=list(range(n)); total=0
    first=idx[0]; rest=idx[1:]
    for jpos,j in enumerate(rest):
        sgn=(-1)**jpos
        sub=[i for i in rest if i!=j]
        Msub=[[M[a][b] for b in sub] for a in sub]
        total+=sgn*M[first][j]*pfaffian(Msub)
    return total

def skew_matrix(beats,n):
    M=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j: M[i][j]=1 if beats[i][j] else -1
    return M

def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in product([0,1],repeat=len(pairs)):
        beats=[[False]*n for _ in range(n)]
        for (i,j),b in zip(pairs,bits):
            if b: beats[i][j]=True
            else: beats[j][i]=True
        yield beats

def iso_key(beats,n):
    best=None
    for p in permutations(range(n)):
        t=tuple(beats[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or t<best: best=t
    return best

def is_square(m):
    if m<0: return False
    r=int(m**0.5)
    return any((r+d)**2==m for d in (-1,0,1))

print("="*70)
print("(A) tournament discriminant det(M): 0 for odd n, a SQUARE Pf^2 for even n")
print("="*70)
print("  n | #iso classes | discriminant values (det) | all squares? | Pfaffian values (even n)")
for n in range(1,7):
    seen=set(); discs=set(); pfs=set()
    for beats in all_tournaments(n):
        k=iso_key(beats,n)
        if k in seen: continue
        seen.add(k)
        M=skew_matrix(beats,n)
        d=det_int(M); discs.add(d)
        if n%2==0: pfs.add(abs(pfaffian(M)))
    allsq=all(is_square(d) for d in discs)
    pfstr = sorted(pfs) if n%2==0 else "(0, odd n)"
    print(f"  {n} | {len(seen):11d} | {sorted(discs)!s:26s} | {str(allsq):>12} | {pfstr}")
print("  -> odd n: det=0 (=0^2). even n: det=Pf^2, a perfect square -- ALWAYS a square,")
print("     matching Aut(T)<=A_n (S643) & the falling-factorial disc model (S644).")

print("\n" + "="*70)
print("(B) RECURSION: n->n-1 flips parity (0 <-> Pf^2); Pfaffian recurses n->n-2 (Mode B)")
print("="*70)
print("  the discriminant alternates: zero (odd) / nonzero-square (even). The square root")
print("  (Pfaffian) of the even-n discriminant is built from (n-2)-Pfaffians by the cofactor")
print("  expansion Pf(M)=sum_j (-1)^j M_1j Pf(M_hat{1,j}) -- the n->n-2 'both-legs' descent")
print("  (CLAUDE.md Mode B). Mode A (n->n-1, add a vertex) is the parity flip.")
# show one explicit even-odd-even chain (transitive + extensions)
def transitive(n):
    beats=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i<j: beats[i][j]=True
    return beats
for n in [2,3,4,5,6]:
    M=skew_matrix(transitive(n),n)
    d=det_int(M); pf=pfaffian(M) if n%2==0 else 0
    print(f"  transitive T_{n}: det={d}, Pf={pf}  ({'square '+str(pf)+'^2' if n%2==0 else 'zero (odd)'})")

print("\n" + "="*70)
print("(C) other invariants and their recursion as n grows")
print("="*70)
print("  n | #iso (A000568) | max|Aut| | #Ham-path values H (odd) | rank(M) values")
def aut_order(beats,n):
    return sum(1 for p in permutations(range(n))
               if all(beats[p[i]][p[j]]==beats[i][j] for i in range(n) for j in range(n)))
def ham_paths(beats,n):
    c=0
    for p in permutations(range(n)):
        if all(beats[p[i]][p[i+1]] for i in range(n-1)): c+=1
    return c
def rank_mod(M,p=10**9+7):
    # rank over a large prime ~ rank over Q
    A=[[x%p for x in row] for row in M]; n=len(A); r=0
    for col in range(n):
        piv=None
        for i in range(r,n):
            if A[i][col]%p!=0: piv=i;break
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        inv=pow(A[r][col],p-2,p)
        A[r]=[(x*inv)%p for x in A[r]]
        for i in range(n):
            if i!=r and A[i][col]%p!=0:
                f=A[i][col]; A[i]=[(A[i][j]-f*A[r][j])%p for j in range(n)]
        r+=1
    return r
for n in range(1,7):
    seen=set(); auts=set(); hs=set(); ranks=set()
    for beats in all_tournaments(n):
        k=iso_key(beats,n)
        if k in seen: continue
        seen.add(k)
        auts.add(aut_order(beats,n)); hs.add(ham_paths(beats,n))
        ranks.add(rank_mod(skew_matrix(beats,n)))
    print(f"  {n} | {len(seen):14d} | {max(auts):8d} | {sorted(hs)!s:24s} | {sorted(ranks)}")
print("  -> A000568 = 1,1,2,4,12,56 (the iso-class count, the master recursion).")
print("  -> H (Rédei) always ODD; rank(M) = n (even, det!=0) or n-1 (odd, det=0): the")
print("     skew rank drops by 1 exactly when n is odd = when the discriminant vanishes.")
print("     So 'discriminant zero' <=> 'rank deficient by 1' <=> n odd -- the parity recursion.")
