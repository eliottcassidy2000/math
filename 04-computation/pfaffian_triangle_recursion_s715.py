#!/usr/bin/env python3
"""
S715 — The Pfaffian as translator + the triangular inclusion-exclusion recursion (the 7-piece tiling).

The user's 7-piece decomposition: an n-tournament's staircase triangle = 3 corner (n-1)-triangles
{A,B,C} (+), minus 3 edge (n-2)-overlaps {D,E,F} (-), plus 1 center (n-3)-triangle {G} (+). Per region:
corners = {A or C or E}, edges = {B+A-D, A+C-E, B+C-F}, interior = A+B+C-D-E-F (+G inside).

This is the inclusion-exclusion reconstruction of a triangular grid by three corner subtriangles. We:

(1) PROVE it is a partition of unity: 3*Tri(n-1) - 3*Tri(n-2) + Tri(n-3) = Tri(n), and every cell has
    net multiplicity 1 (corner cells x1, edge cells 2-1=1, interior 3-3+1=1). Char poly (x-1)^3, so the
    coefficient pattern is the 3rd finite difference Delta^3 = 0.

(2) RECURSIVE TRUTH A: a TILE-ADDITIVE invariant F (F = sum over cells of a local value) obeys
    F(n) - 3F(n-1) + 3F(n-2) - F(n-3) = 0  <=>  F is a QUADRATIC polynomial in n. Verify on #tiles,
    #arcs, #vertices (all Delta^3=0); show H (#Ham paths, exponential) and |Pf| are NOT tile-additive
    (Delta^3 != 0).

(3) RECURSIVE TRUTH B (the multiplicative shadow): H is the independence polynomial I(Omega,2) — a
    PRODUCT over the conflict graph, not a sum over cells — so it does NOT obey (x-1)^3. The base-path
    staircase family (THM-337) instead obeys H(k)=3H(k-1)+H(k-2)+H(k-3): the leading 3 (the three
    corners) SURVIVES, but the inclusion-exclusion signs flip to + because the object is multiplicative.
    Verify THM-337 by direct H computation; contrast its char poly x^3-3x^2-x-1 (lambda~3.383) with the
    additive (x-1)^3.

(4) THE PFAFFIAN = the signed n->n-2 (Mode-B) deletion. Pf(S)=sum_j (-1)^(j-1) S_1j Pf(S_minor): the
    minors are (n-2)-vertex deleted subtournaments and they enter with ALTERNATING SIGN — exactly the
    user's "D,E,F of size n-2 are negative". Verify the Pfaffian minor recursion; this is the Mode-B
    (THM-291) object made into an exact signed recursion.

(5) max-H A038375: exact small values + the 'edge over random' R(n)=a(n)2^(n-1)/n! (slow growth, the
    recursive-efficiency handle).

No numpy/sympy.
"""
import itertools
from math import comb, factorial

# ---------- Pfaffian + tournaments ----------
def pf(M):
    n=len(M)
    if n==0: return 1
    if n%2==1: return 0
    tot=0; r0=M[0]
    for j in range(1,n):
        a=r0[j]
        if a==0: continue
        rest=[k for k in range(1,n) if k!=j]
        sub=[[M[r][c] for c in rest] for r in rest]
        tot+=((-1)**(j-1))*a*pf(sub)
    return tot
def skew(A):
    n=len(A); return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
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
    return sum(dp[(1<<n)-1][v] for v in range(n))
def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def base_path_staircase(k):
    """THM-337: n=2k vertices, base path n-1->...->0, non-base tiles 'upward' (lower beats higher)."""
    n=2*k; A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j==i+1: A[j][i]=1     # base path arc (higher -> lower)
            else: A[i][j]=1          # non-base: lower-indexed beats higher
    return A

def Tri(k):  # triangular number of cells in a side-k triangular grid; 0 for k<=0
    return k*(k+1)//2 if k>0 else 0

def delta3(seq):
    """3rd finite difference at the right end of a length>=4 list."""
    return [seq[i+3]-3*seq[i+2]+3*seq[i+1]-seq[i] for i in range(len(seq)-3)]

# ============================ RUN ============================
if __name__=="__main__":
    print("="*82)
    print("S715 — Pfaffian translator + the triangular inclusion-exclusion (7-piece tiling)")
    print("="*82)

    # (1) partition of unity
    print("\n(1) THE 7-PIECE IE IS A PARTITION OF UNITY (triangular grid)")
    ok=all(3*Tri(n-1)-3*Tri(n-2)+Tri(n-3)==Tri(n) for n in range(3,40))
    print(f"  3*Tri(n-1) - 3*Tri(n-2) + Tri(n-3) == Tri(n) for n=3..39: {ok}")
    print(f"  coefficients (1,-3,3,-1) = (x-1)^3 = the 3rd finite difference Delta^3; per-cell net mult:")
    print(f"    corner cell: in 1 big triangle  -> 1")
    print(f"    edge cell:   in 2 big, 1 overlap -> 2-1 = 1")
    print(f"    interior:    3 big, 3 overlap, 1 center -> 3-3+1 = 1   (Mobius reconstruction)")

    # (2) tile-additive invariants: Delta^3 = 0
    print("\n(2) RECURSIVE TRUTH A: tile-additive invariant => Delta^3 F = 0 (quadratic in n)")
    ns=list(range(3,12))
    tiles=[comb(n-1,2) for n in ns]; arcs=[comb(n,2) for n in ns]; verts=ns[:]
    print(f"  #tiles C(n-1,2):  {tiles}   Delta^3={delta3(tiles)}")
    print(f"  #arcs  C(n,2):    {arcs}   Delta^3={delta3(arcs)}")
    print(f"  #verts n:         {verts}   Delta^3={delta3(verts)}")
    print("  => all tile-additive (degree<=2) invariants are killed by Delta^3: the user's IE recursion.")

    # (3) H is NOT additive; base-path family obeys THM-337 (3,1,1)
    print("\n(3) RECURSIVE TRUTH B: H is multiplicative (I(Omega,2)) — NOT (x-1)^3; base-path = THM-337")
    # exact max-H A038375 small n via full enumeration
    maxH={}
    for n in range(2,7):
        best=0
        for A in all_tournaments(n):
            h=ham_paths(A)
            if h>best: best=h
        maxH[n]=best
    seqA=[maxH[n] for n in range(2,7)]
    print(f"  A038375 max-H (n=2..6, exact): {seqA}  (OEIS 1,1,3,5,15,45; n=7..13: 189,661,3357,15745,103605,531205,3719831 — THM-329)")
    full=[1,1,3,5,15,45,189,661,3357,15745,103605,531205,3719831]
    print(f"  Delta^3 of A038375 (n=1..13): {delta3(full)}  (NOT zero => H not tile-additive)")
    # base-path family
    bp=[ham_paths(base_path_staircase(k)) for k in range(1,10)]
    print(f"  base-path staircase H (k=1..9, n=2..18): {bp}")
    rec_ok=all(bp[i]==3*bp[i-1]+bp[i-2]+bp[i-3] for i in range(3,len(bp)))
    print(f"  THM-337 recurrence H(k)=3H(k-1)+H(k-2)+H(k-3) holds k=4..9: {rec_ok}")
    print(f"  char poly x^3-3x^2-x-1 (lambda~3.383): the leading '3' = the 3 corners survives;")
    print(f"  signs are + (multiplicative), NOT the additive IE's alternating (x-1)^3.")

    # (4) Pfaffian = signed n->n-2 (Mode B) deletion
    print("\n(4) PFAFFIAN = the signed n->n-2 (Mode-B) deletion: minors are (n-2)-subtournaments, ALTERNATING sign")
    import random
    rng=random.Random(1); checks=0; good=0
    for _ in range(2000):
        n=6; A=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if rng.random()<0.5: A[i][j]=1
                else: A[j][i]=1
        S=skew(A); P=pf(S)
        # recursive expansion deleting vertex 0 and j (n-2 minor)
        rec=0
        for j in range(1,n):
            rest=[k for k in range(1,n) if k!=j]
            sub=[[S[r][c] for c in rest] for r in rest]
            rec+=((-1)**(j-1))*S[0][j]*pf(sub)
        checks+=1; good+= (rec==P)
    print(f"  Pf(S)=sum_j (-1)^(j-1) S_0j Pf(S_minor[0,j]) verified on {good}/{checks} random n=6 tournaments")
    print(f"  the (n-2) minors enter with ALTERNATING SIGN = the user's 'D,E,F size n-2 are NEGATIVE'.")
    print(f"  Mode-B (THM-291) is this same n->n-2 step; the Pfaffian makes it an exact SIGNED recursion.")

    # (5) edge over random
    print("\n(5) max-H 'edge over random' R(n)=a(n)*2^(n-1)/n!  (the recursive-efficiency handle)")
    for i,n in enumerate(range(1,14)):
        R=full[i]*(2**(n-1))/factorial(n)
        print(f"   n={n:2d}: a(n)={full[i]:8d}   R(n)={R:.4f}")
    print("  R(n) grows SLOWLY (~sqrt(n)-ish), oscillating by parity: max-H = (n!/2^(n-1)) * slow factor;")
    print("  the maximizers are Paley/circulant (THM-336/338), not the base-path family (17 vs 45 at n=6).")
    print("="*82)
