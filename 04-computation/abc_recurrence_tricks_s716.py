#!/usr/bin/env python3
"""
S716 — Understanding the A+B+C-D-E-F+G recurrence; relations and computational tricks.

The 7 pieces are the 7 = 2^3 - 1 nonempty subsets of a 3-set: {A,B,C}=singletons (+), {D,E,F}=pairs (-),
{G}=triple (+). Signs = (-1)^(|subset|+1) = the Mobius function of the boolean lattice B_3. So
A+B+C-D-E-F+G is INCLUSION-EXCLUSION for the union of the three corner (n-1)-subtriangles, and the
coefficient vector (1,-3,3,-1) is (x-1)^3 = the 3rd finite difference Delta^3.

This session:
(1) the d-SIMPLEX family: corner-IE over a d-simplex gives (x-1)^(d+1) = Delta^(d+1); cell counts are
    degree-d polynomials. Verify d=1,2,3,4.
(2) RELATION + TRICK A (additive/polynomial): a tile-additive (valuative) invariant satisfies the EXACT
    IE A+B+C-D-E-F+G, hence Delta^3=0, hence is QUADRATIC in n: determined by 3 seeds, O(1)/term, closed
    form. Verify on #tiles, #arcs, and the EXACT valuative IE for arc-counts of real sub-tournament unions.
(3) RELATION + TRICK B (multiplicative/exponential): H of a recursive tournament FAMILY is C-finite
    (satisfies a linear recurrence because the boundary width is bounded) => compute H(n) by
    COMPANION-MATRIX EXPONENTIATION in O(log n), vs O(2^n) for direct DP. Auto-discover the recurrence
    from a few terms (minimal-recurrence finder), then matrix-power to huge n. Verify base-path (THM-337).
(4) the boundary: A038375 (max-H) has NO low-order linear recurrence (transcendental growth) — the trick
    works for fixed FAMILIES, not the optimum-over-all-tournaments.

No numpy/sympy. Exact rationals via Fractions; exact integer matrix power.
"""
from fractions import Fraction as Fr
from math import comb
from itertools import combinations

# ---------------- minimal linear recurrence finder ----------------
def find_recurrence(seq, maxord=6):
    """smallest order r and integer-ish coeffs c with seq[i]=sum_j c[j] seq[i-j], verified on all data."""
    n=len(seq)
    for r in range(1,maxord+1):
        if n < 2*r+1: break
        # build system from rows i=r..n-1 ; need >= r equations, use first r for solve, rest to verify
        rows=[[Fr(seq[i-1-j]) for j in range(r)] for i in range(r, 2*r)]
        rhs=[Fr(seq[i]) for i in range(r,2*r)]
        c=gauss_solve(rows,rhs)
        if c is None: continue
        ok=all(sum(c[j]*seq[i-1-j] for j in range(r))==seq[i] for i in range(r,n))
        if ok: return r,c
    return None,None

def gauss_solve(A,b):
    A=[row[:] for row in A]; b=b[:]; n=len(A)
    for k in range(n):
        piv=next((i for i in range(k,n) if A[i][k]!=0),None)
        if piv is None: return None
        A[k],A[piv]=A[piv],A[k]; b[k],b[piv]=b[piv],b[k]
        inv=A[k][k]
        A[k]=[x/inv for x in A[k]]; b[k]=b[k]/inv
        for i in range(n):
            if i!=k and A[i][k]!=0:
                f=A[i][k]; A[i]=[A[i][j]-f*A[k][j] for j in range(n)]; b[i]=b[i]-f*b[k]
    return b

# ---------------- integer matrix power (companion) ----------------
def matmul(A,B):
    n=len(A);m=len(B[0]);k=len(B)
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def matpow(M,e):
    n=len(M); R=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    while e:
        if e&1: R=matmul(R,M)
        M=matmul(M,M); e>>=1
    return R
def companion(coeffs):
    """coeffs c[0..r-1] for s[i]=c0 s[i-1]+...+c_{r-1} s[i-r]; integer companion matrix."""
    r=len(coeffs); M=[[0]*r for _ in range(r)]
    M[0]=[int(c) for c in coeffs]
    for i in range(1,r): M[i][i-1]=1
    return M
def crec(coeffs, seeds, N):
    """N-th term (0-indexed over seeds) via matrix power; seeds = first r terms."""
    r=len(coeffs)
    if N<r: return seeds[N]
    M=companion(coeffs); P=matpow(M,N-(r-1))
    state=[seeds[r-1-i] for i in range(r)]   # [s_{r-1},...,s_0]
    return sum(P[0][j]*state[j] for j in range(r))

# ---------------- tournaments ----------------
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
def base_path_staircase(k):
    n=2*k; A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j==i+1: A[j][i]=1
            else: A[i][j]=1
    return A
def arcs_in(A, S):
    S=list(S); c=0
    for a in range(len(S)):
        for b in range(len(S)):
            if a!=b and A[S[a]][S[b]]: c+=1
    return c

# ============================ RUN ============================
if __name__=="__main__":
    print("="*82)
    print("S716 — the A+B+C-D-E-F+G recurrence: relations & computational tricks")
    print("="*82)

    print("\n(1) A+B+C-D-E-F+G = inclusion-exclusion over B_3; d-simplex => (x-1)^(d+1)")
    print("   pieces = 2^3-1=7 subsets; signs (+,+,+,-,-,-,+) = (-1)^(|S|+1) = Mobius of B_3.")
    for d in (1,2,3,4):
        # d-simplex side-n cell count = C(n+d-1, d) (degree-d polynomial); IE coeffs = C(d+1,k)(-1)^(k+1)
        N=lambda n: comb(n+d-1,d) if n>=0 else 0
        ie_ok=all(sum((-1)**(k+1)*comb(d+1,k)*N(n-k) for k in range(1,d+2))==N(n) for n in range(d+2,30))
        # delta^(d+1) = 0
        seq=[N(n) for n in range(0,12)]
        delta=seq[:]
        for _ in range(d+1):
            delta=[delta[i+1]-delta[i] for i in range(len(delta)-1)]
        print(f"   d={d} ({d}-simplex): corner-IE sum_{{k}}(-1)^(k+1)C({d+1},k)N(n-k)=N(n)? {ie_ok}; "
              f"Delta^{d+1} N == 0? {all(x==0 for x in delta)}  (N(n)=C(n+{d-1},{d}), degree {d})")
    print("   => the user's triangle is the d=2 case: 3 corners, char poly (x-1)^3, quadratic cells.")

    print("\n(2) RELATION+TRICK A (additive/valuative): exact IE => Delta^3=0 => quadratic; 3 seeds suffice")
    tiles=[comb(n-1,2) for n in range(3,14)]
    r,c=find_recurrence(tiles)
    print(f"   #tiles C(n-1,2): minimal recurrence order {r}, coeffs {[str(x) for x in c]}  (expect 3; 3,-3,1)")
    # closed form via 3 seeds + recurrence, O(1)/term
    print(f"   3-seed reconstruction of #tiles from [1,3,6]: "
          f"{[crec([3,-3,1],[1,3,6],N) for N in range(8)]}  (= {tiles[:8]})")
    # EXACT valuative IE on arc-counts of a real tournament's 3 corner sub-tournaments
    import random; rng=random.Random(3); n=8
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<0.5: A[i][j]=1
            else: A[j][i]=1
    # three corners = delete vertex 0 / 1 / 2 ; union = all ; overlaps = delete pairs ; center = delete all 3
    V=set(range(n))
    cA,cB,cC=V-{0},V-{1},V-{2}
    dAB,dAC,dBC=V-{0,1},V-{0,2},V-{1,2}
    g=V-{0,1,2}
    ie = arcs_in(A,cA)+arcs_in(A,cB)+arcs_in(A,cC)-arcs_in(A,dAB)-arcs_in(A,dAC)-arcs_in(A,dBC)+arcs_in(A,g)
    print(f"   EXACT valuative IE on arcs: A+B+C-D-E-F+G = {ie}, full arc-count = {arcs_in(A,V)}  "
          f"(match: {ie==arcs_in(A,V)})  <- additive invariants obey the IE exactly")

    print("\n(3) RELATION+TRICK B (multiplicative): H of a recursive family is C-finite => O(log n) matrix power")
    bp=[ham_paths(base_path_staircase(k)) for k in range(1,10)]   # direct DP, k=1..9
    print(f"   base-path H (direct DP, k=1..9): {bp}")
    r,c=find_recurrence(bp)
    print(f"   AUTO-DISCOVERED recurrence: order {r}, coeffs {[str(x) for x in c]}  (THM-337: 3,1,1)")
    ci=[int(x) for x in c]
    # matrix-power evaluation to huge k (direct DP would be O(2^(2k)) - infeasible past k~13)
    big=[crec(ci, bp[:r], N) for N in (20, 50, 100)]
    print(f"   matrix-power H at k=20,50,100 (O(log k), DP infeasible): {big}")
    # cross-check k=9 by both
    print(f"   cross-check k=9: matrix-power={crec(ci,bp[:r],8)} vs DP={bp[8]}  (match: {crec(ci,bp[:r],8)==bp[8]})")

    print("\n(4) THE BOUNDARY: A038375 (max-H over ALL tournaments) is NOT C-finite (no low-order recurrence)")
    a=[1,1,3,5,15,45,189,661,3357,15745,103605,531205,3719831]
    r,c=find_recurrence(a,maxord=5)
    print(f"   A038375 minimal recurrence (order<=5): {'NONE FOUND' if r is None else (r,c)}  "
          f"=> the trick needs a fixed FAMILY (bounded boundary width), not the optimum.")
    print("="*82)
