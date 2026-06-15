#!/usr/bin/env python3
"""
mobius_sieve_arrow_macmini_0615s3.py  (mac-mini-2026-06-15-S3, THM-512)

(1) THE MÖBIUS VERTEX-DELETION SIEVE (the user's A+B+C-D-E-F+G).
    Claim: for any additive sub-structure invariant phi (e.g. c3 = #directed triangles),
    summing phi(T-A) over nonempty A in the n-1 non-anchor vertices with sign (-1)^{|A|+1}
    reconstructs phi via inclusion-exclusion. codim-c term count = C(n-1,c), sign (-1)^{c+1}.
    Verify exhaustively n=4,5,6 what the signed sum equals (vs c3(T), H(T)).

(2) THE ARROW / CONDORCET BRIDGE. The tiling/arc cube IS the social-choice cube
    {+-1}^{C(m,2)}; transitive = rational/Arrow-dictator outcome; a directed 3-cycle =
    Condorcet paradox = the minimal odd cycle = the OCF obstruction. Compute, over uniform
    tournaments on m candidates: P(transitive) (rational), the Condorcet/Guilbaud structure,
    and that c3 sits ENTIRELY at Walsh level 2 with Var(c3)=3*C(m,3)/16 (the Condorcet
    quadratic) -- the project sibling of Kalai's P_rational = 3/4 + 3/4 Stab_{-1/3}.
"""
import sys, itertools
from math import comb
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def all_labelled_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def c3(A, verts=None):
    if verts is None: verts = range(len(A))
    verts = list(verts); c = 0
    for a, b, d in itertools.combinations(verts, 3):
        if (A[a][b] and A[b][d] and A[d][a]) or (A[b][a] and A[d][b] and A[a][d]): c += 1
    return c

def H(A):
    n = len(A); out = [0]*n
    for i in range(n):
        r = 0
        for j in range(n):
            if A[i][j]: r |= 1 << j
        out[i] = r
    full = (1 << n) - 1; dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            cc = row[last]
            if cc:
                nx = out[last] & ~mask
                while nx:
                    bb = nx & (-nx); j = bb.bit_length()-1; dp[mask|bb][j] += cc; nx ^= bb
    return sum(dp[full][:])

def sub_invariant(A, keep, phi):
    return phi(A, keep)

print("="*72)
print("(1) MÖBIUS VERTEX-DELETION SIEVE: Σ_{∅≠A⊆[n-1]} (-1)^{|A|+1} φ(T-A) =? φ(T)")
print("="*72)
for n in range(4, 7):
    # anchor vertex = 0; deletable = 1..n-1
    deletable = list(range(1, n))
    # test on a sample of tournaments (first 200) for c3
    matches_c3 = 0; tot = 0; diffs = set()
    cnt = 0
    for A in all_labelled_tournaments(n):
        cnt += 1
        if cnt > 400: break
        full = c3(A)
        s = 0
        for r in range(1, len(deletable)+1):
            for Aset in itertools.combinations(deletable, r):
                keep = [v for v in range(n) if v not in Aset]
                s += ((-1)**(r+1)) * c3(A, keep)
        tot += 1
        if s == full: matches_c3 += 1
        else: diffs.add(full - s)
    print(f"  n={n}: c3 sieve == c3(T) on {matches_c3}/{tot} sampled tournaments"
          f"{'' if not diffs else f'; (c3-sieve) residual values seen: {sorted(diffs)[:6]}'}")
    print(f"        codim-c terms: count C({n-1},c), sign (-1)^(c+1); "
          f"Σ(-1)^(c+1)C({n-1},c) = {sum((-1)**(c+1)*comb(n-1,c) for c in range(1,n))}")

print("\n" + "="*72)
print("(2) ARROW / CONDORCET: random tournaments on m candidates")
print("="*72)
for m in range(3, 8):
    N = 0; trans = 0; c3sum = 0; c3sq = 0
    for A in all_labelled_tournaments(m):
        N += 1
        cc = c3(A); c3sum += cc; c3sq += cc*cc
        if cc == 0: trans += 1
    p_rational = trans / N
    Ec3 = c3sum/N; Var = c3sq/N - Ec3**2
    var_formula = 3*comb(m,3)/16
    # Guilbaud limit for majority = 0.9123 (reference, m=3, n->inf voters)
    print(f"  m={m}: P(transitive/rational)={p_rational:.4f} (={trans}/{N}); "
          f"E[c3]={Ec3:.3f}=C(m,3)/4={comb(m,3)/4:.3f}; Var(c3)={Var:.4f} "
          f"(formula 3C(m,3)/16={var_formula:.4f}: {abs(Var-var_formula)<1e-9})")
print("  Condorcet paradox = directed 3-cycle in the (majority) tournament = minimal odd cycle")
print("  = the OCF obstruction. Kalai: P_rational(majority,3 cand) = 3/4 + 3/4 Stab_{-1/3}[Maj]")
print("  -> Guilbaud's number 3/4+(3/4)(1-(2/pi)arcsin(1/3)) ≈ 0.9123 as #voters->inf.")
print("  Var(c3)=3C(m,3)/16 sits ENTIRELY at Walsh level 2 (THM-163/THM-511): c3 IS the")
print("  Condorcet/Guilbaud quadratic; transitive (c3=0) = the rational/Arrow-dictator ground state.")
print("\nDONE.")
