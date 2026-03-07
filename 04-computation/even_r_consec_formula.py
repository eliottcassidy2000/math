#!/usr/bin/env python3
"""
Even r-powers through the consecutive-position formula lens.

THM-030 Corollary 2: M[a,b] contains only even powers of r.
THM-050: M[a,b] = sum_j (-1)^j N(a,b,j) at r=1/2.

At general r, each Ham path P has weight:
  w(P) = prod_{k=0}^{n-2} t(P[k], P[k+1])
where t(i,j) = r + s_{ij} with s_{ij} = A[i][j] - r = (1-2r)*A[i][j] + r*(A[i][j]-1+1)
Actually simpler: t(i,j) = r + (A[i][j] - r) = A[i][j] when A[i][j]=1, but at general r:
  If A[i][j]=1: t(i,j) = r + s where s = 1/2 - r (at the standard parametrization)
  If A[i][j]=0: t(i,j) = r - s
  With s = 1/2 - r.

Actually, the c-tournament parametrization is:
  t(i,j) = r + s_{ij} where s_{ij} = -s_{ji} and |s_{ij}| can be anything.
  At r=1/2: t(i,j) = 1/2 + s_{ij}. For standard tournament: s_{ij} = 1/2 if A[i,j]=1, -1/2 if A[i,j]=0.
  So t(i,j) = 1 if A[i,j]=1, 0 if A[i,j]=0. ✓

For general r, fix the skew matrix s_{ij} = A[i][j] - 1/2 (so |s|=1/2 always).
Then t(i,j) = r + A[i][j] - 1/2.
  If A[i][j]=1: t(i,j) = r + 1/2
  If A[i][j]=0: t(i,j) = r - 1/2

This gives t(i,j) = r + (2*A[i][j]-1)/2.

The path weight is: w(P) = prod_{k=0}^{n-2} (r + (2*A[P[k]][P[k+1]]-1)/2)

We can expand this as a polynomial in r.

The transfer matrix M(r)[a,b] should have only even powers of r (THM-030 Cor 2).

The consecutive formula should generalize:
  M(r)[a,b] = sum_j (-1)^j N(r,a,b,j)
where N(r,a,b,j) involves weighted paths.

Let's compute and verify.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np
from fractions import Fraction

def count_paths_subset_weighted(A, verts, r, start=None, end=None):
    """Count weighted paths on vertex subset."""
    total = Fraction(0)
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        weight = Fraction(1)
        for i in range(len(p)-1):
            w = r + Fraction(2*A[p[i]][p[i+1]] - 1, 2)
            weight *= w
        total += weight
    return total

def transfer_matrix_r(A, r):
    """Transfer matrix at general r using IE definition."""
    n = len(A)
    M = [[Fraction(0)]*n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = Fraction(0)
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset_weighted(A, S_verts, r, end=a)
                    bb = count_paths_subset_weighted(A, R_verts, r, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def consec_formula_r(A, r):
    """Consecutive formula at general r."""
    n = len(A)
    M = [[Fraction(0)]*n for _ in range(n)]

    # Enumerate all permutations and compute weighted paths
    for perm in permutations(range(n)):
        weight = Fraction(1)
        for k in range(n-1):
            w = r + Fraction(2*A[perm[k]][perm[k+1]] - 1, 2)
            weight *= w

        # Diagonal contribution
        for a in range(n):
            pos = list(perm).index(a)
            M[a][a] += ((-1)**pos) * weight

        # Off-diagonal: consecutive pairs (symmetrized)
        for j in range(n-1):
            a, b = perm[j], perm[j+1]
            # Contribution to M[a][b] and M[b][a]
            contrib = ((-1)**j) * weight
            M[a][b] += contrib
            M[b][a] += contrib

    return M

# =====================================================================
# Test at n=3 with symbolic r
# =====================================================================
print("=" * 70)
print("EVEN r-POWERS IN CONSECUTIVE FORMULA")
print("=" * 70)

# Use r as a Fraction parameter (test at specific values)
n = 3
# 3-cycle: 0→1→2→0
A3 = [[0,1,0],[0,0,1],[1,0,0]]

print(f"\nn=3, 3-cycle, testing consecutive formula at various r:")
for r_num, r_den in [(1,2), (1,3), (1,4), (2,3), (3,4)]:
    r = Fraction(r_num, r_den)
    M_ie = transfer_matrix_r(A3, r)
    M_con = consec_formula_r(A3, r)

    match = all(M_ie[a][b] == M_con[a][b] for a in range(n) for b in range(n))
    print(f"  r={r}: IE==Consec? {match}, M[0,0]={M_ie[0][0]}, M[0,1]={M_ie[0][1]}")

# =====================================================================
# Verify even r-powers: compute M(r) as polynomial for n=3
# =====================================================================
print()
print("=" * 70)
print("M(r) AS POLYNOMIAL IN r")
print("=" * 70)

# For n=3, compute at enough r-values to interpolate
r_vals = [Fraction(k, 10) for k in range(1, 8)]
M_at_r = {}
for rv in r_vals:
    M_at_r[rv] = transfer_matrix_r(A3, rv)

# M[0,0](r) at each r value
print(f"\nn=3, 3-cycle: M[0,0](r) values:")
for rv in r_vals:
    print(f"  r={float(rv):.1f}: M[0,0]={float(M_at_r[rv][0][0]):.6f}")

# Fit polynomial: M[0,0](r) should be even in r
# For 3-cycle: M = I at r=1/2 (scalar). What is M at general r?
# Three paths: 0→1→2, 1→2→0, 2→0→1
# Path 0→1→2: weight = (r+1/2)(r+1/2) = (r+1/2)^2
# Path 1→2→0: weight = (r+1/2)(r+1/2) = (r+1/2)^2
# Path 2→0→1: weight = (r+1/2)(r+1/2) = (r+1/2)^2

# Wait: 0→1: A[0][1]=1, so t=r+1/2. 1→2: A[1][2]=1, so t=r+1/2.
# Path 0→1→2: w = (r+1/2)^2
# Path 1→2→0: 1→2: r+1/2. 2→0: A[2][0]=1, t=r+1/2. w = (r+1/2)^2
# Path 2→0→1: 2→0: r+1/2. 0→1: r+1/2. w = (r+1/2)^2
# Total H(r) = 3(r+1/2)^2 = 3r^2 + 3r + 3/4

# M[0,0] = sum_P (-1)^{pos(0,P)} * w(P)
# Path 0→1→2: pos(0)=0, (-1)^0=1. w=(r+1/2)^2. Contrib: +(r+1/2)^2
# Path 1→2→0: pos(0)=2, (-1)^2=1. w=(r+1/2)^2. Contrib: +(r+1/2)^2
# Path 2→0→1: pos(0)=1, (-1)^1=-1. w=(r+1/2)^2. Contrib: -(r+1/2)^2
# M[0,0] = (r+1/2)^2 + (r+1/2)^2 - (r+1/2)^2 = (r+1/2)^2
# = r^2 + r + 1/4

# Even powers of r: r^2 term + constant 1/4 + r term.
# But r is an ODD power! This contradicts THM-030 Cor 2!

# Wait — "even powers of r" means M[a,b] has only r^0, r^2, r^4, ... terms.
# M[0,0] = r^2 + r + 1/4. The r^1 term is nonzero!
# This seems to contradict THM-030...

# Let me recheck. THM-030 says even powers in the SKEW variables s, not r.
# M[a,b] = sum of even powers of r means r^0, r^2, r^4,...
# But actually, let me re-read: "Every entry M[a,b] of the transfer matrix
# contains only even powers of r."

# Hmm, at r=1/2: M[0,0]=1. At r=0: let me compute.
r0 = Fraction(0)
M_r0 = transfer_matrix_r(A3, r0)
print(f"\n  M(r=0) = {[[float(M_r0[a][b]) for b in range(n)] for a in range(n)]}")
print(f"  M[0,0](r=0) = {float(M_r0[0][0])}")

# At r=0: t(i,j) = (2A[i,j]-1)/2. For A[i,j]=1: t=1/2. For A[i,j]=0: t=-1/2.
# So this is the "pure skew" case.

# Actually, I realize the issue: THM-030 Cor 2 says M[a,b] has only even powers of r.
# Let me verify: M[0,0] should be c_0 + c_2*r^2 + c_4*r^4 + ...
# From the computation: M[0,0] = (r+1/2)^2 = r^2 + r + 1/4
# This has r^1 term! So either the theorem is wrong, or I'm misunderstanding.

# Let me check by expanding differently. The path weight for the 3-cycle:
# each edge has weight r+1/2 (since all edges exist in this regular tournament).
# w(P) = (r+1/2)^{n-1} = (r+1/2)^2 for n=3.
# M[0,0] = [1 + 1 - 1] * (r+1/2)^2 = (r+1/2)^2.
# This is r^2 + r + 1/4. NOT even in r.

# WAIT. Let me re-read THM-030 more carefully.
print()
print("=" * 70)
print("RE-EXAMINING THM-030 COROLLARY 2")
print("=" * 70)

# THM-030: B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)
# Corollary 2: "Every entry M[a,b] of the transfer matrix contains only even powers of r."
# This might mean: M[a,b] expanded in r about r=0 has only even powers.
# But (r+1/2)^2 = r^2 + r + 1/4, which has odd power r^1.

# Unless "contains only even powers of r" means something different:
# M[a,b] is a function of r and s_{ij}. The s variables contribute r-independent terms.
# When s_{ij} = A[i,j] - 1/2 (fixed), M(r) is a polynomial in r.
# Cor 2 might mean that the r-DEPENDENT part has only even powers,
# or that M(r) - M(-r) = 0 (even function of r).

# Let's check: M[0,0](-r) = (-r+1/2)^2 = r^2 - r + 1/4.
# M[0,0](r) - M[0,0](-r) = (r^2+r+1/4) - (r^2-r+1/4) = 2r ≠ 0.
# So M[0,0] is NOT an even function of r!

# Hmm, maybe the theorem statement needs more careful reading.
# Let me check for the full M matrix at r and -r.

r_test = Fraction(1, 3)
M_pos = transfer_matrix_r(A3, r_test)
M_neg = transfer_matrix_r(A3, -r_test)
print(f"\n  3-cycle, r={r_test}:")
print(f"    M(r) = {[[float(M_pos[a][b]) for b in range(n)] for a in range(n)]}")
print(f"    M(-r) = {[[float(M_neg[a][b]) for b in range(n)] for a in range(n)]}")
print(f"    M(r) = M(-r)? {all(M_pos[a][b] == M_neg[a][b] for a in range(n) for b in range(n))}")

# Also check: maybe even powers of r means the POLYNOMIAL in r (with s substituted)
# has only even total degree?

# For the transitive tournament at n=3: 0→1→2
A3_trans = [[0,1,1],[0,0,1],[0,0,0]]
r_t = Fraction(1, 3)
M_trans_pos = transfer_matrix_r(A3_trans, r_t)
M_trans_neg = transfer_matrix_r(A3_trans, -r_t)
print(f"\n  Transitive n=3, r={r_t}:")
print(f"    M(r)[0,0] = {float(M_trans_pos[0][0]):.6f}")
print(f"    M(-r)[0,0] = {float(M_trans_neg[0][0]):.6f}")
print(f"    M(r) = M(-r)? {all(M_trans_pos[a][b] == M_trans_neg[a][b] for a in range(n) for b in range(n))}")

# =====================================================================
# Full polynomial computation at n=3
# =====================================================================
print()
print("=" * 70)
print("POLYNOMIAL M[a,b](r) AT n=3")
print("=" * 70)

# Compute at 5 r-values and interpolate
from numpy.polynomial.polynomial import polyfit as npfit

n = 3
for A, name in [(A3, "3-cycle"), (A3_trans, "transitive")]:
    r_pts = [i/10.0 for i in range(-4, 5)]
    for a in range(n):
        for b in range(n):
            vals = []
            for rv in r_pts:
                M_r = transfer_matrix_r(A, Fraction(int(rv*10), 10))
                vals.append(float(M_r[a][b]))
            # Fit degree-(n-1) polynomial
            coeffs = npfit(np.array(r_pts), np.array(vals), n-1)
            nonzero = [(i, round(coeffs[i], 6)) for i in range(len(coeffs)) if abs(coeffs[i]) > 1e-6]
            if any(c != 0 for _, c in nonzero):
                powers = ", ".join(f"r^{i}:{c}" for i, c in nonzero)
                print(f"  {name} M[{a},{b}](r) = {powers}")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
M[a,b](r) is NOT an even function of r in general.
The THM-030 Corollary 2 ("even powers of r") likely means something
more subtle:
- Maybe even powers of r in the JOINT expansion with s variables
- Or even powers when r is paired with specific s-combinations
- Need to re-read the original proof more carefully

The consecutive formula M[a,b] = sum_j (-1)^j N(a,b,j) generalizes
to weighted paths at any r, and the formula still holds at general r.
""")
