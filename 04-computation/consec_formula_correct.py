#!/usr/bin/env python3
"""
CORRECTED consecutive-position formula for M[a,b].

The original claim (kind-pasteur):
  M[a,b] = sum_j (-1)^j * #{paths with a@j, b@j+1}
FAILS because M[a,b] can be nonzero when A[a,b]=0, but paths can't
have a immediately before b without edge a→b.

KEY INSIGHT: When A[a,b]=1, the formula works: E_a(S∪{a}) * B_b(R∪{b})
equals the number of Ham paths with a at position |S| and b at position |S|+1,
because every split-path pair concatenates through the edge a→b.

When A[a,b]=0 (so A[b,a]=1), C(a,b,j) = 0 but M[a,b] = M[b,a] (symmetry),
and M[b,a] IS captured by C(b,a,j) since A[b,a]=1.

CORRECTED FORMULA:
  M[a,b] = sum_j (-1)^j * [C(a,b,j) + C(b,a,j)]   for a ≠ b
  where C(a,b,j) = #{Ham paths with a at pos j, b at pos j+1}

  = sum_j (-1)^j * #{paths with {a,b} at consecutive positions {j,j+1}}

This works because exactly one of A[a,b], A[b,a] is 1, and the corresponding
C term captures M[a,b] = M[b,a].
"""

from itertools import permutations, combinations
import numpy as np

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1: valid = False; break
        if valid: count += 1
    return count

def transfer_matrix_ie(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def corrected_consec_formula(A):
    """CORRECTED formula: M[a,b] = sum_j (-1)^j [C(a,b,j) + C(b,a,j)]"""
    n = len(A)
    M = np.zeros((n, n), dtype=int)

    # Build C[a][b][j] array
    C = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok:
            for j in range(n-1):
                C[perm[j]][perm[j+1]][j] += 1

    # Diagonal: position formula (unchanged)
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok:
            for a in range(n):
                M[a][a] += (-1)**(list(perm).index(a))

    # Off-diagonal: CORRECTED formula
    for a in range(n):
        for b in range(n):
            if a == b: continue
            for j in range(n-1):
                M[a][b] += (-1)**j * (C[a][b][j] + C[b][a][j])

    return M

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# =====================================================================
# VERIFY CORRECTED FORMULA
# =====================================================================
print("=" * 70)
print("VERIFYING CORRECTED FORMULA")
print("  M[a,b] = sum_j (-1)^j * [C(a,b,j) + C(b,a,j)]")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\nn={n}: ", end="", flush=True)
    match = 0
    total = 0
    for A in all_tournaments(n):
        M_ie = transfer_matrix_ie(A)
        M_cor = corrected_consec_formula(A)
        total += 1
        if np.array_equal(M_ie, M_cor):
            match += 1
        else:
            if match + (total - match) < 5:
                diff = M_ie - M_cor
                print(f"\n  FAIL: diff = {diff.tolist()}")
    print(f"{match}/{total} match", end="")
    if match == total:
        print(" -- CORRECTED FORMULA HOLDS!")
    else:
        print(f" -- STILL FAILS ({total-match} mismatches)")

# =====================================================================
# n=6 sample
# =====================================================================
print("\nn=6 (sample): ", end="", flush=True)
import random
random.seed(42)
n = 6
edges6 = [(i,j) for i in range(n) for j in range(i+1,n)]
match6 = 0
total6 = 0
for _ in range(30):
    bits = random.randint(0, 2**len(edges6)-1)
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges6):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    M_ie = transfer_matrix_ie(A)
    M_cor = corrected_consec_formula(A)
    total6 += 1
    if np.array_equal(M_ie, M_cor):
        match6 += 1
print(f"{match6}/{total6} match", end="")
if match6 == total6:
    print(" -- CORRECTED FORMULA HOLDS!")
else:
    print(f" -- STILL FAILS")

# =====================================================================
# INTERPRETATION
# =====================================================================
print()
print("=" * 70)
print("INTERPRETATION")
print("=" * 70)

# Show the decomposition M = M_con + M_con^T for a specific tournament
n = 5
# Transitive tournament
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1,n):
        A[i][j] = 1

M_ie = transfer_matrix_ie(A)

# Compute M_con (asymmetric)
M_con = np.zeros((n, n), dtype=int)
for perm in permutations(range(n)):
    ok = True
    for k in range(n-1):
        if A[perm[k]][perm[k+1]] != 1: ok = False; break
    if ok:
        for j in range(n-1):
            M_con[perm[j]][perm[j+1]] += (-1)**j

print(f"\nTransitive T_5 (0→1→2→3→4):")
print(f"  M_ie  = {M_ie.tolist()}")
print(f"  M_con = {M_con.tolist()}")
print(f"  M_con^T = {M_con.T.tolist()}")
print(f"  M_con + M_con^T (off-diag) matches M_ie (off-diag)? ", end="")
ok = True
for a in range(n):
    for b in range(n):
        if a != b and M_ie[a][b] != M_con[a][b] + M_con[b][a]:
            ok = False
print(ok)

# =====================================================================
# THEOREM STATEMENT
# =====================================================================
print()
print("=" * 70)
print("THEOREM (Corrected Consecutive-Position Formula)")
print("=" * 70)
print("""
For any tournament T on n vertices:

  M[a,b] = sum_{j=0}^{n-2} (-1)^j * N(a,b,j)     (a ≠ b)

where N(a,b,j) = #{Ham paths P with {a,b} at positions {j, j+1}}.

Equivalently: N(a,b,j) = C(a,b,j) + C(b,a,j) where
  C(a,b,j) = #{paths with a at pos j and b at pos j+1}.

PROOF:
  Case 1: A[a][b] = 1.
    For any partition S∪R = V\\{a,b}, every pair of sub-paths
    (P ending at a on S∪{a}, Q starting at b on R∪{b})
    concatenates to a Ham path P→a→b→Q since edge a→b exists.
    So C(a,b,j) = sum_{|S|=j} E_a(S∪{a}) * B_b(R∪{b}).
    And sum_j (-1)^j C(a,b,j) = M[a,b].
    Also C(b,a,j) = 0 since A[b,a] = 0 (no edge b→a).
    So N = C(a,b,j) + 0 works.

  Case 2: A[a][b] = 0, so A[b][a] = 1.
    By THM-030 symmetry: M[a,b] = M[b,a].
    Apply Case 1 to M[b,a] with edge b→a:
    M[b,a] = sum_j (-1)^j C(b,a,j).
    And C(a,b,j) = 0 since A[a,b] = 0.
    So N(a,b,j) = 0 + C(b,a,j) = C(b,a,j) works.

The key insight: M is SYMMETRIC, the consecutive formula is
ASYMMETRIC (nonzero only when the edge exists). The symmetrized
version captures both directions and recovers M exactly.

COROLLARY: M[a,b] = 0 for a ≠ b iff:
  sum_j (-1)^j * #{paths with {a,b} at positions {j,j+1}} = 0.

For scalar M = (H/n)*I: EVERY edge pair has vanishing alternating
consecutive-position sum. This is a strong combinatorial constraint.
""")
