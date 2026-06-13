#!/usr/bin/env python3
"""
CRITICAL CHECK: The consecutive-position formula
  M[a,b] = sum_j (-1)^j * #{paths with a@pos_j, b@pos_{j+1}}
was verified at n=5 by kind-pasteur, but FAILS at n=3,4 (our test).

This script verifies:
1. Does it hold at n=5? (confirm kind-pasteur's result)
2. Does it hold at n=6? (even n)
3. Does it hold at n=7? (next odd n — expensive but sample)
4. WHY does it work at n=5 but not n=3,4?
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

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            if a == b:
                # Diagonal via position formula
                val = 0
                for perm in permutations(range(n)):
                    ok = True
                    for k in range(n-1):
                        if A[perm[k]][perm[k+1]] != 1: ok = False; break
                    if ok:
                        val += (-1)**(list(perm).index(a))
                M[a][a] = val
            else:
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

def consec_formula(A):
    """M via consecutive-position formula (off-diagonal only)."""
    n = len(A)
    M = np.zeros((n, n), dtype=int)

    # Enumerate all Ham paths
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok: paths.append(perm)

    # Diagonal: standard position formula
    for a in range(n):
        for p in paths:
            M[a][a] += (-1)**(list(p).index(a))

    # Off-diagonal: consecutive formula
    for p in paths:
        for j in range(n-1):
            a, b = p[j], p[j+1]
            M[a][b] += (-1)**j

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
# Test at each n
# =====================================================================
for n in [3, 4, 5]:
    print(f"n={n}: ", end="", flush=True)
    match = 0
    total = 0
    first_fail = None
    for A in all_tournaments(n):
        M_def = transfer_matrix(A)
        M_con = consec_formula(A)
        total += 1
        if np.array_equal(M_def, M_con):
            match += 1
        elif first_fail is None:
            first_fail = (A, M_def, M_con)

    print(f"{match}/{total} match", end="")
    if match == total:
        print(" -- FORMULA HOLDS!")
    else:
        print(f" -- FORMULA FAILS ({total-match} mismatches)")
        if first_fail:
            A, M_def, M_con = first_fail
            diff = M_def - M_con
            print(f"  First failure diff (M_def - M_con): {diff.tolist()}")
            # Check: is diff always symmetric?
            print(f"  Diff symmetric? {np.array_equal(diff, diff.T)}")
            # Check: is diff always off-diagonal?
            print(f"  Diff diagonal = {[diff[i][i] for i in range(n)]}")

# =====================================================================
# n=5: also check formula for a!=b only (diagonal always matches)
# =====================================================================
print("\nn=5 detailed: checking off-diagonal only...")
n = 5
match_offdiag = 0
total_offdiag = 0
for A in all_tournaments(n):
    M_def = transfer_matrix(A)
    M_con = consec_formula(A)
    total_offdiag += 1
    offdiag_match = True
    for a in range(n):
        for b in range(n):
            if a != b and M_def[a][b] != M_con[a][b]:
                offdiag_match = False
                break
        if not offdiag_match: break
    if offdiag_match:
        match_offdiag += 1
print(f"  Off-diagonal: {match_offdiag}/{total_offdiag} match")

# =====================================================================
# Analysis: what's the ACTUAL off-diagonal relationship at n=3?
# =====================================================================
print()
print("=" * 70)
print("n=3: WHAT IS THE ACTUAL OFF-DIAGONAL RELATIONSHIP?")
print("=" * 70)

n = 3
for A in all_tournaments(n):
    M_def = transfer_matrix(A)
    M_con = consec_formula(A)
    diff = M_def - M_con

    # Is diff related to something?
    # For a!=b: diff[a][b] = M_def[a][b] - sum_j (-1)^j C(a,b,j)
    # The consecutive formula misses the "non-adjacent" contributions
    # in the inclusion-exclusion.

    H = sum(1 for p in permutations(range(n))
            if all(A[p[k]][p[k+1]] == 1 for k in range(n-1)))

    if H == 1:
        print(f"\n  Transitive (H=1):")
        print(f"    M_def = {M_def.tolist()}")
        print(f"    M_con = {M_con.tolist()}")
        print(f"    diff  = {diff.tolist()}")

        # The difference captures terms where a and b are NOT adjacent
        # in the path. M[a,b] from definition includes S-splits where
        # a and b can be far apart.

        # For n=3, there's only 1 "other" vertex. The S=∅ term gives
        # E_a({a})*B_b({c,b}) and the S={c} term gives E_a({a,c})*B_b({b}).
        # The consecutive formula only captures paths where a→b is an edge.
        # The IE definition captures the product E_a * B_b regardless.
        break

# =====================================================================
# KEY QUESTION: Why does it work at n=5?
# =====================================================================
print()
print("=" * 70)
print("WHY DOES THE FORMULA WORK AT n=5?")
print("=" * 70)
print("""
At n=3,4: M[a,b] ≠ sum_j (-1)^j C(a,b,j)
At n=5:   M[a,b] = sum_j (-1)^j C(a,b,j)  (verified exhaustively)

The IE definition:
  M[a,b] = sum_S (-1)^|S| E_a(S∪{a}) B_b(R∪{b})

The consecutive formula:
  sum_j (-1)^j C(a,b,j) = sum over Ham paths P, for each edge a→b in P
                            at position j, contribute (-1)^j.

These are DIFFERENT combinatorial objects in general!
The IE sum includes terms where a,b are in SEPARATE sub-paths
(not necessarily connected). The consecutive formula only counts
terms where a,b are ADJACENT in a full Ham path.

CONJECTURE: The formula works at n=5 because of a specific
cancellation that happens at this size. This might be related to:
  - The OCF structure (alpha_2 = 0 at n=5)
  - The fact that n-3 = 2 (the "gap" between a,b in split paths)
  - Some algebraic identity specific to 5 vertices

Alternatively: maybe n=5 is NOT the threshold. Let me check n=6.
""")

# =====================================================================
# n=6: Sample test (full enumeration too expensive)
# =====================================================================
print("=" * 70)
print("n=6: SAMPLE TEST")
print("=" * 70)

import random
n = 6

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

tiles6 = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles6.append((a, b))
m6 = len(tiles6)

match6 = 0
fail6 = 0
random.seed(42)
sample = random.sample(range(2**m6), min(50, 2**m6))

for bits in sample:
    A = tiling_to_tournament(bits, n)
    M_def = transfer_matrix(A)
    M_con = consec_formula(A)
    if np.array_equal(M_def, M_con):
        match6 += 1
    else:
        fail6 += 1
        if fail6 <= 2:
            diff = M_def - M_con
            print(f"  FAIL: diff = {diff.tolist()}")

print(f"\nn=6 sample: {match6}/{match6+fail6} match")
if fail6 > 0:
    print("  Formula FAILS at n=6 too!")
else:
    print("  Formula holds for all sampled tournaments!")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
