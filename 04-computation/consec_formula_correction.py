#!/usr/bin/env python3
"""
CORRECTION to kind-pasteur's split_pair_ham_path.py claim.

The consecutive-position formula
  M[a,b] = sum_j (-1)^j * #{paths with a@j, b@j+1}
FAILS at n=3 (6/8 mismatches), n=4 (64/64), n=5 (960/1024), n=6 (50/50).

The diff M_def - M_con has a VERY specific pattern at the transitive tournament:
  diff[i][i+1] = (-1)^{n+i}, all other entries 0.
  (upper bidiagonal with alternating signs)

Let's investigate what the diff IS in general and find the correct formula.
"""

from itertools import permutations, combinations
from collections import defaultdict
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
    """Transfer matrix via inclusion-exclusion definition."""
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

def consec_matrix(A):
    """M via consecutive-position formula."""
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok:
            for j in range(n):
                M[perm[j]][perm[j]] += (-1)**j
            for j in range(n-1):
                M[perm[j]][perm[j+1]] += (-1)**j
    return M

# =====================================================================
# Pattern in the diff D = M_ie - M_con
# =====================================================================
print("=" * 70)
print("DIFF PATTERN: D = M_ie - M_consec")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n--- n={n} ---")
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()

    # Group diffs
    diff_patterns = defaultdict(int)

    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1

        M_ie = transfer_matrix_ie(A)
        M_con = consec_matrix(A)
        diff = M_ie - M_con

        # Canonical form of diff
        diff_tuple = tuple(diff.flatten())
        diff_patterns[diff_tuple] += 1

    print(f"  {len(diff_patterns)} distinct diff patterns from {2**len(edges)} tournaments")

    # Show all patterns
    for diff_flat, count in sorted(diff_patterns.items(), key=lambda x: -x[1]):
        diff = np.array(diff_flat).reshape(n, n)
        if np.all(diff == 0):
            print(f"  count={count}: ZERO (formula matches)")
        else:
            # Describe the diff
            nonzero = [(i,j,int(diff[i][j])) for i in range(n) for j in range(n) if diff[i][j] != 0]
            print(f"  count={count}: {nonzero}")

# =====================================================================
# HYPOTHESIS: The diff D[a][b] might depend on whether there's a
# PATH ARC between a and b (i.e., |a-b|=1 in the natural ordering)
# =====================================================================
print()
print("=" * 70)
print("DIFF vs PATH ARCS (adjacency in natural order)")
print("=" * 70)

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]

# For each tournament, check if diff is always on "path arc" positions
all_path_arc = True
for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    M_ie = transfer_matrix_ie(A)
    M_con = consec_matrix(A)
    diff = M_ie - M_con

    for i in range(n):
        for j in range(n):
            if diff[i][j] != 0 and abs(i-j) != 1:
                all_path_arc = False
                break
        if not all_path_arc: break
    if not all_path_arc: break

print(f"  n=5: diff always on path-arc positions (|i-j|=1)? {all_path_arc}")

# Check n=3, n=4
for test_n in [3, 4]:
    edges_t = [(i,j) for i in range(test_n) for j in range(i+1,test_n)]
    ok = True
    for bits in range(2**len(edges_t)):
        A = [[0]*test_n for _ in range(test_n)]
        for idx, (i,j) in enumerate(edges_t):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        M_ie = transfer_matrix_ie(A)
        M_con = consec_matrix(A)
        diff = M_ie - M_con
        for i in range(test_n):
            for j in range(test_n):
                if diff[i][j] != 0 and abs(i-j) != 1:
                    ok = False
                    break
            if not ok: break
        if not ok: break
    print(f"  n={test_n}: diff always on path-arc positions? {ok}")

# =====================================================================
# The diff is ALWAYS on path-arc positions! So:
# D[i][i+1] and D[i+1][i] can be nonzero but D[i][j]=0 for |i-j|>1
#
# What determines D[i][i+1]?
# =====================================================================
print()
print("=" * 70)
print("WHAT DETERMINES D[i][i+1]?")
print("=" * 70)

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]

# For each tournament, record D[i][i+1] values and the path arc direction
# Path arc (i, i+1): A[i][i+1] = forward, A[i+1][i] = backward (= 1-A[i][i+1])
# The "tiling" model fixes path arcs as A[i+1][i]=1 (backward).
# But general tournaments can have either direction.

patterns = defaultdict(list)
for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    M_ie = transfer_matrix_ie(A)
    M_con = consec_matrix(A)
    diff = M_ie - M_con

    # Path arc directions
    path_arcs = tuple(A[i][i+1] for i in range(n-1))
    # D values on path-arc positions
    d_vals = tuple(int(diff[i][i+1]) for i in range(n-1))
    d_vals_back = tuple(int(diff[i+1][i]) for i in range(n-1))

    patterns[(path_arcs, d_vals, d_vals_back)] = patterns.get(
        (path_arcs, d_vals, d_vals_back), 0) + 1

print(f"  Distinct (path_arcs, D_forward, D_backward) patterns: {len(patterns)}")
for (pa, dv, db), count in sorted(patterns.items()):
    print(f"  path_arcs={pa}, D[i,i+1]={dv}, D[i+1,i]={db}: count={count}")

# =====================================================================
# HYPOTHESIS: D[i][i+1] = (-1)^{i+n} * A[i+1][i] (or similar)
# =====================================================================
print()
print("=" * 70)
print("TESTING D[i][i+1] FORMULAS")
print("=" * 70)

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]

# Test various formulas
formulas_ok = {
    "(-1)^(n+i)": True,
    "(-1)^(n+i) * A[i+1][i]": True,
    "(-1)^(n+i) * (2*A[i+1][i]-1)": True,
    "(-1)^i * A[i][i+1]": True,
}

for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    M_ie = transfer_matrix_ie(A)
    M_con = consec_matrix(A)
    diff = M_ie - M_con

    for i in range(n-1):
        d = int(diff[i][i+1])

        if d != (-1)**(n+i):
            formulas_ok["(-1)^(n+i)"] = False
        if d != (-1)**(n+i) * A[i+1][i]:
            formulas_ok["(-1)^(n+i) * A[i+1][i]"] = False
        if d != (-1)**(n+i) * (2*A[i+1][i]-1):
            formulas_ok["(-1)^(n+i) * (2*A[i+1][i]-1)"] = False
        if d != (-1)**i * A[i][i+1]:
            formulas_ok["(-1)^i * A[i][i+1]"] = False

for f, ok in formulas_ok.items():
    print(f"  D[i][i+1] = {f}: {ok}")

# =====================================================================
# Let me just tabulate D[i][i+1] vs A[i][i+1] for all cases
# =====================================================================
print()
print("  D[i][i+1] vs A[i][i+1] table (sample):")

# Collect unique (i, A[i][i+1], D[i][i+1]) triples
triples = defaultdict(int)
for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    M_ie = transfer_matrix_ie(A)
    M_con = consec_matrix(A)
    diff = M_ie - M_con

    for i in range(n-1):
        triples[(i, A[i][i+1], int(diff[i][i+1]))] = True

print(f"  {'pos i':>5} {'A[i,i+1]':>8} {'D[i,i+1]':>8}")
for (i, a, d) in sorted(triples.keys()):
    print(f"  {i:5d} {a:8d} {d:8d}")

# Also check D[i+1][i]
print(f"\n  {'pos i':>5} {'A[i+1,i]':>8} {'D[i+1,i]':>8}")
triples2 = defaultdict(int)
for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    M_ie = transfer_matrix_ie(A)
    M_con = consec_matrix(A)
    diff = M_ie - M_con

    for i in range(n-1):
        triples2[(i, A[i+1][i], int(diff[i+1][i]))] = True

for (i, a, d) in sorted(triples2.keys()):
    print(f"  {i:5d} {a:8d} {d:8d}")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
