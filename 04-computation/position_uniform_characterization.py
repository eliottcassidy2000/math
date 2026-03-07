#!/usr/bin/env python3
"""
Is "position-uniform" the exact characterization of scalar M at odd n?

Position-uniform: N[v,j] = H/n for all v,j
(every vertex appears equally often in every position across all Ham paths)

This is MUCH stronger than just having M diagonal = H/n.
If it characterizes scalar M, then scalar M implies something very rigid.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np

def tournament_from_adj(A):
    return [row[:] for row in A]

def ham_paths(A):
    n = len(A)
    return [p for p in permutations(range(n))
            if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]

def position_matrix(A):
    """N[v][j] = number of Ham paths where vertex v is at position j."""
    n = len(A)
    N = [[0]*n for _ in range(n)]
    for p in ham_paths(A):
        for j, v in enumerate(p):
            N[v][j] += 1
    return N

def transfer_matrix(A):
    n = len(A)
    M = [[0]*n for _ in range(n)]
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

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        if all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1)):
            count += 1
    return count

def is_scalar(M, n):
    H_over_n = M[0][0]
    for i in range(n):
        for j in range(n):
            if i == j and M[i][j] != H_over_n:
                return False
            if i != j and M[i][j] != 0:
                return False
    return True

def is_position_uniform(N, n):
    if not N:
        return False
    val = N[0][0]
    return all(N[v][j] == val for v in range(n) for j in range(n))

# =====================================================================
print("=" * 70)
print("n=5: POSITION-UNIFORM vs SCALAR M — EXHAUSTIVE CHECK")
print("=" * 70)

n = 5
# Generate all tournaments on n=5 (up to relabeling, but we check all)
# There are C(5,2) = 10 edges, 2^10 = 1024 tournaments.

results = {"PU_and_SM": 0, "PU_not_SM": 0, "SM_not_PU": 0, "neither": 0}
interesting = []

for bits in range(1024):
    A = [[0]*n for _ in range(n)]
    edge_idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> edge_idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            edge_idx += 1

    H = sum(1 for p in permutations(range(n)) if all(A[p[k]][p[k+1]] == 1 for k in range(n-1)))

    if H == 0:
        results["neither"] += 1
        continue

    N = position_matrix(A)
    pu = is_position_uniform(N, n)

    if H % n != 0:
        if pu:
            results["PU_not_SM"] += 1
            interesting.append(("PU_not_SM", bits, A, H, N))
        else:
            results["neither"] += 1
        continue

    M = transfer_matrix(A)
    sm = is_scalar(M, n)

    if pu and sm:
        results["PU_and_SM"] += 1
    elif pu and not sm:
        results["PU_not_SM"] += 1
        interesting.append(("PU_not_SM", bits, A, H, N))
    elif sm and not pu:
        results["SM_not_PU"] += 1
        interesting.append(("SM_not_PU", bits, A, H, N))
    else:
        results["neither"] += 1

print(f"\n  Results (all {1024} tournaments on 5 vertices):")
print(f"    Position-uniform AND scalar M: {results['PU_and_SM']}")
print(f"    Position-uniform but NOT scalar M: {results['PU_not_SM']}")
print(f"    Scalar M but NOT position-uniform: {results['SM_not_PU']}")
print(f"    Neither: {results['neither']}")

if results['PU_not_SM'] == 0 and results['SM_not_PU'] == 0:
    print("\n  CONCLUSION: Position-uniform <=> scalar M at n=5!")
elif results['SM_not_PU'] == 0:
    print("\n  CONCLUSION: Position-uniform => scalar M, but not conversely")
elif results['PU_not_SM'] == 0:
    print("\n  CONCLUSION: Scalar M => position-uniform, but not conversely")
else:
    print("\n  CONCLUSION: Neither implies the other")

for label, bits, A, H, N in interesting[:3]:
    print(f"\n  Example ({label}): bits={bits}, H={H}")
    for row in N:
        print(f"    {row}")

# =====================================================================
print("\n" + "=" * 70)
print("n=3: POSITION-UNIFORM vs SCALAR M")
print("=" * 70)

n = 3
results3 = {"PU_and_SM": 0, "PU_not_SM": 0, "SM_not_PU": 0, "neither": 0}

for bits in range(8):
    A = [[0]*3 for _ in range(3)]
    pairs = [(0,1), (0,2), (1,2)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    H = sum(1 for p in permutations(range(3)) if all(A[p[k]][p[k+1]] == 1 for k in range(2)))
    N = position_matrix(A)
    pu = is_position_uniform(N, 3)

    M = transfer_matrix(A)
    sm = is_scalar(M, 3)

    if pu and sm:
        results3["PU_and_SM"] += 1
    elif pu and not sm:
        results3["PU_not_SM"] += 1
    elif sm and not pu:
        results3["SM_not_PU"] += 1
    else:
        results3["neither"] += 1

print(f"\n  Results (all 8 tournaments on 3 vertices):")
for k, v in results3.items():
    print(f"    {k}: {v}")

# =====================================================================
# Now check: does position-uniform have a nice structural characterization?
# =====================================================================
print("\n" + "=" * 70)
print("STRUCTURE OF POSITION-UNIFORM TOURNAMENTS AT n=5")
print("=" * 70)

print("\n  All position-uniform (= scalar-M) tournaments:")
for bits in range(1024):
    A = [[0]*n for _ in range(n)]
    edge_idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> edge_idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            edge_idx += 1

    H = sum(1 for p in permutations(range(n)) if all(A[p[k]][p[k+1]] == 1 for k in range(n-1)))
    if H == 0 or H % n != 0:
        continue

    N = position_matrix(A)
    if not is_position_uniform(N, n):
        continue

    scores = tuple(sorted(sum(row) for row in A))
    # Check anti-automorphism count
    anti_count = 0
    for perm in permutations(range(n)):
        if all(A[i][j] == A[perm[j]][perm[i]] for i in range(n) for j in range(n)):
            anti_count += 1
    # Check automorphism count
    aut_count = 0
    for perm in permutations(range(n)):
        if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n)):
            aut_count += 1

    print(f"    bits={bits:010b}: H={H}, scores={scores}, |Aut|={aut_count}, |Anti|={anti_count}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
