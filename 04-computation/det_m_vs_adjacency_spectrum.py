#!/usr/bin/env python3
"""
CONNECTION: det(M) vs eigenvalues of the adjacency matrix A.

For tournament T with adjacency matrix A (0/1, A+A^T = J-I):
  - A has eigenvalues {(n-1)/2, lambda_1, ..., lambda_{n-1}} where
    lambda_i have real part -1/2 (for regular tournaments, purely imaginary)
  - The skew-adjacency matrix S = A - A^T has purely imaginary eigenvalues

Our transfer matrix M has real eigenvalues (M is symmetric).

Question: Is det(M) related to det(A), det(S), Pfaffian(S)?

For odd n: det(S) = 0 (skew-symmetric of odd order).
Pfaffian(S) only defined for even order.
But det(A) is nonzero for tournaments.

kind-pasteur-2026-03-06-S25c
"""

import numpy as np
from itertools import permutations

def count_H(T_dict, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T_dict.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def adj_matrix(T_dict, n):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = T_dict.get((i, j), 0)
    return A

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb2 += prod
            val += sign * ea * bb2
        return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_entry(T, n, a, a)
        for b in range(a+1, n):
            val = compute_M_entry(T, n, a, b)
            M[a, b] = val
            M[b, a] = val
    return M

def tournament_from_bits(n, bits_int):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits_int >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


# ============================================================
# n=5: det(M) vs det(A), det(S), eigenvalues of A
# ============================================================
print("=" * 70)
print("n=5: det(M) vs adjacency spectral data")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Process representatives
seen = set()
results = []

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    H = count_H(T, n)
    key = (scores, H)
    if key in seen:
        continue
    seen.add(key)

    A = adj_matrix(T, n)
    S = A - A.T  # Skew-adjacency
    M = compute_M(T, n)

    det_A = round(np.linalg.det(A).real)
    det_M = int(round(np.linalg.det(M.astype(float))))
    det_S = round(np.linalg.det(S).real)  # Should be 0 for odd n

    # Eigenvalues of A
    eig_A = sorted(np.linalg.eigvals(A), key=lambda x: -x.real)

    # per(A) = number of cycle covers
    per_A = 0
    for perm in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= A[i, perm[i]]
        per_A += prod
    per_A = int(round(per_A))

    results.append({
        'scores': scores, 'H': H, 'det_A': det_A, 'det_M': det_M,
        'det_S': det_S, 'per_A': per_A, 'eig_A': eig_A
    })

print(f"\n  {'scores':<20} {'H':>3} {'det(A)':>7} {'det(M)':>7} {'per(A)':>7} {'det(S)':>7}")
print(f"  {'-----':<20} {'---':>3} {'------':>7} {'------':>7} {'------':>7} {'------':>7}")

for r in sorted(results, key=lambda x: x['H']):
    print(f"  {str(r['scores']):<20} {r['H']:>3} {r['det_A']:>7} {r['det_M']:>7} {r['per_A']:>7} {r['det_S']:>7}")

# Check relationships
print("\n  Relationships:")
for r in sorted(results, key=lambda x: x['H']):
    det_A = r['det_A']
    det_M = r['det_M']
    per_A = r['per_A']
    H = r['H']

    # Is det_M related to det_A?
    if det_A != 0:
        ratio = det_M / det_A
        print(f"    H={H}: det(M)/det(A) = {ratio:.4f}, det(M)/H = {det_M/H:.4f}, per(A) = {per_A}")
    else:
        print(f"    H={H}: det(A) = 0!, det(M) = {det_M}, per(A) = {per_A}")


# ============================================================
# Try: det(M) = det(I + A - A^T) or similar?
# ============================================================
print("\n" + "=" * 70)
print("det(M) vs det(I + S) where S = A - A^T")
print("=" * 70)

for bits in [0, 100, 200, 500, 700, 900, 1023]:
    T = tournament_from_bits(n, bits)
    A = adj_matrix(T, n)
    S = A - A.T
    M = compute_M(T, n)

    H = count_H(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

    # Try various matrix expressions
    candidates = {
        'det(I+S)': round(np.linalg.det(np.eye(n) + S).real),
        'det(I-S)': round(np.linalg.det(np.eye(n) - S).real),
        'det(I+2S)': round(np.linalg.det(np.eye(n) + 2*S).real),
        'det(2I+S)': round(np.linalg.det(2*np.eye(n) + S).real),
        'det(A+A^T)': round(np.linalg.det(A + A.T).real),
        'det(A*A^T)': round(np.linalg.det(A @ A.T).real),
        'det(I+A*A^T)': round(np.linalg.det(np.eye(n) + A @ A.T).real),
    }

    print(f"\n  bits={bits}, H={H}, det(M)={det_M}, scores={scores}")
    for name, val in candidates.items():
        match = " <-- MATCH!" if val == det_M else ""
        print(f"    {name:20s} = {val:>10}{match}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
