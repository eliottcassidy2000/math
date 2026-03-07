#!/usr/bin/env python3
"""
Investigate the Pfaffian-Path duality:
  - At odd n: det(S) = 0, tr(M) = H != 0
  - At even n: det(S) = Pf(S)^2, tr(M) = 0

Is there a quantitative relationship between Pf(S) and H(T)?

Also check: det(S) at even n vs H at neighboring odd n.

kind-pasteur-2026-03-06-S25f
"""

import numpy as np
from itertools import permutations, combinations
from math import factorial
import random

def pfaffian(M):
    """Compute Pfaffian of skew-symmetric matrix M (even dimension)."""
    n = M.shape[0]
    assert n % 2 == 0, "Pfaffian requires even dimension"
    if n == 0:
        return 1
    if n == 2:
        return M[0, 1]
    # Use the definition: Pf(M) = sum over perfect matchings
    # For small n, this is feasible
    total = 0
    verts = list(range(n))
    # Generate all perfect matchings via recursive pairing
    def gen_matchings(remaining):
        if len(remaining) == 0:
            yield [], 1
            return
        first = remaining[0]
        rest = remaining[1:]
        for i, partner in enumerate(rest):
            # Pair (first, partner)
            new_remaining = rest[:i] + rest[i+1:]
            # Sign: number of transpositions to bring partner next to first
            sign = (-1) ** i
            for matching, s in gen_matchings(new_remaining):
                yield [(first, partner)] + matching, sign * s

    for matching, sign in gen_matchings(verts):
        prod = sign
        for (i, j) in matching:
            prod *= M[i, j]
        total += prod
    return total

def skew_adj(A):
    """Skew-adjacency matrix: S[i,j] = A[i,j] - A[j,i]."""
    n = len(A)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            S[i, j] = A[i][j] - A[j][i]
    return S

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def tournament_random(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

# =====================================================================
# EVEN n: Pfaffian and H relationship
# =====================================================================
print("=" * 70)
print("PFAFFIAN-PATH DUALITY")
print("=" * 70)

for n in [4, 6]:
    print(f"\n--- n={n} (EVEN) ---")
    print(f"{'seed':>4} {'H':>5} {'det(S)':>8} {'Pf(S)':>6} {'Pf^2':>8} {'H/|Pf|':>8}")

    data = []
    for seed in range(30):
        A = tournament_random(n, seed * 7 + n)
        S = skew_adj(A)
        H = ham_path_count(A)
        d = int(round(np.linalg.det(S)))
        pf = pfaffian(S)
        pf_int = int(round(pf))
        ratio = H / abs(pf_int) if pf_int != 0 else float('inf')
        data.append((seed, H, d, pf_int))
        print(f"{seed:4d} {H:5d} {d:8d} {pf_int:6d} {pf_int**2:8d} {ratio:8.4f}")

    # Collect unique (H, Pf) pairs
    unique = set((H, pf) for (_, H, _, pf) in data)
    print(f"\nUnique (H, Pf) pairs: {len(unique)}")

# =====================================================================
# ODD n: det(S) = 0, but compute adj(S) = cofactor matrix
# =====================================================================
print(f"\n--- n=5 (ODD) ---")
n = 5
for seed in range(10):
    A = tournament_random(n, seed * 11)
    S = skew_adj(A)
    H = ham_path_count(A)
    d = np.linalg.det(S)
    print(f"  seed={seed}: H={H:3d}, det(S)={d:.4f}")

# =====================================================================
# PATH-CYCLE DUALITY: quantitative check
# =====================================================================
print(f"\n--- PATH-CYCLE DUALITY ---")
print(f"At even n: Pf(S) counts CYCLE COVERS (with signs)")
print(f"At odd n: H(T) counts PATHS")
print(f"\nCheck: at n=4, Pf(S) = sum over perfect matchings of S = sum over 2-cycle covers")
print(f"A 2-cycle cover of 4 vertices = two directed 2-cycles = impossible in a tournament!")
print(f"Wait: S[i,j] = A[i,j] - A[j,i] = +-1. So Pf counts SIGNED cycle covers.")

n = 4
print(f"\nExhaustive n=4:")
for bits in range(64):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    scores = tuple(sorted(sum(A[i]) for i in range(4)))
    S = skew_adj(A)
    H = ham_path_count(A)
    pf = int(round(pfaffian(S)))

    # Count directed 3-cycles
    t3 = 0
    for v in combinations(range(4), 3):
        for p in permutations(v):
            if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3)):
                t3 += 1
    t3 //= 3

    if bits < 10 or pf**2 != 1:
        print(f"  bits={bits:2d}: scores={scores}, H={H}, Pf={pf}, Pf^2={pf**2}, t3={t3}")

# =====================================================================
# The key formula: at even n, Pf(S) = (-1)^{n/2} * sum_M prod S[i,j]
# where M ranges over perfect matchings.
# In a tournament: a perfect matching is a set of n/2 disjoint edges.
# Each edge (i,j) contributes S[i,j] = +1 or -1 depending on direction.
# So Pf counts signed matchings, not cycle covers.
#
# But det(S) = Pf(S)^2 does relate to cycle covers (permanent-like).
# Actually det(S) = sum_sigma sgn(sigma) prod S[i,sigma(i)].
# For skew-symmetric S: S[i,sigma(i)] = -S[sigma(i),i].
# A permutation sigma decomposes into cycles.
# Each fixed point contributes S[i,i] = 0, so only fixed-point-free sigma contribute.
# Each 2-cycle (i,j) contributes S[i,j]*S[j,i] = -S[i,j]^2 = -1.
# Each k-cycle (i1,...,ik) contributes prod_{cyclic} S[i_l, i_{l+1}].
# =====================================================================

print(f"\n--- INTERPRETATION ---")
print(f"For n=4 tournaments, Pf(S) = S[0,1]*S[2,3] - S[0,2]*S[1,3] + S[0,3]*S[1,2]")
print(f"Each S[i,j] = A[i,j] - A[j,i] = +-1.")
print(f"So Pf = (+-1)(+-1) - (+-1)(+-1) + (+-1)(+-1)")
print(f"Possible values: Pf in {{-3, -1, 1, 3}}")

n = 4
pf_counts = {}
for bits in range(64):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    S = skew_adj(A)
    pf = int(round(pfaffian(S)))
    H = ham_path_count(A)
    key = (pf, H)
    pf_counts[key] = pf_counts.get(key, 0) + 1

print(f"\n(Pf, H) distribution at n=4:")
for key in sorted(pf_counts.keys()):
    print(f"  Pf={key[0]:+2d}, H={key[1]:2d}: count={pf_counts[key]}")

print("\nDONE")
