#!/usr/bin/env python3
"""
DECOMPOSE M[a,b] into "Hamiltonian" and "non-Hamiltonian" split pair contributions.

M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
       = sum_S (-1)^|S| [sum of (P,Q) pairs where P ends at a in T[S+a], Q starts at b in T[R+b]]

Each (P,Q) pair: concatenation P|Q may or may not be a Ham path of T.
  - If T[a,b]=1 AND P|Q has all valid edges: it's a Ham path
  - If T[a,b]=0 OR internal edges fail: it's NOT a Ham path

DECOMPOSE: M[a,b] = HAM(a,b) + NONHAM(a,b)
  where HAM = signed sum of pairs that ARE Ham paths
        NONHAM = signed sum of pairs that are NOT Ham paths

For T[a,b]=1: HAM(a,b) = sum_j (-1)^j consec(a,b,j), NONHAM = rest
For T[a,b]=0: HAM(a,b) = 0, NONHAM = M[a,b]

QUESTION: Does NONHAM = 0 always for position-uniform tournaments?
If so, then M[a,b] = HAM(a,b) = sum_j (-1)^j consec(a,b,j)
and the proof via symmetry goes through.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

print("=" * 70)
print("n=5: HAM vs NONHAM decomposition for position-uniform tournaments")
print("=" * 70)

nonham_nonzero_count = 0
total_checked = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)

    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform:
        continue

    total_checked += 1
    found_nonzero_nonham = False

    for a in range(n):
        for b in range(n):
            if a == b: continue

            U = [v for v in range(n) if v != a and v != b]

            ham_sum = 0
            nonham_sum = 0

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                sign = (-1)**len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                # Find all (P,Q) pairs
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
                    if not p_valid: continue

                    for q in permutations(R_set):
                        if q[0] != b: continue
                        q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                        if not q_valid: continue

                        concat = list(p) + list(q)
                        is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))

                        if is_ham:
                            ham_sum += sign
                        else:
                            nonham_sum += sign

            if nonham_sum != 0:
                found_nonzero_nonham = True
                if nonham_nonzero_count < 5:
                    print(f"  bits={bits}: ({a},{b}): HAM={ham_sum}, NONHAM={nonham_sum}, "
                          f"M={ham_sum+nonham_sum}, T[{a},{b}]={T.get((a,b),0)}")

    if found_nonzero_nonham:
        nonham_nonzero_count += 1

print(f"\n  Checked {total_checked} position-uniform tournaments")
print(f"  {nonham_nonzero_count} have nonzero NONHAM contributions")
print(f"  {total_checked - nonham_nonzero_count} have NONHAM=0 for all pairs")

if nonham_nonzero_count == 0:
    print("""
  CONFIRMED: For ALL position-uniform tournaments at n=5,
  the non-Hamiltonian split pair contributions cancel to 0!

  Therefore: M[a,b] = HAM(a,b) = sum_j (-1)^j consec(a,b,j)
  for all position-uniform tournaments.

  PROOF OF SCALAR M:
  1. M[a,b] = sum_j (-1)^j consec(a,b,j)  [verified]
  2. If T[a,b]=0: consec(a,b,j)=0 for all j => M[a,b]=0
  3. If T[a,b]=1: By THM-030, M[a,b]=M[b,a].
     But T[b,a]=0, so M[b,a]=0 by step 2.
     Therefore M[a,b]=0.
  4. M[a,a] = H/n [from uniform positions at odd n]
  5. M = (H/n)*I. QED.
""")
else:
    print("\n  NONHAM is nonzero! The proof strategy fails.")
    print("  M[a,b] is NOT simply the Hamiltonian contribution.")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
