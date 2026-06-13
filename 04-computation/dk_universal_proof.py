#!/usr/bin/env python3
"""
dk_universal_proof.py — Prove D_k mod 2^{n-1-k} is universal for ALL k.

STRATEGY: Use the k-fold partition identity (k_fold_partition_identity.py):
For ANY tournament on n vertices and any 2k distinct vertices {v_1,...,v_{2k}}:
  sum over (2k)! perms of A[v_{P_1}][v_{P_2}]*...*A[v_{P_{2k-1}}][v_{P_{2k}}] = (2k)!/2^k

This gives D_S = n!/2^k for non-adjacent position k-tuples S.

For adjacent positions, the contributions depend on the tournament.
The question: do the adjacent contributions have enough factors of 2?

D_k = sum_{|S|=k positions} D_S

D_S = n!/2^k if positions in S are pairwise non-adjacent (gap >= 2).
D_S depends on tournament if some positions are adjacent.

Number of position k-tuples with ALL non-adjacent (no consecutive):
= C(n-k, k) (choosing k from n-1 positions with gaps)

Number with AT LEAST ONE adjacent pair = C(n-1, k) - C(n-k, k)

For the adjacent contributions, we need to show they are divisible
by 2^{n-1-k}. The non-adjacent part is:
  C(n-k, k) * n!/2^k

The tournament-dependent (adjacent) part involves path-counting with
shared vertices, which introduces factors of (n-3)!, (n-4)!, etc.

Let me verify: D_k = C(n-k,k)*n!/2^k + (adjacent contribution).
And check if adjacent contribution mod 2^{n-1-k} is universal.

Author: opus-2026-03-07-S43b
"""
from itertools import permutations, combinations
from collections import Counter
import math
import random

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_DS(A, n, positions):
    """Compute D_S = #{perms with all specified positions forward}."""
    count = 0
    for P in permutations(range(n)):
        if all(A[P[pos]][P[pos+1]] for pos in positions):
            count += 1
    return count

def are_nonadjacent(positions):
    """Check if all positions are pairwise non-adjacent."""
    pos = sorted(positions)
    return all(pos[i+1] - pos[i] >= 2 for i in range(len(pos)-1))

# Verify at n=5
n = 5
m = n*(n-1)//2
print(f"=== n={n}: D_S universality check ===")
print(f"Positions 0..{n-2} = 0..{n-2}")

# For each k, compute D_S for all position k-tuples
for k in range(1, n):
    req = 2**(n-1-k)
    nonadj_count = 0
    adj_count = 0

    nonadj_vals = Counter()
    adj_vals = Counter()

    bits_test = [0, 7, 15, 31]  # Test a few tournaments

    for bits in bits_test[:min(3, len(bits_test))]:
        A = tournament_from_bits(bits, n)

        for S in combinations(range(n-1), k):
            ds = compute_DS(A, n, S)

            if are_nonadjacent(S):
                nonadj_vals[ds] = nonadj_vals.get(ds, 0) + 1
                expected = math.factorial(n) // (2**k)
                if ds != expected:
                    print(f"  NON-ADJ MISMATCH k={k} S={S}: D_S={ds} expected={expected}")
            else:
                adj_vals[(bits, ds % req)] = 1

    # Check if adjacent D_S mod 2^{n-1-k} varies with tournament
    print(f"\n  k={k}: req=2^{n-1-k}={req}")

    # Exhaustive check
    all_ds_by_s = {}
    for S in combinations(range(n-1), k):
        vals = []
        for bits in range(1 << m):
            A = tournament_from_bits(bits, n)
            ds = compute_DS(A, n, S)
            vals.append(ds % req)

        unique_mods = set(vals)
        is_univ = (len(unique_mods) == 1)
        is_nonadj = are_nonadjacent(S)

        if not is_univ:
            print(f"    S={S} adj={not is_nonadj}: D_S mod {req} varies: {sorted(unique_mods)[:5]}...")
        else:
            print(f"    S={S} adj={not is_nonadj}: D_S mod {req} = {unique_mods.pop()} (UNIVERSAL)")

print("\n=== Summary ===")
print("The key question: is D_S mod 2^{n-1-k} universal for ALL position sets S?")
print("If yes, then D_k = sum D_S is universal mod 2^{n-1-k}.")
print("If no, the individual D_S can vary but their SUM might still be universal.")
