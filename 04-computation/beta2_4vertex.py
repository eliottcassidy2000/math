#!/usr/bin/env python3
"""
beta2_4vertex.py - 4-vertex sub-tournament contribution to beta2=0

Every 3-path uses exactly 4 vertices. So 3-paths decompose by
their 4-vertex support. The boundary bd3 maps 4-local to 3-local.

Key question: Does the local Om3 for each 4-vertex subtournament
contribute enough to fill the global Z2?

For a tournament on 4 vertices, there are exactly 2 tournament types:
- Transitive: scores (0,1,2,3) - 24 labelings
- 4-cycle-like: scores (1,1,2,2) - 24 labelings
- WAIT: also scores (0,2,2,2) and (1,1,1,3) with 8 each

Let me compute the exact contribution of each type.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


# ============================================================
# Exhaustive 4-vertex tournament analysis
# ============================================================
print("=" * 70)
print("4-VERTEX TOURNAMENT TYPES")
print("=" * 70)

n = 4
total = 1 << (n*(n-1)//2)

score_types = defaultdict(list)
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, enumerate_allowed_paths(A, n, 1))
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    dd = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

    score_types[scores].append({
        'bits': bits, 'A2': len(a2), 'A3': len(a3),
        'Om2': d_om2, 'Om3': d_om3, 'DD': dd, 'c3': c3
    })

for sc in sorted(score_types.keys()):
    entries = score_types[sc]
    e = entries[0]
    print(f"  {sc}: {len(entries)} tours, c3={e['c3']}, |A2|={e['A2']}, |A3|={e['A3']}, "
          f"DD={e['DD']}, Om2={e['Om2']}, Om3={e['Om3']}")


# ============================================================
# How do 4-vertex subtournaments contribute at n=5?
# ============================================================
print(f"\n{'='*70}")
print("4-VERTEX CONTRIBUTIONS AT n=5")
print("=" * 70)

n = 5

for bits in [0, 10, 341]:  # transitive, c3=3, regular
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))

    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, enumerate_allowed_paths(A, n, 1))
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2, enumerate_allowed_paths(A, n, 1))
    bd3 = build_full_boundary_matrix(a3, a2)

    rk_bd2 = 0
    if d_om2 > 0:
        S = np.linalg.svd(bd2 @ om2, compute_uv=False)
        rk_bd2 = sum(s > 1e-8 for s in S)
    dim_z2 = d_om2 - rk_bd2

    print(f"\nbits={bits}, scores={scores}")
    print(f"  Global: Om2={d_om2}, Z2={dim_z2}, Om3={d_om3}")

    # For each 4-vertex subset, count contribution
    for S_set in combinations(range(n), 4):
        sub_scores = tuple(sorted([sum(A[i][j] for j in S_set if j != i) for i in S_set]))
        local_3 = [i for i, p in enumerate(a3) if set(p).issubset(S_set)]
        local_dd = [i for i in local_3 if A[a3[i][0]][a3[i][2]] and A[a3[i][1]][a3[i][3]]]

        # How many of these are in Om3?
        if d_om3 > 0:
            # Project Om3 onto local 3-paths
            local_om3 = om3[local_3, :]
            rk_local = np.linalg.matrix_rank(local_om3, tol=1e-8)
        else:
            rk_local = 0

        print(f"  {S_set} ({sub_scores}): {len(local_3)} paths, {len(local_dd)} DD, Om3_proj={rk_local}")


# ============================================================
# CRITICAL TEST: beta2=0 for 4-vertex subtournaments
# ============================================================
print(f"\n{'='*70}")
print("BETA2 OF 4-VERTEX SUBTOURNAMENTS")
print("=" * 70)

# At n=4, beta2=0 for ALL tournaments (verified above).
# This is the BASE CASE for an inductive proof!

# Score (0,1,2,3): transitive, A3=4, DD=4, Om3=4, Om2=4, Z2=1
# Score (0,2,2,2): one 3-cycle, A3=2, DD=0, Om3=0, Om2=3, Z2=0
# Score (1,1,1,3): one 3-cycle, A3=2, DD=0, Om3=0, Om2=3, Z2=0
# Score (1,1,2,2): two 3-cycles, A3=3, DD=1, Om3=1, Om2=3, Z2=1

print("\nVerification: beta2=0 for ALL 4-vertex tournaments")
for bits in range(1 << (4*3//2)):
    A = build_adj(4, bits)
    betti = path_betti_numbers(A, 4, max_dim=3)
    if betti[2] != 0:
        print(f"  COUNTEREXAMPLE at bits={bits}!")
        break
else:
    print("  CONFIRMED: beta2=0 for all 64 tournaments on 4 vertices")


# ============================================================
# Can we prove beta2=0 by induction on n?
# ============================================================
print(f"\n{'='*70}")
print("INDUCTION APPROACH")
print("=" * 70)

print("""
INDUCTION STRATEGY:
- Base case: n=4, beta2=0 (verified all 64)
- Inductive step: T on n vertices, assume beta2=0 for all (n-1)-vertex tournaments

For each vertex v in T, T\v is a tournament on n-1 vertices.
By induction, beta2(T\v) = 0.

QUESTION: Can we lift from beta2(T\v)=0 to beta2(T)=0?

The long exact sequence of the pair (T, T\v) gives:
  ... -> H_2(T) -> H_2(T, T\v) -> H_1(T\v) -> H_1(T) -> ...

For beta2(T)=0: need H_2(T) = 0, which means:
  im(H_2(T, T\v) -> H_1(T\v)) = ker(H_1(T\v) -> H_1(T))

This was explored before (HYP-231) and found to be equivalent to beta2=0.

ALTERNATIVE: Use arc-flip induction.
Start from transitive tournament (beta2=0 trivially).
Each flip preserves beta2=0 (HYP-233, verified).
Since any tournament is reachable from transitive by flips,
beta2=0 for all tournaments.

BUT: we need to PROVE that arc-flip preserves beta2=0,
not just verify computationally.
""")

# Check: is the transitive tournament beta2=0?
n = 5
A = build_adj(n, 0)  # bits=0 is transitive
print(f"\nTransitive T_5: betti = {path_betti_numbers(A, n, max_dim=4)}")

n = 6
A = build_adj(n, 0)
print(f"Transitive T_6: betti = {path_betti_numbers(A, n, max_dim=5)}")

n = 7
A = build_adj(n, 0)
print(f"Transitive T_7: betti = {path_betti_numbers(A, n, max_dim=6)}")


# The transitive tournament has all beta_p = 0 for p >= 1!
# It's contractible in path homology.
# So beta2=0 would follow from arc-flip invariance.


print("\n\nDone.")
