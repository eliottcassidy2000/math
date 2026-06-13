#!/usr/bin/env python3
"""
beta2_beta1_check.py - Verify beta1 formula for tournaments

Question: Is beta1 = c3 for all tournaments?
Or is the relationship more complex?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os
import numpy as np
from collections import Counter
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


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    c3 += 1
    return c3


def compute_J2(A, n):
    J2 = 0
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            for b in range(n):
                if b != a and b != c and A[a][b] and A[b][c]:
                    J2 += 1
                    break
    return J2


print("=" * 70)
print("BETA1 AND RANK VERIFICATION")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n--- n={n} ---")
    total = 1 << (n*(n-1)//2)

    b1_vs_c3 = 0
    rk_bd2_data = Counter()

    for bits in range(total):
        A = build_adj(n, bits)
        a0 = enumerate_allowed_paths(A, n, 0)
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        # Compute beta1
        om1 = compute_omega_basis(A, n, 1, a1, a0)
        om2 = compute_omega_basis(A, n, 2, a2, a1)

        d_om1 = om1.shape[1] if om1.ndim == 2 else 0
        d_om2 = om2.shape[1] if om2.ndim == 2 else 0

        bd1 = build_full_boundary_matrix(a1, a0)
        bd2 = build_full_boundary_matrix(a2, a1)

        # rk(bd1|_Om1)
        if d_om1 > 0:
            S1 = np.linalg.svd(bd1 @ om1, compute_uv=False)
            rk_bd1 = sum(s > 1e-8 for s in S1)
        else:
            rk_bd1 = 0

        dim_z1 = d_om1 - rk_bd1

        # rk(bd2|_Om2)
        if d_om2 > 0:
            S2 = np.linalg.svd(bd2 @ om2, compute_uv=False)
            rk_bd2 = sum(s > 1e-8 for s in S2)
        else:
            rk_bd2 = 0

        beta1 = dim_z1 - rk_bd2
        c3 = count_c3(A, n)
        j2 = compute_J2(A, n)

        if beta1 != c3:
            b1_vs_c3 += 1

        scores = tuple(sorted([sum(row) for row in A]))
        rk_bd2_data[(scores, c3, beta1, d_om1, dim_z1, rk_bd2, d_om2, j2)] += 1

    print(f"  beta1 != c3: {b1_vs_c3}/{total}")

    if n <= 5:
        print(f"\n  {'scores':<20} {'c3':>3} {'b1':>3} {'Om1':>4} {'Z1':>3} {'rk2':>4} {'Om2':>4} {'J2':>3} {'cnt':>5}")
        for (sc, c3, b1, om1, z1, rk2, om2, j2), cnt in sorted(rk_bd2_data.items()):
            marker = " " if b1 == c3 else " **"
            print(f"  {str(sc):<20} {c3:>3} {b1:>3} {om1:>4} {z1:>3} {rk2:>4} {om2:>4} {j2:>3} {cnt:>5}{marker}")


# Now verify the correct Z2 formula
print(f"\n{'='*70}")
print("CORRECT Z2 FORMULA")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

print(f"\nn={n}:")
print(f"  Om1 = |A1| = C(n,2) = {n*(n-1)//2}")
print(f"  rk(bd1) = n-1 = {n-1}")
print(f"  Z1 = Om1 - rk(bd1) = {n*(n-1)//2 - (n-1)}")

# The issue: rk(bd2|_Om2) is NOT always C(n-1,2) - c3
# because Om1 might not equal A1 in all cases
# Wait - for tournaments, Om1 = A1 always (every 1-path has allowed 0-faces)
# And rk(bd1) = n-1 always (tournament is strongly connected => connected)

# So Z1 = C(n-1,2) always
# beta1 = Z1 - rk(bd2) = C(n-1,2) - rk(bd2)
# If beta1 = c3, then rk(bd2) = C(n-1,2) - c3

# Let me check directly: is beta1 = c3?

mismatches = 0
for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    c3 = count_c3(A, n)
    if betti[1] != c3:
        mismatches += 1
        scores = tuple(sorted([sum(row) for row in A]))
        if mismatches <= 3:
            print(f"  MISMATCH: bits={bits}, scores={scores}, beta1={betti[1]}, c3={c3}")

print(f"\n  beta1 != c3: {mismatches}/{total}")

if mismatches > 0:
    # What IS beta1?
    print(f"\n  beta1 distribution:")
    b1_dist = Counter()
    for bits in range(total):
        A = build_adj(n, bits)
        betti = path_betti_numbers(A, n, max_dim=2)
        c3 = count_c3(A, n)
        b1_dist[(c3, betti[1])] += 1

    print(f"  {'c3':>3} {'b1':>3} {'count':>5}")
    for (c3, b1), cnt in sorted(b1_dist.items()):
        marker = "" if c3 == b1 else " ***"
        print(f"  {c3:>3} {b1:>3} {cnt:>5}{marker}")


print("\n\nDone.")
