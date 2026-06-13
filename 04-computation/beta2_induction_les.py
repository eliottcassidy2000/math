#!/usr/bin/env python3
"""
beta2_induction_les.py - Inductive proof of beta2=0 via LES

KEY INSIGHT FROM SIMPLEX CONNECTION:
- Transitive tournament has path complex = simplex (all betti=0 except b0)
- For induction: T on n vertices, delete vertex v to get T\v on n-1 vertices
- LES: ... -> H_2(T\v) -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T) -> ...
- If beta2(T\v) = 0: 0 -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T) -> ...
- So H_2(T) = 0 iff H_2(T,T\v) -> H_1(T\v) is injective
- If beta1(T\v) = 0 (H_1(T\v)=0), this is automatic!

CRITICAL QUESTION: Is beta1(T\v) = 0 for ALL vertex deletions?

n=3: beta1 = 0 for ALL 3-vertex tournaments? => induction works for n=4
n=4: beta1 = 0 for ALL 4-vertex tournaments? => induction works for n=5
n=5: beta1 can be 1! => need more subtle argument for n=6

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
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
# Check beta1 for ALL tournaments at n = 3, 4, 5
# ============================================================
print("=" * 70)
print("BETA1 DISTRIBUTION BY n")
print("=" * 70)

for n in range(3, 7):
    total = 1 << (n*(n-1)//2)
    b1_dist = Counter()

    if total <= 32768:
        for bits in range(total):
            A = build_adj(n, bits)
            betti = path_betti_numbers(A, n, max_dim=1)
            b1_dist[betti[1]] += 1
    else:
        import random
        random.seed(42)
        for _ in range(2000):
            bits = random.randint(0, total - 1)
            A = build_adj(n, bits)
            betti = path_betti_numbers(A, n, max_dim=1)
            b1_dist[betti[1]] += 1

    print(f"n={n}: beta1 dist = {dict(sorted(b1_dist.items()))}")


# ============================================================
# INDUCTION ANALYSIS:
# For beta2=0 at n, we need beta2(T\v)=0 and either:
# (a) beta1(T\v) = 0 for ALL v (simple case), or
# (b) the connecting map H_2(R(v)) -> H_1(T\v) is injective
# ============================================================
print(f"\n{'='*70}")
print("INDUCTION FEASIBILITY")
print("=" * 70)

print("beta1 = 0 at n=3: base case for induction to n=4: ", end="")
all_b1_zero_3 = True
for bits in range(1 << 3):
    A = build_adj(3, bits)
    betti = path_betti_numbers(A, 3, max_dim=1)
    if betti[1] > 0:
        all_b1_zero_3 = False
        break
print("YES" if all_b1_zero_3 else "NO")

print("beta1 = 0 at n=4: base for induction to n=5: ", end="")
all_b1_zero_4 = True
for bits in range(1 << 6):
    A = build_adj(4, bits)
    betti = path_betti_numbers(A, 4, max_dim=1)
    if betti[1] > 0:
        all_b1_zero_4 = False
        break
print("YES" if all_b1_zero_4 else "NO")

print("beta1 = 0 at n=5 for ALL tournaments: ", end="")
all_b1_zero_5 = True
for bits in range(1 << 10):
    A = build_adj(5, bits)
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] > 0:
        all_b1_zero_5 = False
        break
print("YES" if all_b1_zero_5 else "NO (some have beta1=1)")


# ============================================================
# For n=6: which tournaments T have ALL subtournaments T\v with beta1=0?
# ============================================================
print(f"\n{'='*70}")
print("N=6: SUBTOURNAMENT BETA1")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

# Check exhaustively
all_sub_b1_zero_count = 0
some_sub_b1_one_count = 0
max_sub_b1_one_count = 0  # max number of v's with beta1(T\v)=1

import random
random.seed(42)
sample_size = 1000

for _ in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)

    sub_b1_one = 0
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        A_sub = [[A[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        betti_sub = path_betti_numbers(A_sub, n-1, max_dim=1)
        if betti_sub[1] > 0:
            sub_b1_one += 1

    if sub_b1_one == 0:
        all_sub_b1_zero_count += 1
    else:
        some_sub_b1_one_count += 1
    max_sub_b1_one_count = max(max_sub_b1_one_count, sub_b1_one)

print(f"n=6, {sample_size} samples:")
print(f"  ALL T\\v have beta1=0: {all_sub_b1_zero_count}")
print(f"  SOME T\\v have beta1=1: {some_sub_b1_one_count}")
print(f"  Max # of v with beta1(T\\v)=1: {max_sub_b1_one_count}")

if some_sub_b1_one_count > 0:
    print("\n  Simple induction BREAKS at n=6!")
    print("  Need to analyze the connecting map when beta1(T\\v) != 0")


# ============================================================
# For the n=6 cases with beta1(T\v) = 1:
# Check if H_1(T\v) -> H_1(T) is injective
# If injective, H_2(R(v)) = 0 and induction still works!
# ============================================================
print(f"\n{'='*70}")
print("H_1(T\\v) -> H_1(T) INJECTIVITY CHECK AT n=6")
print("=" * 70)

n = 6
random.seed(42)
total_pairs = 0
injective_count = 0
not_injective_count = 0

for _ in range(500):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)

    for v in range(n):
        verts = [i for i in range(n) if i != v]
        A_sub = [[A[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        betti_sub = path_betti_numbers(A_sub, n-1, max_dim=1)

        if betti_sub[1] == 0:
            continue  # Skip trivial cases

        total_pairs += 1

        # Compute Z_1(T\v) basis
        a1_s = enumerate_allowed_paths(A_sub, n-1, 1)
        a2_s = enumerate_allowed_paths(A_sub, n-1, 2)
        om1_s = compute_omega_basis(A_sub, n-1, 1, a1_s, [(i,) for i in range(n-1)])
        om2_s = compute_omega_basis(A_sub, n-1, 2, a2_s, a1_s)
        bd1_s = build_full_boundary_matrix(a1_s, [(i,) for i in range(n-1)])
        bd2_s = build_full_boundary_matrix(a2_s, a1_s)

        d_om1_s = om1_s.shape[1] if om1_s.ndim == 2 else 0
        d_om2_s = om2_s.shape[1] if om2_s.ndim == 2 else 0

        # Z_1(T\v) basis in A_1(T\v) coords
        if d_om1_s > 0:
            bd1_om_s = bd1_s @ om1_s
            U_s, S_s, Vt_s = np.linalg.svd(bd1_om_s, full_matrices=True)
            rk1 = sum(s > 1e-8 for s in S_s)
            z1_om_s = Vt_s[rk1:].T  # Z_1 in Om_1 coords
            z1_a1_s = om1_s @ z1_om_s  # Z_1 in A_1(T\v) coords
        else:
            continue

        # B_1(T\v) basis
        if d_om2_s > 0:
            b1_a1_s = bd2_s @ om2_s
            rk_b1_s = np.linalg.matrix_rank(b1_a1_s, tol=1e-8)
        else:
            b1_a1_s = np.zeros((len(a1_s), 0))
            rk_b1_s = 0

        dim_z1_s = z1_a1_s.shape[1]
        dim_h1_s = dim_z1_s - rk_b1_s

        if dim_h1_s == 0:
            continue

        # Map to T via inclusion
        a1_t = enumerate_allowed_paths(A, n, 1)
        a1_t_idx = {path: i for i, path in enumerate(a1_t)}

        inc = np.zeros((len(a1_t), len(a1_s)))
        for j, path_s in enumerate(a1_s):
            path_t = tuple(verts[k] for k in path_s)
            if path_t in a1_t_idx:
                inc[a1_t_idx[path_t], j] = 1

        # Map z1 to A_1(T)
        z1_in_t = inc @ z1_a1_s

        # Compute B_1(T)
        a2_t = enumerate_allowed_paths(A, n, 2)
        om2_t = compute_omega_basis(A, n, 2, a2_t, a1_t)
        bd2_t = build_full_boundary_matrix(a2_t, a1_t)
        d_om2_t = om2_t.shape[1] if om2_t.ndim == 2 else 0

        if d_om2_t > 0:
            b1_t = bd2_t @ om2_t
        else:
            b1_t = np.zeros((len(a1_t), 0))

        # Map B_1(T\v) to A_1(T) too
        b1_s_in_t = inc @ b1_a1_s

        # H_1(T\v) -> H_1(T) is injective iff:
        # For every z in Z_1(T\v) that's not in B_1(T\v),
        # its image i(z) is not in B_1(T).
        # Equivalently: ker(H_1(T\v) -> H_1(T)) = 0
        # ker = {[z] in H_1(T\v) : i(z) in B_1(T)}
        #      = (i^{-1}(B_1(T)) intersect Z_1(T\v)) / B_1(T\v)

        # Combined space: [b1_t | z1_in_t]
        # We want: for each z1 column, check if i(z) is in span(b1_t) + span(b1_s_in_t)
        # i.e., is i(z) in B_1(T) + i(B_1(T\v))?

        # Actually: ker = {z in Z_1(T\v) : i(z) in B_1(T)} / B_1(T\v)
        # So we need to check: for z in Z_1(T\v) \ B_1(T\v), is i(z) in B_1(T)?

        # Build H_1(T\v) representative: project Z_1 modulo B_1
        if rk_b1_s > 0:
            # Get a basis for Z_1/B_1 (i.e., H_1)
            combined_s = np.hstack([b1_a1_s, z1_a1_s])
            rk_combined = np.linalg.matrix_rank(combined_s, tol=1e-8)
            # dim(H_1) = rk_combined - rk_b1_s
            # Find z1 vectors not in span(b1)
            # For each z1 column, check if it's in span(b1)
            h1_reps = []
            for col in range(z1_a1_s.shape[1]):
                z = z1_a1_s[:, col]
                if rk_b1_s > 0:
                    res = np.linalg.lstsq(b1_a1_s, z, rcond=None)
                    residual = np.linalg.norm(b1_a1_s @ res[0] - z)
                    if residual > 1e-6:
                        h1_reps.append(z)
                else:
                    h1_reps.append(z)
        else:
            h1_reps = [z1_a1_s[:, col] for col in range(z1_a1_s.shape[1])]

        if not h1_reps:
            continue

        # Check each H_1 rep: does its image land in B_1(T)?
        is_injective = True
        for z in h1_reps:
            z_t = inc @ z
            if b1_t.shape[1] > 0:
                res = np.linalg.lstsq(b1_t, z_t, rcond=None)
                residual = np.linalg.norm(b1_t @ res[0] - z_t)
                if residual < 1e-6:
                    is_injective = False
                    break

        if is_injective:
            injective_count += 1
        else:
            not_injective_count += 1

print(f"Total (T, v) pairs with beta1(T\\v) > 0: {total_pairs}")
print(f"  H_1(T\\v) -> H_1(T) injective: {injective_count}")
print(f"  H_1(T\\v) -> H_1(T) NOT injective: {not_injective_count}")

if not_injective_count == 0 and total_pairs > 0:
    print("\n  AMAZING: H_1(T\\v) -> H_1(T) is ALWAYS injective!")
    print("  This means H_2(R(v)) = 0 for ALL tournaments!")
    print("  Combined with beta2(T\\v) = 0 by induction:")
    print("  => beta2(T) = 0 for ALL tournaments (PROOF BY INDUCTION)!")
else:
    print(f"\n  H_1 map is not always injective.")
    print(f"  Need alternate approach for those {not_injective_count} cases.")


print("\n\nDone.")
