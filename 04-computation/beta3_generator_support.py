"""
beta3_generator_support.py — Analyze vertex support of H_3 generators

For beta_3=2 tournaments, there are 2 independent H_3 generators.
Key question: do these generators have disjoint vertex support?

If generators need >= 4 vertices each and are vertex-disjoint,
then beta_3 <= floor(n/4). This would explain the observed bound.

Method:
1. Compute Z_3(T) = ker(d_3) in Omega_3(T) — the 3-cycles
2. Compute B_3(T) = im(d_4) in Z_3(T) — the 3-boundaries
3. Extract H_3 generators from Z_3/B_3 quotient
4. Analyze the vertex support of each generator

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import gc
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, enumerate_allowed_paths, build_adj_lists,
    boundary_faces, _gauss_rank_np, _gauss_nullbasis_modp,
    RANK_PRIME
)
from beta3_lean import fast_beta3_lean, _gauss_rank_i32, _gauss_nullbasis_i32

PRIME = 32749


def compute_h3_generators(A, n, prime=PRIME):
    """Compute H_3 generators and their vertex support.

    Returns list of (generator_A3_coords, vertex_support_set) tuples.
    """
    adj = build_adj_lists(A, n)
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p, adj)

    paths_sets = {p: set(ap.get(p, [])) for p in range(6)}
    paths_3 = ap.get(3, [])
    paths_4 = ap.get(4, [])

    if not paths_3:
        return []

    # Build constraint matrix for Omega_3
    non_allowed = {}
    na_count = 0
    entries = []
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if len(set(face)) == len(face) and face not in paths_sets[2]:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
                entries.append((non_allowed[face], j, sign))

    if na_count > 0:
        P3 = np.zeros((na_count, len(paths_3)), dtype=np.int32)
        for row, col, sign in entries:
            P3[row, col] = (P3[row, col] + sign) % prime
        rk_P3, nbasis3 = _gauss_nullbasis_i32(P3, na_count, len(paths_3), prime)
        dim_omega3 = len(paths_3) - rk_P3
        omega3_basis = np.array(nbasis3, dtype=np.int32) if nbasis3 else None
    else:
        dim_omega3 = len(paths_3)
        omega3_basis = np.eye(len(paths_3), dtype=np.int32)

    if dim_omega3 == 0 or omega3_basis is None:
        return []

    # Build boundary d_3: Omega_3 -> A_2
    paths_2 = ap.get(2, [])
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int32)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % prime

    # d_3 restricted to Omega_3
    d3_omega = (bd3.astype(np.int64) @ omega3_basis.T.astype(np.int64) % prime).astype(np.int32)

    # ker(d_3) in Omega_3 = Z_3
    rk_d3, null_d3 = _gauss_nullbasis_i32(d3_omega.copy(), d3_omega.shape[0], d3_omega.shape[1], prime)
    ker_d3_dim = dim_omega3 - rk_d3

    if ker_d3_dim == 0:
        return []

    Z3_omega = np.array(null_d3, dtype=np.int32)
    Z3_A3 = (Z3_omega.astype(np.int64) @ omega3_basis.astype(np.int64) % prime).astype(np.int32)

    # Build boundary d_4: A_4 -> A_3
    if not paths_4:
        # H_3 = Z_3 (no boundaries)
        generators = []
        for i in range(ker_d3_dim):
            gen = Z3_A3[i]
            support = set()
            for j, coeff in enumerate(gen):
                if coeff != 0:
                    support.update(paths_3[j])
            generators.append((gen, support))
        return generators

    idx3 = {p: i for i, p in enumerate(paths_3)}
    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int32)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % prime

    # Omega_4 constraint
    non_allowed4 = {}
    na4_count = 0
    entries4 = []
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if len(set(face)) == len(face) and face not in paths_sets[3]:
                if face not in non_allowed4:
                    non_allowed4[face] = na4_count
                    na4_count += 1
                entries4.append((non_allowed4[face], j, sign))

    if na4_count > 0:
        P4 = np.zeros((na4_count, len(paths_4)), dtype=np.int32)
        for row, col, sign in entries4:
            P4[row, col] = (P4[row, col] + sign) % prime
        rk_P4, nbasis4 = _gauss_nullbasis_i32(P4, na4_count, len(paths_4), prime)
        dim_omega4 = len(paths_4) - rk_P4
        omega4_basis = np.array(nbasis4, dtype=np.int32) if nbasis4 else None
    else:
        dim_omega4 = len(paths_4)
        omega4_basis = np.eye(len(paths_4), dtype=np.int32)

    if dim_omega4 == 0 or omega4_basis is None:
        # H_3 = Z_3
        generators = []
        for i in range(ker_d3_dim):
            gen = Z3_A3[i]
            support = set()
            for j, coeff in enumerate(gen):
                if coeff != 0:
                    support.update(paths_3[j])
            generators.append((gen, support))
        return generators

    # im(d_4) in A_3 coordinates
    im_d4 = (bd4.astype(np.int64) @ omega4_basis.T.astype(np.int64) % prime).astype(np.int32)

    # H_3 = Z_3 / (Z_3 ∩ im(d_4))
    # To get H_3 generators, find which Z_3 vectors are NOT in im(d_4)

    # Combine: [im_d4 | Z_3^T] and find rank
    # Z_3 vectors that are in span(im_d4) are boundaries
    # Others are H_3 generators

    # Method: row-reduce [im_d4 | Z_3^T] and identify which Z_3 columns add rank
    combined = np.concatenate([im_d4, Z3_A3.T], axis=1)
    rk_imd4 = _gauss_rank_i32(im_d4.copy(), prime)
    rk_combined = _gauss_rank_i32(combined.copy(), prime)
    beta3 = rk_combined - rk_imd4

    if beta3 == 0:
        return []

    # To find actual representatives, use extended RREF
    # For each Z_3 vector, check if it's independent of im_d4 + previously added
    generators = []
    current_space = im_d4.copy()
    current_rk = rk_imd4

    for i in range(ker_d3_dim):
        extended = np.concatenate([current_space, Z3_A3[i:i+1].T], axis=1)
        new_rk = _gauss_rank_i32(extended.copy(), prime)
        if new_rk > current_rk:
            gen = Z3_A3[i]
            support = set()
            for j, coeff in enumerate(gen):
                if coeff != 0:
                    support.update(paths_3[j])
            generators.append((gen, support))
            current_space = np.concatenate([current_space, Z3_A3[i:i+1].T], axis=1)
            current_rk = new_rk
            if len(generators) == beta3:
                break

    return generators


def main():
    print("=" * 70)
    print("H_3 GENERATOR VERTEX SUPPORT ANALYSIS")
    print("=" * 70)

    # Find beta_3=2 tournaments at n=8
    n = 8
    rng = np.random.RandomState(12345)
    b3_2_tours = []

    print("\n--- Finding beta_3=2 tournaments at n=8 ---")
    for trial in range(5000):
        A = random_tournament(n, rng)
        gc.collect()
        b3 = fast_beta3_lean(A, n)
        if b3 == 2:
            b3_2_tours.append((trial, A.copy()))
            scores = sorted([int(sum(A[i])) for i in range(n)])
            print(f"  Found b3=2 at trial {trial}, scores={scores}")

    print(f"\nFound {len(b3_2_tours)} beta_3=2 tournaments")

    # Analyze H_3 generators
    print(f"\n{'='*60}")
    print("H_3 GENERATOR ANALYSIS")
    print(f"{'='*60}")

    paths_3_map = {}  # for display

    for trial, A in b3_2_tours:
        scores = sorted([int(sum(A[i])) for i in range(n)])
        print(f"\n  Trial {trial}, scores={scores}")

        generators = compute_h3_generators(A, n)
        print(f"    Found {len(generators)} H_3 generators")

        if len(generators) >= 2:
            supp1 = generators[0][1]
            supp2 = generators[1][1]
            intersection = supp1 & supp2
            union = supp1 | supp2

            print(f"    Generator 1 support: {sorted(supp1)} ({len(supp1)} vertices)")
            print(f"    Generator 2 support: {sorted(supp2)} ({len(supp2)} vertices)")
            print(f"    Intersection: {sorted(intersection)} ({len(intersection)} vertices)")
            print(f"    Union: {sorted(union)} ({len(union)} vertices)")
            print(f"    Disjoint? {len(intersection) == 0}")

            # Count nonzero coefficients
            nz1 = sum(1 for c in generators[0][0] if c != 0)
            nz2 = sum(1 for c in generators[1][0] if c != 0)
            print(f"    Generator 1: {nz1} nonzero paths")
            print(f"    Generator 2: {nz2} nonzero paths")

            # Check if supports are complementary (together cover all vertices)
            missing = set(range(n)) - union
            print(f"    Missing vertices: {sorted(missing)} ({len(missing)} uncovered)")

    # Also check beta_3=1 tournaments for comparison
    print(f"\n{'='*60}")
    print("BETA_3=1 COMPARISON")
    print(f"{'='*60}")

    rng2 = np.random.RandomState(12345)
    b3_1_tours = []
    for trial in range(5000):
        A = random_tournament(n, rng2)
        b3 = fast_beta3_lean(A, n)
        if b3 == 1 and len(b3_1_tours) < 3:
            b3_1_tours.append((trial, A.copy()))
        if len(b3_1_tours) >= 3:
            break

    for trial, A in b3_1_tours[:2]:
        scores = sorted([int(sum(A[i])) for i in range(n)])
        print(f"\n  Trial {trial}, scores={scores}, b3=1")

        generators = compute_h3_generators(A, n)
        if generators:
            supp = generators[0][1]
            nz = sum(1 for c in generators[0][0] if c != 0)
            print(f"    Generator support: {sorted(supp)} ({len(supp)} vertices)")
            print(f"    Nonzero paths: {nz}")
            missing = set(range(n)) - supp
            print(f"    Missing vertices: {sorted(missing)} ({len(missing)})")

    # Floor(n/4) hypothesis check
    print(f"\n{'='*60}")
    print("FLOOR(n/4) HYPOTHESIS")
    print(f"{'='*60}")
    print("If generators can be chosen with disjoint 4-vertex support,")
    print("then beta_3 <= floor(n/4) (pigeonhole).")
    print()
    for trial, A in b3_2_tours:
        scores = sorted([int(sum(A[i])) for i in range(n)])
        generators = compute_h3_generators(A, n)
        if len(generators) >= 2:
            supp1, supp2 = generators[0][1], generators[1][1]
            overlap = len(supp1 & supp2)
            print(f"  Trial {trial} (scores={scores}): "
                  f"|supp1|={len(supp1)}, |supp2|={len(supp2)}, "
                  f"overlap={overlap}")
            if overlap > 0:
                print(f"    OVERLAPPING — BUT generators may be refocusable")
                # Try: for each pair of 4-element subsets that cover different vertices
                # This is more subtle since H_3 lives in a quotient

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
