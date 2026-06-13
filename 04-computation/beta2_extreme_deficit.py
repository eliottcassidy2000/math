#!/usr/bin/env python3
"""
beta2_extreme_deficit.py - Investigate extreme DT deficit cases

At n=7: found gap=6 (!!). What are these tournaments?
Also: what's the structure at n=8?

For extreme deficit cases, analyze:
1. Score sequence
2. Number of transitive triples / 4-tuples
3. dim(Omega_3) vs |DT|
4. How the non-DT Omega_3 elements fill the gap

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
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

def random_tournament(n):
    return random.getrandbits(n*(n-1)//2)

def find_DT_paths(A, n):
    dt = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if not A[a][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c or not A[c][d]: continue
                    if A[b][d]:
                        dt.append((a,b,c,d))
    return dt

def compute_DT_boundary_matrix(dt_paths, allowed_2, n):
    if not dt_paths or not allowed_2:
        return np.zeros((len(allowed_2), 0))
    path_to_idx = {p: i for i, p in enumerate(allowed_2)}
    mat = np.zeros((len(allowed_2), len(dt_paths)))
    for j, (a, b, c, d) in enumerate(dt_paths):
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        signs = [1, -1, 1, -1]
        for face, sign in zip(faces, signs):
            if face in path_to_idx:
                mat[path_to_idx[face], j] += sign
    return mat

def count_transitive_tuples(A, n, k):
    """Count k-tuples where the induced subtournament is transitive."""
    from itertools import combinations
    count = 0
    for subset in combinations(range(n), k):
        # Check if induced tournament on subset is transitive
        # (has a unique topological ordering)
        S = list(subset)
        # A tournament on k vertices is transitive iff it's acyclic
        # iff there's a unique Hamiltonian path
        # Check: for all i < j in the topological order, S[i] -> S[j]
        # Try all orderings... for small k just check acyclicity
        is_trans = True
        for i in range(len(S)):
            for j in range(i+1, len(S)):
                # In a transitive tournament, score sequence is 0,1,...,k-1
                pass
        # Simpler: count 3-cycles
        if k == 3:
            a, b, c = S
            if A[a][b] and A[b][c] and A[c][a]:
                continue
            if A[a][c] and A[c][b] and A[b][a]:
                continue
            count += 1
        elif k == 4:
            # Transitive if no 3-cycle in any triple
            has_3cycle = False
            for triple in combinations(S, 3):
                a, b, c = triple
                if A[a][b] and A[b][c] and A[c][a]:
                    has_3cycle = True; break
                if A[a][c] and A[c][b] and A[b][a]:
                    has_3cycle = True; break
            if not has_3cycle:
                count += 1
    return count

def analyze_tournament(A, n, bits):
    allowed_3 = enumerate_allowed_paths(A, n, 3)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    omega3 = compute_omega_basis(A, n, 3, allowed_3, allowed_2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        return None

    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd2_om = bd2 @ omega2
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    rk2 = int(np.sum(np.abs(S) > 1e-8))
    Z2_dim = dim_O2 - rk2

    if Z2_dim == 0:
        return None

    dt_paths = find_DT_paths(A, n)
    DT_bd = compute_DT_boundary_matrix(dt_paths, allowed_2, n)
    rk_DT = np.linalg.matrix_rank(DT_bd, tol=1e-8)

    Z2_omega = Vt[rk2:].T
    Z2_A2 = omega2 @ Z2_omega
    combined = np.hstack([DT_bd, Z2_A2])
    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
    gap = rk_combined - rk_DT

    # Check full Omega_3 fills
    if dim_O3 > 0:
        bd3 = build_full_boundary_matrix(allowed_3, allowed_2)
        bd3_om = bd3 @ omega3
        rk_full = np.linalg.matrix_rank(bd3_om, tol=1e-8)
        full_combined = np.hstack([bd3_om, Z2_A2])
        rk_full_combined = np.linalg.matrix_rank(full_combined, tol=1e-8)
        beta2 = rk_full_combined - rk_full
    else:
        beta2 = Z2_dim

    scores = tuple(sorted(sum(row) for row in A))
    t3 = count_transitive_tuples(A, n, 3)
    t4 = count_transitive_tuples(A, n, 4)
    c3 = n*(n-1)*(n-2)//6 - t3  # 3-cycles

    return {
        'bits': bits, 'scores': scores,
        'A2': len(allowed_2), 'A3': len(allowed_3),
        'O2': dim_O2, 'O3': dim_O3,
        'Z2': Z2_dim, 'rk_DT': rk_DT,
        'gap': gap, 'DT': len(dt_paths),
        'beta2': beta2,
        't3': t3, 't4': t4, 'c3': c3,
    }


# First, find the gap=6 cases at n=7
print("=" * 70)
print("SEARCHING FOR EXTREME DEFICIT CASES AT n=7")
print("=" * 70)

n = 7
max_gap = 0
extreme_cases = []

random.seed(0)
t0 = time.time()
for trial in range(20000):
    if trial % 5000 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/20000 ({dt:.0f}s), max_gap={max_gap}")

    bits = random_tournament(n)
    A = build_adj(n, bits)
    result = analyze_tournament(A, n, bits)
    if result is None:
        continue

    if result['gap'] > max_gap:
        max_gap = result['gap']

    if result['gap'] >= 3:
        extreme_cases.append(result)

dt = time.time() - t0
print(f"\n  Completed {trial+1} samples in {dt:.0f}s")
print(f"  Max gap found: {max_gap}")
print(f"  Cases with gap >= 3: {len(extreme_cases)}")

for r in extreme_cases[:10]:
    print(f"\n  bits={r['bits']}: scores={r['scores']}")
    print(f"    |A_2|={r['A2']}, |A_3|={r['A3']}, dim(O2)={r['O2']}, dim(O3)={r['O3']}")
    print(f"    Z_2={r['Z2']}, rk(DT)={r['rk_DT']}, gap={r['gap']}, |DT|={r['DT']}")
    print(f"    t3={r['t3']}, t4={r['t4']}, c3={r['c3']}, beta_2={r['beta2']}")

# Now sample n=8
print(f"\n{'='*70}")
print("SAMPLING n=8 (1000 samples)")
print("=" * 70)

n = 8
gap_counts_8 = {}
beta2_fail_8 = 0
max_gap_8 = 0

t0 = time.time()
for trial in range(1000):
    if trial % 250 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/1000 ({dt:.0f}s)")

    bits = random_tournament(n)
    A = build_adj(n, bits)
    result = analyze_tournament(A, n, bits)
    if result is None:
        continue

    gap = result['gap']
    gap_counts_8[gap] = gap_counts_8.get(gap, 0) + 1
    if gap > max_gap_8:
        max_gap_8 = gap
    if result['beta2'] > 0:
        beta2_fail_8 += 1
        print(f"  *** BETA_2 > 0: scores={result['scores']}")

dt = time.time() - t0
print(f"\n  n=8, 1000 samples ({dt:.0f}s)")
print(f"  beta_2 violations: {beta2_fail_8}")
print(f"  Max gap: {max_gap_8}")
print(f"  Gap distribution:")
for gap in sorted(gap_counts_8.keys()):
    total_nontrivial = sum(gap_counts_8.values())
    pct = 100*gap_counts_8[gap]/total_nontrivial
    print(f"    gap={gap}: {gap_counts_8[gap]} ({pct:.1f}%)")

print("\nDone.")
