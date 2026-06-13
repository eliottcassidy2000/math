#!/usr/bin/env python3
"""
beta2_dt_deficit_analysis.py — Analyze the DT-only rank deficit at n=6

At n=6, 960 tournaments have rank deficit exactly 1: DT boundaries span 
a codimension-1 subspace of Z₂. What non-DT Ω₃ element fills the gap?

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def is_dt(A, p):
    """Check if 3-path p=(a,b,c,d) is doubly-transitive."""
    a, b, c, d = p
    return A[a][c] == 1 and A[b][d] == 1

def bad_faces(A, p):
    """Return list of bad face indices for 3-path p=(a,b,c,d).
    Face 0: (b,c,d) — always good
    Face 1: (a,c,d) — bad iff c→a (not a→c)
    Face 2: (a,b,d) — bad iff d→b (not b→d)
    Face 3: (a,b,c) — always good
    """
    a, b, c, d = p
    bads = []
    if A[a][c] == 0:  # c→a, face (a,c,d) is bad
        bads.append(1)
    if A[b][d] == 0:  # d→b, face (a,b,d) is bad
        bads.append(2)
    return bads

print("=" * 70)
print("DT-ONLY RANK DEFICIT ANALYSIS AT n=6")
print("=" * 70)

n = 6
deficit_tournaments = []
all_score_seqs = Counter()
deficit_score_seqs = Counter()

for idx, A in enumerate(all_tournaments(n)):
    if idx % 5000 == 0 and idx > 0:
        print(f"  ... {idx}/{1<<15}", flush=True)
    
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2 or not a3:
        continue
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        continue
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2
    
    if dim_Z2 == 0:
        continue
    
    # Score sequence
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    all_score_seqs[scores] += 1
    
    # DT paths
    dt_paths = [p for p in a3 if is_dt(A, p)]
    
    # DT boundary rank in Z₂
    bd3 = build_full_boundary_matrix(a3, a2)
    
    if not dt_paths:
        deficit_tournaments.append((idx, A, dim_Z2, 0, 0))
        deficit_score_seqs[scores] += 1
        continue
    
    dt_indices = [a3.index(p) for p in dt_paths]
    bd3_dt = bd3[:, dt_indices]
    bd3_dt_om2, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
    
    U2, S2v, Vt2 = np.linalg.svd(bd2 @ om2, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_basis = Vt2[rk2_v:].T
    
    bd3_dt_Z2 = Z2_basis.T @ bd3_dt_om2
    rk_dt = np.linalg.matrix_rank(bd3_dt_Z2, tol=1e-8)
    
    if rk_dt < dim_Z2:
        deficit_tournaments.append((idx, A, dim_Z2, rk_dt, len(dt_paths)))
        deficit_score_seqs[scores] += 1

print(f"\nTotal deficit tournaments: {len(deficit_tournaments)}")
print(f"\nScore sequences of deficit tournaments:")
for ss, count in sorted(deficit_score_seqs.items()):
    total_with_ss = all_score_seqs[ss]
    print(f"  {ss}: {count}/{total_with_ss} ({100*count/total_with_ss:.0f}%)")

# Deep analysis of first few deficit tournaments
print(f"\n{'='*70}")
print("DEEP ANALYSIS OF DEFICIT TOURNAMENTS")
print("=" * 70)

for idx, A, dim_Z2, rk_dt, n_dt in deficit_tournaments[:10]:
    print(f"\n--- Tournament {idx}: dim_Z2={dim_Z2}, rk_dt={rk_dt}, |DT|={n_dt} ---")
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    print(f"  Scores: {scores}")
    
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_Om3 = om3.shape[1] if om3.ndim == 2 else 0
    
    bd3 = build_full_boundary_matrix(a3, a2)
    
    # Classify all 3-paths by bad face pattern
    bad_patterns = Counter()
    for p in a3:
        bf = tuple(bad_faces(A, p))
        bad_patterns[bf] += 1
    print(f"  Bad face patterns: {dict(bad_patterns)}")
    print(f"  DT (0 bad): {bad_patterns[()]}, 1-bad: {bad_patterns.get((1,),0)+bad_patterns.get((2,),0)}, 2-bad: {bad_patterns.get((1,2),0)}")
    print(f"  dim(Ω₃)={dim_Om3}, |A₃|={len(a3)}")
    
    # Find the non-DT Ω₃ elements
    dt_paths = [p for p in a3 if is_dt(A, p)]
    non_dt = [p for p in a3 if not is_dt(A, p)]
    
    # Which non-DT paths contribute to Ω₃?
    # Ω₃ = ker(∂₃ restricted to A₂-compatible elements)
    # Actually om3 gives us the basis. Let's see which a3 paths are in each basis vector.
    
    non_dt_in_om3 = 0
    for col in range(om3.shape[1]):
        vec = om3[:, col]
        active_paths = [a3[i] for i in range(len(a3)) if abs(vec[i]) > 1e-8]
        has_non_dt = any(not is_dt(A, p) for p in active_paths)
        if has_non_dt:
            non_dt_in_om3 += 1
    
    print(f"  Ω₃ basis vectors involving non-DT paths: {non_dt_in_om3}/{dim_Om3}")
    
    # The missing Z₂ direction
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    U2, S2v, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_basis = Vt2[rk2_v:].T
    
    dt_indices = [a3.index(p) for p in dt_paths]
    bd3_dt = bd3[:, dt_indices]
    bd3_dt_om2, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
    bd3_dt_Z2 = Z2_basis.T @ bd3_dt_om2
    
    # SVD of DT boundaries in Z₂
    Udt, Sdt, Vtdt = np.linalg.svd(bd3_dt_Z2, full_matrices=True)
    # The missing direction is the last row of Udt (singular value near 0)
    missing_Z2 = Udt[:, -1]  # in Z₂ coordinates
    missing_Om2 = Z2_basis @ missing_Z2  # in Ω₂ coordinates
    missing_A2 = om2 @ missing_Om2  # in A₂ coordinates
    
    # What 2-paths make up the missing Z₂ direction?
    print(f"  Missing Z₂ direction (in A₂):")
    active_2paths = [(a2[i], missing_A2[i]) for i in range(len(a2)) if abs(missing_A2[i]) > 1e-8]
    for path, coeff in sorted(active_2paths, key=lambda x: -abs(x[1]))[:10]:
        print(f"    {path}: {coeff:.4f}")
    
    # Full ∂₃|Ω₃ in Z₂: does it fill?
    bd3_om = bd3 @ om3
    bd3_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
    bd3_Z2 = Z2_basis.T @ bd3_om2
    rk_full = np.linalg.matrix_rank(bd3_Z2, tol=1e-8)
    print(f"  Full Ω₃ → Z₂ rank: {rk_full} (need {dim_Z2})")
    
    # Which Ω₃ basis elements fill the gap?
    for col in range(om3.shape[1]):
        vec_Z2 = bd3_Z2[:, col]
        # Project onto missing direction
        proj = np.dot(vec_Z2, Udt[:, -1])
        if abs(proj) > 1e-8:
            active_paths = [a3[i] for i in range(len(a3)) if abs(om3[i, col]) > 1e-8]
            dt_status = ['DT' if is_dt(A, p) else f'bad={bad_faces(A, p)}' for p in active_paths]
            print(f"  Gap-filling Ω₃ element {col}: proj={proj:.4f}")
            for p, s in zip(active_paths, dt_status):
                print(f"    {p} [{s}]")
            break  # Just show first one

print("\nDone.")
