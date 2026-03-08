#!/usr/bin/env python3
"""
beta2_section_construction.py — Construct explicit section s: Z₂ → Ω₃

If we can find s such that ∂₃ ∘ s = id on Z₂, then β₂=0.
We look for COMBINATORIAL/FORMULAIC descriptions of such sections.

Strategy: For each tournament, find a right inverse of ∂₃|Ω₃ : Ω₃ → Z₂.
Then analyze the structure of this section map across all tournaments.

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

def compute_section_data(A, n):
    """Compute Z₂, Ω₃, and analyze the section structure."""
    # Get allowed paths
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2:
        return {'dim_Z2': 0, 'dim_Om3': 0, 'section_exists': True, 'trivial': True}
    
    # Omega bases
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        return {'dim_Z2': 0, 'dim_Om3': 0, 'section_exists': True, 'trivial': True}
    
    # ∂₂ restricted to Ω₂
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    
    # Rank of ∂₂|Ω₂
    U2, S2, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2  # ker(∂₂|Ω₂)
    
    if dim_Z2 == 0:
        return {'dim_Z2': 0, 'dim_Om3': 0, 'section_exists': True, 'trivial': True}
    
    # Z₂ basis (kernel of ∂₂|Ω₂ in Ω₂ coordinates)
    Z2_basis = Vt2[rk2:].T  # columns are Z₂ basis vectors in Ω₂ coords
    
    if not a3:
        return {'dim_Z2': dim_Z2, 'dim_Om3': 0, 'section_exists': False, 
                'trivial': False, 'failure': 'no_3paths'}
    
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_Om3 = om3.shape[1] if om3.ndim == 2 else 0
    
    if dim_Om3 == 0:
        return {'dim_Z2': dim_Z2, 'dim_Om3': 0, 'section_exists': False,
                'trivial': False, 'failure': 'no_omega3'}
    
    # ∂₃|Ω₃ in Ω₂ coordinates
    bd3 = build_full_boundary_matrix(a3, a2)
    bd3_om = bd3 @ om3  # in A₂ coordinates
    
    # Project ∂₃|Ω₃ into Ω₂ coordinates
    # bd3_om columns are in span(om2), so express them in Ω₂ coords
    bd3_in_om2, residuals, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
    
    # Now bd3_in_om2 : Ω₃ → Ω₂ (in Ω₂ coordinates)
    # We need: for each z ∈ Z₂, find c ∈ Ω₃ with bd3_in_om2 @ c = z
    
    # Restrict to Z₂ target: project bd3_in_om2 onto Z₂
    # bd3_in_Z2 = Z2_basis.T @ bd3_in_om2
    bd3_in_Z2 = Z2_basis.T @ bd3_in_om2  # dim_Z2 × dim_Om3
    
    # Section exists iff bd3_in_Z2 has rank = dim_Z2
    U3, S3, Vt3 = np.linalg.svd(bd3_in_Z2, full_matrices=True)
    rk3_to_Z2 = int(np.sum(np.abs(S3) > 1e-8))
    
    section_exists = (rk3_to_Z2 == dim_Z2)
    
    if not section_exists:
        return {'dim_Z2': dim_Z2, 'dim_Om3': dim_Om3, 'section_exists': False,
                'rk_to_Z2': rk3_to_Z2}
    
    # Construct the section: right inverse of bd3_in_Z2
    # s = bd3_in_Z2^+ (pseudoinverse) maps Z₂ → Ω₃
    section = np.linalg.pinv(bd3_in_Z2)  # dim_Om3 × dim_Z2
    
    # Verify: bd3_in_Z2 @ section ≈ I_{dim_Z2}
    check = bd3_in_Z2 @ section
    err = np.max(np.abs(check - np.eye(dim_Z2)))
    
    # Analyze sparsity of the section
    section_nonzero = np.sum(np.abs(section) > 1e-8)
    section_total = section.size
    
    # What Ω₃ elements are used by the section?
    # Each column of section (in Ω₃ coords) fills one Z₂ basis element
    used_om3_indices = set()
    for col in range(section.shape[1]):
        for row in range(section.shape[0]):
            if abs(section[row, col]) > 1e-8:
                used_om3_indices.add(row)
    
    return {
        'dim_Z2': dim_Z2,
        'dim_Om2': dim_Om2,
        'dim_Om3': dim_Om3,
        'section_exists': section_exists,
        'section_err': err,
        'section_sparsity': section_nonzero / section_total if section_total > 0 else 0,
        'section_rank': np.linalg.matrix_rank(section, tol=1e-8),
        'used_om3_count': len(used_om3_indices),
        'section': section,
        'bd3_in_Z2': bd3_in_Z2,
        'Z2_basis': Z2_basis,
        'om3': om3,
        'a3': a3,
    }

print("=" * 70)
print("SECTION CONSTRUCTION: s: Z₂ → Ω₃ with ∂₃∘s = id")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n = {n} ---")
    
    section_data = []
    failures = 0
    trivials = 0
    
    for A in all_tournaments(n):
        data = compute_section_data(A, n)
        if data.get('trivial', False):
            trivials += 1
            continue
        if not data['section_exists']:
            failures += 1
            print(f"  FAILURE: dim_Z2={data['dim_Z2']}, dim_Om3={data['dim_Om3']}")
            continue
        section_data.append(data)
    
    total = (1 << (n*(n-1)//2))
    print(f"  Total: {total}, Trivial (Z₂=0): {trivials}, Nontrivial: {len(section_data)}, Failures: {failures}")
    
    if not section_data:
        continue
    
    # Analyze section properties
    dims_Z2 = [d['dim_Z2'] for d in section_data]
    sparsities = [d['section_sparsity'] for d in section_data]
    
    print(f"  dim(Z₂) distribution: {Counter(dims_Z2)}")
    print(f"  Section sparsity: mean={np.mean(sparsities):.3f}, min={min(sparsities):.3f}, max={max(sparsities):.3f}")
    
    # Key question: is there a UNIFORM section formula?
    # For each allowed 3-path type, does it always map the same Z₂ type?
    if n == 5:
        print(f"\n  Analyzing section structure at n={n}...")
        
        # Look at which DT paths are used
        dt_usage = Counter()
        non_dt_usage = Counter()
        
        for data in section_data[:200]:  # sample
            A_paths = data['a3']
            om3 = data['om3']
            section = data['section']
            
            # Get the Ω₃ basis elements used
            for col in range(section.shape[1]):
                for row in range(section.shape[0]):
                    if abs(section[row, col]) > 1e-8:
                        # This Ω₃ basis element is used
                        # What allowed 3-paths does it involve?
                        om3_vec = om3[:, row]
                        nonzero_paths = [A_paths[i] for i in range(len(A_paths)) if abs(om3_vec[i]) > 1e-8]
                        for p in nonzero_paths:
                            a, b, c, d = p
                            # Check if DT
                            # Need adjacency... skip for now
                            pass

# Now the KEY experiment: DT-only section
print(f"\n{'='*70}")
print("DT-ONLY SECTION TEST")
print("Can we build s using ONLY DT 4-paths?")
print("=" * 70)

for n in [5]:
    print(f"\n--- n = {n} ---")
    
    dt_only_success = 0
    dt_only_fail = 0
    total_nontrivial = 0
    
    for idx, A in enumerate(all_tournaments(n)):
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
        
        total_nontrivial += 1
        
        # Find DT 4-paths
        dt_paths = []
        for p in a3:
            a_v, b_v, c_v, d_v = p
            if A[a_v][c_v] == 1 and A[b_v][d_v] == 1:
                dt_paths.append(p)
        
        if not dt_paths:
            dt_only_fail += 1
            continue
        
        # Build boundary of DT paths in Ω₂ coordinates
        bd3 = build_full_boundary_matrix(a3, a2)
        
        # DT path indices in a3
        dt_indices = [a3.index(p) for p in dt_paths]
        
        # Boundary of DT paths
        bd3_dt = bd3[:, dt_indices]  # in A₂ coords
        
        # Project into Ω₂ coords
        bd3_dt_om2, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
        
        # Get Z₂ basis
        U2, S2v, Vt2 = np.linalg.svd(bd2 @ om2, full_matrices=True)
        rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
        Z2_basis = Vt2[rk2_v:].T
        
        # Project DT boundaries onto Z₂
        bd3_dt_Z2 = Z2_basis.T @ bd3_dt_om2  # dim_Z2 × |DT|
        
        # Does DT alone span Z₂?
        rk_dt = np.linalg.matrix_rank(bd3_dt_Z2, tol=1e-8)
        
        if rk_dt == dim_Z2:
            dt_only_success += 1
        else:
            dt_only_fail += 1
            if dt_only_fail <= 5:
                print(f"  DT-only FAIL at tournament {idx}: dim_Z2={dim_Z2}, rk(∂₃|DT→Z₂)={rk_dt}, |DT|={len(dt_paths)}")
    
    print(f"  Nontrivial: {total_nontrivial}")
    print(f"  DT-only success: {dt_only_success}")
    print(f"  DT-only fail: {dt_only_fail}")
    print(f"  DT alone fills Z₂: {'YES' if dt_only_fail == 0 else 'NO'}")

# Now test at n=6
print(f"\n--- n = 6 (DT-only test) ---")
n = 6
dt_only_success = 0
dt_only_fail = 0
total_nontrivial = 0
fail_examples = []

for idx, A in enumerate(all_tournaments(n)):
    if idx % 5000 == 0 and idx > 0:
        print(f"  ... {idx}/{1<<15}")
    
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
    
    total_nontrivial += 1
    
    dt_paths = []
    for p in a3:
        a_v, b_v, c_v, d_v = p
        if A[a_v][c_v] == 1 and A[b_v][d_v] == 1:
            dt_paths.append(p)
    
    if not dt_paths:
        dt_only_fail += 1
        continue
    
    bd3 = build_full_boundary_matrix(a3, a2)
    dt_indices = [a3.index(p) for p in dt_paths]
    bd3_dt = bd3[:, dt_indices]
    bd3_dt_om2, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
    
    U2, S2v, Vt2 = np.linalg.svd(bd2 @ om2, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_basis = Vt2[rk2_v:].T
    
    bd3_dt_Z2 = Z2_basis.T @ bd3_dt_om2
    rk_dt = np.linalg.matrix_rank(bd3_dt_Z2, tol=1e-8)
    
    if rk_dt == dim_Z2:
        dt_only_success += 1
    else:
        dt_only_fail += 1
        if len(fail_examples) < 5:
            fail_examples.append({
                'idx': idx, 'dim_Z2': dim_Z2, 'rk_dt': rk_dt,
                'n_dt': len(dt_paths), 'n_a3': len(a3),
                'dim_Om3': om2.shape[1] if om2.ndim == 2 else 0,
            })

print(f"  Nontrivial: {total_nontrivial}")
print(f"  DT-only success: {dt_only_success}")
print(f"  DT-only fail: {dt_only_fail}")
if fail_examples:
    print(f"  Examples of DT-only failure:")
    for ex in fail_examples:
        print(f"    idx={ex['idx']}: dim_Z2={ex['dim_Z2']}, rk(∂DT→Z₂)={ex['rk_dt']}, |DT|={ex['n_dt']}, |A₃|={ex['n_a3']}")

print("\nDone.")
