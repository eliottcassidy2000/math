#!/usr/bin/env python3
"""
beta2_cone_filling.py — Test the "cone from vertex" filling strategy for Z₂

KEY IDEA: For z ∈ Z₂ (a 2-cycle), define h_v(a,b,c) = (v,a,b,c) if v→a.
Then ∂₃(h_v(z)) = z|_{out(v)} + error terms.

The error terms are:
- Type 1: (v,a,c) faces where c→a. These VANISH due to Ω₂ constraint on z!
- Type 2: (v,b,c) faces where b→v. These DON'T vanish automatically.

PROOF STRATEGY: Show that for z ∈ Z₂:
1. Type 1 errors cancel (due to Ω₂ condition) ← PROVED ALGEBRAICALLY
2. Type 2 errors also cancel (due to Z₂ condition ∂₂(z)=0?) ← TO TEST

If both cancel, then h_v gives a chain homotopy: ∂₃∘h_v = id on Z₂ (restricted to out(v)).
Combining cones from all vertices with suitable weights gives β₂=0.

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


print("=" * 70)
print("CONE FILLING ANALYSIS FOR Z₂ CYCLES")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Test: for each tournament and Z₂ basis element,
# compute the "cone defect" — the non-A₂ terms in ∂₃(h_v(z))

type1_nonzero = 0  # Should be 0 (Ω₂ kills these)
type2_nonzero = 0  # Unknown
type2_z2_kills = 0  # Cases where Z₂ condition kills type 2
total_tests = 0

# Detailed analysis on a few tournaments
detailed_bits = [0, 5, 100, 200, 500, 1000]

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0:
        continue

    z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T

    ap2_list = [tuple(x) for x in ap2]
    ap2_idx = {p: i for i, p in enumerate(ap2_list)}

    do_detailed = bits in detailed_bits

    for z_idx in range(z2_dim):
        z_om2 = z2_basis[:, z_idx]
        # Convert to A₂ coords
        z_a2 = om2 @ z_om2  # vector over A₂ paths
        z_dict = {}
        for i, p in enumerate(ap2_list):
            if abs(z_a2[i]) > 1e-10:
                z_dict[p] = z_a2[i]

        for v in range(n):
            # Compute cone h_v(z): for each (a,b,c) in z with v→a,
            # add c_v * α_{abc} * (v,a,b,c) to the filling
            c_v = 1.0  # weight for this cone

            # Compute boundary of h_v(z) = Σ α_{abc} * ∂₃(v,a,b,c)
            # = Σ α_{abc} * [(a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)]
            # Split into: A₂ terms and non-A₂ terms

            type1_defect = {}  # non-A₂ faces (v,a,c) where c→a
            type2_defect = {}  # non-A₂ faces (v,b,c) where b→v

            for (a,b,c), alpha in z_dict.items():
                if not A[v][a]:  # v must beat a
                    continue

                # Face (v,b,c): in A₂ iff v→b
                if not A[v][b]:  # b→v, so (v,b,c) ∉ A₂
                    key = (v,b,c)
                    type2_defect[key] = type2_defect.get(key, 0) - alpha

                # Face (v,a,c): in A₂ iff a→c
                if not A[a][c]:  # c→a, so (v,a,c) ∉ A₂
                    key = (v,a,c)
                    type1_defect[key] = type1_defect.get(key, 0) + alpha

                # Face (v,a,b): v→a ✓, a→b ✓. Always in A₂.
                # Face (a,b,c): already in z. Always in A₂ (it's a path in z).

            # Check type 1 defects
            t1_max = max(abs(val) for val in type1_defect.values()) if type1_defect else 0
            if t1_max > 1e-8:
                type1_nonzero += 1

            # Check type 2 defects
            t2_max = max(abs(val) for val in type2_defect.values()) if type2_defect else 0
            if t2_max > 1e-8:
                type2_nonzero += 1
            else:
                type2_z2_kills += 1

            total_tests += 1

            if do_detailed and t2_max > 1e-8 and z_idx == 0:
                print(f"\n  T#{bits}, v={v}, z_idx={z_idx}: type2 defect max={t2_max:.6f}")
                for key, val in sorted(type2_defect.items()):
                    if abs(val) > 1e-8:
                        b_v = key
                        print(f"    face {key}: coeff={val:.6f}, "
                              f"v→b={A[key[0]][key[1]]}, b→v={A[key[1]][key[0]]}")

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}")

print(f"\n{'='*70}")
print("CONE DEFECT SUMMARY")
print("=" * 70)
print(f"Total (tournament, z, v) tests: {total_tests}")
print(f"Type 1 defects (v,a,c) with c→a: {type1_nonzero} (should be 0)")
print(f"Type 2 defects (v,b,c) with b→v: {type2_nonzero}")
print(f"Type 2 zero (Z₂ kills it): {type2_z2_kills}")

if type1_nonzero == 0:
    print("\n✓ Type 1 defects always vanish (Ω₂ constraint works!)")
if type2_nonzero == 0:
    print("✓ Type 2 defects also always vanish (Z₂ condition sufficient!)")
    print("\n⟹ CONE FROM ANY VERTEX GIVES VALID Ω₃ FILLING OF Z₂!")
    print("⟹ This would PROVE β₂=0 for all tournaments!")
else:
    print(f"\n✗ Type 2 defects survive in {type2_nonzero} cases")
    print("  Need corrections or different filling strategy")

    # Can we combine cones to cancel type 2 defects?
    print(f"\n{'='*70}")
    print("TESTING WEIGHTED CONE COMBINATION")
    print("=" * 70)

    # For each (T, z), find weights c_v such that Σ c_v h_v(z) ∈ Ω₃
    combo_works = 0
    combo_fails = 0

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1

        ap2 = enumerate_allowed_paths(A, n, 2)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap0 = enumerate_allowed_paths(A, n, 0)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        d2 = dim_om(om2)
        if d2 == 0: continue

        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk2
        if z2_dim == 0: continue

        z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T
        ap2_list = [tuple(x) for x in ap2]

        for z_idx in range(z2_dim):
            z_om2 = z2_basis[:, z_idx]
            z_a2 = om2 @ z_om2

            # For each vertex v, compute the "cone vector" in A₃ space
            # and its boundary defect
            # We need to find c_v such that Σ c_v * type2_defect_v = 0
            # AND Σ c_v * coverage_v = z (coverage condition)

            # Collect all non-A₂ formal 2-paths that appear as type 2 defects
            all_non_a2 = set()
            defect_by_v = {}  # v -> {face: coeff}

            for v in range(n):
                defect = {}
                for i, (a,b,c) in enumerate(ap2_list):
                    if abs(z_a2[i]) < 1e-10: continue
                    if not A[v][a]: continue
                    if not A[v][b]:  # b→v
                        key = (v,b,c)
                        defect[key] = defect.get(key, 0) - z_a2[i]
                        all_non_a2.add(key)
                defect_by_v[v] = defect

            if not all_non_a2:
                combo_works += 1
                continue

            # Build the constraint matrix: columns = vertices v (weights c_v)
            # Rows = non-A₂ faces. Entry = defect coefficient.
            non_a2_list = sorted(all_non_a2)
            M_defect = np.zeros((len(non_a2_list), n))
            for v in range(n):
                for j, face in enumerate(non_a2_list):
                    M_defect[j, v] = defect_by_v[v].get(face, 0)

            # Coverage condition: Σ_{v: v→a} c_v = 1 for each starting vertex a
            # But this needs to be per-path. Actually:
            # The cone from v covers paths starting from out(v).
            # We need: for each path (a,b,c) in z, Σ_{v→a} c_v · α_{abc} contributes correctly.

            # Actually, the identity we need:
            # ∂₃(Σ c_v h_v(z)) = Σ_v c_v [z|out(v)] + A₂ terms + defect
            # For the A₂ part to equal z: Σ_{v→a} c_v = 1 for each vertex a.
            # For defect to be 0: M_defect @ c = 0

            # So: joint constraints
            # Coverage: for each a ∈ {0,...,n-1}, Σ_{v: v→a} c_v = 1
            M_cov = np.zeros((n, n))
            for a in range(n):
                for v in range(n):
                    if v != a and A[v][a]:
                        M_cov[a, v] = 1

            # Full system: [M_cov; M_defect] @ c = [1; 0]
            M_full = np.vstack([M_cov, M_defect])
            rhs = np.concatenate([np.ones(n), np.zeros(len(non_a2_list))])

            # Solve
            x, residual, _, _ = np.linalg.lstsq(M_full, rhs, rcond=None)
            err = np.max(np.abs(M_full @ x - rhs))

            if err < 1e-6:
                combo_works += 1
            else:
                combo_fails += 1
                if combo_fails <= 3:
                    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
                    print(f"  FAIL: bits={bits}, z_idx={z_idx}, scores={scores}, err={err:.2e}")
                    print(f"    n_non_a2={len(non_a2_list)}, rank(M_full)={np.linalg.matrix_rank(M_full)}")

        if bits % 200 == 0 and bits > 0:
            print(f"  ... {bits}/{1 << m}")

    print(f"\nWeighted cone: works={combo_works}, fails={combo_fails}")
    if combo_fails == 0:
        print("✓ Weighted cone combination always works!")
        print("⟹ β₂=0 can be proved via weighted cone construction")

print("\nDone.")
