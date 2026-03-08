#!/usr/bin/env python3
"""
beta2_nonDT_omega3.py - Analyze non-DT elements of Omega_3

At n=6, 960/32768 tournaments have DT deficit (DT boundaries don't span Z_2).
For these, Omega_3 contains non-DT elements whose boundaries fill the gap.

What ARE these non-DT Omega_3 elements?

A 3-path (a,b,c,d) is in Omega_3 iff d_3(a,b,c,d) has all components in A_2.
d_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

For (a,b,c,d) to be DT: need a->c AND b->d (shortcuts).
If NOT DT: either c->a (missing shortcut) or d->b (missing shortcut) or both.

Case 1: c->a, b->d: face (a,c,d) has a->c false (c->a), so NOT in A_2.
  But (a,c,d) might still be "cancellable" if another path shares this bad face.

Case 2: a->c, d->b: face (a,b,d) has b->d false (d->b), so NOT in A_2.

Case 3: c->a, d->b: BOTH shortcuts missing.

A non-DT path (a,b,c,d) is in Omega_3 only if it appears in a LINEAR COMBINATION
where the bad faces cancel. That is, if (a,b,c,d) is not in Omega_3 individually,
it might be part of a cancellation pair p1-p2 where the bad face terms cancel.

Let's find explicit examples of these non-DT Omega_3 elements.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
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

def is_DT(A, a, b, c, d):
    """Check if (a,b,c,d) is a DT 4-path."""
    return A[a][c] == 1 and A[b][d] == 1

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"NON-DT OMEGA_3 ELEMENTS AT n={n}")
print("=" * 70)

# Find a deficit tournament and analyze it in detail
deficit_examples = []

for bits in range(total):
    A = build_adj(n, bits)

    allowed_3 = enumerate_allowed_paths(A, n, 3)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega3 = compute_omega_basis(A, n, 3, allowed_3, allowed_2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    # Count DT among allowed_3
    dt_count = sum(1 for p in allowed_3 if is_DT(A, *p))
    non_dt_count = len(allowed_3) - dt_count

    # Check if Omega_3 has non-DT components
    if dim_O3 > dt_count:
        # Omega_3 has non-DT elements
        scores = tuple(sorted(sum(row) for row in A))
        deficit_examples.append((bits, scores, dim_O3, dt_count, non_dt_count))

        if len(deficit_examples) <= 3:
            print(f"\n  bits={bits}, scores={scores}")
            print(f"    |A_3|={len(allowed_3)}, |DT|={dt_count}, dim(O3)={dim_O3}")

            # Identify which basis vectors have non-DT support
            path_to_idx = {p: i for i, p in enumerate(allowed_3)}
            dt_indices = set(i for i, p in enumerate(allowed_3) if is_DT(A, *p))

            for j in range(min(dim_O3, 5)):
                col = omega3[:, j]
                nonzero = [(i, col[i]) for i in range(len(col)) if abs(col[i]) > 1e-8]
                dt_in_col = sum(1 for i, _ in nonzero if i in dt_indices)
                ndt_in_col = len(nonzero) - dt_in_col
                print(f"    Basis {j}: {len(nonzero)} nonzero ({dt_in_col} DT, {ndt_in_col} non-DT)")

                if ndt_in_col > 0:
                    # Show the non-DT paths
                    for i, coeff in nonzero:
                        p = allowed_3[i]
                        dt_str = "DT" if i in dt_indices else "NTD"
                        missing = []
                        if not A[p[0]][p[2]]: missing.append(f"{p[2]}->{p[0]}")
                        if not A[p[1]][p[3]]: missing.append(f"{p[3]}->{p[1]}")
                        miss_str = f" [missing: {', '.join(missing)}]" if missing else ""
                        print(f"      {p}: {coeff:.4f} ({dt_str}){miss_str}")

    if len(deficit_examples) >= 5:
        break

print(f"\n  Total tournaments with dim(O3) > |DT|: {len(deficit_examples)} (of first {bits+1})")

# For the first deficit example, check the EXACT boundary of the non-DT Omega_3 elements
if deficit_examples:
    bits = deficit_examples[0][0]
    A = build_adj(n, bits)
    print(f"\n{'='*70}")
    print(f"DETAILED BOUNDARY ANALYSIS for bits={bits}")
    print(f"{'='*70}")

    allowed_3 = enumerate_allowed_paths(A, n, 3)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega3 = compute_omega_basis(A, n, 3, allowed_3, allowed_2)
    omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    # Boundary d_3: A_3 -> A_2
    bd3 = build_full_boundary_matrix(allowed_3, allowed_2)

    # d_3 restricted to Omega_3
    bd3_om = bd3 @ omega3  # |A_2| x dim_O3

    # Project boundaries into Omega_2 coordinates
    # The boundary of Omega_3 elements should land in Omega_2
    # (since beta_2 = 0, the boundaries fill Z_2)

    print(f"\n  dim(O3)={dim_O3}, dim(O2)={dim_O2}")
    print(f"  Boundary matrix bd3_om shape: {bd3_om.shape}")

    # Check: which allowed 2-paths appear in the boundaries of non-DT Omega_3 elements?
    dt_indices = set(i for i, p in enumerate(allowed_3) if is_DT(A, *p))

    for j in range(dim_O3):
        col = omega3[:, j]
        has_ndt = any(abs(col[i]) > 1e-8 for i in range(len(col)) if i not in dt_indices)

        if has_ndt:
            # This Omega_3 basis vector has non-DT components
            boundary = bd3_om[:, j]
            nonzero_bd = [(i, boundary[i]) for i in range(len(boundary)) if abs(boundary[i]) > 1e-8]
            print(f"\n  Basis {j} (has non-DT): boundary on {len(nonzero_bd)} 2-paths:")
            for i, coeff in nonzero_bd[:10]:
                p = allowed_2[i]
                tt = "TT" if A[p[0]][p[2]] else "NT"
                print(f"    {p}: {coeff:.4f} ({tt})")

print("\nDone.")
