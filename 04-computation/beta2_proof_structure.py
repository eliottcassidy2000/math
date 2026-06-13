#!/usr/bin/env python3
"""β_2 = 0 proof structure analysis.

Goal: understand the ALGEBRAIC reason why every 2-cycle in a tournament
is a 2-boundary. Focus on explicit cycle-filling constructions.

Key observations to test:
1. For any 2-cycle z ∈ ker(∂_2) ∩ Ω_2, does z decompose into pieces
   each fillable by a SINGLE Ω_3 element?
2. What is the structure of the "bad face" graph — where non-DT 3-paths
   share bad faces and their differences are in Ω_3?
3. Can we always find a spanning tree of the bad-face graph that fills z?

APPROACH: For each tournament at n=5,6, find explicit ker(∂_2) generators,
and for each, construct an explicit filling chain w ∈ Ω_3 with ∂_3(w) = z.
Analyze the STRUCTURE of w.
"""
import numpy as np
from itertools import combinations, permutations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    boundary_coeffs, path_betti_numbers
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def is_transitive_triple(A, a, b, c):
    """Is (a,b,c) a transitive triple: a→b→c and a→c?"""
    return A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1

def is_DT(A, a, b, c, d):
    """Is (a,b,c,d) a doubly-transitive 4-path: a→b→c→d, a→c, b→d?"""
    return (A[a][b] == 1 and A[b][c] == 1 and A[c][d] == 1 and
            A[a][c] == 1 and A[b][d] == 1)

# ===== Part 1: Explicit cycle structure at n=5 =====
print("=" * 70)
print("PART 1: EXPLICIT 2-CYCLE STRUCTURE (n=5)")
print("=" * 70)

n = 5
cycle_types = Counter()
cycle_structure_examples = []
total_cycles = 0

for tournament_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    # Build ∂_2: Ω_2 → A_1
    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2

    # Find kernel of ∂_2|Ω_2
    U, S, Vt = np.linalg.svd(bd2_om)
    rank2 = sum(s > 1e-8 for s in S)
    ker_dim = dim_om2 - rank2

    if ker_dim == 0:
        continue

    # Kernel basis in Ω_2 coordinates
    ker_basis_om = Vt[rank2:].T  # dim_om2 × ker_dim
    # Convert to A_2 coordinates
    ker_basis_a2 = om2 @ ker_basis_om  # |A_2| × ker_dim

    total_cycles += ker_dim

    # Analyze each kernel generator
    for k in range(ker_dim):
        z = ker_basis_a2[:, k]
        z = z / np.max(np.abs(z))  # normalize

        # Which A_2 paths have nonzero coefficient?
        support = [(i, a2_t[i], z[i]) for i in range(len(a2_t)) if abs(z[i]) > 1e-8]

        # Classify: how many are transitive triples?
        n_tt = sum(1 for _, p, _ in support if is_transitive_triple(A, *p))
        n_ntt = len(support) - n_tt

        # What vertices appear?
        verts_used = set()
        for _, p, _ in support:
            verts_used.update(p)

        cycle_types[(len(support), n_tt, n_ntt, len(verts_used))] += 1

        if len(cycle_structure_examples) < 10 and len(support) <= 8:
            cycle_structure_examples.append({
                'tournament': tournament_idx,
                'support_size': len(support),
                'n_tt': n_tt,
                'n_ntt': n_ntt,
                'verts': sorted(verts_used),
                'paths': [(p, round(c, 4)) for _, p, c in support]
            })

print(f"Total 2-cycles (ker generators): {total_cycles}")
print(f"\nCycle types (support_size, n_TT, n_nonTT, n_verts): count")
for k in sorted(cycle_types):
    print(f"  {k}: {cycle_types[k]}")

print(f"\nExamples:")
for ex in cycle_structure_examples[:5]:
    print(f"  T#{ex['tournament']}, verts={ex['verts']}, "
          f"support={ex['support_size']} ({ex['n_tt']} TT + {ex['n_ntt']} nonTT)")
    for p, c in ex['paths']:
        print(f"    {c:+.4f} * {p}")

# ===== Part 2: Explicit filling construction =====
print(f"\n\n{'='*70}")
print("PART 2: FILLING CHAIN CONSTRUCTION (n=5)")
print("=" * 70)

# For each 2-cycle z, find w ∈ Ω_3 with ∂_3(w) = z.
# w is expressed in A_3 coordinates via Ω_3 basis.
# Then analyze the support of w.

n = 5
filling_types = Counter()
filling_examples = []
max_examples = 5

for tournament_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dim_om2 == 0:
        continue

    # ∂_2|Ω_2 kernel
    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    U2, S2, Vt2 = np.linalg.svd(bd2_om)
    rank2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - rank2

    if ker_dim == 0:
        continue

    ker_basis_a2 = om2 @ Vt2[rank2:].T

    # ∂_3|Ω_3 image
    if dim_om3 == 0:
        print(f"  T#{tournament_idx}: ker_dim={ker_dim} but Ω_3=0! β_2>0???")
        continue

    bd3 = build_full_boundary_matrix(a3_t, a2_t)
    bd3_om = bd3 @ om3  # |A_2| × dim_Ω_3

    # For each kernel generator z, find w s.t. bd3_om @ w_om3 = z
    # in A_2 coordinates
    for k in range(ker_dim):
        z = ker_basis_a2[:, k]

        # Solve bd3_om @ w = z (least squares)
        w_om3, res, _, _ = np.linalg.lstsq(bd3_om, z, rcond=None)

        # Check: is residual zero?
        residual = np.max(np.abs(bd3_om @ w_om3 - z))
        if residual > 1e-6:
            print(f"  T#{tournament_idx}: FILLING FAILED, residual={residual}")
            continue

        # Convert w from Ω_3 coordinates to A_3 coordinates
        w_a3 = om3 @ w_om3

        # Analyze support of w
        w_support = [(i, a3_t[i], w_a3[i]) for i in range(len(a3_t))
                     if abs(w_a3[i]) > 1e-8]

        # How many are DT?
        n_dt = sum(1 for _, p, _ in w_support if is_DT(A, *p))
        n_nondt = len(w_support) - n_dt

        filling_types[(len(w_support), n_dt, n_nondt)] += 1

        if len(filling_examples) < max_examples:
            z_support = [(a2_t[i], z[i]) for i in range(len(a2_t)) if abs(z[i]) > 1e-8]
            filling_examples.append({
                'tournament': tournament_idx,
                'cycle_size': len(z_support),
                'filling_size': len(w_support),
                'n_dt': n_dt,
                'n_nondt': n_nondt,
                'cycle': [(p, round(c, 4)) for p, c in z_support],
                'filling': [(p, round(c, 4)) for _, p, c in w_support]
            })

print(f"\nFilling types (support_size, n_DT, n_nonDT): count")
for k in sorted(filling_types):
    print(f"  {k}: {filling_types[k]}")

print(f"\nExamples:")
for ex in filling_examples:
    print(f"\n  T#{ex['tournament']}: cycle({ex['cycle_size']}) → "
          f"filling({ex['filling_size']} = {ex['n_dt']} DT + {ex['n_nondt']} nonDT)")
    print(f"  Cycle z:")
    for p, c in ex['cycle']:
        print(f"    {c:+.4f} * {p}")
    print(f"  Filling w:")
    for p, c in ex['filling']:
        dt = "DT" if is_DT(A, *p) else "  "
        print(f"    {c:+.4f} * {p} [{dt}]")

# ===== Part 3: Source vertex cone =====
print(f"\n\n{'='*70}")
print("PART 3: SOURCE VERTEX CONE FILLING (n=5)")
print("=" * 70)
print("For vertex v with out-degree n-1 (source), can we fill ALL 2-cycles")
print("using 3-chains that 'cone' from v?")

n = 5
source_fill_success = 0
source_fill_fail = 0

for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    U2, S2, Vt2 = np.linalg.svd(bd2_om)
    rank2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    # Find a source vertex (out-deg = n-1)
    source = None
    for v in range(n):
        if sum(A[v]) == n - 1:
            source = v
            break

    if source is None:
        # No source vertex — try highest out-deg
        out_degs = [sum(A[v]) for v in range(n)]
        source = out_degs.index(max(out_degs))

    # Cone from source: for each 2-path (a,b,c) in the cycle,
    # can we prepend source to get (source,a,b,c) ∈ A_3?
    # This requires source→a.
    # Also: (a,source,b,c) requires a→source and source→b.
    # Etc.

    # Actually, let's just check: for each ker generator z,
    # restrict Ω_3 to paths involving source vertex,
    # and see if their images span ker(∂_2).

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    if dim_om3 == 0:
        source_fill_fail += 1
        continue

    bd3 = build_full_boundary_matrix(a3_t, a2_t)
    bd3_om = bd3 @ om3

    # Filter Ω_3 to paths involving source
    source_mask = np.array([source in a3_t[i] for i in range(len(a3_t))])
    # In Ω_3 basis, we need to project
    # Actually, let's use A_3 representation
    # Ω_3 elements in A_3 coords are columns of om3
    # bd3_om_source = bd3 @ (columns of om3 that are supported only on source-paths)
    # This is harder — Ω_3 basis elements may mix source and non-source paths

    # Simpler: among ALL Ω_3 elements (basis columns), compute their
    # images, and check which span ker(∂_2)
    # Then see if we can restrict to source-involving ones

    # Use the full image for now — we already know β_2=0
    # Let's check: what fraction of im(∂_3) is generated by
    # source-related Ω_3 elements?

    # Get DT paths involving source
    dt_source = []
    dt_other = []
    for i, p in enumerate(a3_t):
        if is_DT(A, *p):
            if source in p:
                dt_source.append(i)
            else:
                dt_other.append(i)

    # Image of DT-source paths in A_2
    if dt_source:
        bd3_dt_source = bd3[:, dt_source]
        rank_dt_source = np.linalg.matrix_rank(bd3_dt_source, tol=1e-8)
    else:
        rank_dt_source = 0

    ker_basis_a2 = om2 @ Vt2[rank2:].T

    # Does im(∂_3|DT_source) contain ker(∂_2)?
    # Check by seeing if each ker generator is in the column space of bd3_dt_source
    if dt_source:
        combined = np.hstack([bd3_dt_source, ker_basis_a2])
        rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        fills = (rank_combined == rank_dt_source)
    else:
        fills = (ker_dim == 0)

    if fills:
        source_fill_success += 1
    else:
        source_fill_fail += 1

print(f"  Source-vertex DT fills all 2-cycles: {source_fill_success}")
print(f"  Source-vertex DT insufficient: {source_fill_fail}")

# ===== Part 4: Vertex-pair analysis =====
print(f"\n\n{'='*70}")
print("PART 4: FOR EACH 2-CYCLE, WHICH VERTEX PAIRS PROVIDE DT FILLERS?")
print("=" * 70)

n = 5
all_filled = 0
not_filled = 0

for tournament_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    U2, S2, Vt2 = np.linalg.svd(bd2_om)
    rank2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    ker_basis_a2 = om2 @ Vt2[rank2:].T

    # Get ALL DT paths
    dt_paths = [i for i, p in enumerate(a3_t) if is_DT(A, *p)]

    if not dt_paths:
        not_filled += ker_dim
        continue

    bd3_dt = bd3 = build_full_boundary_matrix(a3_t, a1_t)
    # Wait, I need ∂_3: A_3 → A_2
    bd3 = build_full_boundary_matrix(a3_t, a2_t)

    bd3_dt = bd3[:, dt_paths]

    combined = np.hstack([bd3_dt, ker_basis_a2])
    rank_dt = np.linalg.matrix_rank(bd3_dt, tol=1e-8)
    rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)

    if rank_combined == rank_dt:
        all_filled += ker_dim
    else:
        not_filled += (ker_dim - (rank_combined - rank_dt))
        all_filled += max(0, rank_combined - rank_dt)

print(f"  DT-filled: {all_filled}, NOT DT-filled: {not_filled}")

# ===== Part 5: The key theorem test =====
print(f"\n\n{'='*70}")
print("PART 5: KEY THEOREM TEST")
print("=" * 70)
print("For each non-TT 2-path (a,b,c) with c→a (i.e., 3-cycle a→b→c→a),")
print("find ALL other non-TT 2-paths sharing a non-allowed face with (a,b,c).")
print("The difference of two such paths is in Ω_3 if they share the same bad face.")

n = 5
sharing_stats = Counter()

for A in all_tournaments_gen(n):
    # Find non-TT 2-paths: (a,b,c) with a→b→c but c→a
    non_tt = []
    for p in enumerate_allowed_paths(A, n, 2):
        a, b, c = p
        if not is_transitive_triple(A, a, b, c):
            non_tt.append((a, b, c))

    if not non_tt:
        continue

    # For each non-TT 2-path, find its non-allowed face
    # ∂(a,b,c) = (b,c) - (a,c) + (a,b)
    # Non-allowed face is (a,c) when c→a (since a→c is required for allowed)
    # Actually: face (a,c) is allowed iff a→c. If c→a, then (a,c) is NOT allowed.

    # Group non-TT by their bad face
    bad_face_groups = defaultdict(list)
    for a, b, c in non_tt:
        # Faces: (b,c), (a,c), (a,b)
        # (b,c): allowed iff b→c? YES since it's a 2-path a→b→c
        # (a,b): allowed iff a→b? YES
        # (a,c): allowed iff a→c. NOT allowed iff c→a.
        if A[c][a] == 1:  # c→a, so (a,c) not allowed
            bad_face_groups[('ac', a, c)].append((a, b, c))

    for key, group in bad_face_groups.items():
        sharing_stats[len(group)] += 1

print(f"  Bad-face group sizes (how many non-TT paths share same bad face):")
for k in sorted(sharing_stats):
    print(f"    size {k}: {sharing_stats[k]} groups")

print(f"\nKey: if group size ≥ 2, pairwise differences are in Ω_3")
print(f"(because the non-allowed face cancels).")
print(f"This provides (size-1) independent Ω_3 elements per group.")

print("\nDone.")
