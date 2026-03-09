#!/usr/bin/env python3
"""
saturation_mechanism.py — opus-2026-03-09-S54

The deep analysis revealed:
  GOOD vertex (b3(T\v)=0): delta_ker3 = delta_im4 + 1
  BAD vertex  (b3(T\v)=1): delta_ker3 = delta_im4     (SATURATION)

This script investigates WHY saturation occurs for bad vertices.

Key question: When we add vertex v to T\v to get T, what is the
structure of the "new" chains (those passing through v)?

APPROACH: Decompose Omega_p(T) = Omega_p(T\v) ⊕ Omega_p^v(T)
where Omega_p^v consists of chains that use vertex v.
Then d_p decomposes as a 2x2 block matrix, and the saturation
condition becomes a statement about the rank of the off-diagonal blocks.

Actually, the relative complex R_p = Omega_p(T)/Omega_p(T\v) should
capture exactly the "new" chains through v. The relative homology
H_p(T, T\v) is computed from R_p with induced boundary.

From the LES:
  H_3(T\v) --i*--> H_3(T) --> H_3(T,T\v) --> H_2(T\v) = 0

So H_3(T) --> H_3(T,T\v) is surjective, and:
  beta_3(T) = rank(i*) + dim H_3(T,T\v)

For GOOD: rank(i*)=0, H_3(T,T\v)=1, beta_3=1. CHECK.
For BAD:  rank(i*)=1, H_3(T,T\v)=0, beta_3=1. CHECK.

The saturation for BAD vertices means: adding v creates new ker(d_3)
and new im(d_4) in EQUAL amounts, so there's NO new relative homology.
The existing H_3 generator of T\v "passes through" to T unchanged.

Let's verify this explicitly and examine the relative complex structure.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def get_all_paths(T, n, max_p=7):
    """Return dict p -> list of allowed p-paths (tuples of p+1 vertices)."""
    all_paths = {}
    for p in range(min(n, max_p+1)):
        paths = []
        for verts in combinations(range(n), p+1):
            for perm in permutations(verts):
                if is_allowed_path(T, perm):
                    paths.append(perm)
        all_paths[p] = paths
    return all_paths

def build_omega_basis(all_paths, max_p):
    """Build Omega_p basis via SVD of non-allowed-face relations."""
    tol = 1e-8
    omega = {}
    for p in range(max_p + 1):
        a_p = all_paths.get(p, [])
        if not a_p:
            omega[p] = np.zeros((0, 0))
            continue
        a_pm1_set = set(all_paths.get(p-1, [])) if p > 0 else set()
        na = {}
        for sigma in a_p:
            for i in range(1, len(sigma)-1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in a_pm1_set:
                    if face not in na:
                        na[face] = len(na)
        if not na:
            omega[p] = np.eye(len(a_p))
        else:
            mat = np.zeros((len(na), len(a_p)))
            for j, sigma in enumerate(a_p):
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if face in na:
                        mat[na[face], j] += (-1)**i
            try:
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                if len(a_p) - rank == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
            except:
                omega[p] = np.eye(len(a_p))
    return omega

def compute_projected_boundary(all_paths, omega, max_p):
    """Compute boundary maps and their Omega-projected ranks."""
    tol = 1e-8
    boundary = {}
    for p in range(1, max_p + 1):
        a_p = all_paths.get(p, [])
        a_pm1 = all_paths.get(p-1, [])
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary[p] = mat

    ranks = {}
    for p in range(1, max_p + 1):
        Om_p = omega.get(p, np.zeros((0,0)))
        Om_pm1 = omega.get(p-1, np.zeros((0,0)))
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        try:
            S = np.linalg.svd(dp, compute_uv=False)
            ranks[p] = int(np.sum(S > tol))
        except:
            ranks[p] = 0
    return boundary, ranks

def compute_full_data(T, n):
    """Compute dims, ranks, betti for tournament T on n vertices."""
    all_paths = get_all_paths(T, n)
    max_p = max(all_paths.keys())
    omega = build_omega_basis(all_paths, max_p)
    boundary, ranks = compute_projected_boundary(all_paths, omega, max_p)

    dims = {}
    for p in range(max_p + 1):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0

    betti = {}
    for p in range(max_p + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti[p] = ker - im_next

    return {
        'dims': dims, 'ranks': ranks, 'betti': betti,
        'omega': omega, 'boundary': boundary, 'all_paths': all_paths
    }

def compute_h3_generator(data, n):
    """Compute the H_3 generator as a vector in Omega_3.

    H_3 = ker(d_3) / im(d_4).
    Find a vector in ker(d_3) that is NOT in im(d_4).
    """
    tol = 1e-8
    omega = data['omega']
    boundary = data['boundary']

    Om3 = omega.get(3, np.zeros((0,0)))
    Om2 = omega.get(2, np.zeros((0,0)))
    Om4 = omega.get(4, np.zeros((0,0)))

    if Om3.ndim < 2 or Om3.shape[1] == 0:
        return None
    if Om2.ndim < 2 or Om2.shape[1] == 0:
        return None

    # d_3 in Omega basis
    d3 = Om2.T @ boundary[3] @ Om3

    # ker(d_3)
    U, S, Vt = np.linalg.svd(d3, full_matrices=True)
    rank3 = int(np.sum(S > tol))
    ker_d3 = Vt[rank3:].T  # columns span ker(d_3) in Omega_3 coords

    if ker_d3.shape[1] == 0:
        return None

    # im(d_4) in Omega_3 coords
    if Om4.ndim == 2 and Om4.shape[1] > 0 and 4 in boundary:
        d4 = Om3.T @ boundary[4] @ Om4
        # im(d_4) is column space of d4 projected to ker(d_3)
        # Project d4 columns onto ker(d_3)
        d4_in_ker = ker_d3.T @ (Om3.T @ boundary[4] @ Om4)
    else:
        d4_in_ker = np.zeros((ker_d3.shape[1], 0))

    if d4_in_ker.shape[1] > 0:
        U4, S4, Vt4 = np.linalg.svd(d4_in_ker, full_matrices=True)
        rank4 = int(np.sum(S4 > tol))
        # H_3 = ker(d_3) / im(d_4), represented in ker(d_3) coords
        # The cokernel directions of d4_in_ker
        if rank4 < U4.shape[1]:
            h3_in_ker = U4[:, rank4:]  # columns in ker(d_3) coords
        else:
            return None
    else:
        h3_in_ker = np.eye(ker_d3.shape[1])

    # Convert back to Omega_3 coords
    h3_in_omega = ker_d3 @ h3_in_ker

    # Convert to path coords
    h3_in_paths = Om3 @ h3_in_omega

    return h3_in_paths

def main():
    print("=" * 70)
    print("SATURATION MECHANISM: Why delta_ker3 = delta_im4 for bad vertices")
    print("=" * 70)

    # Part 1: At n=6, compute H_3 generators and check embedding
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    print(f"\nPart 1: H_3 generators at n=6 (Type B)")

    # Find a Type B tournament
    type_b = None
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
        if scores != (2,2,2,3,3,3):
            continue
        data = compute_full_data(T, n)
        if data['betti'].get(3, 0) == 1:
            type_b = (bits, T, data)
            break

    if type_b is None:
        print("  No Type B found!")
        return

    bits, T, data_T = type_b
    print(f"  Found Type B: bits={bits}")
    print(f"  Omega dims: {[data_T['dims'].get(p,0) for p in range(6)]}")
    print(f"  Betti: {[data_T['betti'].get(p,0) for p in range(6)]}")

    # Compute H_3 generator of T
    h3_T = compute_h3_generator(data_T, n)
    if h3_T is not None:
        print(f"  H_3 generator shape: {h3_T.shape}")
        # Show which paths contribute
        paths3 = data_T['all_paths'][3]
        nonzero = [(i, paths3[i], h3_T[i,0]) for i in range(len(paths3)) if abs(h3_T[i,0]) > 1e-8]
        print(f"  H_3 generator has {len(nonzero)} nonzero path coefficients")
        for i, path, coeff in nonzero[:10]:
            print(f"    path {path}: {coeff:.4f}")
        if len(nonzero) > 10:
            print(f"    ... ({len(nonzero)} total)")

    # Part 2: For each vertex v, compute H_3 generator of T\v (if exists)
    # and check how it embeds into T's chain complex
    print(f"\n{'='*70}")
    print("Part 2: Embedding H_3(T\\v) into T's chain complex")
    print("=" * 70)

    # At n=6, ALL deletions of Type B have beta_3=0, so there's nothing to embed.
    # We need n=7 where some deletions have beta_3=1.

    # Let's move to n=7
    print(f"\n  Moving to n=7 for bad-vertex analysis...")
    n7 = 7
    rng = np.random.RandomState(42)
    n7_arcs = n7*(n7-1)//2
    n7_total = 1 << n7_arcs

    found_bad = False
    for trial in range(5000):
        if trial % 1000 == 0:
            print(f"    ... trial {trial}", flush=True)
        bits = rng.randint(0, n7_total)
        T = tournament_from_bits(n7, bits)
        try:
            data_T = compute_full_data(T, n7)
        except:
            continue
        if data_T['betti'].get(3, 0) != 1:
            continue

        # Check each vertex
        for v in range(n7):
            remaining = [i for i in range(n7) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
            try:
                data_sub = compute_full_data(T_sub, n7-1)
            except:
                continue

            if data_sub['betti'].get(3, 0) == 1:
                found_bad = True
                print(f"\n  BAD vertex found: trial={trial}, v={v}")
                score_v = sum(T[v][j] for j in range(n7) if j != v)
                print(f"  score(v)={score_v}")
                print(f"  T betti:   {[data_T['betti'].get(p,0) for p in range(7)]}")
                print(f"  T\\v betti: {[data_sub['betti'].get(p,0) for p in range(6)]}")

                # Compute H_3 generator of T\v
                h3_sub = compute_h3_generator(data_sub, n7-1)
                if h3_sub is not None:
                    paths3_sub = data_sub['all_paths'][3]
                    nonzero_sub = [(i, paths3_sub[i], h3_sub[i,0])
                                   for i in range(len(paths3_sub)) if abs(h3_sub[i,0]) > 1e-8]
                    print(f"  H_3(T\\v) generator: {len(nonzero_sub)} nonzero paths")

                    # Map these paths into T's path numbering
                    # T\v has vertices 0..n7-2 corresponding to remaining[0..n7-2] in T
                    paths3_T = data_T['all_paths'][3]
                    path_idx_T = {p: i for i, p in enumerate(paths3_T)}

                    # Embed h3_sub into T's path space
                    embed = np.zeros(len(paths3_T))
                    mapped = 0
                    for i, path_sub, coeff in nonzero_sub:
                        # Map path from T\v vertex numbering to T vertex numbering
                        path_T = tuple(remaining[k] for k in path_sub)
                        if path_T in path_idx_T:
                            embed[path_idx_T[path_T]] = coeff
                            mapped += 1

                    print(f"  Mapped {mapped}/{len(nonzero_sub)} paths into T")

                    # Check: is embed in ker(d_3^T)?
                    Om3_T = data_T['omega'][3]
                    Om2_T = data_T['omega'][2]
                    if Om3_T.ndim == 2 and Om3_T.shape[1] > 0:
                        embed_omega = Om3_T.T @ embed  # project to Omega_3(T)
                        d3_T = Om2_T.T @ data_T['boundary'][3] @ Om3_T
                        residual_d3 = d3_T @ embed_omega
                        print(f"  |d_3^T(embedded)|: {np.linalg.norm(residual_d3):.2e}")

                        # Check: is embed in im(d_4^T)?
                        Om4_T = data_T['omega'].get(4, np.zeros((0,0)))
                        if Om4_T.ndim == 2 and Om4_T.shape[1] > 0:
                            d4_T = Om3_T.T @ data_T['boundary'][4] @ Om4_T
                            # Check if embed_omega is in column space of d4_T
                            # Solve d4_T @ x = embed_omega
                            try:
                                x, res, _, _ = np.linalg.lstsq(d4_T, embed_omega, rcond=None)
                                residual_im4 = embed_omega - d4_T @ x
                                norm_resid = np.linalg.norm(residual_im4)
                                print(f"  |embed - im(d_4^T)|: {norm_resid:.2e}")
                                if norm_resid < 1e-6:
                                    print(f"  => EMBEDDED H_3 IS A BOUNDARY (killed by d_4)")
                                else:
                                    print(f"  => EMBEDDED H_3 SURVIVES (not in im(d_4))")
                            except:
                                print(f"  lstsq failed")
                        else:
                            print(f"  No Omega_4 in T, embed trivially survives")

                # Also check the RELATIVE complex
                print(f"\n  Relative complex dimensions:")
                for p in range(7):
                    d_T = data_T['dims'].get(p, 0)
                    d_sub = data_sub['dims'].get(p, 0)
                    print(f"    p={p}: dim(T)={d_T}, dim(T\\v)={d_sub}, relative~={d_T - d_sub}")

                delta_ker3 = (data_T['dims'].get(3,0) - data_T['ranks'].get(3,0)) - \
                             (data_sub['dims'].get(3,0) - data_sub['ranks'].get(3,0))
                delta_im4 = data_T['ranks'].get(4,0) - data_sub['ranks'].get(4,0)
                print(f"  delta_ker3={delta_ker3}, delta_im4={delta_im4}, diff={delta_ker3 - delta_im4}")

                break  # one bad vertex per tournament

        if found_bad:
            break

    if not found_bad:
        print("  No bad vertex found in 5000 trials!")
        return

    # Part 3: Check the relative Euler characteristic decomposition
    print(f"\n{'='*70}")
    print("Part 3: Relative complex for ALL vertex types at n=7")
    print("=" * 70)

    # Collect more statistics
    good_stats = defaultdict(list)
    bad_stats = defaultdict(list)
    count = 0

    rng2 = np.random.RandomState(123)
    for trial in range(10000):
        if trial % 2000 == 0:
            print(f"    ... trial {trial}", flush=True)
        bits = rng2.randint(0, n7_total)
        T = tournament_from_bits(n7, bits)
        try:
            data_T = compute_full_data(T, n7)
        except:
            continue
        if data_T['betti'].get(3, 0) != 1:
            continue
        count += 1

        for v in range(n7):
            remaining = [i for i in range(n7) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
            try:
                data_sub = compute_full_data(T_sub, n7-1)
            except:
                continue

            b3_sub = data_sub['betti'].get(3, 0)
            key = "BAD" if b3_sub == 1 else "GOOD"

            # Collect per-level deltas
            deltas = []
            for p in range(7):
                d_T = data_T['dims'].get(p, 0)
                d_sub = data_sub['dims'].get(p, 0)
                r_T = data_T['ranks'].get(p, 0) if p > 0 else 0
                r_sub = data_sub['ranks'].get(p, 0) if p > 0 else 0
                deltas.append((d_T - d_sub, r_T - r_sub))

            if key == "BAD":
                bad_stats['deltas'].append(deltas)
            else:
                good_stats['deltas'].append(deltas)

        if count >= 30:
            break

    print(f"\n  Collected data from {count} beta_3=1 tournaments")
    print(f"  GOOD vertices: {len(good_stats['deltas'])}")
    print(f"  BAD vertices: {len(bad_stats['deltas'])}")

    # Show delta distributions per level
    for name, stats in [("GOOD", good_stats), ("BAD", bad_stats)]:
        print(f"\n  {name} vertices — delta distributions:")
        print(f"    p: (delta_dim, delta_rank) distribution")
        for p in range(7):
            dd = Counter()
            for deltas in stats['deltas']:
                dd[deltas[p]] += 1
            top5 = dd.most_common(5)
            print(f"    p={p}: {dict(top5)}")

    # Part 4: Key algebraic observation
    print(f"\n{'='*70}")
    print("Part 4: Algebraic interpretation")
    print("=" * 70)

    # For BAD vertices: delta_ker3 = delta_im4
    # This means: H_3(T,T\v) = 0
    # From LES: 0 = H_2(T\v) -> H_3(T\v) -i*-> H_3(T) -> H_3(T,T\v)=0 -> 0
    # So i_* is SURJECTIVE from H_3(T\v) to H_3(T)!
    # Since both are 1-dimensional, i_* is an ISOMORPHISM.

    # For GOOD vertices: delta_ker3 = delta_im4 + 1
    # This means: H_3(T,T\v) = 1
    # From LES: 0 -> H_3(T) -> H_3(T,T\v)=F -> 0
    # So H_3(T) ≅ H_3(T,T\v) = F.

    print("""
  SATURATION THEOREM (computed):

  For beta_3(T) = 1 tournaments and any vertex v:

  GOOD vertex (beta_3(T\\v) = 0):
    delta_ker3 = delta_im4 + 1
    H_3(T, T\\v) ≅ F  (1-dimensional)
    H_3(T) ≅ H_3(T, T\\v)  (isomorphism from LES)

  BAD vertex (beta_3(T\\v) = 1):
    delta_ker3 = delta_im4  (SATURATION)
    H_3(T, T\\v) = 0
    i_*: H_3(T\\v) -> H_3(T) is an ISOMORPHISM

  The saturation condition means: adding vertex v creates new
  kernel in d_3 AND new image in d_4 in PERFECTLY EQUAL amounts.
  No new homology is created — the H_3 class passes through.

  This is equivalent to: the H_3 generator of T\\v, when embedded
  into T's chain complex, is NOT killed by the new d_4 boundaries.

  ALGEBRAIC INTERPRETATION:
  The new Omega_4 chains through v have their d_4-images in the
  "new" part of Omega_3 (chains through v), NOT in the "old" part
  (chains avoiding v). The boundary of a path through v always
  involves at least one face through v (specifically, the face
  obtained by removing v itself from the path).

  Wait — this is exactly the RELATIVE COMPLEX structure!
  The boundary map d_4: Omega_4 -> Omega_3 restricted to "new"
  chains maps INTO the "new" Omega_3 space (plus the old space).
  But the crucial point is that im(d_4) on new chains does NOT
  hit the old H_3 generator.
""")

    # Part 5: Check if d_4(new Omega_4) intersects old Omega_3
    print(f"{'='*70}")
    print("Part 5: Does d_4(new chains) hit old Omega_3?")
    print("=" * 70)

    # Go back to the bad-vertex case from Part 2
    rng3 = np.random.RandomState(42)
    for trial in range(5000):
        bits = rng3.randint(0, n7_total)
        T = tournament_from_bits(n7, bits)
        try:
            data_T = compute_full_data(T, n7)
        except:
            continue
        if data_T['betti'].get(3, 0) != 1:
            continue

        for v in range(n7):
            remaining = [i for i in range(n7) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
            try:
                data_sub = compute_full_data(T_sub, n7-1)
            except:
                continue

            if data_sub['betti'].get(3, 0) != 1:
                continue

            # BAD vertex found!
            # Identify "old" vs "new" paths in T's Omega_4
            paths4_T = data_T['all_paths'].get(4, [])
            paths3_T = data_T['all_paths'].get(3, [])

            old_p4 = [i for i, p in enumerate(paths4_T) if v not in p]
            new_p4 = [i for i, p in enumerate(paths4_T) if v in p]
            old_p3 = [i for i, p in enumerate(paths3_T) if v not in p]
            new_p3 = [i for i, p in enumerate(paths3_T) if v in p]

            print(f"\n  BAD vertex v={v} (trial={trial}):")
            print(f"  Omega_4: {len(paths4_T)} paths ({len(old_p4)} old, {len(new_p4)} new)")
            print(f"  Omega_3: {len(paths3_T)} paths ({len(old_p3)} old, {len(new_p3)} new)")

            # d_4 boundary matrix (in raw path coords)
            bd4 = data_T['boundary'].get(4, np.zeros((0,0)))
            if bd4.shape[0] == 0:
                print(f"  No boundary d_4")
                continue

            # d_4 applied to NEW Omega_4 paths -> project onto OLD Omega_3 paths
            if new_p4 and old_p3:
                d4_new_to_old = bd4[np.ix_(old_p3, new_p4)]
                rank_cross = np.linalg.matrix_rank(d4_new_to_old, tol=1e-8)
                print(f"  d_4(new Omega_4 -> old Omega_3) block: {d4_new_to_old.shape}, rank={rank_cross}")

                # This rank tells us: do the NEW 5-paths have boundaries
                # that overlap with the OLD 4-paths?
                if rank_cross > 0:
                    print(f"  => YES, d_4(new) has rank-{rank_cross} overlap with old Omega_3")
                else:
                    print(f"  => NO, d_4(new) is entirely in new Omega_3")

            # Now check specifically against the H_3 generator of T\v
            h3_sub = compute_h3_generator(data_sub, n7-1)
            if h3_sub is not None:
                paths3_sub = data_sub['all_paths'][3]
                path_idx_T = {p: i for i, p in enumerate(paths3_T)}

                # Embed h3_sub into T
                embed = np.zeros(len(paths3_T))
                for i in range(len(paths3_sub)):
                    if abs(h3_sub[i,0]) > 1e-8:
                        path_T = tuple(remaining[k] for k in paths3_sub[i])
                        if path_T in path_idx_T:
                            embed[path_idx_T[path_T]] = h3_sub[i,0]

                # Check: is embed in column space of d4_new (the new part of d_4)?
                d4_new = bd4[:, new_p4] if new_p4 else np.zeros((bd4.shape[0], 0))
                if d4_new.shape[1] > 0:
                    # Project embed: does d4_new @ x ≈ embed?
                    x, res, _, _ = np.linalg.lstsq(d4_new, embed, rcond=None)
                    resid = np.linalg.norm(embed - d4_new @ x)
                    print(f"  |H_3(T\\v) - im(d_4^new)|: {resid:.2e}")

                # Also check against FULL d_4
                d4_full = bd4
                if d4_full.shape[1] > 0:
                    x2, res2, _, _ = np.linalg.lstsq(d4_full, embed, rcond=None)
                    resid2 = np.linalg.norm(embed - d4_full @ x2)
                    print(f"  |H_3(T\\v) - im(d_4^full)|: {resid2:.2e}")
                    if resid2 < 1e-6:
                        print(f"  WARNING: H_3(T\\v) IS in im(d_4^T)!")
                    else:
                        print(f"  CONFIRMED: H_3(T\\v) survives in H_3(T)")

            break  # one example
        break  # one tournament

if __name__ == '__main__':
    main()
