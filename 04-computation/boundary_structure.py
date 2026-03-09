#!/usr/bin/env python3
"""
boundary_structure.py — opus-2026-03-09-S54

KEY STRUCTURAL INSIGHT:

A new 5-path σ = (a₀,...,aᵢ=v,...,a₄) through v has boundary:
  d(σ) = Σⱼ (-1)ʲ d_j(σ)

where d_j removes the j-th vertex. If v is at position i:
  - d_i(σ) = (a₀,...,â_i,...,a₄) — the ONLY old face (avoids v)
  - d_j(σ) for j≠i — new faces (still contain v)

So the "old projection" of d(σ) is just (-1)^i times the single old path.

QUESTION: Can a linear combination of these single-old-path contributions
equal the H_3 generator z of T\v?

If not, then i(z) ∉ im(d_4^T) even at the "old projection" level,
proving i_*-injectivity without needing Omega_4 analysis.

But from Part 5 of saturation_mechanism.py, the rank of the cross-block
was 63 (= all of old Omega_3). So the old faces DO span all of old Omega_3.
This means the constraint must come from the OMEGA QUOTIENT, not the
raw face structure.

Let's investigate: what happens after Omega projection?
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

def build_omega(all_paths, max_p):
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

def main():
    print("=" * 70)
    print("BOUNDARY STRUCTURE: Old-face projection of new d_4")
    print("=" * 70)

    n = 7
    n_arcs = n*(n-1)//2
    n_total = 1 << n_arcs
    tol = 1e-8
    rng = np.random.RandomState(42)

    # Find a beta_3=1 tournament with a BAD vertex
    for trial in range(5000):
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        # Quick beta check
        all_paths = {}
        for p in range(min(n, 7)):
            paths = []
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
            all_paths[p] = paths
        max_p = max(all_paths.keys())
        omega = build_omega(all_paths, max_p)

        # Boundary maps
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
                for k in range(len(sigma)):
                    face = sigma[:k] + sigma[k+1:]
                    if face in idx:
                        mat[idx[face], j] += (-1)**k
            boundary[p] = mat

        # Betti
        dims = {p: omega[p].shape[1] if omega[p].ndim == 2 else 0 for p in range(max_p+1)}
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

        betti = {}
        for p in range(max_p + 1):
            ker = dims[p] - ranks.get(p, 0)
            im_next = ranks.get(p+1, 0)
            betti[p] = ker - im_next

        if betti.get(3, 0) != 1:
            continue

        # Check for BAD vertex
        bad_v = None
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]

            sub_paths = {}
            for p in range(min(n-1, 6)):
                paths = []
                for verts in combinations(range(n-1), p+1):
                    for perm in permutations(verts):
                        if is_allowed_path(T_sub, perm):
                            paths.append(perm)
                sub_paths[p] = paths

            sub_mp = max(sub_paths.keys())
            sub_omega = build_omega(sub_paths, sub_mp)

            sub_boundary = {}
            for p in range(1, sub_mp + 1):
                a_p = sub_paths.get(p, [])
                a_pm1 = sub_paths.get(p-1, [])
                if not a_p or not a_pm1:
                    sub_boundary[p] = np.zeros((0, 0))
                    continue
                idx_s = {path: i for i, path in enumerate(a_pm1)}
                mat = np.zeros((len(a_pm1), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for k in range(len(sigma)):
                        face = sigma[:k] + sigma[k+1:]
                        if face in idx_s:
                            mat[idx_s[face], j] += (-1)**k
                sub_boundary[p] = mat

            sub_dims = {p: sub_omega[p].shape[1] if sub_omega[p].ndim == 2 else 0 for p in range(sub_mp+1)}
            sub_ranks = {}
            for p in range(1, sub_mp + 1):
                Om_p = sub_omega.get(p, np.zeros((0,0)))
                Om_pm1 = sub_omega.get(p-1, np.zeros((0,0)))
                if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
                    sub_ranks[p] = 0
                    continue
                dp = Om_pm1.T @ sub_boundary[p] @ Om_p
                try:
                    S = np.linalg.svd(dp, compute_uv=False)
                    sub_ranks[p] = int(np.sum(S > tol))
                except:
                    sub_ranks[p] = 0
            sub_betti = {}
            for p in range(sub_mp + 1):
                ker = sub_dims[p] - sub_ranks.get(p, 0)
                im_next = sub_ranks.get(p+1, 0)
                sub_betti[p] = ker - im_next

            if sub_betti.get(3, 0) == 1:
                bad_v = v
                break

        if bad_v is None:
            continue

        v = bad_v
        remaining = [i for i in range(n) if i != v]
        print(f"\n  Found: trial={trial}, v={v}")
        print(f"  T betti: {[betti.get(p,0) for p in range(7)]}")
        print(f"  T\\v betti: {[sub_betti.get(p,0) for p in range(6)]}")

        # Compute H_3 generator of T\v
        Om3_sub = sub_omega[3]
        Om2_sub = sub_omega[2]
        d3_sub_proj = Om2_sub.T @ sub_boundary[3] @ Om3_sub
        U3, S3, Vt3 = np.linalg.svd(d3_sub_proj, full_matrices=True)
        r3_sub = int(np.sum(S3 > tol))
        ker_d3_sub = Vt3[r3_sub:].T  # in Omega_3(T\v) coords

        if 4 in sub_boundary and sub_omega.get(4, np.zeros((0,0))).ndim == 2 and sub_omega[4].shape[1] > 0:
            d4_sub_proj = Om3_sub.T @ sub_boundary[4] @ sub_omega[4]
            d4_in_ker = ker_d3_sub.T @ (Om3_sub.T @ sub_boundary[4] @ sub_omega[4])
            if d4_in_ker.shape[1] > 0:
                U4s, S4s, _ = np.linalg.svd(d4_in_ker, full_matrices=True)
                r4_sub = int(np.sum(S4s > tol))
                h3_sub_ker = U4s[:, r4_sub:]
            else:
                h3_sub_ker = np.eye(ker_d3_sub.shape[1])
        else:
            h3_sub_ker = np.eye(ker_d3_sub.shape[1])

        # h3_sub_ker is the H_3 generator in ker(d_3) coordinates
        # Convert to Omega_3(T\v) coordinates
        h3_sub_omega = ker_d3_sub @ h3_sub_ker  # shape (dim_Omega3_sub, 1)

        # Convert to T\v path coordinates
        h3_sub_paths = Om3_sub @ h3_sub_omega  # shape (num_paths_3_sub, 1)

        # Embed into T's path coordinates
        paths3_T = all_paths[3]
        paths3_sub = sub_paths[3]
        path_idx_T = {p: i for i, p in enumerate(paths3_T)}

        embed_paths = np.zeros(len(paths3_T))
        for j in range(len(paths3_sub)):
            if abs(h3_sub_paths[j, 0]) > tol:
                path_T = tuple(remaining[k] for k in paths3_sub[j])
                if path_T in path_idx_T:
                    embed_paths[path_idx_T[path_T]] = h3_sub_paths[j, 0]

        # Project to Omega_3(T)
        Om3_T = omega[3]
        embed_omega = Om3_T.T @ embed_paths

        # Verify it's in ker(d_3^T)
        Om2_T = omega[2]
        d3_T_proj = Om2_T.T @ boundary[3] @ Om3_T
        resid_d3 = d3_T_proj @ embed_omega
        print(f"\n  |d_3^T(embedded z)|: {np.linalg.norm(resid_d3):.2e}")

        # Now: is embed_omega in im(d_4^T)?
        # d_4^T: Omega_4(T) -> Omega_3(T) in Omega coords
        Om4_T = omega.get(4, np.zeros((0,0)))
        if Om4_T.ndim == 2 and Om4_T.shape[1] > 0:
            d4_T_proj = Om3_T.T @ boundary[4] @ Om4_T
            # Is embed_omega in column space of d4_T_proj?
            x, _, _, _ = np.linalg.lstsq(d4_T_proj, embed_omega, rcond=None)
            resid_im4 = np.linalg.norm(embed_omega - d4_T_proj @ x)
            print(f"  |z - im(d_4^T)| in Omega coords: {resid_im4:.2e}")
            if resid_im4 > tol:
                print(f"  => z survives: i_* is INJECTIVE")
            else:
                print(f"  => z killed: i_* NOT injective")

        # DECOMPOSE d_4^T into old + new parts
        paths4_T = all_paths.get(4, [])
        old_p4 = [i for i, p in enumerate(paths4_T) if v not in p]
        new_p4 = [i for i, p in enumerate(paths4_T) if v in p]

        print(f"\n  Omega_4(T) decomposition:")
        print(f"    Total 4-paths: {len(paths4_T)} ({len(old_p4)} old, {len(new_p4)} new)")

        # d_4 in Omega coordinates, split by old/new 4-paths
        # First: project 4-paths to Omega_4(T) coordinates
        # old_p4 and new_p4 are indices in the raw path list
        # The Omega_4 basis is Om4_T: shape (n_paths, dim_omega4)

        # Project d_4 onto old Omega_4 part vs new part
        # d4_T_proj = Om3_T.T @ boundary[4] @ Om4_T
        # This maps all of Omega_4(T) to Omega_3(T).
        # But we want to separate by whether the underlying path uses v.

        # The raw boundary matrix is boundary[4]: shape (n_paths_3, n_paths_4)
        # Restrict to new paths:
        bd4_new = boundary[4][:, new_p4]  # raw boundary of new 4-paths
        bd4_old = boundary[4][:, old_p4]  # raw boundary of old 4-paths

        # Project new 4-paths to Omega_4 space
        # This is subtle: the Omega_4 basis Om4_T mixes old and new paths
        # Let's check: how much of Omega_4 is "new" vs "old"?

        # Indicator: which Omega_4 basis vectors use vertex v?
        om4_new_weight = np.zeros(Om4_T.shape[1])
        for j in range(Om4_T.shape[1]):
            col = Om4_T[:, j]
            new_mass = sum(col[i]**2 for i in new_p4)
            total_mass = sum(col[i]**2 for i in range(len(paths4_T)))
            om4_new_weight[j] = new_mass / max(total_mass, 1e-15)

        print(f"    Omega_4 dim: {Om4_T.shape[1]}")
        print(f"    Omega_4 'new' weight per basis: {[f'{w:.2f}' for w in om4_new_weight]}")

        # Key analysis: decompose d_4^T image in Omega_3(T)
        # Specifically, project embed_omega onto im(d_4^T) and check components

        # Method: find the orthogonal projection of embed_omega onto im(d_4^T)
        U_d4, S_d4, Vt_d4 = np.linalg.svd(d4_T_proj, full_matrices=True)
        r4 = int(np.sum(S_d4 > tol))
        proj_onto_im = U_d4[:, :r4] @ U_d4[:, :r4].T @ embed_omega
        proj_perp = embed_omega - proj_onto_im

        print(f"\n  Projection of embedded z onto im(d_4^T):")
        print(f"    |projection|: {np.linalg.norm(proj_onto_im):.6f}")
        print(f"    |perpendicular|: {np.linalg.norm(proj_perp):.6f}")
        print(f"    |z|: {np.linalg.norm(embed_omega):.6f}")

        # Now compute im(d_4^{T\v}) projected into Omega_3(T)
        # The old d_4 image (from T\v) should NOT hit z (since z is a generator of H_3(T\v))
        # The NEW d_4 image should also not hit z (for i_* injectivity)

        # Old im(d_4): from T\v's Omega_4 to T\v's Omega_3, embedded in T's Omega_3
        Om4_sub = sub_omega.get(4, np.zeros((0,0)))
        if Om4_sub.ndim == 2 and Om4_sub.shape[1] > 0:
            # d_4 of T\v in T\v's Omega coords
            d4_sub_in_omega = Om3_sub.T @ sub_boundary[4] @ Om4_sub
            # Convert each column to T's path coords, then to T's Omega_3
            d4_sub_paths = Om3_sub @ d4_sub_in_omega
            # Embed each column
            n_cols = d4_sub_paths.shape[1]
            d4_sub_in_T_omega = np.zeros((Om3_T.shape[1], n_cols))
            for c in range(n_cols):
                col_embed = np.zeros(len(paths3_T))
                for j in range(len(paths3_sub)):
                    if abs(d4_sub_paths[j, c]) > tol:
                        path_T = tuple(remaining[k] for k in paths3_sub[j])
                        if path_T in path_idx_T:
                            col_embed[path_idx_T[path_T]] = d4_sub_paths[j, c]
                d4_sub_in_T_omega[:, c] = Om3_T.T @ col_embed

            # Check: is z in the span of d4_sub_in_T_omega?
            x_sub, _, _, _ = np.linalg.lstsq(d4_sub_in_T_omega, embed_omega, rcond=None)
            resid_sub = np.linalg.norm(embed_omega - d4_sub_in_T_omega @ x_sub)
            print(f"\n  |z - im(d_4^{{T\\v}} embedded)|: {resid_sub:.2e}")
            print(f"  (should be large — z is NOT a boundary in T\\v)")

        # NOW: What's the NEW im(d_4) that gets added?
        # new_im = im(d_4^T) but NOT in im(d_4^{T\v})
        # For this, compute the combined image

        print(f"\n  Image dimensions:")
        print(f"    rank(d_4^T): {r4}")
        print(f"    rank(d_4^{{T\\v}}): {sub_ranks.get(4, 0)}")

        # Key test: span of old d_4 + new d_4 vs just new d_4
        if Om4_T.shape[1] > 0:
            # Direct computation: project z onto NEW part of im(d_4^T)
            # New part = im(d_4^T) ∩ (complement of im(d_4^{T\v}))
            # Actually, let's just check: does adding new Omega_4 increase
            # the rank of projection onto z's direction?

            # z direction in Omega_3(T)
            z_hat = embed_omega / np.linalg.norm(embed_omega)

            # Correlation of each d_4^T column with z
            correlations = d4_T_proj.T @ z_hat
            print(f"\n  Correlation of d_4^T columns with z:")
            print(f"    max |corr|: {np.max(np.abs(correlations)):.6f}")
            print(f"    mean |corr|: {np.mean(np.abs(correlations)):.6f}")
            print(f"    # columns with |corr| > 0.01: {np.sum(np.abs(correlations) > 0.01)}")

        # FINAL KEY TEST: What is z in the H_3(T) generator coordinates?
        # If i_* is an isomorphism, then z should be (proportional to) the H_3(T) generator.
        # Compute H_3(T) generator
        d3_T_full = Om2_T.T @ boundary[3] @ Om3_T
        U3_T, S3_T, Vt3_T = np.linalg.svd(d3_T_full, full_matrices=True)
        r3_T = int(np.sum(S3_T > tol))
        ker_d3_T = Vt3_T[r3_T:].T

        if Om4_T.shape[1] > 0:
            d4_in_ker_T = ker_d3_T.T @ (Om3_T.T @ boundary[4] @ Om4_T)
            if d4_in_ker_T.shape[1] > 0:
                U4_T, S4_T, _ = np.linalg.svd(d4_in_ker_T, full_matrices=True)
                r4_T = int(np.sum(S4_T > tol))
                h3_T_ker = U4_T[:, r4_T:]  # H_3(T) generator in ker(d_3) coords
            else:
                h3_T_ker = np.eye(ker_d3_T.shape[1])
        else:
            h3_T_ker = np.eye(ker_d3_T.shape[1])

        h3_T_omega = ker_d3_T @ h3_T_ker
        h3_T_hat = h3_T_omega[:, 0] / np.linalg.norm(h3_T_omega[:, 0])

        # Express z in terms of h3_T_hat
        z_proj = np.dot(embed_omega, h3_T_hat)
        z_perp = np.linalg.norm(embed_omega - z_proj * h3_T_hat)

        print(f"\n  H_3(T) generator test:")
        print(f"    |z projected onto H_3(T)|: {abs(z_proj):.6f}")
        print(f"    |z perpendicular to H_3(T)|: {z_perp:.2e}")
        print(f"    |z|: {np.linalg.norm(embed_omega):.6f}")
        if z_perp < tol:
            print(f"    => z IS proportional to H_3(T) generator (ratio: {z_proj / np.linalg.norm(embed_omega):.4f})")
            print(f"    => i_* is an ISOMORPHISM: H_3(T\\v) ≅ H_3(T)")

        break

    # Part 2: Check the face structure of new 5-paths
    print(f"\n{'='*70}")
    print("Part 2: Face structure of new 5-paths")
    print("=" * 70)

    paths4_T = all_paths.get(4, [])
    new_5paths = [paths4_T[i] for i in new_p4]

    # For each new 5-path, identify position of v and the old face
    face_stats = Counter()
    for sigma in new_5paths:
        v_pos = sigma.index(v)
        face_stats[v_pos] += 1

    print(f"  Position of v in new 5-paths:")
    for pos, cnt in sorted(face_stats.items()):
        print(f"    position {pos}: {cnt} paths")

    # Check: for old faces (removing v), what 4-paths appear?
    old_faces_used = Counter()
    for sigma in new_5paths:
        v_pos = sigma.index(v)
        old_face = sigma[:v_pos] + sigma[v_pos+1:]
        old_faces_used[old_face] += 1

    print(f"\n  Old faces (4-paths avoiding v) from new 5-paths:")
    print(f"    {len(old_faces_used)} distinct old faces generated")
    print(f"    Total old 4-paths in T: {len(old_p4)} (as raw paths)")

    # How many of these old faces are in the support of z?
    nonzero_z = set()
    for j in range(len(paths3_T)):
        if abs(embed_paths[j]) > tol:
            nonzero_z.add(paths3_T[j])

    old_faces_in_z = sum(1 for f in old_faces_used if f in nonzero_z)
    print(f"    Old faces in support of z: {old_faces_in_z}/{len(old_faces_used)}")
    print(f"    z support size: {len(nonzero_z)} paths")

if __name__ == '__main__':
    main()
