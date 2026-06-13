#!/usr/bin/env python3
"""
beta2_arcflip_chainmap.py - Explicit chain-level analysis of arc-flip invariance

For each arc flip u->v to v->u, explicitly compute:
1. How Omega_2 basis changes (which basis vectors are lost/gained)
2. How the boundary d_2 restricted to Omega_2 changes
3. How d_3 restricted to Omega_3 changes
4. Whether there's a natural chain map Omega_*(T) -> Omega_*(T') preserving exactness

Key insight to explore: the flip changes the "junk face" structure.
- In T: (u,v) is an arc, so (u,v) can be a face of TT triples involving u->v
- In T': (v,u) is an arc instead, so (v,u) is now the face

Can we construct an explicit isomorphism Z_2(T) -> Z_2(T') that also maps B_2(T) -> B_2(T')?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from itertools import combinations
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

def flip_arc(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits

def get_full_chain_data(A, n, max_p=4):
    """Get full chain complex data including explicit bases."""
    paths = {}
    omega = {}
    for p in range(max_p+1):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(max_p+1)}

    bds = {}
    bd_om = {}
    rks = {}

    for p in range(1, max_p+1):
        if len(paths.get(p, [])) > 0 and len(paths.get(p-1, [])) > 0:
            bds[p] = build_full_boundary_matrix(paths[p], paths[p-1])
        else:
            bds[p] = None

    for p in range(2, max_p+1):
        if dims[p] > 0 and bds.get(p) is not None:
            raw = bds[p] @ omega[p]
            if dims[p-1] > 0:
                coords, _, _, _ = np.linalg.lstsq(omega[p-1], raw, rcond=None)
                bd_om[p] = coords
                rks[p] = np.linalg.matrix_rank(coords, tol=1e-8)
            else:
                bd_om[p] = raw
                rks[p] = np.linalg.matrix_rank(raw, tol=1e-8)
        else:
            bd_om[p] = None
            rks[p] = 0

    # d_1 rank
    if bds.get(1) is not None and dims[1] > 0:
        raw1 = bds[1] @ omega[1]
        rks[1] = np.linalg.matrix_rank(raw1, tol=1e-8)
    else:
        rks[1] = 0

    betas = {}
    for p in range(max_p):
        z_p = dims[p] - rks.get(p, 0)
        b_p = z_p - rks.get(p+1, 0)
        betas[p] = b_p

    return paths, omega, dims, bds, bd_om, rks, betas


def classify_omega2_basis(omega2, paths2, A, n):
    """Classify each Omega_2 basis vector as TT-type or NT-cancellation-type."""
    if omega2.shape[1] == 0:
        return [], []

    tt_indices = []
    nt_indices = []

    for col in range(omega2.shape[1]):
        vec = omega2[:, col]
        support = np.where(np.abs(vec) > 1e-10)[0]

        # Check if single TT triple
        if len(support) == 1:
            idx = support[0]
            a, b, c = paths2[idx]
            if A[a][b] and A[b][c] and A[a][c]:
                tt_indices.append(col)
                continue

        # Check if all support paths are TT
        all_tt = True
        for idx in support:
            a, b, c = paths2[idx]
            if not (A[a][b] and A[b][c] and A[a][c]):
                all_tt = False
                break

        if all_tt:
            tt_indices.append(col)
        else:
            nt_indices.append(col)

    return tt_indices, nt_indices


print("=" * 70)
print("CHAIN-LEVEL ANALYSIS OF ARC-FLIP INVARIANCE")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

# Find a rich example: flip where both dOm2 and db1 are nonzero
print(f"\nSearching for instructive examples at n={n}...")

# Precompute all chain data
print("Precomputing all tournaments...")
t0 = time.time()
all_chain = {}
for bits in range(total):
    A = build_adj(n, bits)
    all_chain[bits] = get_full_chain_data(A, n)
dt = time.time() - t0
print(f"  Done in {dt:.1f}s")

# ANALYSIS 1: Track exactly which Omega_2 basis vectors are affected
print(f"\n{'='*70}")
print("ANALYSIS 1: How does Omega_2 relate between T and T'?")
print("=" * 70)

# For a specific flip, look at the overlap between Omega_2(T) and Omega_2(T')
# They live in different ambient spaces (different allowed paths) but we can
# compare by embedding in the space of ALL 2-paths

def get_omega2_in_full_space(omega2, paths2, n):
    """Embed Omega_2 basis in the full C(n,3)*P(3,3) space of ordered triples."""
    # Index all ordered triples (a,b,c) with distinct a,b,c
    all_triples = [(a,b,c) for a in range(n) for b in range(n) for c in range(n)
                   if a != b and b != c and a != c]
    triple_idx = {t: i for i, t in enumerate(all_triples)}

    if omega2.shape[1] == 0:
        return np.zeros((len(all_triples), 0)), all_triples

    full = np.zeros((len(all_triples), omega2.shape[1]))
    for col in range(omega2.shape[1]):
        for row in range(omega2.shape[0]):
            if abs(omega2[row, col]) > 1e-10:
                t = paths2[row]
                full[triple_idx[t], col] = omega2[row, col]

    return full, all_triples

# Detailed analysis of one example per delta type
examples_found = {}
for bits in range(total):
    A = build_adj(n, bits)
    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

            dOm2 = d2[2] - d1[2]
            db1 = b2[1] - b1[1]
            key = (dOm2, db1)

            if key not in examples_found:
                if A[i][j] == 1:
                    u, v = i, j
                else:
                    u, v = j, i
                examples_found[key] = (bits, bits2, u, v)

print(f"\nFound {len(examples_found)} delta types. Analyzing subspace overlap:")

for key in sorted(examples_found.keys()):
    dOm2, db1 = key
    bits, bits2, u, v = examples_found[key]
    A = build_adj(n, bits)
    A2 = build_adj(n, bits2)

    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]
    p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

    # Embed both Omega_2 in full triple space
    full1, triples = get_omega2_in_full_space(o1[2], p1[2], n)
    full2, _ = get_omega2_in_full_space(o2[2], p2[2], n)

    # Compute subspace intersection dimension
    if full1.shape[1] > 0 and full2.shape[1] > 0:
        combined = np.hstack([full1, full2])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        rk1 = np.linalg.matrix_rank(full1, tol=1e-8)
        rk2 = np.linalg.matrix_rank(full2, tol=1e-8)
        intersection_dim = rk1 + rk2 - rk_combined
    else:
        intersection_dim = 0
        rk1 = full1.shape[1]
        rk2 = full2.shape[1]

    drk3 = r2[3] - r1[3]
    print(f"\n  (dOm2={dOm2:+d}, db1={db1:+d}): flip {u}->{v}")
    print(f"    dim(Om2): {d1[2]} -> {d2[2]}")
    print(f"    Om2 intersection dim: {intersection_dim} (of {rk1} and {rk2})")
    print(f"    Om2 exclusive to T: {rk1 - intersection_dim}")
    print(f"    Om2 exclusive to T': {rk2 - intersection_dim}")
    print(f"    rk(d3): {r1[3]} -> {r2[3]} (d={drk3:+d})")
    print(f"    beta: {[b1[p] for p in range(4)]} -> {[b2[p] for p in range(4)]}")


# ANALYSIS 2: The exact boundary map change
print(f"\n{'='*70}")
print("ANALYSIS 2: Structure of d_3 change under arc flip")
print("=" * 70)
print("Looking at how im(d_3) changes in Z_2")

# For a specific non-trivial example, compute Z_2(T) and im(d_3) in T
# Then Z_2(T') and im(d_3) in T'
# Map both to the common ambient space and compare

key = (0, 1)  # dOm2=0, db1=+1 is the most interesting
if key in examples_found:
    bits, bits2, u, v = examples_found[key]
    A = build_adj(n, bits)
    A2 = build_adj(n, bits2)

    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]
    p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

    print(f"\nExample: bits={bits}, flip {u}->{v} to {v}->{u}")
    print(f"  beta(T): {[b1[p] for p in range(4)]}")
    print(f"  beta(T'): {[b2[p] for p in range(4)]}")

    # Z_2(T) = ker(d_2) in Omega_2(T)
    # d_2: Omega_2 -> Omega_1
    # In Omega coords: bdo1[2] maps Omega_2 coords to Omega_1 coords
    if bdo1.get(2) is not None and bdo1[2].shape[1] > 0:
        U, S, Vt = np.linalg.svd(bdo1[2], full_matrices=True)
        rk = int(np.sum(np.abs(S) > 1e-8))
        Z2_coords_T = Vt[rk:].T  # Omega_2 coords of Z_2(T)
        # Map to full triple space
        Z2_full_T = o1[2][:, :] @ Vt[rk:].T if Vt.shape[0] > rk else np.zeros((len(p1[2]), 0))
        # Actually, we want the vectors in the full space
        full1, triples = get_omega2_in_full_space(o1[2], p1[2], n)
        Z2_in_full_T = full1 @ Vt[rk:].T if Vt.shape[0] > rk else np.zeros((len(triples), 0))
        print(f"  dim Z_2(T): {Z2_in_full_T.shape[1]}")
    else:
        Z2_in_full_T = np.zeros((60, 0))
        print(f"  dim Z_2(T): 0")

    if bdo2.get(2) is not None and bdo2[2].shape[1] > 0:
        U2, S2, Vt2 = np.linalg.svd(bdo2[2], full_matrices=True)
        rk2_val = int(np.sum(np.abs(S2) > 1e-8))
        full2, _ = get_omega2_in_full_space(o2[2], p2[2], n)
        Z2_in_full_T2 = full2 @ Vt2[rk2_val:].T if Vt2.shape[0] > rk2_val else np.zeros((len(triples), 0))
        print(f"  dim Z_2(T'): {Z2_in_full_T2.shape[1]}")
    else:
        Z2_in_full_T2 = np.zeros((60, 0))
        print(f"  dim Z_2(T'): 0")

    # Intersection of Z_2 subspaces
    if Z2_in_full_T.shape[1] > 0 and Z2_in_full_T2.shape[1] > 0:
        combined_Z = np.hstack([Z2_in_full_T, Z2_in_full_T2])
        rk_comb = np.linalg.matrix_rank(combined_Z, tol=1e-8)
        rk_z1 = np.linalg.matrix_rank(Z2_in_full_T, tol=1e-8)
        rk_z2 = np.linalg.matrix_rank(Z2_in_full_T2, tol=1e-8)
        z_inter = rk_z1 + rk_z2 - rk_comb
        print(f"  Z_2 intersection dim: {z_inter}")
        print(f"  Z_2 exclusive to T: {rk_z1 - z_inter}")
        print(f"  Z_2 exclusive to T': {rk_z2 - z_inter}")


# ANALYSIS 3: Rank formula via delta(|A_2|) and delta(|A_3|)
print(f"\n{'='*70}")
print("ANALYSIS 3: |A_p| (allowed path counts) under arc flip")
print("=" * 70)

# For each p, track how |allowed_p| changes
delta_A_counts = {}
for bits in range(total):
    A = build_adj(n, bits)
    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

            if A[i][j] == 1:
                uu, vv = i, j
            else:
                uu, vv = j, i

            scores = [sum(row) for row in A]
            du, dv = scores[uu], scores[vv]

            dA2 = len(p2[2]) - len(p1[2])
            dA3 = len(p2[3]) - len(p1[3])
            dOm2 = d2[2] - d1[2]
            dOm3 = d2[3] - d1[3]

            key = (du - dv, dA2, dA3, dOm2, dOm3)
            delta_A_counts[key] = delta_A_counts.get(key, 0) + 1

print("(du-dv, dA2, dA3, dOm2, dOm3) distribution:")
for key in sorted(delta_A_counts.keys()):
    du_dv, dA2, dA3, dOm2, dOm3 = key
    dNT2 = dOm2 - dA2  # This is wrong; need to separate TT from NT
    print(f"  du-dv={du_dv:+d}: dA2={dA2:+d}, dA3={dA3:+d}, "
          f"dOm2={dOm2:+d}, dOm3={dOm3:+d}, count={delta_A_counts[key]}")


# ANALYSIS 4: Is delta(|A_2|) determined by (du, dv)?
print(f"\n{'='*70}")
print("ANALYSIS 4: delta(|A_2|) formula")
print("=" * 70)

dA2_by_score = {}
for bits in range(total):
    A = build_adj(n, bits)
    p1 = all_chain[bits][0]
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2 = all_chain[bits2][0]

            if A[i][j] == 1:
                uu, vv = i, j
            else:
                uu, vv = j, i

            du, dv = scores[uu], scores[vv]
            dA2 = len(p2[2]) - len(p1[2])

            key = (du, dv)
            if key not in dA2_by_score:
                dA2_by_score[key] = set()
            dA2_by_score[key].add(dA2)

print("delta(|A_2|) by (du, dv):")
deterministic = 0
non_det = 0
for key in sorted(dA2_by_score.keys()):
    du, dv = key
    vals = dA2_by_score[key]
    status = "DETERMINED" if len(vals) == 1 else "NOT determined"
    if len(vals) == 1:
        deterministic += 1
    else:
        non_det += 1
    if len(vals) <= 5:
        print(f"  (du={du}, dv={dv}): {sorted(vals)} -- {status}")
    else:
        print(f"  (du={du}, dv={dv}): {len(vals)} values -- {status}")
print(f"  Determined: {deterministic}, Not: {non_det}")


# ANALYSIS 5: HYP-227/228 verification - known formulas
print(f"\n{'='*70}")
print("ANALYSIS 5: Verify delta(|A_2|) = 2*(du - dv - 1) (HYP-228)")
print("=" * 70)

errors_228 = 0
for bits in range(total):
    A = build_adj(n, bits)
    p1 = all_chain[bits][0]
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2 = all_chain[bits2][0]

            if A[i][j] == 1:
                uu, vv = i, j
            else:
                uu, vv = j, i

            du, dv = scores[uu], scores[vv]
            dA2 = len(p2[2]) - len(p1[2])
            predicted = 2*(du - dv - 1)

            if dA2 != predicted:
                errors_228 += 1

print(f"  Errors: {errors_228}")
if errors_228 == 0:
    print("  CONFIRMED: delta(|A_2|) = 2*(du - dv - 1)")
    print("  Note: This is for allowed 2-paths, not TT triples")


# ANALYSIS 6: Nullity-rank decomposition
print(f"\n{'='*70}")
print("ANALYSIS 6: Rank-nullity tracking")
print("=" * 70)
print("For d_2: Om_2 -> Om_1: rk + nullity = dim(Om_2)")
print("nullity = dim(Z_2) = dim(Om_2) - rk(d_2)")
print("d_3: Om_3 -> Om_2: rk + nullity = dim(Om_3)")
print("beta_2 = dim(Z_2) - rk(d_3)")
print()
print("We need: delta(beta_2) = 0")
print("  <=> delta(dim Z_2) = delta(rk d_3)")
print("  <=> delta(dim Om_2) - delta(rk d_2) = delta(rk d_3)")
print("  <=> delta(dim Om_2) + delta(beta_1) = delta(rk d_3)  [since drk2=-db1]")
print()
print("Let's track delta(rk d_3) vs delta(dim Om_3) - delta(nullity d_3)")

# nullity(d_3) = dim(Z_3) = dim(Om_3) - rk(d_3)
# So delta(rk d_3) = delta(dim Om_3) - delta(dim Z_3)
# And beta_2 = 0 <=> delta(dim Om_3) - delta(dim Z_3) = delta(dim Om_2) + delta(beta_1)

null_tracking = {}
for bits in range(total):
    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

            dOm2 = d2[2] - d1[2]
            dOm3 = d2[3] - d1[3]
            drk2 = r2[2] - r1[2]
            drk3 = r2[3] - r1[3]
            dZ2 = dOm2 - drk2
            dZ3 = dOm3 - drk3  # delta(nullity d_3) = delta(dim Z_3)
            db1 = b1[1] - b2[1]  # Note: drk2 = -db1, so db1_actual = b2[1]-b1[1] = -drk2
            db1_actual = b2[1] - b1[1]

            key = (dOm2, dOm3, drk2, drk3, db1_actual)
            null_tracking[key] = null_tracking.get(key, 0) + 1

print("(dOm2, dOm3, drk2, drk3, db1) distribution:")
for key in sorted(null_tracking.keys()):
    dOm2, dOm3, drk2, drk3, db1 = key
    dZ2 = dOm2 - drk2
    dZ3 = dOm3 - drk3
    check = "OK" if dZ2 == drk3 else "FAIL"
    print(f"  dOm2={dOm2:+d}, dOm3={dOm3:+d}, drk2={drk2:+d}, drk3={drk3:+d}, "
          f"db1={db1:+d} | dZ2={dZ2:+d}, dZ3={dZ3:+d} | {check} | count={null_tracking[key]}")


# ANALYSIS 7: Is drk3 = dOm3 - dZ3 where dZ3 follows a formula?
print(f"\n{'='*70}")
print("ANALYSIS 7: Can we express drk3 as dOm3 - dZ3?")
print("=" * 70)
print("Since drk3 = dOm2 + db1 (our target identity),")
print("this means dOm3 - dZ3 = dOm2 + db1")
print("=> dZ3 = dOm3 - dOm2 - db1")
print("=> delta(dim Z_3) = delta(dim Om_3) - delta(dim Om_2) - delta(beta_1)")
print()

# Verify this identity
errors_z3 = 0
for bits in range(total):
    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

            dOm2 = d2[2] - d1[2]
            dOm3 = d2[3] - d1[3]
            drk3 = r2[3] - r1[3]
            dZ3 = dOm3 - drk3
            db1 = b2[1] - b1[1]

            predicted_dZ3 = dOm3 - dOm2 - db1
            if dZ3 != predicted_dZ3:
                errors_z3 += 1

print(f"  Verification errors: {errors_z3}")
if errors_z3 == 0:
    print("  CONFIRMED: delta(Z_3) = delta(Om_3) - delta(Om_2) - delta(beta_1)")
    print("  This is EQUIVALENT to beta_2=0 invariance, not an independent fact.")


# ANALYSIS 8: Key observation -- dim(Om_3) and dim(Om_2) formulas?
print(f"\n{'='*70}")
print("ANALYSIS 8: Are dim(Om_2) and dim(Om_3) determined by local data?")
print("=" * 70)

from collections import defaultdict

dOm2_by_local = defaultdict(set)
dOm3_by_local = defaultdict(set)

for bits in range(total):
    A = build_adj(n, bits)
    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2, o2, d2, bd2, bdo2, r2, b2 = all_chain[bits2]

            if A[i][j] == 1:
                uu, vv = i, j
            else:
                uu, vv = j, i

            du, dv = scores[uu], scores[vv]

            # Local data: common out-neighbors, common in-neighbors, etc.
            cout = sum(1 for w in range(n) if w != uu and w != vv and A[uu][w] and A[vv][w])
            cin = sum(1 for w in range(n) if w != uu and w != vv and A[w][uu] and A[w][vv])
            puv = sum(1 for w in range(n) if w != uu and w != vv and A[uu][w] and A[w][vv])
            pvu = sum(1 for w in range(n) if w != uu and w != vv and A[vv][w] and A[w][uu])

            local_key = (du, dv, cout, cin, puv, pvu)

            dOm2 = d2[2] - d1[2]
            dOm3 = d2[3] - d1[3]

            dOm2_by_local[local_key].add(dOm2)
            dOm3_by_local[local_key].add(dOm3)

det2 = sum(1 for v in dOm2_by_local.values() if len(v) == 1)
ndet2 = sum(1 for v in dOm2_by_local.values() if len(v) > 1)
det3 = sum(1 for v in dOm3_by_local.values() if len(v) == 1)
ndet3 = sum(1 for v in dOm3_by_local.values() if len(v) > 1)

print(f"  delta(Om_2) determined by local data: {det2}/{det2+ndet2}")
print(f"  delta(Om_3) determined by local data: {det3}/{det3+ndet3}")

if ndet2 > 0:
    print("  NON-DETERMINED delta(Om_2) examples:")
    for key in sorted(dOm2_by_local.keys()):
        if len(dOm2_by_local[key]) > 1:
            print(f"    (du={key[0]}, dv={key[1]}, cout={key[2]}, cin={key[3]}, "
                  f"puv={key[4]}, pvu={key[5]}): {sorted(dOm2_by_local[key])}")


# ANALYSIS 9: The key - maybe dim(Om_2) has a formula involving |TT| + correction?
print(f"\n{'='*70}")
print("ANALYSIS 9: dim(Om_2) vs |TT| (number of transitive triples)")
print("=" * 70)

om2_vs_tt = defaultdict(set)
for bits in range(total):
    A = build_adj(n, bits)
    p1, o1, d1, bd1, bdo1, r1, b1 = all_chain[bits]

    # Count TT triples
    tt = 0
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if a != b and b != c and a != c:
                    if A[a][b] and A[b][c] and A[a][c]:
                        tt += 1

    om2_vs_tt[tt].add(d1[2])

print("TT count -> dim(Om_2) values:")
for tt in sorted(om2_vs_tt.keys()):
    vals = sorted(om2_vs_tt[tt])
    if len(vals) == 1:
        print(f"  |TT|={tt}: dim(Om_2)={vals[0]} (DETERMINED)")
    else:
        print(f"  |TT|={tt}: dim(Om_2) in {vals}")

# Check if dim(Om_2) = |TT|/6 + correction
print("\n  Check dim(Om_2) = |TT|/6 + something:")
for tt in sorted(om2_vs_tt.keys()):
    vals = sorted(om2_vs_tt[tt])
    print(f"  |TT|={tt}, |TT|/6={tt/6:.1f}, dim(Om_2) in {vals}, diff={[v-tt//6 for v in vals]}")


print("\nDone.")
