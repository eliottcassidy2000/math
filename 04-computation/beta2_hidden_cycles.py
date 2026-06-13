#!/usr/bin/env python3
"""beta2_hidden_cycles.py - Analyze the "hidden cycles" at bad vertices

When b1(T)=0 and b1(T\v)=1: there exists a 1-cycle z_v in Z_1(T\v)
that's NOT in im(d_2(T\v)) but IS in im(d_2(T)).

The filling w of z_v in T necessarily uses 2-paths through v.

KEY QUESTIONS:
1. What do the hidden cycles z_v look like?
2. How do hidden cycles at different bad vertices interact?
3. Can we derive a contradiction from having n bad vertices?

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def compute_homology_data(A, n):
    """Compute full homology data: Z1, im(d2), b1, and H1 generator if b1>0."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)
    path1_list = [tuple(p) for p in paths1]
    path2_list = [tuple(p) for p in paths2] if paths2 else []

    if not paths1:
        return {'b1': 0, 'paths1': [], 'Z1': None, 'im_d2': None, 'h1_gen': None}

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return {'b1': 0, 'paths1': path1_list, 'Z1': None, 'im_d2': None, 'h1_gen': None}

    D1 = build_full_boundary_matrix(path1_list, paths0)
    D1_om = D1 @ omega1

    U1, S1, Vt1 = np.linalg.svd(D1_om, full_matrices=True)
    rk_d1 = int(sum(s > 1e-8 for s in S1))

    # Z1 basis in path coordinates
    Z1_omega = Vt1[rk_d1:].T  # columns are Z1 basis vectors in Omega1 coords
    Z1 = omega1 @ Z1_omega  # columns are Z1 basis vectors in path coords
    dim_Z1 = Z1.shape[1]

    # im(d2) in path coordinates
    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        if dim_O2 > 0:
            D2 = build_full_boundary_matrix(path2_list, path1_list)
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))
            im_d2 = D2_om  # columns are d2 images
        else:
            rk_d2 = 0
            im_d2 = np.zeros((len(path1_list), 0))
    else:
        rk_d2 = 0
        im_d2 = np.zeros((len(path1_list), 0))

    b1 = dim_Z1 - rk_d2

    # H1 generator if b1 > 0
    h1_gen = None
    if b1 > 0 and rk_d2 > 0:
        # Project Z1 onto orthogonal complement of im_d2
        U2, S2, _ = np.linalg.svd(im_d2, full_matrices=True)
        proj = np.eye(len(path1_list)) - U2[:, :rk_d2] @ U2[:, :rk_d2].T
        Z1_proj = proj @ Z1
        norms = np.linalg.norm(Z1_proj, axis=0)
        idx = np.argmax(norms)
        if norms[idx] > 1e-8:
            h1_gen = Z1_proj[:, idx]
            h1_gen = h1_gen / h1_gen[np.argmax(np.abs(h1_gen))]
    elif b1 > 0:
        h1_gen = Z1[:, 0]
        h1_gen = h1_gen / h1_gen[np.argmax(np.abs(h1_gen))]

    return {
        'b1': b1, 'paths1': path1_list, 'paths2': path2_list,
        'Z1': Z1, 'im_d2': im_d2, 'h1_gen': h1_gen,
        'rk_d2': rk_d2, 'dim_Z1': dim_Z1
    }


# ============================================================
# Part 1: Analyze hidden cycles at bad vertices, n=5
# ============================================================
print("=" * 70)
print("HIDDEN CYCLES AT BAD VERTICES, n=5")
print("=" * 70)

n = 5
total = 2 ** (n * (n - 1) // 2)
examples_found = 0

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    # Quick b1 check
    data_T = compute_homology_data(A, n)
    if data_T['b1'] != 0:
        continue

    # Check each vertex
    bad_vs = []
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        data_Tv = compute_homology_data(B_sub, n - 1)
        if data_Tv['b1'] == 1:
            bad_vs.append((v, data_Tv, vlist))

    if len(bad_vs) >= 2 and examples_found < 3:
        examples_found += 1
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        print(f"\nExample {examples_found}: bits={bits}, scores={scores}, #bad={len(bad_vs)}")

        for v, data_Tv, vlist in bad_vs:
            print(f"\n  Bad vertex v={v} (d_out={sum(A[v])})")
            gen = data_Tv['h1_gen']
            if gen is not None:
                # Convert to original vertex labels
                path1_sub = data_Tv['paths1']
                print(f"    H1 generator of T\\{v} (in T\\{v} labels):")
                for i, p in enumerate(path1_sub):
                    if abs(gen[i]) > 1e-8:
                        # Convert to original labels
                        orig_p = tuple(vlist[x] for x in p)
                        print(f"      {orig_p}: {gen[i]:.4f}")

        # Check: do the hidden cycles at different bad vertices interact?
        if len(bad_vs) >= 2:
            print(f"\n  Hidden cycle interaction:")
            v1, data1, vlist1 = bad_vs[0]
            v2, data2, vlist2 = bad_vs[1]
            gen1 = data1['h1_gen']
            gen2 = data2['h1_gen']

            # Embed both into full arc space and compute overlap
            full_paths1 = data_T['paths1']
            full_gen1 = np.zeros(len(full_paths1))
            full_gen2 = np.zeros(len(full_paths1))

            path1_sub1 = data1['paths1']
            path1_sub2 = data2['paths1']

            for i, p in enumerate(path1_sub1):
                orig_p = tuple(vlist1[x] for x in p)
                if orig_p in full_paths1:
                    j = full_paths1.index(orig_p)
                    full_gen1[j] = gen1[i]

            for i, p in enumerate(path1_sub2):
                orig_p = tuple(vlist2[x] for x in p)
                if orig_p in full_paths1:
                    j = full_paths1.index(orig_p)
                    full_gen2[j] = gen2[i]

            # Inner product
            ip = np.dot(full_gen1, full_gen2)
            # Are they linearly independent?
            mat = np.column_stack([full_gen1, full_gen2])
            sv = np.linalg.svd(mat, compute_uv=False)
            rank = int(sum(s > 1e-8 for s in sv))
            print(f"    Inner product: {ip:.4f}")
            print(f"    Rank of [z_{v1}, z_{v2}]: {rank}")

    if examples_found >= 3:
        break

# ============================================================
# Part 2: Statistics on hidden cycles
# ============================================================
print(f"\n\n{'=' * 70}")
print("HIDDEN CYCLE STATISTICS, n=5")
print("=" * 70)

n = 5
total = 2 ** (n * (n - 1) // 2)
pair_ranks = Counter()
pair_ip = []
num_bad_dist = Counter()

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    data_T = compute_homology_data(A, n)
    if data_T['b1'] != 0:
        continue

    bad_data = []
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        data_Tv = compute_homology_data(B_sub, n - 1)
        if data_Tv['b1'] == 1:
            # Embed hidden cycle into full space
            gen = data_Tv['h1_gen']
            if gen is None:
                continue
            full_gen = np.zeros(len(data_T['paths1']))
            path1_sub = data_Tv['paths1']
            for i, p in enumerate(path1_sub):
                orig_p = tuple(vlist[x] for x in p)
                if orig_p in data_T['paths1']:
                    j = data_T['paths1'].index(orig_p)
                    full_gen[j] = gen[i]
            bad_data.append((v, full_gen))

    num_bad = len(bad_data)
    num_bad_dist[num_bad] += 1

    if num_bad >= 2:
        # Check pairwise rank
        for i in range(len(bad_data)):
            for j in range(i + 1, len(bad_data)):
                mat = np.column_stack([bad_data[i][1], bad_data[j][1]])
                sv = np.linalg.svd(mat, compute_uv=False)
                rank = int(sum(s > 1e-8 for s in sv))
                pair_ranks[rank] += 1
                ip = np.dot(bad_data[i][1], bad_data[j][1])
                pair_ip.append(ip)

        # Check: rank of ALL hidden cycles together
        if num_bad >= 2:
            mat_all = np.column_stack([bd[1] for bd in bad_data])
            sv_all = np.linalg.svd(mat_all, compute_uv=False)
            rank_all = int(sum(s > 1e-8 for s in sv_all))

print(f"\n#bad vertices distribution:")
for k in sorted(num_bad_dist.keys()):
    print(f"  {k} bad: {num_bad_dist[k]}")

print(f"\nPairwise rank of hidden cycles:")
for r in sorted(pair_ranks.keys()):
    print(f"  rank {r}: {pair_ranks[r]}")

if pair_ip:
    print(f"\nInner product distribution: min={min(pair_ip):.4f} max={max(pair_ip):.4f}")
    print(f"  mean={np.mean(pair_ip):.4f}")


# ============================================================
# Part 3: Hidden cycles at n=6, check independence
# ============================================================
print(f"\n\n{'=' * 70}")
print("HIDDEN CYCLE RANK vs #BAD VERTICES, n=6")
print("=" * 70)

n = 6
total = 2 ** (n * (n - 1) // 2)
rank_vs_bad = Counter()
num_checked = 0

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    data_T = compute_homology_data(A, n)
    if data_T['b1'] != 0:
        continue

    bad_gens = []
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        data_Tv = compute_homology_data(B_sub, n - 1)
        if data_Tv['b1'] == 1:
            gen = data_Tv['h1_gen']
            if gen is None:
                continue
            full_gen = np.zeros(len(data_T['paths1']))
            path1_sub = data_Tv['paths1']
            for i, p in enumerate(path1_sub):
                orig_p = tuple(vlist[x] for x in p)
                if orig_p in data_T['paths1']:
                    j = data_T['paths1'].index(orig_p)
                    full_gen[j] = gen[i]
            bad_gens.append(full_gen)

    num_bad = len(bad_gens)
    if num_bad >= 2:
        mat = np.column_stack(bad_gens)
        sv = np.linalg.svd(mat, compute_uv=False)
        rank = int(sum(s > 1e-8 for s in sv))
        rank_vs_bad[(num_bad, rank)] += 1
    elif num_bad == 1:
        rank_vs_bad[(1, 1)] += 1

    num_checked += 1
    if num_checked % 5000 == 0:
        print(f"  checked {num_checked}...")

print(f"\n(#bad, rank_of_hidden_cycles):")
for (nb, rk), count in sorted(rank_vs_bad.items()):
    print(f"  #bad={nb}, rank={rk}: {count}")


# ============================================================
# Part 4: Critical check - does rank(hidden cycles) < #bad?
# ============================================================
print(f"\n{'=' * 70}")
print("RANK DEFICIENCY OF HIDDEN CYCLES")
print("=" * 70)

# If the hidden cycles are NOT independent, this gives a structural
# constraint that might lead to a proof.

for (nb, rk), count in sorted(rank_vs_bad.items()):
    deficiency = nb - rk
    if deficiency > 0:
        print(f"  DEFICIENT: #bad={nb}, rank={rk}, deficiency={deficiency}: {count} cases")

no_deficiency = all(nb == rk for (nb, rk) in rank_vs_bad.keys())
if no_deficiency:
    print("  No deficiency found - hidden cycles are always independent!")
else:
    print("  Hidden cycles are sometimes linearly dependent!")


print("\n\nDone.")
