"""
beta3_proportionality_proof.py — Algebraic proof strategy for beta_3 <= 1

Key observation from beta3_homology_structure.py:
ALL ker(d_3) basis vectors project to PROPORTIONAL elements in H_3.
This directly implies dim(H_3) <= 1.

The mechanism: ker(d_3) has a large dimension (e.g., 7 at n=6 Type B),
but im(d_4) captures all but one direction of it.

WHY would this proportionality hold? Investigate:
1. Is there an algebraic identity forcing this?
2. Does the tournament structure constrain the cokernel?
3. What's special about d_4 that makes its rank = ker(d_3) - 1 (not less)?

Also: investigate the "3-cycle analogy" — for beta_1 <= 1, the proof uses
"all directed 3-cycles are homologous" (THM-103). Is there an analogous
statement for beta_3?

Author: kind-pasteur-S46 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def compute_omega_basis(ap, p):
    if not ap.get(p, []):
        return np.zeros((0, 0)), 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(len(ap[p])), len(ap[p])

    P = np.zeros((na_count, len(ap[p])))
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    if ns.shape[1] > 0:
        return ns, ns.shape[1]
    else:
        return np.zeros((len(ap[p]), 0)), 0


def analyze_h3_structure(A, n):
    """Full analysis of H_3 structure: ker(d_3), im(d_4), quotient."""
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(min(6, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    dim3 = omega_dims.get(3, 0)
    if dim3 == 0:
        return None

    # Build d_3 in Omega basis
    bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_coeffs(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3  # boundary in Omega_2 coords x Omega_3 dim
    U3, S3, V3t = np.linalg.svd(d3_om, full_matrices=True)
    rank_d3 = int(sum(s > 1e-8 for s in S3))
    ker_d3 = dim3 - rank_d3

    if ker_d3 == 0:
        return None

    # Kernel of d_3: last ker_d3 rows of V3t
    ker_basis = V3t[rank_d3:, :].T  # dim3 x ker_d3

    # Build d_4 in Omega basis
    dim4 = omega_dims.get(4, 0)
    if dim4 == 0:
        return {
            'ker_d3': ker_d3, 'rank_d4_in_ker': 0, 'beta3': ker_d3,
            'omega_dims': omega_dims, 'ker_basis': ker_basis
        }

    bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap.get(4, [])):
        for sign, face in boundary_coeffs(path):
            if face in idx3:
                bd4[idx3[face], j] += sign

    O4 = omega_bases[4]
    d4_om = bd4 @ O4  # A_3 coords x Omega_4 dim

    # Project d4 image into Omega_3 basis
    O3pinv = np.linalg.pinv(O3)
    d4_omega3 = O3pinv @ d4_om  # Omega_3 dim x Omega_4 dim

    # Now project into ker(d_3)
    # d4_in_ker = ker_basis^T @ d4_omega3
    d4_in_ker = ker_basis.T @ d4_omega3  # ker_d3 x dim4

    sv_in_ker = np.linalg.svd(d4_in_ker, compute_uv=False)
    rank_d4_in_ker = int(sum(s > 1e-8 for s in sv_in_ker))

    beta3 = ker_d3 - rank_d4_in_ker

    return {
        'ker_d3': ker_d3,
        'rank_d4_in_ker': rank_d4_in_ker,
        'beta3': beta3,
        'omega_dims': omega_dims,
        'ker_basis': ker_basis,
        'd4_in_ker': d4_in_ker,
        'dim4': dim4,
    }


def main():
    print("=" * 70)
    print("BETA_3 PROPORTIONALITY — Algebraic proof investigation")
    print("=" * 70)

    # Part 1: At n=6, analyze the d_4-in-ker matrix structure
    print("\n--- Part 1: d_4-in-ker structure at n=6 (Type B, scores 2,2,2,3,3,3) ---")
    n = 6
    type_b_count = 0
    type_a_count = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        result = analyze_h3_structure(A, n)
        if result is None:
            continue

        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores == (2, 2, 2, 3, 3, 3):
            type_b_count += 1
            if type_b_count <= 3:
                print(f"\n  Type B #{type_b_count} (bits={bits}):")
                print(f"    ker(d_3) = {result['ker_d3']}")
                print(f"    rank(d_4 in ker) = {result['rank_d4_in_ker']}")
                print(f"    beta_3 = {result['beta3']}")
                print(f"    Omega dims: {result['omega_dims']}")
                print(f"    d_4-in-ker matrix shape: {result['d4_in_ker'].shape}")

                # Check: are all rows of d4_in_ker proportional?
                d4k = result['d4_in_ker']
                if d4k.shape[0] > 1 and d4k.shape[1] > 0:
                    # Find a nonzero row
                    norms = np.linalg.norm(d4k, axis=1)
                    ref_idx = np.argmax(norms)
                    ref = d4k[ref_idx]
                    if norms[ref_idx] > 1e-10:
                        ratios = []
                        for i in range(d4k.shape[0]):
                            if i == ref_idx:
                                continue
                            # Check proportionality
                            cross = np.abs(np.cross(ref[:2] if len(ref) >= 2 else [ref[0], 0],
                                                    d4k[i, :2] if len(d4k[i]) >= 2 else [d4k[i, 0], 0]))
                            # Better: use SVD rank check of the 2-row matrix
                            mat = np.vstack([ref, d4k[i]])
                            sv = np.linalg.svd(mat, compute_uv=False)
                            rk = int(sum(s > 1e-8 for s in sv))
                            ratios.append(rk)
                        print(f"    Row-pair ranks: {ratios} (all 1 = proportional)")

        elif scores == (1, 1, 1, 4, 4, 4):
            type_a_count += 1

    print(f"\n  Total Type B (2,2,2,3,3,3): {type_b_count}")
    print(f"  Total Type A (1,1,1,4,4,4): {type_a_count}")

    # Part 2: Look at the cokernel direction at n=6 Type B
    print("\n--- Part 2: Cokernel direction analysis ---")
    cokernel_vecs = []
    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        result = analyze_h3_structure(A, n)
        if result is None or result['beta3'] == 0:
            continue

        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores != (2, 2, 2, 3, 3, 3):
            continue

        # The cokernel element: kernel of d4_in_ker^T
        d4k = result['d4_in_ker']
        ker_basis = result['ker_basis']

        U, S, Vt = np.linalg.svd(d4k.T, full_matrices=True)
        rank = int(sum(s > 1e-8 for s in S))
        # Null space of d4k = right null of d4k^T
        # Actually cokernel = ker(d_3) / im(d_4 in ker)
        # Cokernel basis = null space of d4k^T acting on ker(d_3)
        # d4k: ker_d3 x dim4. Its row space spans im(d_4) in ker(d_3).
        # Cokernel = ker(d_3) / row_space(d4k).
        # The cokernel direction is the left null space of d4k.
        U_left, S_left, Vt_left = np.linalg.svd(d4k, full_matrices=True)
        rank_left = int(sum(s > 1e-8 for s in S_left))
        null_left = U_left[:, rank_left:]  # ker_d3 x (ker_d3 - rank)
        if null_left.shape[1] == 1:
            cokernel_vecs.append(null_left[:, 0])

    if cokernel_vecs:
        print(f"  Number of cokernel vectors collected: {len(cokernel_vecs)}")
        # Check if all are the same (up to sign)
        ref = cokernel_vecs[0]
        ref_norm = ref / np.linalg.norm(ref)
        all_same = True
        for v in cokernel_vecs[1:]:
            v_norm = v / np.linalg.norm(v)
            dot = abs(np.dot(ref_norm, v_norm))
            if abs(dot - 1.0) > 1e-6:
                all_same = False
                break
        if all_same:
            print("  ALL cokernel directions are the same (up to sign)!")
            print(f"  Cokernel vector (in ker(d_3) basis): {ref_norm}")
        else:
            print("  Cokernel directions vary across tournaments")
            # Check pairwise dot products
            dots = []
            for i in range(min(10, len(cokernel_vecs))):
                for j in range(i+1, min(10, len(cokernel_vecs))):
                    v1 = cokernel_vecs[i] / np.linalg.norm(cokernel_vecs[i])
                    v2 = cokernel_vecs[j] / np.linalg.norm(cokernel_vecs[j])
                    dots.append(abs(np.dot(v1, v2)))
            print(f"  Pairwise |dot products| (first 10): min={min(dots):.4f}, max={max(dots):.4f}, mean={np.mean(dots):.4f}")
            # Print a few cokernel dimensions
            dims = [len(v) for v in cokernel_vecs]
            print(f"  Cokernel vector dimensions: {Counter(dims)}")

    # Part 3: Analogy with beta_1: "all X are homologous"
    # For beta_1: all directed 3-cycles are homologous in Omega_1.
    # For beta_3: what are the "elementary" 3-cycles in Omega_3?
    # A 3-cycle in Omega_3 means a cycle: sum of Omega_3 elements with zero boundary.
    print("\n--- Part 3: Elementary H_3 generators at n=6 ---")
    n = 6
    gen_count = 0
    gen_path_counts = Counter()
    gen_vertex_counts = Counter()

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        result = analyze_h3_structure(A, n)
        if result is None or result['beta3'] == 0:
            continue

        gen_count += 1
        # The H_3 generator is the cokernel direction in Omega_3 coordinates
        ker_basis = result['ker_basis']
        d4k = result.get('d4_in_ker', None)

        if d4k is not None and d4k.size > 0:
            U_left, S_left, Vt_left = np.linalg.svd(d4k, full_matrices=True)
            rank_left = int(sum(s > 1e-8 for s in S_left))
            null_left = U_left[:, rank_left:]
            if null_left.shape[1] != 1:
                continue
            h3_gen_omega = ker_basis @ null_left[:, 0]
        else:
            # No d_4 => H_3 = ker(d_3) entirely. If beta_3=1, ker_d3=1.
            if result['ker_d3'] != 1:
                continue
            h3_gen_omega = ker_basis[:, 0]

        # Back to allowed 3-paths
        ap3 = enumerate_allowed_paths(A, n, 3)
        O3, d3 = compute_omega_basis({p: enumerate_allowed_paths(A, n, p) for p in range(min(6, n))}, 3)
        h3_gen_paths = O3 @ h3_gen_omega

        # Vertex sets used
        vsets = set()
        nonzero_paths = 0
        for i, coeff in enumerate(h3_gen_paths):
            if abs(coeff) > 1e-10:
                nonzero_paths += 1
                vsets.add(frozenset(ap3[i]))

        gen_path_counts[nonzero_paths] += 1
        gen_vertex_counts[len(vsets)] += 1

    print(f"  Total H_3 generators: {gen_count}")
    print(f"  #nonzero allowed 3-paths in generator: {dict(sorted(gen_path_counts.items()))}")
    print(f"  #vertex sets in generator: {dict(sorted(gen_vertex_counts.items()))}")

    # Part 4: n=7 — ker(d_3) dimensions and d4-rank structure
    print("\n--- Part 4: ker(d_3) and d_4-rank structure at n=7 ---")
    n = 7
    rng7 = np.random.RandomState(77)
    N7 = 200
    ker_rank_pairs = Counter()
    b3_count7 = 0

    for trial in range(N7):
        A = random_tournament(n, rng7)
        result = analyze_h3_structure(A, n)
        if result is None:
            continue
        b3_count7 += 1

        ker = result['ker_d3']
        rk = result['rank_d4_in_ker']
        ker_rank_pairs[(ker, rk, result['beta3'])] += 1

        if (trial + 1) % 50 == 0:
            print(f"  n=7: {trial+1}/{N7} done, {b3_count7} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count7}/{N7}")
    print(f"  (ker_d3, rank_d4_in_ker, beta_3) distribution:")
    for (k, r, b), cnt in sorted(ker_rank_pairs.items()):
        print(f"    ker={k}, rank_d4={r}, beta3={b}: {cnt}")

    # Part 5: At n=7, is there any case with beta_3 > 1?
    print("\n--- Part 5: Exhaustive beta_3 check at n=7 (large sample) ---")
    n = 7
    rng_big = np.random.RandomState(12345)
    N_big = 2000
    beta3_dist = Counter()
    max_beta3 = 0

    for trial in range(N_big):
        A = random_tournament(n, rng_big)
        result = analyze_h3_structure(A, n)
        if result is None:
            beta3_dist[0] += 1
            continue

        beta3_dist[result['beta3']] += 1
        if result['beta3'] > max_beta3:
            max_beta3 = result['beta3']
            if max_beta3 > 1:
                print(f"  FOUND beta_3 > 1! trial={trial}, beta_3={max_beta3}")

        if (trial + 1) % 500 == 0:
            print(f"  n=7: {trial+1}/{N_big} done", flush=True)

    print(f"\n  beta_3 distribution (n=7, {N_big} samples):")
    for b, cnt in sorted(beta3_dist.items()):
        pct = 100*cnt/N_big
        print(f"    beta_3={b}: {cnt} ({pct:.1f}%)")
    print(f"  Max beta_3: {max_beta3}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
