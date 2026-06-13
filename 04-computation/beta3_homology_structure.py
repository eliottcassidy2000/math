"""
beta3_homology_structure.py — Why is H_3 at most 1-dimensional?

For beta_1: THM-103 proves b1 <= 1 by showing all 3-cycles are homologous.
Key: any two 3-cycles sharing a directed edge differ by a boundary.

For beta_3: can we show all elements of ker(d_3) are "homologous mod im(d_4)"?
I.e., ker(d_3) / im(d_4) is at most 1-dimensional?

Strategy:
1. Find a basis for ker(d_3) at n=6 (Type B: 7-dimensional)
2. Check if any two basis elements differ by a boundary from Omega_4
3. Look for a "common generator" that all ker(d_3) elements are multiples of (mod im(d_4))

Also investigate:
- Does the H_3 generator have a nice combinatorial description?
- Is there a "3-cycle analogue" of the shared-edge homology argument?

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
    if not ap[p]:
        return np.zeros((0, 0)), 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

    apm1_set = set(ap[p-1])
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


def analyze_tournament(bits, n):
    """Deep analysis of H_3 structure for a given tournament."""
    A = bits_to_adj(bits, n)
    scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))

    ap = {}
    for p in range(n):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(n):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    # Build d_3: Omega_3 -> A_2
    dim3 = omega_dims.get(3, 0)
    if dim3 == 0:
        return None

    bd3 = np.zeros((len(ap[2]), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap[2])}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_coeffs(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3

    U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
    rank_d3 = int(sum(s > 1e-8 for s in S3))
    ker_d3_basis = Vt3[rank_d3:].T  # Columns in Omega_3 coordinates
    ker_dim = ker_d3_basis.shape[1] if ker_d3_basis.ndim == 2 else 0

    if ker_dim == 0:
        return None

    # Build d_4: Omega_4 -> A_3
    dim4 = omega_dims.get(4, 0)
    if dim4 == 0:
        # im(d_4) = 0, beta_3 = ker_dim
        return {
            'scores': scores, 'ker_dim': ker_dim, 'im_dim': 0,
            'beta3': ker_dim, 'om4_zero': True
        }

    bd4 = np.zeros((len(ap[3]), len(ap[4])))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap[4]):
        for sign, face in boundary_coeffs(path):
            if face in idx3:
                bd4[idx3[face], j] += sign

    O4 = omega_bases[4]
    d4_om = bd4 @ O4

    # Convert d4 image to Omega_3 coordinates
    O3pinv = np.linalg.pinv(O3)
    d4_omega3 = O3pinv @ d4_om

    U4, S4, Vt4 = np.linalg.svd(d4_omega3, full_matrices=True)
    rank_d4 = int(sum(s > 1e-8 for s in S4))
    im_d4_basis = U4[:, :rank_d4]

    beta3 = ker_dim - rank_d4
    if beta3 <= 0:
        return None

    # Now analyze the quotient ker(d_3) / im(d_4)
    # Project each ker(d_3) basis vector onto the quotient
    quotient_reps = []
    for i in range(ker_dim):
        v = ker_d3_basis[:, i]
        if rank_d4 > 0:
            proj = im_d4_basis @ (im_d4_basis.T @ v)
            residual = v - proj
        else:
            residual = v.copy()
        quotient_reps.append(residual)

    # Check: are all quotient representatives parallel?
    # (This would mean H_3 is 1-dimensional)
    # Compare via cross products / ratios
    nonzero_reps = [r for r in quotient_reps if np.linalg.norm(r) > 1e-6]

    if len(nonzero_reps) >= 2:
        # Check if they're proportional
        ref = nonzero_reps[0] / np.linalg.norm(nonzero_reps[0])
        all_proportional = True
        for r in nonzero_reps[1:]:
            r_normed = r / np.linalg.norm(r)
            # Check if r_normed is +/- ref
            if abs(abs(np.dot(ref, r_normed)) - 1) > 1e-6:
                all_proportional = False
                break

        print(f"  ker(d_3) basis has {ker_dim} vectors, {len(nonzero_reps)} project to nonzero quotient elements")
        print(f"  All quotient projections proportional: {all_proportional}")

    # Express the H_3 generator in A_3 coordinates
    gen_omega = nonzero_reps[0] / np.linalg.norm(nonzero_reps[0])
    gen_A3 = O3 @ gen_omega

    nonzero_paths = [(ap[3][j], gen_A3[j]) for j in range(len(ap[3])) if abs(gen_A3[j]) > 1e-6]
    vsets = set(frozenset(p) for p, c in nonzero_paths)

    # Analyze coefficient structure
    coeffs = [c for p, c in nonzero_paths]
    abs_coeffs = sorted(set(abs(c) for c in coeffs))

    # What are the vertex sets?
    vset_list = sorted(vsets, key=lambda s: sorted(s))

    return {
        'scores': scores,
        'ker_dim': ker_dim,
        'im_dim': rank_d4,
        'beta3': beta3,
        'num_paths': len(nonzero_paths),
        'num_vsets': len(vsets),
        'vset_sizes': Counter(len(vs) for vs in vsets),
        'all_proportional': all_proportional if len(nonzero_reps) >= 2 else True,
        'abs_coeffs': abs_coeffs,
        'signs': Counter(int(np.sign(c)) for p, c in nonzero_paths),
        'om4_zero': False,
        'dim3': dim3,
        'dim4': dim4,
        'vertex_coverage': len(set().union(*vsets))
    }


def main():
    print("=" * 70)
    print("H_3 HOMOLOGY STRUCTURE — WHY dim(H_3) <= 1?")
    print("=" * 70)

    # Part 1: Detailed analysis at n=6, Type B (scores 2,2,2,3,3,3)
    print("\n--- Part 1: n=6 Type B analysis ---")
    n = 6
    type_b_count = 0
    type_a_count = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores != (2, 2, 2, 3, 3, 3):
            continue

        result = analyze_tournament(bits, n)
        if result is None or result['beta3'] == 0:
            continue

        type_b_count += 1
        if type_b_count <= 5:
            print(f"\n  bits={bits}:")
            print(f"    ker(d_3)={result['ker_dim']}, im(d_4)={result['im_dim']}, beta_3={result['beta3']}")
            print(f"    dim(Om3)={result['dim3']}, dim(Om4)={result['dim4']}")
            print(f"    Generator: {result['num_paths']} paths on {result['num_vsets']} vertex sets")
            print(f"    Signs: {dict(result['signs'])}")
            print(f"    All quotient projections proportional: {result['all_proportional']}")
            print(f"    Abs coefficients: {[f'{c:.6f}' for c in result['abs_coeffs']]}")
            print(f"    Vertex coverage: {result['vertex_coverage']}/{n}")

    print(f"\n  Total Type B with beta_3=1: {type_b_count}")

    # Part 2: Type A analysis
    print("\n--- Part 2: n=6 Type A analysis ---")
    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores != (1, 1, 1, 4, 4, 4):
            continue

        result = analyze_tournament(bits, n)
        if result is None or result['beta3'] == 0:
            continue

        type_a_count += 1
        if type_a_count <= 3:
            print(f"\n  bits={bits}:")
            if result.get('om4_zero'):
                print(f"    Om4=0, ker(d_3)={result['ker_dim']}")
            else:
                print(f"    ker(d_3)={result['ker_dim']}, im(d_4)={result['im_dim']}")
                print(f"    Generator: {result['num_paths']} paths on {result['num_vsets']} vertex sets")
                print(f"    All proportional: {result.get('all_proportional', 'N/A')}")

    print(f"\n  Total Type A with beta_3=1: {type_a_count}")

    # Part 3: n=7 analysis — focus on the quotient proportionality
    print("\n--- Part 3: n=7 sampled ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 500
    proportional_count = 0
    non_prop_count = 0
    b3_count = 0

    for trial in range(N):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        # Quick beta_3 check
        ap = {}
        for p in range(n):
            ap[p] = enumerate_allowed_paths(A, n, p)

        omega_bases = {}
        omega_dims = {}
        for p in range(n):
            if not ap.get(p, []):
                omega_bases[p] = np.zeros((0, 0))
                omega_dims[p] = 0
            else:
                omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

        dim3 = omega_dims.get(3, 0)
        dim4 = omega_dims.get(4, 0)
        if dim3 == 0:
            continue

        # Build d_3
        bd3 = np.zeros((len(ap[2]), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap[2])}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign

        O3 = omega_bases[3]
        d3_om = bd3 @ O3
        sv3 = np.linalg.svd(d3_om, compute_uv=False)
        rank_d3 = int(sum(s > 1e-8 for s in sv3))
        ker_dim = dim3 - rank_d3

        if ker_dim == 0:
            continue

        if dim4 == 0:
            b3_count += 1
            if ker_dim == 1:
                proportional_count += 1
            continue

        # Build d_4
        bd4 = np.zeros((len(ap[3]), len(ap[4])))
        idx3 = {path: i for i, path in enumerate(ap[3])}
        for j, path in enumerate(ap[4]):
            for sign, face in boundary_coeffs(path):
                if face in idx3:
                    bd4[idx3[face], j] += sign

        O4 = omega_bases[4]
        d4_om = bd4 @ O4
        O3pinv = np.linalg.pinv(O3)
        d4_omega3 = O3pinv @ d4_om

        sv4 = np.linalg.svd(d4_omega3, compute_uv=False)
        rank_d4 = int(sum(s > 1e-8 for s in sv4))
        beta3 = ker_dim - rank_d4

        if beta3 > 0:
            b3_count += 1
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            if beta3 == 1:
                proportional_count += 1
            else:
                non_prop_count += 1
                print(f"  beta_3={beta3} > 1! trial={trial}, scores={scores}")

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N} done, {b3_count} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count}")
    print(f"  beta_3=1: {proportional_count}")
    print(f"  beta_3>1: {non_prop_count}")

    # Part 4: Specific investigation — does dim(ker(d_3)) correlate with something simple?
    print("\n--- Part 4: ker(d_3) dimensions at n=6 (exhaustive) ---")
    n = 6
    ker_d3_counter = Counter()
    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        ap3 = enumerate_allowed_paths(A, n, 3)
        ap2 = enumerate_allowed_paths(A, n, 2)
        if not ap3:
            ker_d3_counter[0] += 1
            continue

        # Compute Omega_3 basis
        O3, dim3 = compute_omega_basis({2: ap2, 3: ap3}, 3)
        if dim3 == 0:
            ker_d3_counter[0] += 1
            continue

        # Build d_3
        bd3 = np.zeros((len(ap2), len(ap3)))
        idx2 = {path: i for i, path in enumerate(ap2)}
        for j, path in enumerate(ap3):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign

        d3_om = bd3 @ O3
        sv = np.linalg.svd(d3_om, compute_uv=False)
        rank_d3 = int(sum(s > 1e-8 for s in sv))
        ker_dim = dim3 - rank_d3
        ker_d3_counter[ker_dim] += 1

    print(f"  dim(ker(d_3)) distribution:")
    for k in sorted(ker_d3_counter.keys()):
        cnt = ker_d3_counter[k]
        pct = 100*cnt/2**(n*(n-1)//2)
        print(f"    ker(d_3)={k}: {cnt} ({pct:.1f}%)")

    # Part 5: For Type B, analyze the 36 paths in the generator
    print("\n--- Part 5: Generator path analysis for first Type B tournament ---")
    n = 6
    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores != (2, 2, 2, 3, 3, 3):
            continue

        ap = {}
        for p in range(6):
            ap[p] = enumerate_allowed_paths(A, n, p)

        omega_bases = {}
        omega_dims = {}
        for p in range(6):
            if not ap.get(p, []):
                omega_bases[p] = np.zeros((0, 0))
                omega_dims[p] = 0
            else:
                omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

        dim3 = omega_dims.get(3, 0)
        dim4 = omega_dims.get(4, 0)
        if dim3 == 0 or dim4 == 0:
            continue

        # Full computation
        bd3 = np.zeros((len(ap[2]), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap[2])}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign

        O3 = omega_bases[3]
        d3_om = bd3 @ O3
        U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
        rank_d3 = int(sum(s > 1e-8 for s in S3))
        ker_basis = Vt3[rank_d3:].T
        ker_dim = ker_basis.shape[1]

        if ker_dim == 0:
            continue

        bd4 = np.zeros((len(ap[3]), len(ap[4])))
        idx3 = {path: i for i, path in enumerate(ap[3])}
        for j, path in enumerate(ap[4]):
            for sign, face in boundary_coeffs(path):
                if face in idx3:
                    bd4[idx3[face], j] += sign

        O4 = omega_bases[4]
        d4_om = bd4 @ O4
        O3pinv = np.linalg.pinv(O3)
        d4_omega3 = O3pinv @ d4_om

        U4, S4, Vt4 = np.linalg.svd(d4_omega3, full_matrices=True)
        rank_d4 = int(sum(s > 1e-8 for s in S4))
        beta3 = ker_dim - rank_d4

        if beta3 != 1:
            continue

        # Found a Type B tournament!
        im_basis = U4[:, :rank_d4]

        # Get the generator
        for i in range(ker_dim):
            v = ker_basis[:, i]
            proj = im_basis @ (im_basis.T @ v)
            residual = v - proj
            if np.linalg.norm(residual) > 1e-6:
                gen_omega = residual / np.linalg.norm(residual)
                gen_A3 = O3 @ gen_omega

                print(f"  bits={bits}, scores={scores}")
                print(f"  Generator in A_3 coordinates ({len(ap[3])} total 3-paths):")
                paths_by_vset = defaultdict(list)
                for j, path in enumerate(ap[3]):
                    c = gen_A3[j]
                    if abs(c) > 1e-6:
                        vs = frozenset(path)
                        paths_by_vset[vs].append((path, c))

                print(f"  Organized by vertex set ({len(paths_by_vset)} sets):")
                for vs in sorted(paths_by_vset.keys(), key=lambda s: sorted(s)):
                    paths = paths_by_vset[vs]
                    path_strs = [f"  {'->'.join(map(str,p))}:{c:.4f}" for p, c in paths]
                    print(f"    {sorted(vs)}: {', '.join(path_strs)}")

                # Check: what is the subtournament type on each vertex set?
                print(f"\n  Subtournament analysis:")
                for vs in sorted(paths_by_vset.keys(), key=lambda s: sorted(s)):
                    vlist = sorted(vs)
                    # Count directed 3-cycles in the subtournament
                    c3_sub = 0
                    for i, j, k in combinations(vlist, 3):
                        if A[i][j] and A[j][k] and A[k][i]: c3_sub += 1
                        if A[i][k] and A[k][j] and A[j][i]: c3_sub += 1
                    sub_scores = sorted([sum(A[vi][vj] for vj in vlist if vj != vi) for vi in vlist])
                    paths = paths_by_vset[vs]
                    sum_coeffs = sum(c for p, c in paths)
                    print(f"    {sorted(vs)}: sub_scores={sub_scores}, c3={c3_sub}, #paths={len(paths)}, sum_coeffs={sum_coeffs:.4f}")

                break
        break

    print("\nDONE.")


if __name__ == '__main__':
    main()
