"""
h3_generator.py — When beta_3 = 1, what does the H_3 generator look like?

If beta_3 = 1, there's exactly one independent 3-cycle (in Omega_3) that is
in ker(d_3) but not in im(d_4). What structure does this cycle have?

Understanding the generator could explain:
1. Why beta_3 in {0,1} (unique generator type)
2. Why beta_1*beta_3 = 0 (generator incompatible with beta_1=1)
3. How to prove the Boolean property

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import defaultdict, Counter
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

def compute_h3_generator(A, n):
    """If beta_3 > 0, find and return the H_3 generator."""
    ap = {}
    for p in range(6):  # Need up to p=5 for im(d_5)
        if p <= n - 1:
            ap[p] = enumerate_allowed_paths(A, n, p)
        else:
            ap[p] = []

    # Compute Omega bases for p=2,3,4
    omega_bases = {}
    omega_dims = {}
    for p in range(6):
        if p > n - 1 or len(ap[p]) == 0:
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
            continue
        if p == 0:
            omega_bases[p] = np.eye(len(ap[p]))
            omega_dims[p] = len(ap[p])
            continue

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
            omega_bases[p] = np.eye(len(ap[p]))
            omega_dims[p] = len(ap[p])
        else:
            P = np.zeros((na_count, len(ap[p])))
            for j, path in enumerate(ap[p]):
                for sign, face in boundary_coeffs(path):
                    if face in non_allowed:
                        P[non_allowed[face], j] += sign
            U, S, Vt = np.linalg.svd(P, full_matrices=True)
            rank = int(sum(s > 1e-10 for s in S))
            ns = Vt[rank:].T
            if ns.shape[1] > 0:
                omega_bases[p] = ns
                omega_dims[p] = ns.shape[1]
            else:
                omega_bases[p] = np.zeros((len(ap[p]), 0))
                omega_dims[p] = 0

    dim_3 = omega_dims.get(3, 0)
    if dim_3 == 0:
        return None, 0, None

    # Build d_3: Omega_3 -> Omega_2
    bd3 = np.zeros((len(ap[2]), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap[2])}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_coeffs(path):
            if face in idx2:
                bd3[idx2[face], j] += sign
    d3_om = bd3 @ omega_bases[3]
    U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
    rank_d3 = int(sum(s > 1e-8 for s in S3))
    ker_d3_basis = Vt3[rank_d3:].T  # Columns are ker(d_3) basis vectors
    ker_d3_dim = ker_d3_basis.shape[1] if ker_d3_basis.ndim == 2 else 0

    if ker_d3_dim == 0:
        return None, 0, None

    # Build d_4: Omega_4 -> Omega_3
    dim_4 = omega_dims.get(4, 0)
    if dim_4 == 0:
        # im(d_4) = 0, so beta_3 = ker_d3_dim
        # Express generator in A_3 basis
        gen_omega = ker_d3_basis[:, 0]  # First kernel vector (in Omega_3 basis)
        gen_A3 = omega_bases[3] @ gen_omega  # In A_3 basis
        return gen_A3, ker_d3_dim, ap[3]

    bd4 = np.zeros((len(ap[3]), len(ap[4])))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap[4]):
        for sign, face in boundary_coeffs(path):
            if face in idx3:
                bd4[idx3[face], j] += sign
    d4_om = bd4 @ omega_bases[4]

    # im(d_4) in A_3 coordinates: columns of d4_om span im(d_4) in A_3 space
    # Convert to Omega_3 coordinates using pseudoinverse of omega_bases[3]
    # omega_bases[3] maps Omega_3 coords -> A_3 coords
    # To go back: c = (O3^T O3)^{-1} O3^T v = O3^+ v
    O3 = omega_bases[3]  # shape: (len(ap[3]), dim_3)
    O3pinv = np.linalg.pinv(O3)  # shape: (dim_3, len(ap[3]))

    # d4 in Omega_3 coords: O3pinv @ d4_om maps Omega_4 coords -> Omega_3 coords
    d4_omega3 = O3pinv @ d4_om
    sv4 = np.linalg.svd(d4_omega3, compute_uv=False)
    rank_d4 = int(sum(s > 1e-8 for s in sv4))

    beta_3 = ker_d3_dim - rank_d4
    if beta_3 <= 0:
        return None, 0, None

    # Find the H_3 generator: a vector in ker(d_3) not in im(d_4)
    # Both in Omega_3 coordinates
    U4, S4, Vt4 = np.linalg.svd(d4_omega3, full_matrices=True)
    im_d4_basis = U4[:, :rank_d4]  # Image basis in Omega_3 coordinates

    # Simple approach: take first ker(d_3) basis vector that's not in im(d_4)
    for i in range(ker_d3_dim):
        v = ker_d3_basis[:, i]
        if rank_d4 > 0:
            # Project v onto im(d_4) and check if residual is nonzero
            proj = im_d4_basis @ (im_d4_basis.T @ v)
            residual = v - proj
            if np.linalg.norm(residual) > 1e-6:
                gen_omega = residual / np.linalg.norm(residual)
                gen_A3 = omega_bases[3] @ gen_omega
                return gen_A3, beta_3, ap[3]
        else:
            gen_omega = v
            gen_A3 = omega_bases[3] @ gen_omega
            return gen_A3, beta_3, ap[3]

    return None, beta_3, ap[3]


def main():
    print("=" * 70)
    print("H_3 GENERATOR ANALYSIS")
    print("=" * 70)

    # Part 1: n=6 exhaustive, find all beta_3=1 tournaments and their generators
    print("\n--- Part 1: n=6 exhaustive ---")
    n = 6
    N = 2**(n*(n-1)//2)
    generators = []

    for bits in range(N):
        A = bits_to_adj(bits, n)
        gen, beta3, a3 = compute_h3_generator(A, n)
        if beta3 > 0 and gen is not None:
            # Find which 3-paths have nonzero coefficient
            nonzero_paths = []
            for i, coeff in enumerate(gen):
                if abs(coeff) > 1e-6:
                    nonzero_paths.append((a3[i], coeff))

            # What vertex sets are involved?
            vertex_sets = set()
            for path, coeff in nonzero_paths:
                vertex_sets.add(frozenset(path))

            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            generators.append({
                'bits': bits,
                'scores': scores,
                'num_paths': len(nonzero_paths),
                'num_vsets': len(vertex_sets),
                'vertex_sets': vertex_sets,
                'paths': nonzero_paths[:10],
                'beta3': beta3
            })

    print(f"  Found {len(generators)} tournaments with beta_3 > 0")

    # Analyze generator structure
    path_counts = Counter(g['num_paths'] for g in generators)
    vset_counts = Counter(g['num_vsets'] for g in generators)
    score_counts = Counter(g['scores'] for g in generators)

    print(f"\n  Number of nonzero 3-paths in generator:")
    for k in sorted(path_counts.keys()):
        print(f"    {k} paths: {path_counts[k]} tournaments")

    print(f"\n  Number of vertex sets in generator:")
    for k in sorted(vset_counts.keys()):
        print(f"    {k} sets: {vset_counts[k]} tournaments")

    print(f"\n  Score sequences:")
    for sc, cnt in score_counts.most_common():
        print(f"    {sc}: {cnt}")

    # Print a few generators in detail
    print(f"\n  Detailed generators (first 5):")
    for g in generators[:5]:
        print(f"\n    bits={g['bits']}, scores={g['scores']}, "
              f"{g['num_paths']} paths, {g['num_vsets']} vertex sets")
        for path, coeff in g['paths']:
            vset = frozenset(path)
            print(f"      {path}: coeff={coeff:.4f}, vertices={sorted(vset)}")

    # Part 2: Do all generators use the SAME vertex sets?
    print("\n--- Part 2: Vertex set structure of generators ---")
    vset_sizes = Counter()
    for g in generators:
        for vs in g['vertex_sets']:
            vset_sizes[len(vs)] += 1

    print(f"  Vertex set sizes in generators:")
    for size in sorted(vset_sizes.keys()):
        print(f"    size {size}: {vset_sizes[size]} occurrences")

    # Part 3: n=7 sampled
    print("\n--- Part 3: n=7 sampled ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 500
    gen7_data = []

    for trial in range(N):
        A = random_tournament(n, rng)
        gen, beta3, a3 = compute_h3_generator(A, n)
        if beta3 > 0 and gen is not None:
            nonzero_paths = []
            for i, coeff in enumerate(gen):
                if abs(coeff) > 1e-6:
                    nonzero_paths.append((a3[i], coeff))
            vertex_sets = set()
            for path, coeff in nonzero_paths:
                vertex_sets.add(frozenset(path))
            gen7_data.append({
                'trial': trial,
                'num_paths': len(nonzero_paths),
                'num_vsets': len(vertex_sets),
                'beta3': beta3
            })

        if (trial + 1) % 100 == 0:
            print(f"  {trial+1}/{N}: {len(gen7_data)} beta_3>0 found", flush=True)

    if gen7_data:
        path_counts7 = Counter(g['num_paths'] for g in gen7_data)
        vset_counts7 = Counter(g['num_vsets'] for g in gen7_data)
        print(f"\n  n=7: {len(gen7_data)} generators found")
        print(f"  Nonzero path counts: {dict(sorted(path_counts7.items()))}")
        print(f"  Vertex set counts: {dict(sorted(vset_counts7.items()))}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
