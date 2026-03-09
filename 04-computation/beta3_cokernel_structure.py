"""
beta3_cokernel_structure.py — Algebraic structure of the beta_3 cokernel

When beta_3(T) = 1, we have rank(d_4) = ker(d_3) - 1. This means im(d_4)
has codimension 1 in ker(d_3). What is the "missing" direction?

Questions investigated:
1. What is the cokernel generator (the 1-dim complement of im(d_4) in ker(d_3))?
2. Is it the same for all beta_3=1 tournaments? (S46: NO, varies)
3. What is its support on allowed 3-paths? Is there a pattern?
4. Is there a linear functional on ker(d_3) that characterizes the quotient?
5. Does the cokernel generator have a simple combinatorial description?

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament, enumerate_all_allowed,
    compute_omega_basis_numpy, boundary_faces, build_adj_lists,
    full_chain_complex, deletion_adj
)


def get_cokernel_generator(A, n, data=None):
    """For a tournament with beta_3=1, compute the cokernel generator.

    Returns (cokernel_vec, ker_d3_basis, im_d4_vecs, omega3_paths)
    or None if beta_3 != 1.
    """
    if data is None:
        data = full_chain_complex(A, n, max_p=5)

    b3 = data['bettis'].get(3, 0)
    if b3 != 1:
        return None

    ap = data['ap']
    omega_bases = data['omega_bases']
    omega_dims = data['omega_dims']

    # Get d_3 and d_4 in Omega coordinates
    # d_3: Omega_3 -> Omega_2
    bd3 = np.zeros((len(ap[2]), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap[2])}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3  # len(ap[2]) x dim_Omega_3

    # ker(d_3) in Omega_3 coordinates
    U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
    rank_d3 = int(sum(s > 1e-8 for s in S3))
    ker_d3_basis = Vt3[rank_d3:]  # rows are kernel basis vectors in Omega_3 coords
    ker_dim = omega_dims[3] - rank_d3

    # d_4: Omega_4 -> Omega_3
    if omega_dims.get(4, 0) > 0:
        bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
        idx3 = {path: i for i, path in enumerate(ap[3])}
        for j, path in enumerate(ap.get(4, [])):
            for sign, face in boundary_faces(path):
                if face in idx3:
                    bd4[idx3[face], j] += sign

        O4 = omega_bases[4]
        d4_om = bd4 @ O4  # in A_3 coords
        O3pinv = np.linalg.pinv(O3)
        d4_omega3 = O3pinv @ d4_om  # in Omega_3 coords

        # im(d_4) in ker(d_3) coords
        im_d4_in_ker = ker_d3_basis @ d4_omega3

        # SVD to find what im(d_4) spans in ker(d_3) space
        U4, S4, Vt4 = np.linalg.svd(im_d4_in_ker, full_matrices=True)
        rank_d4_in_ker = int(sum(s > 1e-8 for s in S4))

        # The cokernel = left null space of im_d4_in_ker
        coker_basis = U4[:, rank_d4_in_ker:]  # ker_dim x (ker_dim - rank_d4_in_ker)
    else:
        # No Omega_4 means d_4 = 0, cokernel = all of ker(d_3)
        coker_basis = np.eye(ker_dim)
        d4_omega3 = None

    # cokernel generator in Omega_3 coordinates
    coker_omega3 = coker_basis.T @ ker_d3_basis  # should be 1 x dim_Omega_3

    # Convert to A_3 coordinates (actual 3-path coefficients)
    coker_A3 = coker_omega3 @ O3.T  # 1 x len(ap[3])
    coker_A3 = coker_A3.flatten()

    # Normalize
    norm = np.max(np.abs(coker_A3))
    if norm > 1e-10:
        coker_A3 /= norm

    return {
        'coker_A3': coker_A3,
        'coker_omega3': coker_omega3.flatten(),
        'ker_dim': ker_dim,
        'ap3': ap[3],
        'omega3_dim': omega_dims[3],
        'rank_d3': rank_d3,
        'rank_d4': data['ranks'].get(4, 0),
    }


def analyze_support(coker_A3, paths):
    """Analyze the support of the cokernel generator on 3-paths."""
    nonzero = [(i, coker_A3[i]) for i in range(len(coker_A3)) if abs(coker_A3[i]) > 1e-8]

    # Group by vertex set
    vset_groups = defaultdict(list)
    for i, c in nonzero:
        vset = tuple(sorted(paths[i]))
        vset_groups[vset].append((paths[i], c))

    return nonzero, vset_groups


def main():
    print("=" * 70)
    print("BETA_3 COKERNEL STRUCTURE ANALYSIS")
    print("=" * 70)

    # Part 1: Exhaustive at n=6
    print("\n--- Part 1: Cokernel analysis at n=6 ---")
    n = 6
    total = 2 ** (n*(n-1)//2)

    type_a_cokers = []
    type_b_cokers = []
    type_a_supports = []
    type_b_supports = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue

        result = get_cokernel_generator(A, n, data)
        if result is None:
            continue

        # Classify by Omega_4 dimension (Type A: dim_O4=0, Type B: dim_O4>0)
        if data['omega_dims'].get(4, 0) == 0:
            type_a_cokers.append(result)
            _, vsg = analyze_support(result['coker_A3'], result['ap3'])
            type_a_supports.append(vsg)
        else:
            type_b_cokers.append(result)
            _, vsg = analyze_support(result['coker_A3'], result['ap3'])
            type_b_supports.append(vsg)

        if len(type_a_cokers) + len(type_b_cokers) <= 5:
            print(f"\n  bits={bits}, ker_dim={result['ker_dim']}, omega3_dim={result['omega3_dim']}")
            print(f"  rank(d3)={result['rank_d3']}, rank(d4)={result['rank_d4']}")
            nonzero, vsg = analyze_support(result['coker_A3'], result['ap3'])
            print(f"  Cokernel support: {len(nonzero)} paths, {len(vsg)} vertex sets")
            for vset, pcs in sorted(vsg.items()):
                coeffs = [f"{c:.3f}" for _, c in pcs]
                print(f"    {vset}: {coeffs}")

    print(f"\n  Total beta_3=1: {len(type_a_cokers)+len(type_b_cokers)} ({len(type_a_cokers)} Type A, {len(type_b_cokers)} Type B)")

    # Part 2: Analyze whether cokernel generators are structurally similar
    print("\n--- Part 2: Cokernel similarity analysis ---")

    # For Type B, check coefficient patterns
    print(f"\n  Type B ({len(type_b_cokers)} tournaments):")
    if type_b_cokers:
        # How many distinct coefficient absolute-value patterns?
        patterns = set()
        for r in type_b_cokers:
            c = np.abs(r['coker_A3'])
            c = np.sort(c[c > 1e-8])
            patterns.add(tuple(np.round(c, 4)))
        print(f"  Distinct |coefficient| patterns: {len(patterns)}")

        # Number of nonzero paths in cokernel
        nz_counts = [sum(abs(r['coker_A3']) > 1e-8) for r in type_b_cokers]
        print(f"  Nonzero paths: min={min(nz_counts)}, max={max(nz_counts)}, mean={np.mean(nz_counts):.1f}")

        # Number of vertex sets
        vs_counts = [len(s) for s in type_b_supports]
        print(f"  Vertex sets: min={min(vs_counts)}, max={max(vs_counts)}, mean={np.mean(vs_counts):.1f}")

        # Are coefficients always rational? Check if they're close to simple fractions
        all_coeffs = set()
        for r in type_b_cokers:
            for c in r['coker_A3']:
                if abs(c) > 1e-8:
                    all_coeffs.add(round(abs(c), 6))
        print(f"  Distinct |coefficient| values (rounded): {sorted(all_coeffs)[:20]}")

    # For Type A
    print(f"\n  Type A ({len(type_a_cokers)} tournaments):")
    if type_a_cokers:
        nz_counts = [sum(abs(r['coker_A3']) > 1e-8) for r in type_a_cokers]
        print(f"  Nonzero paths: min={min(nz_counts)}, max={max(nz_counts)}, mean={np.mean(nz_counts):.1f}")
        vs_counts = [len(s) for s in type_a_supports]
        print(f"  Vertex sets: min={min(vs_counts)}, max={max(vs_counts)}, mean={np.mean(vs_counts):.1f}")

    # Part 3: Check if cokernel is related to a known combinatorial object
    print("\n--- Part 3: Cokernel combinatorial description ---")
    # For first few Type B tournaments, show full cokernel
    for idx in range(min(3, len(type_b_cokers))):
        r = type_b_cokers[idx]
        print(f"\n  Type B example #{idx+1}: ker_dim={r['ker_dim']}")
        nonzero, vsg = analyze_support(r['coker_A3'], r['ap3'])
        print(f"  {len(nonzero)} nonzero 3-paths, {len(vsg)} vertex sets:")

        # Check if coefficients are +/-1 after rescaling
        if nonzero:
            max_c = max(abs(c) for _, c in nonzero)
            rescaled = [(i, c/max_c) for i, c in nonzero]
            close_to_pm1 = all(abs(abs(c) - 1.0) < 0.01 for _, c in rescaled)
            print(f"  All +/-1 after rescale? {close_to_pm1}")

            if not close_to_pm1:
                # Show all distinct ratios
                ratios = sorted(set(round(c/max_c, 4) for _, c in rescaled))
                print(f"  Distinct coefficient ratios: {ratios}")

    # Part 4: vertex-set level structure
    print("\n--- Part 4: Vertex-set decomposition ---")
    print("  For each beta_3=1 tournament, how does cokernel decompose by 4-vertex sets?")

    # Each 3-path spans a 4-vertex set. Group coefficients by vertex set.
    for idx in range(min(5, len(type_b_cokers))):
        r = type_b_cokers[idx]
        nonzero, vsg = analyze_support(r['coker_A3'], r['ap3'])
        print(f"\n  Type B #{idx+1}: {len(vsg)} vertex sets")
        for vset in sorted(vsg.keys()):
            pcs = vsg[vset]
            coeffs_sum = sum(c for _, c in pcs)
            print(f"    {vset}: {len(pcs)} paths, sum={coeffs_sum:.4f}")

    # Part 5: Is cokernel a "uniform sum" over something?
    print("\n--- Part 5: Testing if cokernel = uniform sum over 3-cycles ---")
    # The cokernel might be sum_C sign(C) * [cycle_chain_C] for 3-cycles
    for idx in range(min(5, len(type_b_cokers))):
        r = type_b_cokers[idx]
        # Reconstruct the adjacency matrix for this tournament
        # (We need to re-find it... let's just use the first few)
        pass  # Will implement if pattern emerges from Part 3-4

    # Part 6: n=7 cokernel analysis
    print("\n--- Part 6: Cokernel at n=7 (sampled) ---")
    n = 7
    rng = np.random.RandomState(42)
    n7_cokers = []
    checked = 0
    for trial in range(2000):
        A = random_tournament(n, rng)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue
        checked += 1

        result = get_cokernel_generator(A, n, data)
        if result:
            n7_cokers.append(result)

        if checked >= 30:
            break

    print(f"  Found {len(n7_cokers)} beta_3=1 tournaments from {trial+1} total")
    if n7_cokers:
        nz_counts = [sum(abs(r['coker_A3']) > 1e-8) for r in n7_cokers]
        print(f"  Nonzero paths: min={min(nz_counts)}, max={max(nz_counts)}, mean={np.mean(nz_counts):.1f}")

        # Check coefficient structure
        for idx in range(min(3, len(n7_cokers))):
            r = n7_cokers[idx]
            nonzero, vsg = analyze_support(r['coker_A3'], r['ap3'])
            print(f"\n  n=7 example #{idx+1}: ker_dim={r['ker_dim']}, {len(nonzero)} nonzero, {len(vsg)} vertex sets")
            if nonzero:
                max_c = max(abs(c) for _, c in nonzero)
                rescaled = [(i, c/max_c) for i, c in nonzero]
                ratios = sorted(set(round(abs(c), 4) for _, c in rescaled))
                print(f"  |coeff| ratios: {ratios[:15]}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
