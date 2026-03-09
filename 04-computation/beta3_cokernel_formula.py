"""
beta3_cokernel_formula.py — Find a formula for the H_3 generator

From beta3_cokernel_deep.py: 3 * cokernel has integer coefficients in {-3,...,3}.
The factor of 3 suggests involvement of 3-cycles.

Strategy: For each Type B tournament at n=6, extract 3*cokernel (integer vector),
then try to express it as a linear combination of known combinatorial objects:
1. "Cycle chains" — boundaries of 3-cycle homology classes
2. Path reversal signatures
3. Alternating sums over vertex orderings
4. Connection to the Ihara zeta function / L-function

Key test: express 3*coker[path] as a function of local arc structure.

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import permutations
from collections import defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, enumerate_all_allowed, compute_omega_basis_numpy,
    boundary_faces, full_chain_complex
)


def get_integer_cokernel(A, n):
    """Get 3*cokernel as integer vector. Returns None if beta_3 != 1."""
    data = full_chain_complex(A, n, max_p=5)
    if data['bettis'].get(3, 0) != 1:
        return None, None

    ap = data['ap']
    omega_bases = data['omega_bases']
    omega_dims = data['omega_dims']

    bd3 = np.zeros((len(ap[2]), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap[2])}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3
    U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
    rank_d3 = int(sum(s > 1e-8 for s in S3))
    ker_d3_basis = Vt3[rank_d3:]

    if omega_dims.get(4, 0) > 0:
        bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
        idx3 = {path: i for i, path in enumerate(ap[3])}
        for j, path in enumerate(ap.get(4, [])):
            for sign, face in boundary_faces(path):
                if face in idx3:
                    bd4[idx3[face], j] += sign
        O4 = omega_bases[4]
        d4_om = bd4 @ O4
        O3pinv = np.linalg.pinv(O3)
        d4_omega3 = O3pinv @ d4_om
        im_d4_in_ker = ker_d3_basis @ d4_omega3
        U4, S4, Vt4 = np.linalg.svd(im_d4_in_ker, full_matrices=True)
        rank_d4_in_ker = int(sum(s > 1e-8 for s in S4))
        coker_basis = U4[:, rank_d4_in_ker:]
    else:
        coker_basis = np.eye(ker_d3_basis.shape[0])

    coker_omega3 = coker_basis.T @ ker_d3_basis
    coker_A3 = (coker_omega3 @ O3.T).flatten()

    # Normalize so max |coefficient| = 1
    norm = np.max(np.abs(coker_A3))
    if norm > 1e-10:
        coker_A3 /= norm

    # Multiply by 3 and round to integers
    int_coker = np.round(3 * coker_A3).astype(int)
    return int_coker, data['ap'][3]


def path_local_features(path, A, n):
    """Extract local features of a 3-path (a,b,c,d).
    Returns dict of features."""
    a, b, c, d = path
    features = {}

    # Which edges are "forward" (score-increasing)?
    scores = [int(sum(A[i])) for i in range(n)]
    features['score_diff'] = sum(scores[path[i]] - scores[path[i+1]] for i in range(3))

    # How many 3-cycles does each consecutive triple participate in?
    # Triple (a,b,c): is it a 3-cycle? a->b (yes, part of path), b->c (yes), c->a?
    features['abc_cycle'] = int(A[c][a] == 1)  # c->a makes a 3-cycle
    features['bcd_cycle'] = int(A[d][b] == 1)  # d->b makes a 3-cycle

    # "Backward edges" count: edges from later to earlier in path
    features['back_edges'] = sum([
        int(A[b][a] == 1),  # b->a (always 0 since a->b in path)
        int(A[c][a] == 1),
        int(A[c][b] == 1),  # always 0 since b->c
        int(A[d][a] == 1),
        int(A[d][b] == 1),
        int(A[d][c] == 1),  # always 0 since c->d
    ])
    # Since a->b, b->c, c->d are all 1, back_edges counts among {c->a, d->a, d->b}

    # Position of "source" (vertex with max out-degree in subtournament)
    sub_v = [a, b, c, d]
    sub_scores = [sum(A[v][w] for w in sub_v if w != v) for v in sub_v]
    features['max_score_pos'] = sub_scores.index(max(sub_scores))
    features['min_score_pos'] = sub_scores.index(min(sub_scores))

    # Number of directed 3-cycles in the subtournament
    c3 = 0
    for i in range(4):
        for j in range(i+1, 4):
            for k in range(j+1, 4):
                vi, vj, vk = sub_v[i], sub_v[j], sub_v[k]
                if (A[vi][vj] and A[vj][vk] and A[vk][vi]) or \
                   (A[vi][vk] and A[vk][vj] and A[vj][vi]):
                    c3 += 1
    features['c3'] = c3

    return features


def main():
    print("=" * 70)
    print("BETA_3 COKERNEL FORMULA SEARCH")
    print("=" * 70)

    n = 6
    total = 2 ** (n*(n-1)//2)

    # Part 1: Collect all Type B cokernel data
    print("\n--- Part 1: Integer cokernel vs local features ---")

    all_path_data = []  # (path, coeff, features, bits)
    count_b = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        int_coker, ap3 = get_integer_cokernel(A, n)
        if int_coker is None:
            continue

        # Check if Type B (has Omega_4)
        data = full_chain_complex(A, n, max_p=5)
        if data['omega_dims'].get(4, 0) == 0:
            continue  # Skip Type A
        count_b += 1

        for i, path in enumerate(ap3):
            coeff = int_coker[i]
            feat = path_local_features(path, A, n)
            all_path_data.append((path, coeff, feat, bits))

    print(f"  Collected {len(all_path_data)} path-coefficient pairs from {count_b} Type B tournaments")

    # Part 2: Check which feature combinations determine the coefficient
    print("\n--- Part 2: Coefficient vs features ---")

    # Group by (c3, abc_cycle, bcd_cycle, back_edges)
    groups = defaultdict(list)
    for path, coeff, feat, bits in all_path_data:
        key = (feat['c3'], feat['abc_cycle'], feat['bcd_cycle'], feat['back_edges'])
        groups[key].append(coeff)

    print("  (c3, abc_cycle, bcd_cycle, back_edges) -> coefficient distribution:")
    for key in sorted(groups.keys()):
        from collections import Counter
        dist = Counter(groups[key])
        print(f"    {key}: {dict(sorted(dist.items()))}")

    # Part 3: Try finer grouping
    print("\n--- Part 3: Finer features ---")
    groups2 = defaultdict(list)
    for path, coeff, feat, bits in all_path_data:
        key = (feat['c3'], feat['abc_cycle'], feat['bcd_cycle'])
        groups2[key].append(coeff)

    print("  (c3, abc_cycle, bcd_cycle) -> coefficient distribution:")
    for key in sorted(groups2.keys()):
        from collections import Counter
        dist = Counter(groups2[key])
        total = sum(dist.values())
        print(f"    {key}: {dict(sorted(dist.items()))} (n={total})")

    # Part 4: Test specific formulas
    print("\n--- Part 4: Testing specific formulas ---")

    # Formula 1: coeff = abc_cycle + bcd_cycle - 1 (shifted)
    # Formula 2: coeff = 2*abc_cycle + bcd_cycle - 1
    # Formula 3: coeff related to back_edges
    # Formula 4: coeff = (-1)^{position of sink} * (...)

    formulas = {
        'abc+bcd-1': lambda f: f['abc_cycle'] + f['bcd_cycle'] - 1,
        '2abc+bcd-1': lambda f: 2*f['abc_cycle'] + f['bcd_cycle'] - 1,
        'abc-bcd': lambda f: f['abc_cycle'] - f['bcd_cycle'],
        'back-2': lambda f: f['back_edges'] - 2,
    }

    for name, formula in formulas.items():
        matches = 0
        total_count = 0
        for path, coeff, feat, bits in all_path_data:
            pred = formula(feat)
            if pred == coeff:
                matches += 1
            total_count += 1
        print(f"  {name}: {matches}/{total_count} = {100*matches/total_count:.1f}%")

    # Part 5: For the nonzero coefficients, is there a clean relationship?
    print("\n--- Part 5: Nonzero coefficient analysis ---")
    nonzero_data = [(p, c, f, b) for p, c, f, b in all_path_data if c != 0]
    print(f"  {len(nonzero_data)} nonzero out of {len(all_path_data)} total")

    # For nonzero paths: what is the relationship between coeff and features?
    groups3 = defaultdict(list)
    for path, coeff, feat, bits in nonzero_data:
        key = (feat['c3'], feat['abc_cycle'], feat['bcd_cycle'], feat['back_edges'])
        groups3[key].append(coeff)

    print("  Nonzero (c3, abc_cycle, bcd_cycle, back_edges) -> coefficient:")
    for key in sorted(groups3.keys()):
        from collections import Counter
        dist = Counter(groups3[key])
        print(f"    {key}: {dict(sorted(dist.items()))}")

    # Part 6: Sign analysis - does the global tournament structure determine signs?
    print("\n--- Part 6: Sign determination ---")
    # For each tournament, check if sign pattern depends on tournament orientation
    # relative to the score ordering
    sign_patterns = defaultdict(int)
    for bits in range(min(total, 10000)):
        A = bits_to_adj(bits, n)
        int_coker, ap3 = get_integer_cokernel(A, n)
        if int_coker is None:
            continue
        data = full_chain_complex(A, n, max_p=5)
        if data['omega_dims'].get(4, 0) == 0:
            continue

        # Pattern: (number of +3, number of +2, number of +1, ...)
        vals = tuple(sorted(int_coker[int_coker != 0]))
        sign_patterns[vals] += 1

    print("  Distinct sign patterns (sorted nonzero coefficients):")
    for pattern, count in sorted(sign_patterns.items(), key=lambda x: -x[1])[:20]:
        print(f"    {pattern}: {count}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
