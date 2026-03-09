"""
beta3_cokernel_deep.py — Deep structure of the universal cokernel pattern

Key finding from beta3_cokernel_structure.py:
- All 240 Type B beta_3=1 tournaments at n=6 have the SAME |coefficient| pattern
- Coefficients are ±{1/3, 2/3, 1}
- Vertex sets contribute 2 or 3 paths each (6 with 3, 9 with 2; 6*3+9*2=36)

Questions:
1. What distinguishes 3-path vs 2-path vertex sets? (Subtournament type? 3-cycle content?)
2. Are the coefficients 1/3, 2/3, 1 related to the number of 3-cycles?
3. What determines the sign? Can we express it via the tournament adjacency?
4. Is the cokernel = sum of elementary chains weighted by some combinatorial invariant?
5. If we multiply by 3, do we get integer coefficients with a clean formula?

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament, enumerate_all_allowed,
    compute_omega_basis_numpy, boundary_faces,
    full_chain_complex, count_3cycles
)


def get_cokernel_A3(A, n, data=None):
    """Get cokernel generator in A_3 (3-path) coordinates. Returns None or vector."""
    if data is None:
        data = full_chain_complex(A, n, max_p=5)
    if data['bettis'].get(3, 0) != 1:
        return None

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

    norm = np.max(np.abs(coker_A3))
    if norm > 1e-10:
        coker_A3 /= norm
    return coker_A3


def subtournament_type(A, vset):
    """Classify 4-vertex subtournament by c3 count and score."""
    vlist = sorted(vset)
    c3 = 0
    for i in range(4):
        for j in range(i+1, 4):
            for k in range(j+1, 4):
                a, b, c = vlist[i], vlist[j], vlist[k]
                if (A[a][b] == 1 and A[b][c] == 1 and A[c][a] == 1) or \
                   (A[a][c] == 1 and A[c][b] == 1 and A[b][a] == 1):
                    c3 += 1
    scores = sorted([sum(A[vlist[i]][vlist[j]] for j in range(4) if i != j) for i in range(4)])
    return c3, tuple(scores)


def count_ham_paths(A, vset):
    """Count Hamiltonian paths in induced subtournament."""
    vlist = sorted(vset)
    count = 0
    for perm in permutations(vlist):
        ok = True
        for i in range(len(perm)-1):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count


def main():
    print("=" * 70)
    print("BETA_3 COKERNEL DEEP STRUCTURE")
    print("=" * 70)

    n = 6
    total = 2 ** (n*(n-1)//2)

    # Part 1: Relate vertex set properties to cokernel coefficients
    print("\n--- Part 1: Vertex set properties vs cokernel at n=6 (Type B) ---")

    type_b_data = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue
        if data['omega_dims'].get(4, 0) == 0:
            continue  # Skip Type A

        coker = get_cokernel_A3(A, n, data)
        if coker is None:
            continue

        ap3 = data['ap'][3]

        # For each 4-vertex set, get subtournament type and cokernel info
        vset_info = {}
        for i, path in enumerate(ap3):
            vset = tuple(sorted(path))
            if vset not in vset_info:
                c3, scores = subtournament_type(A, vset)
                nhp = count_ham_paths(A, vset)
                vset_info[vset] = {
                    'c3': c3, 'scores': scores, 'nhp': nhp,
                    'paths': [], 'coeffs': []
                }
            vset_info[vset]['paths'].append(path)
            vset_info[vset]['coeffs'].append(coker[i])

        type_b_data.append({
            'bits': bits, 'A': A, 'vset_info': vset_info, 'coker': coker
        })

        if len(type_b_data) == 1:
            print(f"\n  First Type B example (bits={bits}):")
            for vset in sorted(vset_info.keys()):
                info = vset_info[vset]
                print(f"    {vset}: c3={info['c3']}, scores={info['scores']}, "
                      f"HP={info['nhp']}, #paths={len(info['paths'])}")
                for p, c in zip(info['paths'], info['coeffs']):
                    print(f"      {p}: coeff={c:.4f} (x3={3*c:.4f})")

    print(f"\n  Total Type B: {len(type_b_data)}")

    # Part 2: Is the number of paths per vset determined by subtournament type?
    print("\n--- Part 2: Paths per vertex set vs subtournament type ---")
    paths_by_type = defaultdict(list)
    for d in type_b_data:
        for vset, info in d['vset_info'].items():
            key = (info['c3'], info['scores'])
            paths_by_type[key].append(len(info['paths']))

    for key in sorted(paths_by_type.keys()):
        vals = Counter(paths_by_type[key])
        print(f"  c3={key[0]}, scores={key[1]}: paths per vset = {dict(vals)}")

    # Part 3: Integer structure after multiplying by 3
    print("\n--- Part 3: 3x cokernel — integer structure ---")
    for idx, d in enumerate(type_b_data[:5]):
        coker3 = 3 * d['coker']
        rounded = np.round(coker3)
        is_int = np.allclose(coker3, rounded, atol=0.01)
        vals = sorted(set(int(round(x)) for x in coker3 if abs(x) > 0.01))
        print(f"  Example {idx}: 3*coker integer? {is_int}, values: {vals}")

    # Part 4: Sign pattern — what determines the sign?
    print("\n--- Part 4: Sign pattern analysis ---")
    # For each Type B tournament, check if sign correlates with path direction
    # relative to score ordering
    for idx in range(min(3, len(type_b_data))):
        d = type_b_data[idx]
        A = d['A']
        scores = [int(sum(A[i])) for i in range(n)]
        print(f"\n  Example {idx} (bits={d['bits']}), scores={scores}:")

        coker = d['coker']
        ap3 = list(d['vset_info'].values())[0]['paths']  # just grab the paths from data

        # For each path, classify by:
        # - number of "forward" edges (i->j where score(i) > score(j))
        # - whether path visits vertices in score-increasing order
        pos_paths = []
        neg_paths = []
        for vset in sorted(d['vset_info'].keys()):
            info = d['vset_info'][vset]
            for p, c in zip(info['paths'], info['coeffs']):
                fwd = sum(1 for k in range(3) if scores[p[k]] > scores[p[k+1]])
                if c > 0.01:
                    pos_paths.append((p, c, fwd, info['c3']))
                elif c < -0.01:
                    neg_paths.append((p, c, fwd, info['c3']))

        print(f"  Positive coefficients ({len(pos_paths)}):")
        fwd_dist = Counter(p[2] for p in pos_paths)
        c3_dist = Counter(p[3] for p in pos_paths)
        print(f"    Forward-edge dist: {dict(fwd_dist)}, c3 dist: {dict(c3_dist)}")

        print(f"  Negative coefficients ({len(neg_paths)}):")
        fwd_dist = Counter(p[2] for p in neg_paths)
        c3_dist = Counter(p[3] for p in neg_paths)
        print(f"    Forward-edge dist: {dict(fwd_dist)}, c3 dist: {dict(c3_dist)}")

    # Part 5: Is cokernel related to "cycle chain"?
    print("\n--- Part 5: Testing cycle chain hypothesis ---")
    # Hypothesis: the cokernel generator is the sum over 3-cycles C of
    # alpha(C) * (boundary chain of C), where alpha depends on C
    # A 3-cycle C = (a,b,c) with a->b->c->a has boundary chain:
    # (b,c) - (a,c) + (a,b) in dimension 1
    # But we need dimension 3, so this doesn't directly apply.

    # Alternative: cokernel is related to the Euler chain or some global invariant
    # Let's test: is cokernel in the kernel of ALL boundary maps?
    # (It should be, since it's in H_3 = ker(d_3)/im(d_4))
    for d in type_b_data[:1]:
        A = d['A']
        data = full_chain_complex(A, n, max_p=5)
        ap = data['ap']
        coker = d['coker']

        # Check d_3(coker) = 0 in A_2 coordinates
        bd3 = np.zeros((len(ap[2]), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap[2])}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_faces(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign

        d3_coker = bd3 @ coker
        print(f"  |d_3(coker)| = {np.max(np.abs(d3_coker)):.2e}")

    # Part 6: Count how many 3-paths in Omega_3 per vertex set
    print("\n--- Part 6: Omega_3 paths per vertex set (Type B) ---")
    d = type_b_data[0]
    A = d['A']
    data = full_chain_complex(A, n, max_p=5)
    ap3 = data['ap'][3]

    # How many allowed 3-paths per vertex set?
    ap3_by_vset = defaultdict(list)
    for p in ap3:
        ap3_by_vset[tuple(sorted(p))].append(p)

    # But Omega_3 is a SUBSPACE, not all of A_3. How many Omega basis vectors
    # have support on each vertex set?
    O3 = data['omega_bases'][3]
    print(f"  dim(Omega_3) = {data['omega_dims'][3]}, dim(A_3) = {len(ap3)}")
    print(f"  Number of allowed 3-paths per vertex set:")
    for vset in sorted(ap3_by_vset.keys()):
        c3, scores = subtournament_type(A, vset)
        nhp = count_ham_paths(A, vset)
        print(f"    {vset}: {len(ap3_by_vset[vset])} paths, c3={c3}, scores={scores}, HP={nhp}")

    # Part 7: Universal sign structure? Check if all Type B tournaments have same
    # sign pattern after identifying vertex sets by subtournament type
    print("\n--- Part 7: Sign universality check ---")
    # For each Type B tournament, compute the sum of cokernel over each vertex set
    sum_by_c3 = defaultdict(list)
    for d in type_b_data:
        for vset, info in d['vset_info'].items():
            s = sum(info['coeffs'])
            sum_by_c3[info['c3']].append(abs(round(s, 4)))

    for c3_val in sorted(sum_by_c3.keys()):
        vals = Counter(sum_by_c3[c3_val])
        print(f"  c3={c3_val}: |sum| distribution = {dict(vals)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
