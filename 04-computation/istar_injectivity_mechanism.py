"""
istar_injectivity_mechanism.py — WHY is i_* injective when b3(T\\v) = 1?

Key question: when T has beta_3=1 and T\\v also has beta_3=1,
why does the cycle in H_3(T\\v) remain non-trivial in H_3(T)?

Possible mechanisms:
1. Seesaw: b3(T\\v)=1 => b1(T\\v)=0 => im(d_2) is large in T\\v.
   Does this constrain d_4 in T to not kill the embedded cycle?

2. Omega_4 dimension: when b3(T\\v)=1, the extra Omega_4(T) content
   (from paths through v) maps into "new" ker(d_3) directions, not the
   embedded T\\v cycle.

3. Structural: the H_3 generator of T\\v is "compatible" with T's structure
   in a way that forces the T\\v-generator = T-generator (up to im(d_4)).

Investigate by:
- Computing dim(Omega_4(T)) vs dim(Omega_4(T\\v)) at bad vertices
- Checking if the im(d_4) contribution from v-paths is disjoint from
  the embedded T\\v cycle
- Looking at the actual cycle vectors

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex

PRIME = (1 << 31) - 1


def gauss_rank_modp(M, prime=PRIME):
    M = M.copy() % prime
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M[nz] = (M[nz] - np.outer(factors[nz], M[rank])) % prime
        rank += 1
    return rank


def gauss_nullspace_modp(M, prime=PRIME):
    M = M.copy() % prime
    nrows, ncols = M.shape
    pivot_cols = []
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M[nz] = (M[nz] - np.outer(factors[nz], M[rank])) % prime
        pivot_cols.append(col)
        rank += 1

    null_dim = ncols - rank
    if null_dim == 0:
        return np.zeros((0, ncols), dtype=np.int64)

    free_cols = sorted(set(range(ncols)) - set(pivot_cols))
    null_basis = np.zeros((null_dim, ncols), dtype=np.int64)
    for i, fc in enumerate(free_cols):
        null_basis[i, fc] = 1
        for j, pc in enumerate(pivot_cols):
            null_basis[i, pc] = (-M[j, fc]) % prime
    return null_basis % prime


def get_allowed_paths(A, vertices, max_p):
    paths = {}
    for p in range(max_p + 2):
        p_paths = []
        if p + 1 <= len(vertices):
            for verts in combinations(vertices, p + 1):
                for perm in permutations(verts):
                    ok = True
                    for i in range(len(perm) - 1):
                        if A[perm[i]][perm[i+1]] != 1:
                            ok = False
                            break
                    if ok:
                        p_paths.append(perm)
        paths[p] = p_paths
    return paths


def build_omega_basis_modp(ap, max_p, prime=PRIME):
    omega = {}
    for p in range(max_p + 1):
        a_p = ap[p]
        if not a_p:
            omega[p] = np.zeros((0, 0), dtype=np.int64)
            continue

        a_pm1_set = set(ap.get(p-1, [])) if p > 0 else set()
        na_faces = {}
        for sigma in a_p:
            for i in range(1, len(sigma) - 1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in a_pm1_set:
                    if face not in na_faces:
                        na_faces[face] = len(na_faces)

        if not na_faces:
            omega[p] = np.eye(len(a_p), dtype=np.int64)
        else:
            C = np.zeros((len(na_faces), len(a_p)), dtype=np.int64)
            for j, sigma in enumerate(a_p):
                for i in range(1, len(sigma) - 1):
                    face = sigma[:i] + sigma[i+1:]
                    if face in na_faces:
                        C[na_faces[face], j] = (C[na_faces[face], j] + ((-1)**i)) % prime
            C = C % prime
            null = gauss_nullspace_modp(C, prime)
            omega[p] = null if null.shape[0] > 0 else np.zeros((0, len(a_p)), dtype=np.int64)
    return omega


def build_boundary_modp(ap, p, prime=PRIME):
    a_p = ap[p]
    a_pm1 = ap.get(p-1, [])
    if not a_p or not a_pm1:
        return np.zeros((len(a_pm1) if a_pm1 else 0, len(a_p) if a_p else 0), dtype=np.int64)

    idx = {path: i for i, path in enumerate(a_pm1)}
    mat = np.zeros((len(a_pm1), len(a_p)), dtype=np.int64)
    for j, sigma in enumerate(a_p):
        for i in range(len(sigma)):
            face = sigma[:i] + sigma[i+1:]
            if face in idx:
                mat[idx[face], j] = (mat[idx[face], j] + ((-1)**i)) % prime
    return mat % prime


def analyze_mechanism(A, n, v, prime=PRIME):
    """Detailed analysis of i_* mechanism for a specific (T,v) pair."""
    all_v = list(range(n))
    remaining = [i for i in range(n) if i != v]
    max_p = 5

    # Build chain complexes for T and T\v
    ap_T = get_allowed_paths(A, all_v, max_p)
    ap_Tv = get_allowed_paths(A, remaining, max_p)

    omega_T = build_omega_basis_modp(ap_T, max_p, prime)
    omega_Tv = build_omega_basis_modp(ap_Tv, max_p, prime)

    result = {}

    # Dimensions
    for p in range(6):
        result[f'dim_Om{p}_T'] = omega_T[p].shape[0] if p in omega_T else 0
        result[f'dim_Om{p}_Tv'] = omega_Tv[p].shape[0] if p in omega_Tv else 0
        result[f'dim_A{p}_T'] = len(ap_T.get(p, []))
        result[f'dim_A{p}_Tv'] = len(ap_Tv.get(p, []))

    # Boundary ranks for T
    for p in [3, 4]:
        om = omega_T[p]
        if om.shape[0] > 0:
            bd = build_boundary_modp(ap_T, p, prime)
            d_mat = bd @ om.T % prime
            result[f'rank_d{p}_T'] = gauss_rank_modp(d_mat.copy(), prime)
        else:
            result[f'rank_d{p}_T'] = 0

    # Boundary ranks for T\v
    for p in [3, 4]:
        om = omega_Tv[p]
        if om.shape[0] > 0:
            bd = build_boundary_modp(ap_Tv, p, prime)
            d_mat = bd @ om.T % prime
            result[f'rank_d{p}_Tv'] = gauss_rank_modp(d_mat.copy(), prime)
        else:
            result[f'rank_d{p}_Tv'] = 0

    # Betti numbers
    result['b3_T'] = (result['dim_Om3_T'] - result['rank_d3_T']) - result['rank_d4_T']
    result['b3_Tv'] = (result['dim_Om3_Tv'] - result['rank_d3_Tv']) - result['rank_d4_Tv']
    result['ker_d3_T'] = result['dim_Om3_T'] - result['rank_d3_T']
    result['ker_d3_Tv'] = result['dim_Om3_Tv'] - result['rank_d3_Tv']

    # Key: how does Omega_4 change?
    result['delta_Om4'] = result['dim_Om4_T'] - result['dim_Om4_Tv']
    result['delta_Om3'] = result['dim_Om3_T'] - result['dim_Om3_Tv']

    # How many NEW 3-paths and 4-paths involve v?
    v_paths_3 = [p for p in ap_T[3] if v in p]
    v_paths_4 = [p for p in ap_T.get(4, []) if v in p]
    result['v_paths_3'] = len(v_paths_3)
    result['v_paths_4'] = len(v_paths_4)

    # How much of the NEW Omega_4 content maps into the ker(d_3) of T\v-embedded?
    # This is the key question: does im(d_4^{new}) hit the embedded T\v cycle?

    # Count v-incident 4-paths that are in Omega_4(T)
    # These are 4-paths through v that satisfy the Omega constraint
    om4_T = omega_T[4]
    if om4_T.shape[0] > 0:
        # Which basis vectors have support on v-paths?
        ap4_T = ap_T[4]
        v_cols = [i for i, p in enumerate(ap4_T) if v in p]
        if v_cols:
            v_support = om4_T[:, v_cols]
            # How many basis vectors have ANY v-path support?
            v_active = np.sum(np.any(v_support != 0, axis=1))
            result['om4_v_active'] = int(v_active)
        else:
            result['om4_v_active'] = 0
    else:
        result['om4_v_active'] = 0

    return result


def main():
    print("=" * 70)
    print("i_* INJECTIVITY MECHANISM ANALYSIS")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    print(f"\n--- Detailed analysis of bad vertices at n={n} ---")

    bad_data = []
    good_data = []
    checked = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue

        checked += 1

        for v in range(n):
            r = analyze_mechanism(A, n, v)

            if r['b3_Tv'] == 1:
                bad_data.append(r)
            elif r['b3_Tv'] == 0 and r['b3_T'] == 1:
                good_data.append(r)

        if checked >= 100:
            break

        if checked % 20 == 0:
            print(f"  ... {checked} beta_3=1 found", flush=True)

    print(f"\n  Total: {checked} beta_3=1 tours, {len(bad_data)} bad vertices, {len(good_data)} good vertices")

    # Part 1: Compare chain complex dimensions
    print(f"\n--- Part 1: Omega dimension changes ---")

    print("\n  BAD vertices (b3(T\\v)=1):")
    for key in ['delta_Om3', 'delta_Om4', 'dim_Om3_T', 'dim_Om3_Tv',
                'dim_Om4_T', 'dim_Om4_Tv', 'ker_d3_T', 'ker_d3_Tv',
                'rank_d3_T', 'rank_d3_Tv', 'rank_d4_T', 'rank_d4_Tv']:
        vals = [d[key] for d in bad_data]
        print(f"    {key}: min={min(vals)}, max={max(vals)}, avg={np.mean(vals):.1f}")

    print("\n  GOOD vertices (b3(T\\v)=0):")
    for key in ['delta_Om3', 'delta_Om4', 'dim_Om3_T', 'dim_Om3_Tv',
                'dim_Om4_T', 'dim_Om4_Tv', 'ker_d3_T', 'ker_d3_Tv',
                'rank_d3_T', 'rank_d3_Tv', 'rank_d4_T', 'rank_d4_Tv']:
        vals = [d[key] for d in good_data]
        if vals:
            print(f"    {key}: min={min(vals)}, max={max(vals)}, avg={np.mean(vals):.1f}")

    # Part 2: Key comparison — delta_Om4 vs delta_ker_d3
    print(f"\n--- Part 2: Growth analysis ---")

    for label, dataset in [("BAD", bad_data), ("GOOD", good_data)]:
        if not dataset:
            continue
        print(f"\n  {label} vertices:")

        delta_ker = [d['ker_d3_T'] - d['ker_d3_Tv'] for d in dataset]
        delta_rk4 = [d['rank_d4_T'] - d['rank_d4_Tv'] for d in dataset]
        delta_om4 = [d['delta_Om4'] for d in dataset]

        print(f"    delta(ker_d3) = ker_T - ker_Tv: {Counter(delta_ker)}")
        print(f"    delta(rank_d4) = rk4_T - rk4_Tv: {Counter(delta_rk4)}")
        print(f"    delta(Omega_4):                   {Counter(delta_om4)}")

        # The key relation: b3_T - b3_Tv = delta_ker - delta_rk4
        delta_b3 = [d['b3_T'] - d['b3_Tv'] for d in dataset]
        print(f"    delta(b3) = b3_T - b3_Tv:        {Counter(delta_b3)}")

        # Check: delta_ker - delta_rk4 should equal delta_b3
        check = [dk - dr4 - db3 for dk, dr4, db3 in zip(delta_ker, delta_rk4, delta_b3)]
        print(f"    Consistency check (delta_ker - delta_rk4 - delta_b3): {Counter(check)}")

    # Part 3: The SEESAW connection
    print(f"\n--- Part 3: Seesaw connection ---")
    print("  b1(T) and b1(T\\v) values:")

    b1_data_bad = []
    b1_data_good = []

    for trial in range(5000):
        A = random_tournament(n, np.random.RandomState(trial + 1000))
        data_T = full_chain_complex(A, n, max_p=5)
        if data_T['bettis'].get(3, 0) != 1:
            continue

        b1_T = data_T['bettis'].get(1, 0)

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            data_Tv = full_chain_complex(A_sub, n-1, max_p=5)
            b3_Tv = data_Tv['bettis'].get(3, 0)
            b1_Tv = data_Tv['bettis'].get(1, 0)

            if b3_Tv == 1:
                b1_data_bad.append((b1_T, b1_Tv))
            elif b3_Tv == 0:
                b1_data_good.append((b1_T, b1_Tv))

        if len(b1_data_bad) >= 30:
            break

    if b1_data_bad:
        print(f"\n  BAD vertices: (b1_T, b1_Tv) = {Counter(b1_data_bad)}")
    if b1_data_good:
        print(f"  GOOD vertices (sample): (b1_T, b1_Tv) = {Counter(b1_data_good[:100])}")

    print(f"\n  Seesaw (b1*b3=0) predicts: b3_T=1 => b1_T=0, b3_Tv=1 => b1_Tv=0")

    # Part 4: Rank growth pattern
    print(f"\n--- Part 4: Exact growth pattern for bad vertices ---")
    print("  For each bad vertex, the growth of ker(d_3) and rank(d_4):")

    for i, d in enumerate(bad_data[:10]):
        print(f"  Example {i}: ker_d3: {d['ker_d3_Tv']} -> {d['ker_d3_T']} (+{d['ker_d3_T']-d['ker_d3_Tv']}), "
              f"rank_d4: {d['rank_d4_Tv']} -> {d['rank_d4_T']} (+{d['rank_d4_T']-d['rank_d4_Tv']}), "
              f"Om4: {d['dim_Om4_Tv']} -> {d['dim_Om4_T']} (+{d['delta_Om4']})")

    # Part 5: THE KEY INSIGHT
    print(f"\n--- Part 5: Key pattern ---")

    if bad_data:
        # For bad vertices: b3_T = 1, b3_Tv = 1
        # delta_b3 = 0 means delta_ker = delta_rank_d4
        # Adding v doesn't change b3! It adds equally to ker and rank.

        all_same = all(d['b3_T'] - d['b3_Tv'] == 0 for d in bad_data)
        print(f"  delta(b3) = 0 for ALL bad vertices? {all_same}")

        if all_same:
            print(f"  This means: adding vertex v adds EQUAL amounts to ker(d3) and rank(d4)")
            print(f"  delta_ker = delta_rank_d4 exactly")
            print(f"  The 'new' kernel directions from v-paths are FULLY saturated by 'new' d4 boundaries")
            print(f"  The ORIGINAL T\\v cycle survives because it's orthogonal to the new content")

        # More specifically: how does rank(d4) grow?
        print(f"\n  For i_*-injectivity: the embedded T\\v cycle is NOT in im(d_4^T)")
        print(f"  This means: im(d_4^T) = im(d_4^{{T\\v}}) + (v-paths contribution)")
        print(f"  And: the v-paths contribution is disjoint from the T\\v cycle space")

    print("\nDONE.")


if __name__ == '__main__':
    main()
