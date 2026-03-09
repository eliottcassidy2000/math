#!/usr/bin/env python3
"""
istar_mechanism_deep.py — opus-2026-03-09-S54

WHY is i_*: H_3(T\v) → H_3(T) injective when b3(T\v)=1 and b3(T)=1?

The H_3 generator of T\v is a specific element z ∈ ker(d_3) / im(d_4)
in the (n-1)-vertex complex. When we include T\v ↪ T, z maps to
i_*(z) ∈ ker(d_3^T) / im(d_4^T).

i_* injective means: i_*(z) ∉ im(d_4^T), i.e., z doesn't become a boundary
when we add vertex v and its incident arcs.

QUESTION: What is the relationship between im(d_4^T) and im(d_4^{T\v})?
Adding v creates new Omega_4 elements (4-paths through v). Their boundaries
add to im(d_4^T). The question is whether these new boundaries "hit" z.

From the seesaw: when b3(T)=1, b1(T)=0. And b1(T\v)=0 (hereditary seesaw).
This means both T and T\v have im(d_2) at its maximum rank.

Let me compute: for beta_3=1 tournaments, what are the dimensions of
the key spaces in the chain complex?
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

def compute_chain_complex_data(T, n):
    """Compute full chain complex dimensions and ranks for a tournament."""
    tol = 1e-8
    all_paths = {}
    for p in range(0, min(n, 8)):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    # Build Omega (allowed paths modulo non-allowed faces)
    max_p = max(all_paths.keys())
    def build_omega(pd, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = pd.get(p, [])
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(pd.get(p-1, [])) if p > 0 else set()
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

    omega = build_omega(all_paths, max_p)
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
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary[p] = mat

    # Compute Omega-projected boundary
    ranks = {}
    dims = {}
    for p in range(max_p + 1):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
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

    # Compute Betti numbers
    betti = {}
    for p in range(max_p + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti[p] = ker - im_next

    return {
        'dims': dims,
        'ranks': ranks,
        'betti': betti,
        'paths_count': {p: len(all_paths.get(p, [])) for p in range(max_p + 1)},
    }

def main():
    print("=" * 70)
    print("i_*-INJECTIVITY MECHANISM: Deep analysis")
    print("=" * 70)

    # At n=6: exhaustive search for beta_3=1 tournaments
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    print(f"\n  Part 1: Chain complex data for beta_3=1 tournaments at n=6")

    type_a_data = []  # score (1,1,1,4,4,4)
    type_b_data = []  # score (2,2,2,3,3,3)

    for bits in range(n_total):
        if bits % 5000 == 0:
            print(f"    ... {bits}/{n_total}", flush=True)
        T = tournament_from_bits(n, bits)
        data = compute_chain_complex_data(T, n)
        if data['betti'].get(3, 0) == 1:
            scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
            entry = {
                'bits': bits,
                'scores': scores,
                **data,
            }
            if scores == (1,1,1,4,4,4):
                type_a_data.append(entry)
            else:
                type_b_data.append(entry)

    print(f"\n  Type A (1,1,1,4,4,4): {len(type_a_data)} tournaments")
    if type_a_data:
        d = type_a_data[0]
        print(f"    Omega dims: {[d['dims'].get(p,0) for p in range(6)]}")
        print(f"    d_p ranks:  {[d['ranks'].get(p,0) for p in range(1,6)]}")
        print(f"    Betti:      {[d['betti'].get(p,0) for p in range(6)]}")
        print(f"    Path counts: {[d['paths_count'].get(p,0) for p in range(6)]}")

    print(f"\n  Type B (2,2,2,3,3,3): {len(type_b_data)} tournaments")
    if type_b_data:
        d = type_b_data[0]
        print(f"    Omega dims: {[d['dims'].get(p,0) for p in range(6)]}")
        print(f"    d_p ranks:  {[d['ranks'].get(p,0) for p in range(1,6)]}")
        print(f"    Betti:      {[d['betti'].get(p,0) for p in range(6)]}")
        print(f"    Path counts: {[d['paths_count'].get(p,0) for p in range(6)]}")

    # Check: are all dims the same within each type?
    for name, data_list in [("Type A", type_a_data), ("Type B", type_b_data)]:
        dim_set = set()
        rank_set = set()
        for d in data_list:
            dim_set.add(tuple(d['dims'].get(p,0) for p in range(6)))
            rank_set.add(tuple(d['ranks'].get(p,0) for p in range(1,6)))
        print(f"\n  {name}: {len(dim_set)} distinct dim profiles, {len(rank_set)} distinct rank profiles")

    # Part 2: For each beta_3=1 tournament at n=6, check the DELETION chain complex
    # and the INCLUSION map i_*
    print(f"\n{'='*70}")
    print("Part 2: Deletion complex comparison")
    print("=" * 70)

    # For a Type B tournament, compare T and T\v chain complex data
    if type_b_data:
        d = type_b_data[0]
        T = tournament_from_bits(n, d['bits'])
        print(f"\n  Type B tournament (bits={d['bits']}):")
        print(f"    T dims:  {[d['dims'].get(p,0) for p in range(6)]}")
        print(f"    T ranks: {[d['ranks'].get(p,0) for p in range(1,6)]}")
        print(f"    T betti: {[d['betti'].get(p,0) for p in range(6)]}")

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            sub_data = compute_chain_complex_data(T_sub, n-1)
            score_v = sum(T[v][j] for j in range(n) if j != v)
            print(f"    T\\{v} (score_v={score_v}): dims={[sub_data['dims'].get(p,0) for p in range(5)]}, "
                  f"ranks={[sub_data['ranks'].get(p,0) for p in range(1,5)]}, "
                  f"betti={[sub_data['betti'].get(p,0) for p in range(5)]}")

    # Part 3: The key question — what changes when we add v?
    # For Type B: T has Omega_4 content (dim > 0), T\v has Omega_3 = dim(ker d_3 T\v)
    # When we add v, new Omega_4 paths appear (paths through v).
    # These create new im(d_4) content.
    # The question: does the NEW im(d_4) hit the H_3 generator of T\v?
    #
    # For i_* to be injective: the H_3 generator of T\v, when embedded in T's Omega_3,
    # must NOT be in im(d_4^T).
    #
    # The delta (change) in Omega_p dimensions when adding v:
    print(f"\n{'='*70}")
    print("Part 3: Delta dimensions (T vs T\\v)")
    print("=" * 70)

    for name, data_list in [("Type A", type_a_data[:1]), ("Type B", type_b_data[:1])]:
        if not data_list:
            continue
        d = data_list[0]
        T = tournament_from_bits(n, d['bits'])
        print(f"\n  {name} (bits={d['bits']}):")

        for v in range(min(3, n)):  # check first 3 vertices
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            sub_data = compute_chain_complex_data(T_sub, n-1)

            score_v = sum(T[v][j] for j in range(n) if j != v)
            print(f"\n    Vertex v={v} (score={score_v}):")
            for p in range(6):
                d_T = d['dims'].get(p, 0)
                d_Tv = sub_data['dims'].get(p, 0)
                r_T = d['ranks'].get(p, 0) if p > 0 else 0
                r_Tv = sub_data['ranks'].get(p, 0) if p > 0 else 0
                tv_label = "T\\v"
                print(f"      p={p}: dim(Omega_p^T)={d_T}, dim(Omega_p^{tv_label})={d_Tv}, "
                      f"delta_dim={d_T - d_Tv}, "
                      f"rank(d_p^T)={r_T}, rank(d_p^{tv_label})={r_Tv}, "
                      f"delta_rank={r_T - r_Tv}")

    # Part 4: Check the SATURATION property
    # HYP-383: bad vertices have δ(β₃)=0: adding v adds EQUAL kernel and im(d_4)
    print(f"\n{'='*70}")
    print("Part 4: Saturation check (delta ker(d_3) vs delta im(d_4))")
    print("=" * 70)

    n7 = 7
    n7_arcs = n7*(n7-1)//2
    n7_total = 1 << n7_arcs
    rng = np.random.RandomState(42)

    sat_data = defaultdict(Counter)
    count_b3_1 = 0

    for trial in range(1500):
        if trial % 300 == 0:
            print(f"    ... trial {trial}", flush=True)
        bits = rng.randint(0, n7_total)
        T = tournament_from_bits(n7, bits)
        try:
            T_data = compute_chain_complex_data(T, n7)
        except:
            continue
        if T_data['betti'].get(3, 0) != 1:
            continue
        count_b3_1 += 1

        for v in range(n7):
            remaining = [i for i in range(n7) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
            try:
                sub_data = compute_chain_complex_data(T_sub, n7-1)
            except:
                continue

            # Compute deltas
            delta_dim3 = T_data['dims'].get(3, 0) - sub_data['dims'].get(3, 0)
            delta_dim4 = T_data['dims'].get(4, 0) - sub_data['dims'].get(4, 0)
            delta_rank3 = T_data['ranks'].get(3, 0) - sub_data['ranks'].get(3, 0)
            delta_rank4 = T_data['ranks'].get(4, 0) - sub_data['ranks'].get(4, 0)
            delta_ker3 = (T_data['dims'].get(3,0) - T_data['ranks'].get(3,0)) - \
                         (sub_data['dims'].get(3,0) - sub_data['ranks'].get(3,0))
            delta_im4 = T_data['ranks'].get(4, 0) - sub_data['ranks'].get(4, 0)

            b3_sub = sub_data['betti'].get(3, 0)
            is_bad = b3_sub == 1

            key = "BAD" if is_bad else "GOOD"
            sat_data[key][(delta_ker3, delta_im4)] += 1

        if count_b3_1 >= 50:
            break

    print(f"\n  Found {count_b3_1} beta_3=1 tournaments at n=7")
    for key in ["GOOD", "BAD"]:
        print(f"\n  {key} vertices:")
        print(f"    (delta_ker3, delta_im4): count")
        for vals, cnt in sorted(sat_data[key].items()):
            delta_b3 = vals[0] - vals[1]
            print(f"    {vals}: {cnt} (delta_beta3 = {delta_b3})")

if __name__ == '__main__':
    main()
