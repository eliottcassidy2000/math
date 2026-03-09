"""
relative_h3_structure.py — Understanding WHY dim H_3(T, T\\v) <= 1

The relative chain complex C_*^rel consists of Omega chains through vertex v.
We analyze the structure of:
  - C_3^rel (through-v 3-cycles in Omega_3)
  - C_2^rel (through-v 2-chains in Omega_2)
  - C_4^rel (through-v 4-chains in Omega_4)
  - d_3^rel: C_3^rel -> C_2^rel
  - d_4^rel: C_4^rel -> C_3^rel

Key question: why is ker(d_3^rel) - rank(d_4^rel) <= 1?

Approach: investigate dimensional structure of relative complex.

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


def compute_relative_complex_dims(A, n, v):
    """
    Compute dimensions of relative chain complex C_*(T, T\\v).

    Returns dict with:
      rel_omega_dim[p] = dim(Omega_p^rel)
      rel_ker_d[p] = dim(ker(d_p^rel))
      rel_rank_d[p] = rank(d_p^rel)
      rel_h[p] = beta_p^rel = dim(H_p^rel)
    """
    # Full complex
    ap_full = {}
    for p in range(min(6, n)):
        ap_full[p] = enumerate_allowed_paths(A, n, p)

    # Compute through-v paths
    through_v = {}
    for p in range(min(6, n)):
        through_v[p] = [path for path in ap_full[p] if v in path]

    # Compute Omega bases for full complex
    ob_full = {}
    od_full = {}
    for p in range(min(6, n)):
        if not ap_full.get(p, []):
            ob_full[p] = np.zeros((0, 0))
            od_full[p] = 0
        else:
            ob_full[p], od_full[p] = compute_omega_basis(ap_full, p)

    # Compute relative Omega dimensions
    tv_idx = {}
    for p in range(min(6, n)):
        tv_idx[p] = {path: i for i, path in enumerate(through_v[p])}

    rel_bases = {}
    rel_dims = {}
    for p in range(min(6, n)):
        if od_full.get(p, 0) == 0 or not through_v[p]:
            rel_bases[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue

        num_full = len(ap_full[p])
        num_tv = len(through_v[p])
        proj = np.zeros((num_tv, num_full))
        full_idx = {path: i for i, path in enumerate(ap_full[p])}
        for path, ti in tv_idx[p].items():
            proj[ti, full_idx[path]] = 1.0
        proj_omega = proj @ ob_full[p]
        if proj_omega.size == 0:
            rel_bases[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue
        sv = np.linalg.svd(proj_omega, compute_uv=False)
        rank = int(sum(s > 1e-10 for s in sv))
        rel_dims[p] = rank
        if rank == 0:
            rel_bases[p] = np.zeros((num_tv, 0))
            continue
        U, S, Vt = np.linalg.svd(proj_omega, full_matrices=False)
        rel_bases[p] = U[:, :rank]

    # Compute boundary maps in relative complex
    rel_ranks = {}
    rel_kers = {}
    rel_bettis = {}
    for p in range(1, min(6, n)):
        if rel_dims.get(p, 0) == 0:
            rel_ranks[p] = 0
            rel_kers[p] = 0
            continue
        # d_p^rel: C_p^rel -> C_{p-1}^rel
        bd_tv = np.zeros((len(through_v.get(p-1, [])), len(through_v[p])))
        for j, path in enumerate(through_v[p]):
            for sign, face in boundary_coeffs(path):
                if face in tv_idx.get(p-1, {}):
                    bd_tv[tv_idx[p-1][face], j] += sign

        if rel_dims.get(p-1, 0) > 0 and rel_dims[p] > 0:
            d_rel = rel_bases[p-1].T @ bd_tv @ rel_bases[p]
        elif rel_dims[p] > 0:
            d_rel = np.zeros((0, rel_dims[p]))
        else:
            d_rel = np.zeros((rel_dims.get(p-1, 0), 0))

        sv = np.linalg.svd(d_rel, compute_uv=False) if d_rel.size > 0 else np.array([])
        rel_ranks[p] = int(sum(s > 1e-8 for s in sv))
        rel_kers[p] = rel_dims[p] - rel_ranks[p]

    # Compute Betti numbers
    for p in range(min(6, n)):
        ker = rel_kers.get(p, rel_dims.get(p, 0))
        im_next = rel_ranks.get(p+1, 0)
        rel_bettis[p] = max(0, ker - im_next)

    return {
        'tv_counts': {p: len(through_v[p]) for p in range(min(6, n))},
        'rel_dims': rel_dims,
        'rel_ranks': rel_ranks,
        'rel_kers': rel_kers,
        'rel_bettis': rel_bettis,
    }


def main():
    print("=" * 70)
    print("RELATIVE H_3 STRUCTURE — Why dim H_3(T, T\\v) <= 1")
    print("=" * 70)

    # Part 1: Detailed relative complex dimensions at n=6 (beta_3>0)
    print("\n--- Part 1: Relative complex at n=6, beta_3>0 tournaments ---")
    n = 6
    dim_profiles = Counter()
    rank_profiles = Counter()

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        # Quick beta_3 check via scores
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores != (1, 1, 1, 4, 4, 4) and scores != (2, 2, 2, 3, 3, 3):
            continue

        # Only process beta_3>0 tournaments (rough filter by score, then verify)
        ap = {}
        for p in range(min(6, n)):
            ap[p] = enumerate_allowed_paths(A, n, p)
        ob3, od3 = compute_omega_basis(ap, 3)
        if od3 == 0:
            continue

        bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign
        d3 = bd3 @ ob3
        sv = np.linalg.svd(d3, compute_uv=False)
        rk3 = int(sum(s > 1e-8 for s in sv))
        if od3 - rk3 == 0:
            continue

        # This tournament has beta_3 > 0. Compute relative complex for each v.
        for v in range(n):
            result = compute_relative_complex_dims(A, n, v)
            dims = tuple(result['rel_dims'].get(p, 0) for p in range(6))
            ranks = tuple(result['rel_ranks'].get(p, 0) for p in range(1, 6))
            dim_profiles[dims] += 1
            rank_profiles[(dims, ranks)] += 1

    print(f"  Relative Omega dimension profiles (Omega_0^rel,...,Omega_5^rel):")
    for dims, cnt in dim_profiles.most_common():
        print(f"    {dims}: {cnt}")
    print(f"\n  (dims, ranks) profiles:")
    for (dims, ranks), cnt in rank_profiles.most_common():
        ker3 = dims[3] - ranks[2]  # rank of d_3^rel is ranks[2] (d_3 goes C_3->C_2)
        im4 = ranks[3]  # rank of d_4^rel
        h3 = ker3 - im4
        print(f"    dims={dims}, ranks={ranks}, ker_d3={ker3}, im_d4={im4}, H_3={h3}: {cnt}")

    # Part 2: Position of v in through-v paths
    print("\n--- Part 2: Position of v in through-v allowed 3-paths at n=6 ---")
    n = 6
    pos_dist = Counter()
    count = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores != (2, 2, 2, 3, 3, 3):
            continue

        ap = {}
        for p in range(min(6, n)):
            ap[p] = enumerate_allowed_paths(A, n, p)
        ob3, od3 = compute_omega_basis(ap, 3)
        if od3 == 0:
            continue

        bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign
        d3 = bd3 @ ob3
        sv = np.linalg.svd(d3, compute_uv=False)
        rk3 = int(sum(s > 1e-8 for s in sv))
        if od3 - rk3 == 0:
            continue
        count += 1
        if count > 10:
            break

        for v in range(n):
            tv_paths = [p for p in ap[3] if v in p]
            for path in tv_paths:
                pos = path.index(v)
                pos_dist[pos] += 1

    print(f"  Position of v in through-v 3-paths (first 10 Type B tournaments):")
    for pos in sorted(pos_dist.keys()):
        print(f"    position {pos}: {pos_dist[pos]}")

    # Part 3: n=7 relative complex dimensions
    print("\n--- Part 3: Relative complex at n=7, beta_3>0 tournaments ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 300
    dim_profiles7 = Counter()
    b3_count = 0

    for trial in range(N):
        A = random_tournament(n, rng)

        ap = {}
        for p in range(min(6, n)):
            ap[p] = enumerate_allowed_paths(A, n, p)
        ob3, od3 = compute_omega_basis(ap, 3)
        if od3 == 0:
            continue

        bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign
        d3 = bd3 @ ob3
        sv = np.linalg.svd(d3, compute_uv=False)
        rk3 = int(sum(s > 1e-8 for s in sv))
        if od3 - rk3 == 0:
            continue
        b3_count += 1

        # Pick one representative vertex
        v = 0
        result = compute_relative_complex_dims(A, n, v)
        dims = tuple(result['rel_dims'].get(p, 0) for p in range(6))
        dim_profiles7[dims] += 1

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N} done, {b3_count} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count}/{N}")
    print(f"  Relative Omega dims (v=0):")
    for dims, cnt in dim_profiles7.most_common():
        print(f"    {dims}: {cnt}")

    # Part 4: Key ratio: dim(C_3^rel) vs dim(C_2^rel) and dim(C_4^rel)
    print("\n--- Part 4: Key dimensional ratios at n=6 ---")
    n = 6
    ratio_data = []

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        if scores not in [(1, 1, 1, 4, 4, 4), (2, 2, 2, 3, 3, 3)]:
            continue

        ap = {}
        for p in range(min(6, n)):
            ap[p] = enumerate_allowed_paths(A, n, p)
        ob3, od3 = compute_omega_basis(ap, 3)
        if od3 == 0:
            continue
        bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
        idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
        for j, path in enumerate(ap[3]):
            for sign, face in boundary_coeffs(path):
                if face in idx2:
                    bd3[idx2[face], j] += sign
        d3 = bd3 @ ob3
        sv = np.linalg.svd(d3, compute_uv=False)
        rk3 = int(sum(s > 1e-8 for s in sv))
        if od3 - rk3 == 0:
            continue

        for v in range(n):
            result = compute_relative_complex_dims(A, n, v)
            d = result['rel_dims']
            ratio_data.append((scores, d.get(2,0), d.get(3,0), d.get(4,0),
                             result['rel_bettis'].get(3, 0)))

    # Summarize
    type_summary = defaultdict(list)
    for scores, d2, d3, d4, h3 in ratio_data:
        type_summary[(scores, d2, d3, d4)].append(h3)

    for (sc, d2, d3, d4), h3_vals in sorted(type_summary.items()):
        h3_dist = Counter(h3_vals)
        print(f"  scores={sc}, (d2,d3,d4)=({d2},{d3},{d4}): H_3 dist={dict(h3_dist)}")

    # Part 5: Relative complex for transitive tournaments (should be trivial)
    print("\n--- Part 5: Sanity check — relative complex for transitive tournament ---")
    n = 6
    # Transitive: A[i][j] = 1 iff i < j
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1

    for v in range(n):
        result = compute_relative_complex_dims(A, n, v)
        d = result['rel_dims']
        b = result['rel_bettis']
        print(f"  v={v}: dims={d}, bettis={b}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
