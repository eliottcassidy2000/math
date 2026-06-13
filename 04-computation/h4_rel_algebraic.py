"""
h4_rel_algebraic.py — Investigate WHY H_4(T,T\\v) = 0 for beta_3=1 tournaments

Strategy: The relative chain complex R_* = Omega_*(T)/Omega_*(T\\v) consists of
paths that USE vertex v. We want to understand the algebraic structure that
forces H_4(R) = 0.

Key observations from computation:
- R_4 = "new" 5-paths through v
- d_4^rel: R_4 -> R_3 (faces of new 5-paths that are new 4-paths)
- d_5^rel: R_5 -> R_4 (new 6-paths through v)
- H_4(R) = ker(d_4^rel) / im(d_5^rel)

We want: rank(d_4^rel|_{R_4}) + rank(d_5^rel|_{R_5}) = dim(R_4)
i.e., d_4^rel is as surjective as possible and d_5^rel fills the rest.

Actually, H_4(R) = 0 means ker(d_4) = im(d_5) on R_4.

This script investigates:
1. Dimensions of R_p for p=3,4,5
2. Ranks of d_4^rel and d_5^rel
3. Whether there's a pattern linking R_4 dimension to beta_3

Author: opus-2026-03-09-S55
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed, boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp,
    RANK_PRIME
)

PRIME = RANK_PRIME


def compute_relative_complex_detailed(A, n, v, max_p=6):
    """Compute relative chain complex R_* with detailed boundary rank info.

    Returns: dict with R_dims, R_ranks, R_kers, R_bettis for each degree.
    """
    remaining = [i for i in range(n) if i != v]
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n-1, max_p)

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    # For each degree, find "new" paths (in T but not in T\v when embedded)
    R_bases = {}
    for p in range(max_p + 1):
        t_paths = ap_T.get(p, [])
        tv_paths_embedded = set(embed_path(q) for q in ap_Tv.get(p, []))
        new_paths = [q for q in t_paths if tuple(q) not in tv_paths_embedded]
        R_bases[p] = new_paths

    R_dims = {p: len(R_bases[p]) for p in range(max_p + 1)}

    # Build boundary matrices for the relative complex
    # d_p^rel: R_p -> R_{p-1} ∪ R_{p-1}
    # Actually d_p maps new p-paths to a mix of new and old (p-1)-paths.
    # In the quotient, we only keep the R_{p-1} components.

    # But wait — the relative complex is C_*(T)/C_*(T\v).
    # The boundary d_p in the quotient maps [sigma] in R_p to [d_p(sigma)] in R_{p-1}.
    # Since d_p(sigma) might have faces in C_{p-1}(T\v), those project to 0.
    # So d_p^rel[sigma] = sum of NEW faces of sigma.

    R_ranks = {}
    R_kers = {}

    for p in range(1, max_p + 1):
        if R_dims.get(p, 0) == 0 or R_dims.get(p-1, 0) == 0:
            R_ranks[p] = 0
            R_kers[p] = R_dims.get(p, 0)
            continue

        # Build boundary matrix: R_p -> R_{p-1}
        new_paths_p = R_bases[p]
        new_paths_pm1 = R_bases[p-1]
        pm1_index = {tuple(q): i for i, q in enumerate(new_paths_pm1)}

        nrows = len(new_paths_pm1)
        ncols = len(new_paths_p)
        bd = np.zeros((nrows, ncols), dtype=np.int64)

        for j, sigma in enumerate(new_paths_p):
            faces = boundary_faces(sigma)
            for sign, face in faces:
                face_t = tuple(face)
                if face_t in pm1_index:
                    bd[pm1_index[face_t], j] = (bd[pm1_index[face_t], j] + sign) % PRIME

        R_ranks[p] = _gauss_rank_np(bd % PRIME, PRIME)
        R_kers[p] = ncols - R_ranks[p]

    R_kers[0] = R_dims.get(0, 0)
    R_ranks[0] = 0

    # Betti numbers
    R_bettis = {}
    for p in range(max_p + 1):
        ker_p = R_kers.get(p, 0)
        im_from_above = R_ranks.get(p + 1, 0)
        R_bettis[p] = ker_p - im_from_above

    return {
        'R_dims': R_dims,
        'R_ranks': R_ranks,
        'R_kers': R_kers,
        'R_bettis': R_bettis,
    }


def analyze_beta3_tournaments(n=7, num_samples=200):
    """Analyze relative complex structure for beta_3=1 tournaments."""
    print(f"\n--- Part 1: Relative complex analysis, n={n} ---")

    found = 0
    t0 = time.time()

    # Collect data
    good_data = []  # (v, R_dims, R_ranks, R_bettis)
    bad_data = []

    dim_profiles_good = Counter()
    dim_profiles_bad = Counter()
    rank_profiles_good = Counter()
    rank_profiles_bad = Counter()

    while found < num_samples:
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p=n-1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            cc_sub = full_chain_complex_modp(A_sub, n-1, max_p=n-2)
            b3_sub = cc_sub['bettis'].get(3, 0)

            rel = compute_relative_complex_detailed(A, n, v, max_p=n-1)

            dims_key = tuple(rel['R_dims'].get(p, 0) for p in range(n))
            ranks_key = tuple(rel['R_ranks'].get(p, 0) for p in range(1, n))

            if b3_sub == 1:
                bad_data.append((v, rel))
                dim_profiles_bad[dims_key] += 1
                rank_profiles_bad[ranks_key] += 1
            else:
                good_data.append((v, rel))
                dim_profiles_good[dims_key] += 1
                rank_profiles_good[ranks_key] += 1

        if found % 20 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{num_samples} beta_3=1 found, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tournaments, {len(good_data)} good verts, {len(bad_data)} bad verts, {elapsed:.1f}s")

    # Check H_4^rel = 0
    h4_violations = 0
    for v, rel in good_data + bad_data:
        if rel['R_bettis'].get(4, 0) != 0:
            h4_violations += 1
    print(f"\n  H_4^rel = 0 violations: {h4_violations}")

    # Print dimension profiles
    print(f"\n  GOOD vertex R_dim profiles (top 10):")
    for dims, count in dim_profiles_good.most_common(10):
        print(f"    {dims}: {count}")

    print(f"\n  BAD vertex R_dim profiles (top 10):")
    for dims, count in dim_profiles_bad.most_common(10):
        print(f"    {dims}: {count}")

    # Key question: what makes d_4^rel + d_5^rel cover all of R_4?
    print(f"\n  GOOD vertex rank profiles d_1..d_{n-1} (top 10):")
    for ranks, count in rank_profiles_good.most_common(10):
        print(f"    {ranks}: {count}")

    print(f"\n  BAD vertex rank profiles d_1..d_{n-1} (top 10):")
    for ranks, count in rank_profiles_bad.most_common(10):
        print(f"    {ranks}: {count}")

    return good_data, bad_data


def analyze_h4_mechanism(n=7, num_samples=50):
    """Deep analysis: WHY is H_4(R) = 0?

    For H_4(R) = 0, need: dim(ker d_4^rel) = dim(im d_5^rel)
    i.e., rank(d_5^rel) = dim(R_4) - rank(d_4^rel)

    Let's check if there's a universal relationship.
    """
    print(f"\n--- Part 2: H_4=0 mechanism analysis, n={n} ---")

    found = 0
    t0 = time.time()

    # For each beta_3=1 tournament, for each vertex v, record:
    # (R_4 dim, rank d_4, rank d_5, ker d_4, im d_5)
    records = []

    while found < num_samples:
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p=n-1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            rel = compute_relative_complex_detailed(A, n, v, max_p=n-1)

            dim_R4 = rel['R_dims'].get(4, 0)
            rank_d4 = rel['R_ranks'].get(4, 0)
            rank_d5 = rel['R_ranks'].get(5, 0)
            ker_d4 = rel['R_kers'].get(4, 0)
            h4 = rel['R_bettis'].get(4, 0)

            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            cc_sub = full_chain_complex_modp(A_sub, n-1, max_p=n-2)
            b3_sub = cc_sub['bettis'].get(3, 0)

            vtype = "BAD" if b3_sub == 1 else "GOOD"
            records.append({
                'type': vtype,
                'dim_R4': dim_R4,
                'rank_d4': rank_d4,
                'rank_d5': rank_d5,
                'ker_d4': ker_d4,
                'h4': h4,
                'dim_R3': rel['R_dims'].get(3, 0),
                'dim_R5': rel['R_dims'].get(5, 0),
            })

    elapsed = time.time() - t0
    print(f"  Done: {found} tournaments, {len(records)} vertex records, {elapsed:.1f}s")

    # Analyze the mechanism
    # Group by type
    for vtype in ["GOOD", "BAD"]:
        recs = [r for r in records if r['type'] == vtype]
        print(f"\n  {vtype} vertices ({len(recs)}):")

        # Key relationship: rank_d5 = ker_d4 when H_4=0
        rank_d5_eq_ker_d4 = sum(1 for r in recs if r['rank_d5'] == r['ker_d4'])
        print(f"    rank(d5) = ker(d4): {rank_d5_eq_ker_d4}/{len(recs)}")

        # What fraction of R_4 is killed by d_4?
        if recs:
            avg_d4_ratio = np.mean([r['rank_d4'] / max(r['dim_R4'], 1) for r in recs])
            avg_d5_ratio = np.mean([r['rank_d5'] / max(r['dim_R4'], 1) for r in recs])
            print(f"    avg rank(d4)/dim(R4) = {avg_d4_ratio:.3f}")
            print(f"    avg rank(d5)/dim(R4) = {avg_d5_ratio:.3f}")
            print(f"    avg (rank_d4 + rank_d5)/dim(R4) = {avg_d4_ratio + avg_d5_ratio:.3f}")

        # Distribution of (dim_R4, rank_d4, rank_d5)
        triple_dist = Counter()
        for r in recs:
            triple_dist[(r['dim_R4'], r['rank_d4'], r['rank_d5'])] += 1
        print(f"    (dim_R4, rank_d4, rank_d5) distribution (top 10):")
        for triple, count in triple_dist.most_common(10):
            d, r4, r5 = triple
            h4 = d - r4 - r5
            print(f"      {triple}: {count}  (h4={h4})")


def analyze_vertex_role(n=7, num_samples=30):
    """Analyze how vertex v sits in the new paths.

    Key question: Do new paths always go THROUGH v (not just start/end)?
    What is the position distribution of v in new paths?
    """
    print(f"\n--- Part 3: Vertex position analysis in new paths ---")

    found = 0
    t0 = time.time()

    pos_counts = defaultdict(Counter)  # degree -> position -> count

    while found < num_samples:
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p=n-1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        v = 0  # Just analyze vertex 0 for speed
        remaining = [i for i in range(n) if i != v]
        A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]

        ap_T = enumerate_all_allowed(A, n, n-1)
        ap_Tv = enumerate_all_allowed(A_sub, n-1, n-2)

        def embed_path(path_tv):
            return tuple(remaining[i] for i in path_tv)

        for p in range(2, n):
            t_paths = ap_T.get(p, [])
            tv_paths_embedded = set(embed_path(q) for q in ap_Tv.get(p, []))
            new_paths = [q for q in t_paths if tuple(q) not in tv_paths_embedded]

            for path in new_paths:
                path_t = tuple(path)
                if v in path_t:
                    pos = path_t.index(v)
                    pos_counts[p][pos] += 1
                else:
                    pos_counts[p]['absent'] += 1

    elapsed = time.time() - t0
    print(f"  Done: {found} tournaments, {elapsed:.1f}s")

    for p in sorted(pos_counts.keys()):
        total = sum(pos_counts[p].values())
        absent = pos_counts[p].get('absent', 0)
        print(f"\n  Degree {p}: {total} new paths, {absent} DON'T contain v ({100*absent/max(total,1):.1f}%)")
        for pos in range(p + 1):
            c = pos_counts[p].get(pos, 0)
            if c > 0:
                print(f"    position {pos}: {c} ({100*c/max(total,1):.1f}%)")


if __name__ == '__main__':
    print("=" * 70)
    print("H_4(T,T\\v) = 0 MECHANISM ANALYSIS")
    print("Why is relative H_4 always zero for beta_3=1 tournaments?")
    print("=" * 70)

    good_data, bad_data = analyze_beta3_tournaments(n=7, num_samples=100)
    analyze_h4_mechanism(n=7, num_samples=50)
    analyze_vertex_role(n=7, num_samples=30)

    print("\nDONE.")
