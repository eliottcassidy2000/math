"""
quasi_iso_algebraic.py — Investigate the quasi-isomorphism T\\v -> T for BAD vertices

When b3(T)=1 and b3(T\\v)=1, ALL relative homology H_p(T,T\\v) = 0.
This means the inclusion is a quasi-isomorphism.

WHY? Key observations:
1. At n=7, T\\v has Betti (1,0,0,1,0,0) — exactly like T=(1,0,0,1,0,0,0)
2. Both have chi = 0
3. The "new" paths through v are exact at every degree

This script investigates:
- The chain complex dimension differences (Omega_p(T) vs Omega_p(T\\v))
- Whether the dimension deficit follows a pattern
- Whether the Euler characteristic of the relative complex is always 0

Author: opus-2026-03-09-S55
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    full_chain_complex_modp,
)


def analyze_chi_universal(n=7, num_samples=200):
    """Check chi(T) for ALL tournaments, not just beta_3=1."""
    print(f"\n--- Part 1: Euler characteristic analysis, n={n} ---")
    max_p = n - 1
    rng = np.random.RandomState(42)

    chi_by_betti = Counter()  # (betti_profile) -> set of chi values
    chi_dist = Counter()

    t0 = time.time()
    for i in range(num_samples):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p)
        bettis = cc['bettis']
        profile = tuple(bettis.get(p, 0) for p in range(max_p + 1))
        chi = sum((-1)**p * bettis.get(p, 0) for p in range(max_p + 1))

        chi_dist[chi] += 1
        if profile not in chi_by_betti:
            chi_by_betti[profile] = set()
        chi_by_betti[profile].add(chi)

    elapsed = time.time() - t0
    print(f"  {num_samples} tournaments in {elapsed:.1f}s")

    print(f"\n  chi distribution:")
    for chi, count in sorted(chi_dist.items()):
        print(f"    chi={chi}: {count}")

    print(f"\n  chi by Betti profile:")
    for profile, chis in sorted(chi_by_betti.items()):
        print(f"    {profile}: chi={chis}")


def analyze_omega_deficits(n=7, num_samples=100):
    """For beta_3=1 tournaments, compare Omega dims of T and T\\v."""
    print(f"\n--- Part 2: Omega dimension deficits for beta_3=1, n={n} ---")
    max_p = n - 1
    rng = np.random.RandomState(42)

    deficit_profiles = Counter()  # BAD vertex deficit profiles
    deficit_profiles_good = Counter()  # GOOD vertex deficit profiles

    found = 0
    t0 = time.time()

    for trial in range(5000):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1
        if found > num_samples:
            break

        omega_T = cc['omega_dims']
        ranks_T = cc['ranks']

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            cc_sub = full_chain_complex_modp(A_sub, n-1, min(max_p, n-2))
            omega_Tv = cc_sub['omega_dims']
            b3_sub = cc_sub['bettis'].get(3, 0)

            deficit = tuple(omega_T.get(p, 0) - omega_Tv.get(p, 0) for p in range(max_p + 1))

            # Euler char of deficit
            chi_deficit = sum((-1)**p * d for p, d in enumerate(deficit))

            if b3_sub == 1:
                deficit_profiles[deficit] += 1
            else:
                deficit_profiles_good[deficit] += 1

    elapsed = time.time() - t0
    print(f"  {found} beta_3=1 tournaments, {elapsed:.1f}s")

    print(f"\n  BAD vertex Omega deficits (top 10):")
    for deficit, count in deficit_profiles.most_common(10):
        chi = sum((-1)**p * d for p, d in enumerate(deficit))
        print(f"    {deficit}: {count}, chi_deficit={chi}")

    print(f"\n  GOOD vertex Omega deficits (top 10):")
    for deficit, count in deficit_profiles_good.most_common(10):
        chi = sum((-1)**p * d for p, d in enumerate(deficit))
        print(f"    {deficit}: {count}, chi_deficit={chi}")


def analyze_rank_deficits(n=7, num_samples=50):
    """Compare boundary ranks of T and T\\v for bad vertices.

    If the inclusion is a quasi-iso, the rank deficits must satisfy
    a specific relationship with the Omega deficits.

    For BAD vertices: bettis match at all p, so:
    omega_deficit[p] = rank_deficit[p] + rank_deficit[p+1] (for each p>=1)
    where rank_deficit[p] = rank(d_p^T) - rank(d_p^{T\\v})
    """
    print(f"\n--- Part 3: Boundary rank deficits for beta_3=1, n={n} ---")
    max_p = n - 1
    rng = np.random.RandomState(42)

    found = 0
    t0 = time.time()

    verified = 0
    total_bad = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1
        if found > num_samples:
            break

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            cc_sub = full_chain_complex_modp(A_sub, n-1, min(max_p, n-2))
            b3_sub = cc_sub['bettis'].get(3, 0)

            if b3_sub != 1:
                continue
            total_bad += 1

            # Check: omega_deficit = rank_deficit[p] + rank_deficit[p+1]
            ok = True
            for p in range(max_p + 1):
                od = cc['omega_dims'].get(p, 0) - cc_sub['omega_dims'].get(p, 0)
                rd_p = cc['ranks'].get(p, 0) - cc_sub['ranks'].get(p, 0)
                rd_pp1 = cc['ranks'].get(p + 1, 0) - cc_sub['ranks'].get(p + 1, 0)
                betti_T_p = cc['bettis'].get(p, 0)
                betti_Tv_p = cc_sub['bettis'].get(p, 0)

                # Betti deficit should be 0 for quasi-iso
                betti_deficit = betti_T_p - betti_Tv_p
                # omega_deficit = rank_deficit[p] + rank_deficit[p+1] + betti_deficit
                expected = rd_p + rd_pp1 + betti_deficit
                if od != expected:
                    ok = False
                    if total_bad <= 3:
                        print(f"  MISMATCH at p={p}: omega_deficit={od}, "
                              f"rank_deficit[p]+[p+1]+betti_deficit={expected}")

            if ok:
                verified += 1

    elapsed = time.time() - t0
    print(f"  {found} tours, {total_bad} bad verts, {verified}/{total_bad} verified, {elapsed:.1f}s")


def investigate_degree_1_rank(n=7, num_samples=50):
    """Check how in-/out-degree of v relates to rank changes.

    For vertex v in tournament T, let d_out(v) = number of edges v->w.
    Does d_out(v) determine the Omega/rank deficits?
    """
    print(f"\n--- Part 4: Vertex degree vs rank deficit, n={n} ---")
    max_p = n - 1
    rng = np.random.RandomState(42)

    # For each (d_out, is_bad), collect omega deficits
    deficit_by_degree = {}

    found = 0
    t0 = time.time()

    for trial in range(5000):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1
        if found > num_samples:
            break

        for v in range(n):
            d_out = int(sum(A[v]))
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            cc_sub = full_chain_complex_modp(A_sub, n-1, min(max_p, n-2))
            b3_sub = cc_sub['bettis'].get(3, 0)
            is_bad = b3_sub == 1

            key = (d_out, is_bad)
            deficit = tuple(cc['omega_dims'].get(p, 0) - cc_sub['omega_dims'].get(p, 0) for p in range(max_p + 1))

            if key not in deficit_by_degree:
                deficit_by_degree[key] = Counter()
            deficit_by_degree[key][deficit] += 1

    elapsed = time.time() - t0
    print(f"  {found} tours, {elapsed:.1f}s")

    for (d_out, is_bad), deficits in sorted(deficit_by_degree.items()):
        vtype = "BAD" if is_bad else "GOOD"
        total = sum(deficits.values())
        unique = len(deficits)
        print(f"\n  d_out={d_out}, {vtype} ({total} verts, {unique} unique deficits):")
        for deficit, count in deficits.most_common(3):
            chi = sum((-1)**p * d for p, d in enumerate(deficit))
            print(f"    {deficit}: {count} (chi={chi})")


if __name__ == '__main__':
    print("=" * 70)
    print("QUASI-ISOMORPHISM ANALYSIS FOR BAD VERTICES")
    print("=" * 70)

    analyze_chi_universal(n=7, num_samples=500)
    analyze_omega_deficits(n=7, num_samples=80)
    analyze_rank_deficits(n=7, num_samples=50)
    investigate_degree_1_rank(n=7, num_samples=50)

    print("\nDONE.")
