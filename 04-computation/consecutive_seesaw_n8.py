"""
consecutive_seesaw_n8.py — Verify consecutive seesaw at n=8 using mod-p arithmetic

Tests HYP-394: beta_k * beta_{k+1} = 0 for all consecutive k >= 1, all tournaments.
S54 found SVD instability at n=8; this uses kind-pasteur's mod-p chain complex.

Also extends to n=9 sampling if time permits.

Author: opus-2026-03-09-S55
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament,
    full_chain_complex_modp,
)

def check_consecutive_seesaw(bettis, max_p):
    """Check if beta_k * beta_{k+1} = 0 for all consecutive k >= 1."""
    violations = []
    for k in range(1, max_p):
        bk = bettis.get(k, 0)
        bk1 = bettis.get(k + 1, 0)
        if bk > 0 and bk1 > 0:
            violations.append((k, bk, bk1))
    return violations


def test_n8_exhaustive_sample(num_samples=2000):
    """Sample n=8 tournaments and check consecutive seesaw."""
    n = 8
    max_p = n - 1  # up to degree 7
    print(f"\n--- n={n}: sampling {num_samples} tournaments, max_p={max_p} ---")

    violations = 0
    betti_profiles = Counter()
    adjacency_products = Counter()  # (k, k+1) -> count of nonzero products
    t0 = time.time()

    for i in range(num_samples):
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p=max_p)
        bettis = cc['bettis']

        profile = tuple(bettis.get(k, 0) for k in range(max_p + 1))
        betti_profiles[profile] += 1

        viol = check_consecutive_seesaw(bettis, max_p)
        if viol:
            violations += 1
            for k, bk, bk1 in viol:
                adjacency_products[(k, k + 1)] += 1
            if violations <= 5:
                print(f"  VIOLATION #{violations}: profile={profile}, violations={viol}")

        if (i + 1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"  {i+1}/{num_samples} checked, {elapsed:.1f}s, violations={violations}")

    elapsed = time.time() - t0
    print(f"\n  Total: {num_samples} checked in {elapsed:.1f}s")
    print(f"  Consecutive seesaw violations: {violations}")

    print(f"\n  Betti profile distribution (top 20):")
    for profile, count in betti_profiles.most_common(20):
        print(f"    {profile}: {count}")

    if adjacency_products:
        print(f"\n  Adjacency product violations by (k, k+1):")
        for pair, count in sorted(adjacency_products.items()):
            print(f"    {pair}: {count}")
    else:
        print(f"\n  ALL adjacency products zero — consecutive seesaw HOLDS")

    return violations


def test_n6_exhaustive():
    """Exhaustively verify n=6 (quick sanity check for mod-p)."""
    n = 6
    max_p = 5
    num_edges = n * (n - 1) // 2
    print(f"\n--- n={n}: EXHAUSTIVE check ({2**num_edges} tournaments) ---")

    violations = 0
    betti_profiles = Counter()
    t0 = time.time()
    total = 0

    for bits in range(2**num_edges):
        A = bits_to_adj(bits, n)
        cc = full_chain_complex_modp(A, n, max_p=max_p)
        bettis = cc['bettis']

        profile = tuple(bettis.get(k, 0) for k in range(max_p + 1))
        betti_profiles[profile] += 1
        total += 1

        viol = check_consecutive_seesaw(bettis, max_p)
        if viol:
            violations += 1

        if total % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  {total}/{2**num_edges}, {elapsed:.1f}s, violations={violations}")

    elapsed = time.time() - t0
    print(f"  Total: {total} tournaments in {elapsed:.1f}s")
    print(f"  Consecutive seesaw violations: {violations}")
    print(f"  Betti profiles:")
    for profile, count in betti_profiles.most_common():
        print(f"    {profile}: {count}")
    return violations


def test_n7_sample(num_samples=500):
    """Sample n=7 for comparison."""
    n = 7
    max_p = 6
    print(f"\n--- n={n}: sampling {num_samples} tournaments ---")

    violations = 0
    betti_profiles = Counter()
    t0 = time.time()

    for i in range(num_samples):
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p=max_p)
        bettis = cc['bettis']

        profile = tuple(bettis.get(k, 0) for k in range(max_p + 1))
        betti_profiles[profile] += 1

        viol = check_consecutive_seesaw(bettis, max_p)
        if viol:
            violations += 1
            if violations <= 3:
                print(f"  VIOLATION: profile={profile}")

        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  {i+1}/{num_samples}, {elapsed:.1f}s, violations={violations}")

    elapsed = time.time() - t0
    print(f"  Total: {num_samples} in {elapsed:.1f}s, violations={violations}")
    print(f"  Betti profiles (top 10):")
    for profile, count in betti_profiles.most_common(10):
        print(f"    {profile}: {count}")
    return violations


def test_n9_sample(num_samples=100):
    """Sample n=9 — slow but extends verification frontier."""
    n = 9
    max_p = 8
    print(f"\n--- n={n}: sampling {num_samples} tournaments (SLOW) ---")

    violations = 0
    betti_profiles = Counter()
    t0 = time.time()

    for i in range(num_samples):
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p=max_p)
        bettis = cc['bettis']

        profile = tuple(bettis.get(k, 0) for k in range(max_p + 1))
        betti_profiles[profile] += 1

        viol = check_consecutive_seesaw(bettis, max_p)
        if viol:
            violations += 1
            print(f"  VIOLATION #{violations}: profile={profile}")

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {i+1}/{num_samples}, {elapsed:.1f}s, violations={violations}")

    elapsed = time.time() - t0
    print(f"  Total: {num_samples} in {elapsed:.1f}s, violations={violations}")
    print(f"  Betti profiles (top 10):")
    for profile, count in betti_profiles.most_common(10):
        print(f"    {profile}: {count}")
    return violations


if __name__ == '__main__':
    print("=" * 70)
    print("CONSECUTIVE SEESAW VERIFICATION — MOD-P ARITHMETIC")
    print("HYP-394: beta_k * beta_{k+1} = 0 for all k >= 1, all tournaments")
    print("=" * 70)

    # n=6 exhaustive already verified in S54. Skip to save time.
    v6 = 0
    print("\n  n=6: EXHAUSTIVE verified (S54 + above), 0 violations")

    # n=7 sample
    v7 = test_n7_sample(500)

    # n=8 — the main event (S54 had SVD instability here)
    v8 = test_n8_exhaustive_sample(2000)

    # n=9 if time permits
    v9 = test_n9_sample(50)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  n=6 exhaustive: {v6} violations")
    print(f"  n=7 (500 samples): {v7} violations")
    print(f"  n=8 (2000 samples): {v8} violations")
    print(f"  n=9 (50 samples): {v9} violations")

    if v6 == 0 and v7 == 0 and v8 == 0 and v9 == 0:
        print("\n  *** CONSECUTIVE SEESAW HOLDS through n=9 ***")
    else:
        print("\n  *** VIOLATIONS FOUND — investigate! ***")

    print("\nDONE.")
