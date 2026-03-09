"""
higher_betti_survey.py — Survey of higher Betti numbers for beta_3=1 tournaments

Key questions:
1. Does beta_3=1 imply beta_4=0? (needed for chi(T)=0)
2. Does chi(T)=0 hold for ALL beta_3=1 tournaments?
3. When does beta_5 first appear?
4. Does the LES dichotomy extend to beta_5?

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament, full_chain_complex_modp
)


def main():
    print("=" * 70)
    print("HIGHER BETTI NUMBER SURVEY")
    print("=" * 70)

    # Part 1: n=7 exhaustive — beta_3 vs beta_4 mutual exclusion
    print("\n--- Part 1: n=7 exhaustive beta_3 vs beta_4 ---")
    n = 7
    total = 2 ** (n*(n-1)//2)

    b3_b4_dist = Counter()
    chi_dist_b3_1 = Counter()
    all_betti_b3_1 = Counter()
    count = 0

    t0 = time.time()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        res = full_chain_complex_modp(A, n, max_p=6)
        b = res['bettis']
        b3 = b.get(3, 0)
        b4 = b.get(4, 0)
        b3_b4_dist[(b3, b4)] += 1

        if b3 == 1:
            betti_vec = tuple(b.get(p, 0) for p in range(7))
            all_betti_b3_1[betti_vec] += 1
            chi = sum((-1)**p * b.get(p, 0) for p in range(7))
            chi_dist_b3_1[chi] += 1

        count += 1
        if count % 100000 == 0:
            elapsed = time.time() - t0
            rate = count / elapsed
            remaining = (total - count) / rate
            print(f"  {count}/{total} ({count*100/total:.1f}%), {rate:.0f}/s, ETA {remaining:.0f}s",
                  flush=True)

    t1 = time.time()
    print(f"\n  n=7: {total} tournaments in {t1-t0:.1f}s")
    print(f"\n  (beta_3, beta_4) distribution:")
    for key, cnt in sorted(b3_b4_dist.items()):
        print(f"    {key}: {cnt}")

    print(f"\n  chi(T) distribution for beta_3=1:")
    for key, cnt in sorted(chi_dist_b3_1.items()):
        print(f"    chi={key}: {cnt}")

    print(f"\n  Full Betti vectors for beta_3=1:")
    for key, cnt in sorted(all_betti_b3_1.items()):
        print(f"    {key}: {cnt}")

    # Part 2: n=8 sampled — beta_3 vs beta_4
    print("\n--- Part 2: n=8 sampled beta_3 vs beta_4 ---")
    n = 8
    rng = np.random.RandomState(42)
    n_samples = 500

    b3_b4_dist_n8 = Counter()
    chi_dist_b3_1_n8 = Counter()
    all_betti_b3_1_n8 = Counter()

    t0 = time.time()
    for trial in range(n_samples):
        A = random_tournament(n, rng)
        res = full_chain_complex_modp(A, n, max_p=7)
        b = res['bettis']
        b3 = b.get(3, 0)
        b4 = b.get(4, 0)
        b3_b4_dist_n8[(b3, b4)] += 1

        if b3 == 1:
            betti_vec = tuple(b.get(p, 0) for p in range(8))
            all_betti_b3_1_n8[betti_vec] += 1
            chi = sum((-1)**p * b.get(p, 0) for p in range(8))
            chi_dist_b3_1_n8[chi] += 1

        if (trial+1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/{n_samples}, {elapsed:.1f}s", flush=True)

    t1 = time.time()
    print(f"\n  n=8: {n_samples} tournaments in {t1-t0:.1f}s")
    print(f"\n  (beta_3, beta_4) distribution:")
    for key, cnt in sorted(b3_b4_dist_n8.items()):
        print(f"    {key}: {cnt}")

    print(f"\n  chi(T) for beta_3=1:")
    for key, cnt in sorted(chi_dist_b3_1_n8.items()):
        print(f"    chi={key}: {cnt}")

    print(f"\n  Full Betti vectors for beta_3=1:")
    for key, cnt in sorted(all_betti_b3_1_n8.items()):
        print(f"    {key}: {cnt}")

    # Part 3: beta_5 search at n=9 (small sample)
    print("\n--- Part 3: beta_5 search at n=9 ---")
    n = 9
    rng = np.random.RandomState(42)
    n_samples = 50  # n=9 is slow

    b5_found = 0
    betti_dist_n9 = Counter()
    interesting_n9 = []

    t0 = time.time()
    for trial in range(n_samples):
        A = random_tournament(n, rng)
        res = full_chain_complex_modp(A, n, max_p=8)
        b = res['bettis']
        betti_vec = tuple(b.get(p, 0) for p in range(9))
        betti_dist_n9[betti_vec] += 1

        if b.get(5, 0) > 0:
            b5_found += 1
            interesting_n9.append((trial, betti_vec))

        if any(b.get(p, 0) > 0 for p in [3, 4, 5, 6, 7]):
            if len(interesting_n9) < 20:
                interesting_n9.append((trial, betti_vec))

        if (trial+1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/{n_samples}, {elapsed:.1f}s, b5_found={b5_found}", flush=True)

    t1 = time.time()
    print(f"\n  n=9: {n_samples} tournaments in {t1-t0:.1f}s")
    print(f"  beta_5 > 0 found: {b5_found}")
    print(f"\n  Distinct Betti vectors:")
    for key, cnt in sorted(betti_dist_n9.items()):
        print(f"    {key}: {cnt}")

    if interesting_n9:
        print(f"\n  Interesting cases (b3+ or b4+ or b5+):")
        for trial, bv in interesting_n9[:10]:
            print(f"    trial={trial}: {bv}")

    # Part 4: Conclusions
    print("\n--- Part 4: Conclusions ---")
    if all(chi == 0 for chi in chi_dist_b3_1.keys()):
        print("  CONFIRMED at n=7: chi(T) = 0 for ALL beta_3=1 tournaments")
        print("  This means chi_rel = -chi(T\\v) for all vertices")
        print("  => chi_rel dichotomy follows from chi(T\\v) values alone")
    else:
        print(f"  chi(T) for beta_3=1 at n=7: {dict(chi_dist_b3_1)}")

    b34_coexist = sum(cnt for (b3, b4), cnt in b3_b4_dist.items()
                      if b3 > 0 and b4 > 0)
    if b34_coexist == 0:
        print("  CONFIRMED at n=7: beta_3 > 0 and beta_4 > 0 NEVER coexist")
    else:
        print(f"  beta_3 and beta_4 coexist in {b34_coexist} tournaments at n=7")

    print("\nDONE.")


if __name__ == '__main__':
    main()
