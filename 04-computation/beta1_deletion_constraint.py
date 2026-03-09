"""
beta1_deletion_constraint.py — Is beta_1(T\\v) = 0 for ALL vertices of beta_3=1 tournaments?

Prediction from relative Euler characteristic:
  chi_rel = chi(T) - chi(T\\v) = -1 for good vertices, 0 for bad vertices
  chi(T) = 0 at n=6 (beta_0=1, beta_1=0, beta_2=0, beta_3=1)
  => chi(T\\v) = 1 for good vertices
  chi(T\\v) = 1 - beta_1(T\\v) at n=5
  => beta_1(T\\v) = 0 for ALL good vertices at n=6

This would mean: beta_3=1 implies ALL deletions are beta_1-free too.
Combined with the seesaw: beta_3(T)=1 => beta_1(T)=0, this gives a
"global beta_1-free" property for beta_3=1 tournaments.

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament, full_chain_complex
)


def main():
    print("=" * 70)
    print("BETA_1 DELETION CONSTRAINT FOR BETA_3=1 TOURNAMENTS")
    print("=" * 70)

    # Part 1: n=6 exhaustive
    print("\n--- Part 1: n=6 exhaustive ---")
    n = 6
    total = 2 ** (n*(n-1)//2)

    b1_del_dist_b3_1 = Counter()
    b1_del_dist_b3_0_sample = Counter()
    chi_T_b3_1 = Counter()
    chi_Tv_b3_1 = Counter()
    count_b3_1 = 0
    violations = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        data = full_chain_complex(A, n, max_p=5)
        b3 = data['bettis'].get(3, 0)

        if b3 == 1:
            count_b3_1 += 1
            bettis_T = data['bettis']
            chi_T = sum((-1)**p * bettis_T.get(p, 0) for p in range(6))
            chi_T_b3_1[chi_T] += 1

            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
                data_sub = full_chain_complex(A_sub, n-1, max_p=5)

                b1_Tv = data_sub['bettis'].get(1, 0)
                b3_Tv = data_sub['bettis'].get(3, 0)
                b1_del_dist_b3_1[b1_Tv] += 1

                bettis_Tv = data_sub['bettis']
                chi_Tv = sum((-1)**p * bettis_Tv.get(p, 0) for p in range(6))
                chi_Tv_b3_1[(b1_Tv, b3_Tv, chi_Tv)] += 1

                if b1_Tv > 0:
                    violations += 1
                    if violations <= 3:
                        print(f"  VIOLATION: bits={bits}, v={v}, b1(T\\v)={b1_Tv}, b3(T\\v)={b3_Tv}")

        elif b3 == 0 and len(b1_del_dist_b3_0_sample) < 50:
            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
                data_sub = full_chain_complex(A_sub, n-1, max_p=3)
                b1_Tv = data_sub['bettis'].get(1, 0)
                b1_del_dist_b3_0_sample[b1_Tv] += 1

        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{total}", flush=True)

    print(f"\n  n=6: {count_b3_1} tournaments with beta_3=1")
    print(f"  beta_1(T\\v) distribution for beta_3=1 tours: {dict(sorted(b1_del_dist_b3_1.items()))}")
    print(f"  Violations (b1(T\\v) > 0): {violations}")
    print(f"  chi(T) distribution: {dict(sorted(chi_T_b3_1.items()))}")
    print(f"\n  (b1_Tv, b3_Tv, chi_Tv) distribution:")
    for key, cnt in sorted(chi_Tv_b3_1.items()):
        print(f"    {key}: {cnt}")

    print(f"\n  For comparison, beta_3=0 tours (sample):")
    print(f"  beta_1(T\\v) distribution: {dict(sorted(b1_del_dist_b3_0_sample.items()))}")

    # Part 2: n=7
    print("\n--- Part 2: n=7 sampled ---")
    n = 7
    rng = np.random.RandomState(42)

    b1_del_b3_1_n7 = Counter()
    chi_data_n7 = []
    checked = 0
    violations_n7 = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        data = full_chain_complex(A, n, max_p=5)
        if data['bettis'].get(3, 0) != 1:
            continue

        checked += 1
        bettis_T = data['bettis']
        chi_T = sum((-1)**p * bettis_T.get(p, 0) for p in range(7))

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            data_sub = full_chain_complex(A_sub, n-1, max_p=5)

            b1_Tv = data_sub['bettis'].get(1, 0)
            b3_Tv = data_sub['bettis'].get(3, 0)
            b1_del_b3_1_n7[(b1_Tv, b3_Tv)] += 1

            bettis_Tv = data_sub['bettis']
            chi_Tv = sum((-1)**p * bettis_Tv.get(p, 0) for p in range(7))

            if b1_Tv > 0:
                violations_n7 += 1
                if violations_n7 <= 3:
                    score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
                    print(f"  VIOLATION: trial={trial}, v={v}, b1(T\\v)={b1_Tv}, b3(T\\v)={b3_Tv}, "
                          f"chi_T={chi_T}, chi_Tv={chi_Tv}, score={score}")

        if checked >= 100:
            break

        if checked % 20 == 0:
            print(f"  ... {checked} beta_3=1 found", flush=True)

    print(f"\n  n=7: {checked} tournaments with beta_3=1")
    print(f"  (b1_Tv, b3_Tv) distribution:")
    for key, cnt in sorted(b1_del_b3_1_n7.items()):
        print(f"    {key}: {cnt}")
    print(f"  Violations (b1(T\\v) > 0): {violations_n7}")

    # Part 3: Conclusion
    print("\n--- Part 3: Summary ---")
    if violations == 0 and violations_n7 == 0:
        print("  CONFIRMED: beta_1(T\\v) = 0 for ALL vertices of beta_3=1 tournaments!")
        print("  This is STRONGER than the seesaw (which only says b1*b3=0 per tournament)")
        print("  It says: if b3(T)=1, then ALL deletions have b1=0")
        print("  This forces H_2(T,T\\v) = 0 for ALL vertices (from LES)")
        print("  Combined with H_0^rel = 0, H_1^rel = 0, the relative complex simplifies to:")
        print("    H_p^rel = 0 for p != 3, and H_3^rel in {0,1}")
    elif violations == 0:
        print(f"  At n=6: CONFIRMED (0 violations)")
        print(f"  At n=7: {violations_n7} violations found")
    else:
        print(f"  REFUTED: {violations} violations at n=6, {violations_n7} at n=7")

    print("\nDONE.")


if __name__ == '__main__':
    main()
