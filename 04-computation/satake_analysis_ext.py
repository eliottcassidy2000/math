"""
satake_analysis_ext.py — Extended Satake NDR tournament analysis (S56).

Investigations:
  1. Full H distribution at n=13 — all 64 values
  2. Structure of maximizer S={7,...,12} — cyclic interval hypothesis
  3. Cyclic interval tournament vs Satake at all n=5,9,13 and primes q≡5(8)
  4. Satake primes q≡5(8): 5,13,29,37 — H gap from maximum
  5. H distribution at n=7 (Paley) vs n=13 (non-Paley) comparison
  6. Does cyclic interval ALWAYS maximize H for circulant tournaments?

Author: kind-pasteur-2026-03-12-S56
"""
import sys
import time
sys.path.insert(0, '04-computation')

from satake_ndrt_h import (
    satake_connection_set, quartic_cosets,
    ham_count_circulant, ham_paths_from_0,
    count_3cycles, all_circulant_H
)


def cyclic_interval_S(n):
    """S = {ceil(n/2), ..., n-1} — the 'upper half' cyclic interval."""
    k = (n - 1) // 2
    start = n - k
    return list(range(start, n))


def all_maximizers(results):
    """Return all S achieving maximum H."""
    H_max = results[0][0]
    return [(H, S) for H, S in results if H == H_max]


def H_rank(results, S_target):
    """Return (rank, H) of S_target in results."""
    S_t = sorted(S_target)
    for i, (H, S) in enumerate(results):
        if sorted(S) == S_t:
            return i + 1, H
    return None, None


def satake_primes_small():
    """Primes q ≡ 5 mod 8, q < 50."""
    from sympy import isprime
    return [q for q in range(5, 50, 8) if isprime(q)]


def main():
    print("=" * 70)
    print("SATAKE EXTENDED ANALYSIS (S56)")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Section 1: Full H distribution at n=13
    # -----------------------------------------------------------------------
    print("\n--- Full H distribution at n=13 ---")
    t0 = time.time()
    results13 = all_circulant_H(13)
    print(f"[{time.time()-t0:.1f}s, {len(results13)} total circulant tournaments]")

    from collections import Counter
    H_counts = Counter(H for H, _ in results13)
    print(f"Distinct H values: {len(H_counts)}")
    print("Full distribution (H: count):")
    for H in sorted(H_counts, reverse=True):
        print(f"  H={H}: {H_counts[H]} tournament(s)")

    # -----------------------------------------------------------------------
    # Section 2: Cyclic interval hypothesis
    # -----------------------------------------------------------------------
    print("\n--- Cyclic interval hypothesis ---")
    print("Is S = {ceil(n/2),...,n-1} always the maximizer?")

    for n in [5, 7, 9, 11, 13]:
        results = all_circulant_H(n)
        S_ci = cyclic_interval_S(n)
        rank, H = H_rank(results, S_ci)
        H_max = results[0][0]
        maxers = all_maximizers(results)
        print(f"  n={n}: cyclic_interval S={S_ci}, H={H}, rank={rank}/{len(results)}, "
              f"max={H_max}, is_max={'YES' if H==H_max else f'NO (gap={H_max-H})'}")
        if len(maxers) <= 4:
            for Hm, Sm in maxers:
                print(f"    maximizer: S={Sm}")

    # -----------------------------------------------------------------------
    # Section 3: Satake primes ≡ 5 mod 8
    # -----------------------------------------------------------------------
    print("\n--- Satake primes q == 5 (mod 8) ---")
    primes_5mod8 = satake_primes_small()
    print(f"Primes q == 5 mod 8 up to 50: {primes_5mod8}")

    for q in primes_5mod8:
        S_sat = satake_connection_set(q)
        S_ci = cyclic_interval_S(q)
        H_sat, paths_sat = ham_count_circulant(q, S_sat)
        H_ci, paths_ci = ham_count_circulant(q, S_ci)
        # For small q, get exact rank
        if q <= 13:
            results = all_circulant_H(q)
            rank_sat, _ = H_rank(results, S_sat)
            rank_ci, _ = H_rank(results, S_ci)
            H_max = results[0][0]
            print(f"  q={q}: H_sat={H_sat} (rank {rank_sat}/{len(results)}), "
                  f"H_ci={H_ci} (rank {rank_ci}/{len(results)}), H_max={H_max}")
        else:
            gap_sat_ci = H_ci - H_sat
            print(f"  q={q}: H_sat={H_sat}, H_ci={H_ci}, gap(ci-sat)={gap_sat_ci}")

    # -----------------------------------------------------------------------
    # Section 4: Structure of n=13 maximizer tournaments
    # -----------------------------------------------------------------------
    print("\n--- Structure of n=13 maximizer tournaments ---")
    results13 = all_circulant_H(13)
    maxers13 = all_maximizers(results13)
    print(f"{len(maxers13)} maximizer(s) at H_max={results13[0][0]}:")
    for H, S in maxers13:
        S_neg = sorted((13 - s) % 13 for s in S)
        is_ci = sorted(S) == sorted(cyclic_interval_S(13))
        print(f"  S={S}  -S={S_neg}  cyclic_interval={is_ci}")

    # Also check second-best
    H_2nd = sorted(set(H for H, _ in results13), reverse=True)[1]
    maxers13_2nd = [(H, S) for H, S in results13 if H == H_2nd]
    print(f"\n{len(maxers13_2nd)} tournament(s) at H_2nd={H_2nd}:")
    for H, S in maxers13_2nd:
        S_neg = sorted((13 - s) % 13 for s in S)
        print(f"  S={S}  -S={S_neg}")

    # -----------------------------------------------------------------------
    # Section 5: Coset structure — which cosets form the maximizer?
    # -----------------------------------------------------------------------
    print("\n--- Coset structure at n=13 ---")
    cosets13, g13 = quartic_cosets(13)
    print(f"Quartic cosets mod 13 (g={g13}):")
    for i, c in enumerate(cosets13):
        print(f"  C_{i} = {c}")

    S_max = sorted(results13[0][1])
    S_sat = sorted(satake_connection_set(13))
    print(f"\nMaximizer S_max = {S_max}")
    print(f"Satake S_sat   = {S_sat}")

    for i, c in enumerate(cosets13):
        overlap_max = set(c) & set(S_max)
        overlap_sat = set(c) & set(S_sat)
        print(f"  C_{i}={c}: overlap_max={sorted(overlap_max)}, overlap_sat={sorted(overlap_sat)}")

    # -----------------------------------------------------------------------
    # Section 6: Compare n=7 (Paley) vs n=13 H distributions
    # -----------------------------------------------------------------------
    print("\n--- n=7 Paley comparison ---")
    results7 = all_circulant_H(7)
    H_counts7 = Counter(H for H, _ in results7)
    print(f"n=7: {len(results7)} tournaments, {len(H_counts7)} distinct H values")
    for H in sorted(H_counts7, reverse=True):
        S_reps = [S for hh, S in results7 if hh == H][:2]
        print(f"  H={H}: {H_counts7[H]} tournaments, e.g. {S_reps[0]}")

    S_paley7 = [1, 2, 4]
    rank_p7, H_p7 = H_rank(results7, S_paley7)
    S_ci7 = cyclic_interval_S(7)
    rank_ci7, H_ci7 = H_rank(results7, S_ci7)
    print(f"  Paley S={S_paley7}: H={H_p7}, rank {rank_p7}/{len(results7)}")
    print(f"  CycInt S={S_ci7}: H={H_ci7}, rank {rank_ci7}/{len(results7)}")


if __name__ == '__main__':
    main()
