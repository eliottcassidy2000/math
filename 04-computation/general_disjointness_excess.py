#!/usr/bin/env python3
"""
general_disjointness_excess.py -- Extend THM-142 to all cycle lengths k

THM-143 gives: co_occ_k(d) = a_k + b_k * d  for Interval
                co_occ_k(d) = c_k (constant)  for Paley

The disjointness excess between Interval and Paley comes from the
#{ov >= 2} pairs term in inclusion-exclusion.

For k=3 (triangles):
  Two k-sets share >= 2 vertices iff they share exactly 2 (can't share 3 w/o being equal).
  #{ov>=2} = sum over (vertex, distance d) of C(co_occ(d), 2)
           = p * sum_{d=1}^m C(co_occ(d), 2)

For k=5,7,...:
  Two k-sets can share 2, 3, ..., k-1 vertices.
  #{ov>=2} = sum_{j=2}^{k-1} #{pairs sharing exactly j vertices} * 1
  But the co_occ gives only the PAIRWISE vertex overlap count.
  We need the FULL inclusion-exclusion for disjointness:

  disjoint = C(n,2) - sum_v C(n_v, 2) + sum_{v<w} C(n_{v,w}, 2) - ...

  where n_{v,w,...} = #{k-sets containing all of {v,w,...}}.

  For the DIFFERENCE (Interval - Paley), many terms cancel (same n, same n_v).
  The first non-cancelling term involves the co_occ = n_{0,d}.

This script:
1. Computes exact disjointness excess by brute force for small primes
2. Checks if the co_occ-level I-E accounts for the excess
3. Derives the general formula using THM-143

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations
from math import comb
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def has_ham_cycle(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a] + A[a][c] * A[c][b] * A[b][a]) > 0
    dp = set()
    dp.add((1 << 0, 0))
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    dp.add((mask | (1 << w), w))
    full = (1 << k) - 1
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            return True
    return False


def get_cycle_sets(A, p, k):
    """Get all k-vertex sets supporting a Ham cycle."""
    result = []
    for subset in combinations(range(p), k):
        if has_ham_cycle(A, list(subset)):
            result.append(frozenset(subset))
    return result


def count_disjoint_pairs(cycle_sets):
    n = len(cycle_sets)
    disjoint = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not (cycle_sets[i] & cycle_sets[j]):
                disjoint += 1
    return disjoint


def overlap_distribution(cycle_sets):
    """Count pairs by overlap size."""
    n = len(cycle_sets)
    dist = defaultdict(int)
    for i in range(n):
        for j in range(i + 1, n):
            ov = len(cycle_sets[i] & cycle_sets[j])
            dist[ov] += 1
    return dict(dist)


def co_occ_formula(p, k, d):
    """THM-143 exact formula."""
    m = (p - 1) // 2
    d_eff = min(d, p - d)
    return comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - (m - d_eff) * comb(m - 2, k - 3)


def main():
    print("=" * 70)
    print("GENERAL DISJOINTNESS EXCESS ANALYSIS")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        paley_exists = (p % 4 == 3)
        if not paley_exists:
            continue

        A_pal = build_adj(p, S_qr)
        A_int = build_adj(p, S_int)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"  Paley S={S_qr}")
        print(f"  Interval S={S_int}")
        print(f"{'='*70}")

        for k in range(3, min(p, 10), 2):
            if 2 * k > p:
                continue

            print(f"\n--- k={k} ---")

            t0 = time.time()
            cs_pal = get_cycle_sets(A_pal, p, k)
            cs_int = get_cycle_sets(A_int, p, k)
            t1 = time.time()

            n_pal = len(cs_pal)
            n_int = len(cs_int)

            # Disjoint pairs
            disj_pal = count_disjoint_pairs(cs_pal)
            disj_int = count_disjoint_pairs(cs_int)
            excess = disj_int - disj_pal
            t2 = time.time()

            print(f"  Paley: {n_pal} vertex sets, {disj_pal} disjoint pairs")
            print(f"  Interval: {n_int} vertex sets, {disj_int} disjoint pairs")
            print(f"  Excess: {excess} ({t2-t0:.1f}s)")

            # Overlap distribution
            ov_pal = overlap_distribution(cs_pal)
            ov_int = overlap_distribution(cs_int)

            print(f"\n  Overlap distribution:")
            all_sizes = sorted(set(list(ov_pal.keys()) + list(ov_int.keys())))
            for s in all_sizes:
                print(f"    ov={s}: Paley={ov_pal.get(s,0):>6}, "
                      f"Interval={ov_int.get(s,0):>6}, "
                      f"diff={ov_int.get(s,0)-ov_pal.get(s,0):>6}")

            # Inclusion-exclusion analysis
            print(f"\n  Inclusion-Exclusion decomposition:")

            # Level 0: C(n, 2) [total pairs]
            total_pal = comb(n_pal, 2)
            total_int = comb(n_int, 2)
            print(f"    Level 0: C(n,2): Pal={total_pal}, Int={total_int}, "
                  f"diff={total_int-total_pal}")

            # Level 1: sum_v C(n_v, 2) [pairs sharing vertex v]
            # For circulant: n_v = n_0 for all v, so sum = p * C(n_0, 2)
            n0_pal = n_pal * k // p
            n0_int = n_int * k // p
            lev1_pal = p * comb(n0_pal, 2)
            lev1_int = p * comb(n0_int, 2)
            print(f"    Level 1: p*C(n_0,2): Pal={lev1_pal}, Int={lev1_int}, "
                  f"diff={lev1_int-lev1_pal}")
            print(f"      n_0: Pal={n0_pal}, Int={n0_int}")

            # Level 2: sum_{v<w} C(n_{vw}, 2) [pairs sharing vertex pair {v,w}]
            # For circulant: this is p * sum_{d=1}^m C(co_occ(d), 2)
            lev2_pal_co = n0_pal * (k - 1) // (p - 1)  # constant co_occ for Paley
            lev2_pal = p * m * comb(lev2_pal_co, 2)

            lev2_int = 0
            for d in range(1, m + 1):
                c = co_occ_formula(p, k, d)
                lev2_int += comb(c, 2)
            lev2_int *= p

            print(f"    Level 2: p*sum C(co_occ(d),2): Pal={lev2_pal}, Int={lev2_int}, "
                  f"diff={lev2_int-lev2_pal}")

            # Check: does I-E truncated at level 2 give disjoint pairs?
            ie2_pal = total_pal - lev1_pal + lev2_pal
            ie2_int = total_int - lev1_int + lev2_int
            ie2_excess = ie2_int - ie2_pal

            print(f"\n    I-E level 2 prediction:")
            print(f"      Paley: {ie2_pal} (actual {disj_pal}, err={ie2_pal-disj_pal})")
            print(f"      Interval: {ie2_int} (actual {disj_int}, err={ie2_int-disj_int})")
            print(f"      Excess: {ie2_excess} (actual {excess}, err={ie2_excess-excess})")

            # For k=3, level 2 should be EXACT
            if k == 3:
                assert ie2_pal == disj_pal, "k=3 Paley I-E should be exact!"
                assert ie2_int == disj_int, "k=3 Interval I-E should be exact!"
                print(f"      k=3: I-E level 2 is EXACT (confirmed)")

            # Level 3 (if k >= 5): sum over vertex triples
            if k >= 5:
                # For circulant: p * sum_{d1<d2} C(co_occ_triple(d1,d2), 2)
                # co_occ_triple(d1,d2) = #{k-sets through {0, d1, d2}}
                # We compute this brute-force for validation
                lev3_val = 0
                for tri in combinations(range(p), 3):
                    # Count k-sets containing all 3 vertices
                    count = 0
                    for fs in cs_int:
                        if all(v in fs for v in tri):
                            count += 1
                    lev3_val += comb(count, 2)

                lev3_pal_val = 0
                for tri in combinations(range(p), 3):
                    count = 0
                    for fs in cs_pal:
                        if all(v in fs for v in tri):
                            count += 1
                    lev3_pal_val += comb(count, 2)

                ie3_int = total_int - lev1_int + lev2_int - lev3_val
                ie3_pal = total_pal - lev1_pal + lev2_pal - lev3_pal_val

                print(f"\n    Level 3: sum_triples C(n_triple, 2):")
                print(f"      Paley={lev3_pal_val}, Interval={lev3_val}")
                print(f"    I-E level 3 prediction:")
                print(f"      Paley: {ie3_pal} (actual {disj_pal}, err={ie3_pal-disj_pal})")
                print(f"      Interval: {ie3_int} (actual {disj_int}, err={ie3_int-disj_int})")
                print(f"      Excess: {ie3_int-ie3_pal} (actual {excess})")

    # ====== ANALYTICAL EXCESS FROM THM-143 ======
    print(f"\n{'='*70}")
    print("ANALYTICAL DISJOINTNESS EXCESS FROM THM-143")
    print("=" * 70)

    print("\nFor k=3, THM-142 excess = p(p-1)(p+1)(p-3)/192:")
    for p in [7, 11, 19, 23]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2

        # Interval: sum C(d, 2) for d=1..m
        s_int = sum(comb(d, 2) for d in range(1, m + 1))

        # Paley: constant co_occ = (p-1)/4
        c_pal = (p - 1) // 4
        s_pal = m * comb(c_pal, 2)

        # Level-2 difference only (for k=3, this IS the excess):
        # But we also need to account for Level-0 and Level-1 differences.
        # For k=3 with same c3: n_pal = n_int = C(p,3)/something... NOT necessarily equal.
        # Actually, for Interval vs Paley, the NUMBER of 3-cycles can differ!

        # Total 3-cycles (vertex sets)
        n_int = p * m * (m + 1) // 6  # from sum co_occ: n_0 = sum d = m(m+1)/2, n = n_0*p/3
        # Actually n_0 = sum_{d=1}^{m} co_occ(d) / 2 = (m)(m+1)/2 / 2... no.
        # n_0 = sum_{d=1}^{p-1} co_occ(d) / (k-1) = 2*sum_{d=1}^m d / 2 = sum_{d=1}^m d = m(m+1)/2
        n0_int = m * (m + 1) // 2
        n_int = n0_int * p // 3

        # Paley: n_0 = (p-1)/4 * (p-1) / 2 = (p-1)^2/8? No.
        # For Paley, co_occ = (p-1)/4 = m/2. n_0 = m/2 * (p-1) / (k-1) = m*(p-1)/4
        # Wait: sum co_occ(d) for d=1..p-1 = c_pal * (p-1) = (p-1)^2/4.
        # n_0 = sum / (k-1) = (p-1)^2 / 8.
        # n = n_0 * p / 3 = p(p-1)^2 / 24.
        n0_pal = (p - 1) ** 2 // 8
        n_pal_from_formula = n0_pal * p // 3

        # Or just use: c3 = p*(p-1)/6 for both when they have same c3.
        # But Paley and Interval DON'T have the same c3 in general!
        c3_int = n_int
        c3_pal = n_pal_from_formula

        print(f"\n  p={p}: c3_int={c3_int}, c3_pal={c3_pal}")

    # ====== SLOPE DETERMINES EXCESS SCALING ======
    print(f"\n{'='*70}")
    print("EXCESS SCALING FROM SLOPE")
    print("=" * 70)

    print("\nThe key quantity: sum_{d=1}^m C(a+bd, 2) - m*C(c, 2)")
    print("where a,b from THM-143 (Interval) and c = constant co_occ (Paley)")
    print()

    for k in [3, 5, 7]:
        print(f"k={k}:")
        for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
            if p % 4 != 3:
                continue
            m = (p - 1) // 2
            b = comb(m - 2, k - 3)
            a = comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - m * b

            # Interval sum
            s_int = sum(comb(a + b * d, 2) for d in range(1, m + 1))

            # Paley co_occ (constant) — need to know the value
            # For Paley vertex-transitive + edge-transitive:
            # co_occ_k_Paley = n_0_Paley * (k-1) / (p-1)
            # But we don't have n_0_Paley without computing...
            # For the EXCESS, what matters is the variance contribution from the slope b.
            # The extra #{ov>=2} from Interval's gradient is:
            # p * [sum C(a+bd,2) - m*C(a + b*(m+1)/2, 2)]  approximately
            # = p * b^2 * sum d^2 terms / 2 approximately

            # Actually, for exact comparison we need Paley's co_occ.
            # Let's compute the Interval's VARIANCE contribution:
            # V = sum (co_occ(d) - mean)^2 / m
            mean = a + b * (m + 1) / 2
            var = sum((a + b*d - mean)**2 for d in range(1, m + 1)) / m

            print(f"  p={p}: slope={b}, var_co_occ={var:.1f}, "
                  f"p*var={p*var:.1f}")


if __name__ == '__main__':
    main()
