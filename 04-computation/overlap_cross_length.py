#!/usr/bin/env python3
"""
overlap_cross_length.py -- Cross-length disjointness and overlap structure

Extends the overlap weight analysis to study:
1. Cross-length disjointness: (3,5), (3,7), (5,5), (5,7) cycle pairs
2. How co_occ_k(d) from THM-143 controls cross-length disjointness
3. The alpha_j decomposition BY CYCLE LENGTH COMPOSITION
4. Overlap weight distribution as a function of cycle lengths

The key insight to test: THM-143 gives co_occ_k(d) for EACH k.
The co-occurrence CROSS-PRODUCT co_occ_{k1,k2}(d) = expected overlap
between a k1-cycle and a k2-cycle depends on BOTH co-occurrence profiles.

If co_occ_k(d) is linear with slope b_k = C(m-2, k-3), then the
slope INCREASES with k. This means larger cycles have STEEPER gradients,
amplifying the Interval advantage at higher alpha_j levels.

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from math import comb
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def has_ham_cycle(A, verts):
    """Check if vertex subset has a directed Hamiltonian cycle."""
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


def enumerate_cycle_vertex_sets(A, p, max_k=None):
    """Enumerate all vertex sets that support a directed Hamiltonian cycle, by length."""
    if max_k is None:
        max_k = p
    by_k = {}
    for k in range(3, max_k + 1, 2):
        sets_k = []
        for subset in combinations(range(p), k):
            if has_ham_cycle(A, list(subset)):
                sets_k.append(frozenset(subset))
        by_k[k] = sets_k
    return by_k


def co_occ_formula(p, k, d):
    """THM-143 formula for co_occ_k(d) for Interval tournament."""
    m = (p - 1) // 2
    if d > m:
        d = p - d  # symmetry
    return comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - (m - d) * comb(m - 2, k - 3)


def cross_length_disjointness(by_k, k1, k2):
    """Count disjoint pairs between cycles of length k1 and k2.

    Returns: (disjoint_count, total_count, overlap_distribution)
    """
    sets1 = by_k.get(k1, [])
    sets2 = by_k.get(k2, [])

    if k1 == k2:
        # Same length: only count unordered pairs
        total = len(sets1) * (len(sets1) - 1) // 2
        overlap_dist = defaultdict(int)
        disjoint = 0
        for i in range(len(sets1)):
            for j in range(i + 1, len(sets1)):
                ov = len(sets1[i] & sets1[j])
                overlap_dist[ov] += 1
                if ov == 0:
                    disjoint += 1
        return disjoint, total, dict(sorted(overlap_dist.items()))
    else:
        # Different lengths: all pairs
        total = len(sets1) * len(sets2)
        overlap_dist = defaultdict(int)
        disjoint = 0
        for s1 in sets1:
            for s2 in sets2:
                ov = len(s1 & s2)
                overlap_dist[ov] += 1
                if ov == 0:
                    disjoint += 1
        return disjoint, total, dict(sorted(overlap_dist.items()))


def cross_co_occurrence(by_k, k1, k2, p):
    """Compute co_occ_{k1,k2}(d) = #{(C1,C2) : |C1|=k1, |C2|=k2, both contain 0 and d}.

    For circulant: co_occ_{k1,k2}(d) depends only on min(d, p-d).
    """
    sets1 = by_k.get(k1, [])
    sets2 = by_k.get(k2, [])

    co_occ = [0] * p
    for d in range(1, p):
        count = 0
        for s1 in sets1:
            if 0 not in s1 or d not in s1:
                continue
            for s2 in sets2:
                if k1 == k2 and s2 <= s1:  # avoid double counting for same k
                    continue
                if 0 not in s2 or d not in s2:
                    continue
                # This is a pair of cycles both containing 0 and d
                count += 1
        co_occ[d] = count
    return co_occ


def alpha_by_composition(by_k, p):
    """Decompose alpha_j (independent set count of size j) by cycle length composition.

    alpha_2 = sum_{(k1,k2) with k1+k2<=p} #{disjoint (C1,C2) pairs with |C1|=k1, |C2|=k2}

    Returns: {(k1,k2): count} for each cycle length pair
    """
    all_ks = sorted(by_k.keys())

    # alpha_1 by length
    alpha_1_by_k = {k: len(by_k[k]) for k in all_ks}

    # alpha_2 by (k1, k2)
    alpha_2_by_kk = {}
    for i, k1 in enumerate(all_ks):
        for k2 in all_ks[i:]:
            if k1 + k2 > p:
                continue
            disj, total, _ = cross_length_disjointness(by_k, k1, k2)
            if disj > 0 or k1 + k2 <= p:
                alpha_2_by_kk[(k1, k2)] = disj

    # alpha_3: triples of mutually disjoint cycles
    # Only compute for small enough cycle collections
    all_cycles = []
    for k in all_ks:
        for fs in by_k[k]:
            all_cycles.append((fs, k))

    n_cyc = len(all_cycles)
    alpha_3_by_kkk = {}
    if n_cyc <= 200:
        for i in range(n_cyc):
            for j in range(i + 1, n_cyc):
                if all_cycles[i][0] & all_cycles[j][0]:
                    continue
                for l in range(j + 1, n_cyc):
                    if all_cycles[i][0] & all_cycles[l][0]:
                        continue
                    if all_cycles[j][0] & all_cycles[l][0]:
                        continue
                    # Mutually disjoint triple
                    ks = tuple(sorted([all_cycles[i][1], all_cycles[j][1], all_cycles[l][1]]))
                    alpha_3_by_kkk[ks] = alpha_3_by_kkk.get(ks, 0) + 1

    return alpha_1_by_k, alpha_2_by_kk, alpha_3_by_kkk


def overlap_weight_matrix_analysis(by_k, p):
    """Analyze the full overlap weight matrix restricted to 3-cycles.

    For each pair of 3-cycles, compute:
    - overlap = |V(C1) ∩ V(C2)|
    - gap structure: the gaps between the cycles on Z_p
    """
    c3 = by_k.get(3, [])
    n3 = len(c3)

    # Compute the overlap histogram by gap
    # For two 3-cycles {0,a,b} and {0,c,d}, the overlap depends on
    # whether {a,b} ∩ {c,d} is nonempty
    # But for general cycles, we need to check all rotations

    # Key idea: for circulant tournament, the overlap between C_i and C_j
    # depends on the "relative position" of their vertex sets on Z_p

    # Let's compute: for each pair of gap-type triples, what's the
    # expected overlap when placed randomly vs on Z_p?

    # First: gap types for 3-cycles
    # A 3-cycle uses 3 vertices spanning gaps (d1, d2, d3) with d1+d2+d3 = p
    # For circulant, we can normalize so one vertex is 0

    gap_types = defaultdict(list)
    for idx, fs in enumerate(c3):
        verts = sorted(fs)
        # Translate so min vertex is 0
        shifted = [v - verts[0] for v in verts]
        gaps = tuple(sorted([shifted[1], shifted[2] - shifted[1], p - shifted[2]]))
        gap_types[gaps].append(idx)

    return gap_types


def interval_slope_amplification(p):
    """Compute the co-occurrence slope b_k = C(m-2, k-3) for each k.

    The slope measures how fast co-occurrence drops with distance d.
    If b_k increases with k, then LARGER cycles have STEEPER gradients,
    meaning the Interval advantage is AMPLIFIED at higher alpha_j levels.
    """
    m = (p - 1) // 2
    slopes = {}
    for k in range(3, p + 1, 2):
        slopes[k] = comb(m - 2, k - 3)
    return slopes


def co_occ_variance_by_k(p):
    """Compute Var(co_occ_k) over d for each k using THM-143 formula.

    co_occ_k(d) = a_k - b_k * d  for 1 <= d <= m
    Var over d = b_k^2 * Var(d) = b_k^2 * (m-1)(m+1)/12

    This variance drives disjointness excess.
    """
    m = (p - 1) // 2
    var_d = (m - 1) * (m + 1) / 12  # variance of uniform on {1,...,m}

    result = {}
    for k in range(3, p + 1, 2):
        b_k = comb(m - 2, k - 3)
        # co_occ values
        co_occ_vals = [co_occ_formula(p, k, d) for d in range(1, m + 1)]
        mean = sum(co_occ_vals) / len(co_occ_vals)
        var = sum((v - mean)**2 for v in co_occ_vals) / len(co_occ_vals)

        # Theoretical variance = b_k^2 * Var(uniform on {1,...,m})
        var_theory = b_k**2 * var_d

        result[k] = {
            'slope': b_k,
            'mean': mean,
            'var': var,
            'var_theory': var_theory,
            'min': min(co_occ_vals),
            'max': max(co_occ_vals)
        }
    return result


def disjointness_excess_by_k(p):
    """Predict the disjointness excess for (k,k)-pairs using the co-occurrence structure.

    THM-142 shows that for k=3: disjointness excess = p(p-1)(p+1)(p-3)/192
    Can we generalize this to arbitrary k?

    The disjointness excess comes from the co-occurrence variance:
    More variance in co_occ => more bimodal overlap distribution => more disjoint pairs.

    For level-2 inclusion-exclusion:
    alpha_2 approx = C(n_k, 2) - sum_d p * co_occ_k(d) * (co_occ_k(d) - 1) / 2 * correction

    Actually, let's think about it differently.
    For two randomly chosen k-cycle vertex sets through vertex 0:
    P(disjoint | both through 0) depends on co_occ_k.

    Total disjoint k-k pairs = (total k-k pairs) * P(disjoint)
    P(disjoint) = 1 - P(share at least one vertex)

    For Paley (constant co_occ): P(share vertex d) = co_occ^2 / n_k^2
    For Interval (varying co_occ): P(share vertex d) = (co_occ(d))^2 / n_k^2

    By Jensen's inequality on convex function x^2:
    E[co_occ^2] >= (E[co_occ])^2
    So Interval has MORE co-occurrence on some pairs (overlap=2) and LESS on others (overlap=0)

    This is exactly the bimodal redistribution seen in general_disjointness_excess.py!
    """
    m = (p - 1) // 2

    result = {}
    for k in range(3, min(p, 14), 2):  # limit to prevent huge computation
        # Count n_k (number of k-cycle vertex sets) using THM-143
        # n_k = p * co_occ_k(d) / k for any d... no, that's not right
        # n_k = p * (sum of co_occ_k(d) for d=1..p-1) / (k*(k-1))
        # Because each k-cycle contains C(k,2) vertex pairs, each contributing to co_occ

        co_occ_vals = [co_occ_formula(p, k, d) for d in range(1, m + 1)]
        # Total co-occurrence = sum over d of co_occ_k(d)
        # By circulant symmetry, co_occ(d) = co_occ(p-d), so total = 2 * sum_{d=1}^m co_occ(d)
        total_co_occ = 2 * sum(co_occ_vals)

        # Each k-cycle has k vertices, and each pair (0, d) is counted once
        # So n_k = total_co_occ / (k-1) + 1 ??? No...
        # Actually n_k * (k-1) = total_co_occ (each cycle through 0 has k-1 other vertices)
        # Wait: co_occ_k(d) = #{k-cycles through BOTH 0 and d}
        # sum_d co_occ_k(d) = #{(C, d) : C is k-cycle through 0, d in C, d != 0}
        # = sum_{C through 0} (|C| - 1) = n_{k,0} * (k-1)
        # where n_{k,0} = #{k-cycles through vertex 0} = n_k * k / p by circulant
        # So total_co_occ = n_{k,0} * (k-1) = n_k * k * (k-1) / p
        # Therefore n_k = p * total_co_occ / (k * (k-1))

        n_k = p * total_co_occ // (k * (k - 1))
        # Verify this is an integer
        n_k_exact = p * total_co_occ / (k * (k - 1))

        # For Paley, co_occ is constant, so:
        co_occ_paley = co_occ_vals[0]  # would be constant for Paley
        # Actually we need to compute Paley co_occ separately
        # For now, use the MEAN co_occ as Paley approximation
        mean_co_occ = sum(co_occ_vals) / len(co_occ_vals)

        # The co_occ variance
        var_co_occ = sum((v - mean_co_occ)**2 for v in co_occ_vals) / len(co_occ_vals)

        # Disjointness excess prediction (first-order):
        # Excess disjoint pairs ~ p * var_co_occ * (something)
        # From THM-142 for k=3: excess = p * (p-1) * var_co_occ * 2 / ???

        b_k = comb(m - 2, k - 3)

        result[k] = {
            'n_k': n_k,
            'n_k_exact': n_k_exact,
            'mean_co_occ': mean_co_occ,
            'var_co_occ': var_co_occ,
            'slope': b_k,
            'co_occ_range': (min(co_occ_vals), max(co_occ_vals)),
        }
    return result


def main():
    print("=" * 70)
    print("CROSS-LENGTH DISJOINTNESS AND OVERLAP STRUCTURE")
    print("=" * 70)

    # ====== SLOPE AMPLIFICATION ======
    print("\n" + "=" * 70)
    print("1. CO-OCCURRENCE SLOPE BY CYCLE LENGTH (THM-143)")
    print("=" * 70)

    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        slopes = interval_slope_amplification(p)
        print(f"\n  p={p}, m={m}:")
        print(f"    {'k':>4} {'b_k=C(m-2,k-3)':>15} {'C(2m-1,k-2)':>12} {'a_k':>12} {'co_occ(1)':>10} {'co_occ(m)':>10}")
        for k in range(3, min(p + 1, 18), 2):
            b_k = slopes[k]
            a_k = comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - m * b_k
            co_1 = co_occ_formula(p, k, 1)
            co_m = co_occ_formula(p, k, m)
            print(f"    {k:>4} {b_k:>15} {comb(2*m-1,k-2):>12} {a_k:>12} {co_1:>10} {co_m:>10}")

        # Slope growth rate
        if m >= 3:
            ratios = []
            for k in range(5, min(p + 1, 18), 2):
                prev = slopes[k - 2]
                if prev > 0:
                    ratios.append(slopes[k] / prev)
            if ratios:
                print(f"    Slope ratio b_{{k+2}}/b_k: {[f'{r:.2f}' for r in ratios]}")

    # ====== CO-OCCURRENCE VARIANCE BY k ======
    print("\n" + "=" * 70)
    print("2. CO-OCCURRENCE VARIANCE BY CYCLE LENGTH")
    print("=" * 70)

    for p in [7, 11, 13, 17, 19]:
        m = (p - 1) // 2
        var_data = co_occ_variance_by_k(p)
        print(f"\n  p={p}, m={m}:")
        print(f"    {'k':>4} {'slope':>8} {'mean':>10} {'var':>12} {'var/mean^2':>12} {'range':>20}")
        for k in sorted(var_data):
            d = var_data[k]
            vm2 = d['var'] / d['mean']**2 if d['mean'] > 0 else 0
            print(f"    {k:>4} {d['slope']:>8} {d['mean']:>10.1f} {d['var']:>12.1f} {vm2:>12.4f} [{d['min']}-{d['max']}]")

    # ====== COMPUTATIONAL VERIFICATION ======
    print("\n" + "=" * 70)
    print("3. COMPUTATIONAL CROSS-LENGTH DISJOINTNESS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)

            print(f"\n  p={p}, {name}, S={S}")
            print(f"  {'-'*60}")

            t0 = time.time()
            by_k = enumerate_cycle_vertex_sets(A, p)
            t1 = time.time()

            total_cycles = sum(len(v) for v in by_k.values())
            print(f"  Total cycle vertex sets: {total_cycles} ({t1-t0:.1f}s)")
            for k in sorted(by_k):
                print(f"    c_{k} = {len(by_k[k])}")

            # Cross-length disjointness
            ks = sorted(by_k.keys())
            print(f"\n  Cross-length disjointness:")
            print(f"    {'(k1,k2)':>10} {'disj':>8} {'total':>10} {'ratio':>8} {'max_k1+k2':>12}")
            for i, k1 in enumerate(ks):
                for k2 in ks[i:]:
                    if k1 + k2 > p:
                        print(f"    ({k1},{k2}):>10 {'N/A':>8} {'N/A':>10} {'N/A':>8} {k1+k2:>12} > p")
                        continue
                    disj, total, ov_dist = cross_length_disjointness(by_k, k1, k2)
                    ratio = disj / total if total > 0 else 0
                    print(f"    ({k1},{k2}){' ':>{7-len(f'({k1},{k2})')}} {disj:>8} {total:>10} {ratio:>8.4f} {k1+k2:>12}")
                    if len(ov_dist) <= 6:
                        print(f"              overlap dist: {ov_dist}")

            # Alpha decomposition by composition
            print(f"\n  Alpha decomposition by cycle composition:")
            a1, a2, a3 = alpha_by_composition(by_k, p)

            print(f"    alpha_1 = {sum(a1.values())}:")
            for k, cnt in sorted(a1.items()):
                print(f"      k={k}: {cnt}")

            a2_total = sum(a2.values())
            print(f"    alpha_2 = {a2_total}:")
            for (k1, k2), cnt in sorted(a2.items()):
                if cnt > 0:
                    print(f"      ({k1},{k2}): {cnt}")

            if a3:
                a3_total = sum(a3.values())
                print(f"    alpha_3 = {a3_total}:")
                for ks_tuple, cnt in sorted(a3.items()):
                    if cnt > 0:
                        print(f"      {ks_tuple}: {cnt}")

            # Cross-length co-occurrence
            if total_cycles <= 200:
                print(f"\n  Cross-length co-occurrence co_occ_{{k1,k2}}(d):")
                for k1 in ks:
                    for k2 in ks:
                        if k2 < k1:
                            continue
                        co_occ = cross_co_occurrence(by_k, k1, k2, p)
                        vals = [co_occ[d] for d in range(1, m + 1)]
                        print(f"    ({k1},{k2}): {vals}")

    # ====== DISJOINTNESS EXCESS PREDICTIONS ======
    print("\n" + "=" * 70)
    print("4. DISJOINTNESS EXCESS PREDICTIONS BY CYCLE LENGTH")
    print("=" * 70)

    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")
        excess_data = disjointness_excess_by_k(p)

        print(f"    {'k':>4} {'n_k':>10} {'mean_co':>10} {'var_co':>12} {'slope':>8} {'range':>20}")
        for k in sorted(excess_data):
            d = excess_data[k]
            print(f"    {k:>4} {d['n_k']:>10.0f} {d['mean_co_occ']:>10.1f} {d['var_co_occ']:>12.2f} {d['slope']:>8} [{d['co_occ_range'][0]}-{d['co_occ_range'][1]}]")

        # Key ratio: var_co_occ(k) / mean_co_occ(k)^2
        # This is the "coefficient of variation squared" which controls the
        # magnitude of the disjointness excess relative to the mean
        print(f"\n    Coefficient of variation CV^2 = var/mean^2:")
        for k in sorted(excess_data):
            d = excess_data[k]
            if d['mean_co_occ'] > 0:
                cv2 = d['var_co_occ'] / d['mean_co_occ']**2
                print(f"      k={k}: CV^2 = {cv2:.6f}")

    # ====== INTERVAL vs PALEY: PREDICTED vs ACTUAL CYCLE COUNTS ======
    print("\n" + "=" * 70)
    print("5. INTERVAL vs PALEY CYCLE COUNTS AND ALPHA CONTRIBUTIONS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        print(f"\n  p={p}:")

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)
            by_k = enumerate_cycle_vertex_sets(A, p)

            all_cycles = []
            for k in sorted(by_k):
                for fs in by_k[k]:
                    all_cycles.append(fs)

            n = len(all_cycles)

            # Compute alpha_1, alpha_2
            alpha_1 = n
            alpha_2 = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if not (all_cycles[i] & all_cycles[j]):
                        alpha_2 += 1

            # Compute alpha_3 if feasible
            alpha_3 = 0
            if n <= 200:
                for i in range(n):
                    for j in range(i + 1, n):
                        if all_cycles[i] & all_cycles[j]:
                            continue
                        for l in range(j + 1, n):
                            if not (all_cycles[i] & all_cycles[l]) and \
                               not (all_cycles[j] & all_cycles[l]):
                                alpha_3 += 1

            H_partial = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
            print(f"\n    {name}: alpha_1={alpha_1}, alpha_2={alpha_2}, alpha_3={alpha_3}")
            print(f"      2*alpha_1 = {2*alpha_1}")
            print(f"      4*alpha_2 = {4*alpha_2}")
            print(f"      8*alpha_3 = {8*alpha_3}")
            print(f"      H_3 (lower bound) = {H_partial}")

    # ====== KEY INSIGHT: SLOPE GROWTH RATE ======
    print("\n" + "=" * 70)
    print("6. KEY INSIGHT: SLOPE GROWTH IMPLIES AMPLIFIED ADVANTAGE")
    print("=" * 70)

    for p in [7, 11, 13, 17, 19, 23, 29, 31]:
        m = (p - 1) // 2
        slopes = interval_slope_amplification(p)

        # The slope b_k = C(m-2, k-3)
        # For k=3: b_3 = 1
        # For k=5: b_5 = C(m-2, 2) = (m-2)(m-3)/2
        # For k=7: b_7 = C(m-2, 4) = (m-2)(m-3)(m-4)(m-5)/24
        # Ratio b_5/b_3 = (m-2)(m-3)/2
        # Ratio b_7/b_5 = (m-4)(m-5)/12

        b3 = slopes.get(3, 0)
        b5 = slopes.get(5, 0)
        b7 = slopes.get(7, 0)

        r53 = b5 / b3 if b3 > 0 else 0
        r75 = b7 / b5 if b5 > 0 else 0

        # The slope controls variance which controls disjointness excess
        # Higher slope => more variance => more disjoint pairs => higher alpha_2 contribution
        # If slope grows super-linearly in k, then higher-order alpha dominate at large p

        print(f"  p={p:>3}, m={m:>3}: b_3={b3:>6}, b_5={b5:>8}, b_7={b7:>10}, "
              f"b_5/b_3={r53:>8.1f}, b_7/b_5={r75:>8.2f}")

    print("\n  Key: b_{k+2}/b_k grows with m, meaning the co-occurrence gradient")
    print("  steepens FASTER for larger cycles at larger primes.")
    print("  This explains why alpha_2+ (disjoint pairs from larger cycles)")
    print("  eventually overwhelm alpha_1 (total cycles) as p grows.")


if __name__ == '__main__':
    main()
