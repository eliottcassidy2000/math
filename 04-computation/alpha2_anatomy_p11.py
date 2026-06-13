#!/usr/bin/env python3
"""
alpha2_anatomy_p11.py -- OCF cancellation anatomy at p=11

Extends alpha2_walsh_anatomy.py to p=11.
At p=11: m=5, 32 orientations, max 3 disjoint 3-cycles.
H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3

Uses Held-Karp DP for CORRECT simple cycle enumeration (no Fourier overcounting).
Only enumerates 3-cycles and 5-cycles for speed (7,9,11-cycles are expensive).
For full alpha_j we need ALL odd cycles -- try a staged approach.

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations
from collections import defaultdict
import time
import sys


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def count_ham_cycles(A, verts):
    """Count directed Hamiltonian cycles on verts using Held-Karp DP."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def get_all_odd_cycles(A, p, max_k=None):
    """Get all directed odd cycles as (vertex_set, multiplicity) pairs.

    Returns list of frozensets, with each frozenset repeated by its
    Hamiltonian cycle count (i.e., if a 5-vertex set has 2 distinct
    directed 5-cycles, it appears twice).
    """
    if max_k is None:
        max_k = p
    all_cycles = []
    by_len = defaultdict(int)

    for k in range(3, max_k + 1, 2):
        count_k = 0
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            for _ in range(n_cyc):
                all_cycles.append(frozenset(subset))
            count_k += n_cyc
        by_len[k] = count_k

    return all_cycles, dict(by_len)


def count_independent_sets(cycles, max_j=None):
    """Count independent sets of size j in the conflict graph.

    alpha_j = number of j-element subsets of cycles that are pairwise vertex-disjoint.
    """
    n = len(cycles)
    if max_j is None:
        max_j = n

    # For small cycle counts, brute force
    # alpha_1 = n (every single cycle is an independent set of size 1)
    alpha = {0: 1, 1: n}

    # alpha_2: pairs of disjoint cycles
    alpha_2 = 0
    for i in range(n):
        for j in range(i+1, n):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    alpha[2] = alpha_2

    # alpha_3: triples of pairwise disjoint cycles
    # At p=11, max 3 disjoint 3-cycles (using 9 of 11 vertices)
    alpha_3 = 0
    for i in range(n):
        for j in range(i+1, n):
            if cycles[i] & cycles[j]:
                continue
            ij = cycles[i] | cycles[j]
            for k in range(j+1, n):
                if not (ij & cycles[k]):
                    alpha_3 += 1
    alpha[3] = alpha_3

    return alpha


def main():
    p = 11
    m = (p - 1) // 2  # 5
    pairs = [(s, p - s) for s in range(1, m + 1)]
    half = 1 << m  # 32

    print("=" * 70)
    print(f"ALPHA WALSH ANATOMY at p={p}")
    print(f"m={m}, {half} orientations")
    print("=" * 70)

    # First pass: enumerate all cycles for each orientation
    # At p=11, we need k=3,5,7,9,11
    # k=3: C(11,3)=165 subsets, fast
    # k=5: C(11,5)=462 subsets, each DP is 2^5=32 states, fast
    # k=7: C(11,7)=330 subsets, each DP is 2^7=128 states, moderate
    # k=9: C(11,9)=55 subsets, each DP is 2^9=512 states, moderate
    # k=11: C(11,11)=1 subset, DP is 2^11=2048 states, fast
    # Total per orientation: ~1013 subsets with small DP => should be fast

    data = {}
    t0 = time.time()

    for bits in range(half):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]

        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1

        all_cycles, by_len = get_all_odd_cycles(A, p)
        alpha = count_independent_sets(all_cycles)

        H = sum((1 << (j+1)) * alpha.get(j, 0) for j in range(1, max(alpha.keys())+1 if alpha else 1))
        H += 1  # alpha_0 = 1 contributes 2^0 = 1... wait
        # H = I(Omega, 2) = sum_j alpha_j * 2^j
        H = sum(alpha.get(j, 0) * (2**j) for j in range(max(alpha.keys())+1 if alpha else 1))

        data[bits] = {
            'sigma': sigma,
            'S': S,
            'all_cycles': all_cycles,
            'by_len': by_len,
            'alpha': alpha,
            'H': H,
        }

        elapsed = time.time() - t0
        print(f"\n  [{bits+1}/{half}, {elapsed:.1f}s] bits={bits:05b}, sigma={sigma}, S={S}")
        print(f"    by_len: {dict(sorted(by_len.items()))}")
        print(f"    alpha: {dict(sorted(alpha.items()))}")
        print(f"    H = {H}")
        sys.stdout.flush()

    # Summary
    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)

    H_vals = [data[bits]['H'] for bits in range(half)]
    print(f"  H values: {H_vals}")
    print(f"  H unique: {sorted(set(H_vals))}")
    print(f"  H mean: {sum(H_vals)/len(H_vals):.4f}")

    alpha1_vals = [data[bits]['alpha'][1] for bits in range(half)]
    alpha2_vals = [data[bits]['alpha'][2] for bits in range(half)]
    alpha3_vals = [data[bits]['alpha'].get(3, 0) for bits in range(half)]
    print(f"\n  alpha_1 values: {alpha1_vals}")
    print(f"  alpha_1 unique: {sorted(set(alpha1_vals))}")
    print(f"  alpha_2 values: {alpha2_vals}")
    print(f"  alpha_2 unique: {sorted(set(alpha2_vals))}")
    print(f"  alpha_3 values: {alpha3_vals}")
    print(f"  alpha_3 unique: {sorted(set(alpha3_vals))}")

    # Cycle counts by length
    for k in [3, 5, 7, 9, 11]:
        vals = [data[bits]['by_len'].get(k, 0) for bits in range(half)]
        print(f"\n  c_{k} values: {vals}")
        print(f"  c_{k} unique: {sorted(set(vals))}")

    # Walsh decomposition
    print(f"\n\n{'='*70}")
    print("WALSH DECOMPOSITION (degree 2)")
    print("=" * 70)

    for name, vals in [("H", H_vals),
                       ("alpha_1", alpha1_vals),
                       ("alpha_2", alpha2_vals),
                       ("alpha_3", alpha3_vals)]:
        print(f"\n  {name}:")
        mean = sum(vals) / half
        print(f"    Mean: {mean:.4f}")

        for a in range(m):
            for b in range(a + 1, m):
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += vals[bits] * sa * sb
                h = total / half
                if abs(h) > 0.001:
                    gap_a, gap_b = a + 1, b + 1
                    chi_ab = legendre(gap_a * gap_b, p)
                    sign_ok = (h > 0) == (chi_ab > 0) if h != 0 else True
                    print(f"    h_hat[{{{a},{b}}}] = {h:>12.4f}, "
                          f"gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}, "
                          f"sign {'OK' if sign_ok else 'FAIL'}")

    # OCF decomposition
    print(f"\n\n{'='*70}")
    print("OCF DECOMPOSITION: h_hat_H = 2*h_hat_a1 + 4*h_hat_a2 + 8*h_hat_a3")
    print("=" * 70)

    for a in range(m):
        for b in range(a + 1, m):
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)

            h_vals_ab = {}
            for name, vals in [("alpha_1", alpha1_vals),
                               ("alpha_2", alpha2_vals),
                               ("alpha_3", alpha3_vals),
                               ("H", H_vals)]:
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += vals[bits] * sa * sb
                h_vals_ab[name] = total / half

            h_a1, h_a2, h_a3, h_H = (h_vals_ab["alpha_1"], h_vals_ab["alpha_2"],
                                       h_vals_ab["alpha_3"], h_vals_ab["H"])

            recon = 2*h_a1 + 4*h_a2 + 8*h_a3

            print(f"\n  Pair ({a},{b}), gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}:")
            print(f"    h_hat_alpha1 = {h_a1:>12.4f} (sign {'+' if h_a1>0 else '-'})")
            print(f"    h_hat_alpha2 = {h_a2:>12.4f} (sign {'+' if h_a2>0 else '-'})")
            print(f"    h_hat_alpha3 = {h_a3:>12.4f} (sign {'+' if h_a3>0 else '-'})")
            print(f"    2*h_hat_a1   = {2*h_a1:>12.4f}")
            print(f"    4*h_hat_a2   = {4*h_a2:>12.4f}")
            print(f"    8*h_hat_a3   = {8*h_a3:>12.4f}")
            print(f"    h_hat_H      = {h_H:>12.4f} (reconstructed = {recon:>12.4f})")

            total_contrib = abs(2*h_a1) + abs(4*h_a2) + abs(8*h_a3)
            if total_contrib > 0:
                print(f"    Contribution: a1={abs(2*h_a1)/total_contrib*100:.1f}% "
                      f"a2={abs(4*h_a2)/total_contrib*100:.1f}% "
                      f"a3={abs(8*h_a3)/total_contrib*100:.1f}%")

    # Degree-4 Walsh analysis
    print(f"\n\n{'='*70}")
    print("WALSH DECOMPOSITION (degree 4)")
    print("=" * 70)

    for name, vals in [("H", H_vals),
                       ("alpha_1", alpha1_vals),
                       ("alpha_2", alpha2_vals),
                       ("alpha_3", alpha3_vals)]:
        print(f"\n  {name} degree-4:")
        for indices in combinations(range(m), 4):
            total = 0
            for bits in range(half):
                prod = 1
                for idx in indices:
                    prod *= (1 if bits & (1 << idx) else -1)
                total += vals[bits] * prod
            h = total / half
            if abs(h) > 0.001:
                gaps = tuple(i+1 for i in indices)
                prod_gaps = 1
                for g in gaps:
                    prod_gaps *= g
                chi_prod = legendre(prod_gaps, p)
                sign_ok = (h > 0) == (chi_prod > 0) if h != 0 else True
                print(f"    h_hat[{set(indices)}] = {h:>12.4f}, "
                      f"gaps={gaps}, chi(prod)={chi_prod:+d}, "
                      f"sign {'OK' if sign_ok else 'FAIL'}")

    # Cycle-length Walsh at degree 2
    print(f"\n\n{'='*70}")
    print("CYCLE-LENGTH WALSH (degree 2)")
    print("=" * 70)

    for k in [3, 5, 7, 9, 11]:
        vals = [data[bits]['by_len'].get(k, 0) for bits in range(half)]
        if max(vals) == min(vals):
            print(f"\n  c_{k}: CONSTANT = {vals[0]}")
            continue

        print(f"\n  c_{k}:")
        print(f"    Mean: {sum(vals)/half:.4f}")

        for a in range(m):
            for b in range(a + 1, m):
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += vals[bits] * sa * sb
                h = total / half
                if abs(h) > 0.001:
                    gap_a, gap_b = a + 1, b + 1
                    chi_ab = legendre(gap_a * gap_b, p)
                    sign_ok = (h > 0) == (chi_ab > 0) if h != 0 else True
                    print(f"    h_hat[{{{a},{b}}}] = {h:>12.4f}, "
                          f"gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}, "
                          f"sign {'OK' if sign_ok else 'FAIL'}")

    # alpha_j sign pattern summary
    print(f"\n\n{'='*70}")
    print("SIGN PATTERN SUMMARY")
    print("=" * 70)

    for a in range(m):
        for b in range(a + 1, m):
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)

            signs = {}
            for name, vals in [("alpha_1", alpha1_vals),
                               ("alpha_2", alpha2_vals),
                               ("alpha_3", alpha3_vals),
                               ("H", H_vals)]:
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += vals[bits] * sa * sb
                h = total / half
                signs[name] = '+' if h > 0 else ('-' if h < 0 else '0')

            print(f"  ({a},{b}) chi={chi_ab:+d}: "
                  f"a1={signs['alpha_1']} a2={signs['alpha_2']} "
                  f"a3={signs['alpha_3']} H={signs['H']}")

    print(f"\n  Total time: {time.time()-t0:.1f}s")
    print("\nDONE.")


if __name__ == '__main__':
    main()
