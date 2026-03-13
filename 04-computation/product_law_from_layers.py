#!/usr/bin/env python3
"""
product_law_from_layers.py -- Derive the degree-4 product law from Gauss layering

The PRODUCT LAW observed at degree-4 Walsh:
  h_hat_4[{a,b,c,d}] proportional to eps_3(a,b)*eps_3(c,d) for some pairing

can be explained by the Gauss sum layering theorem:
  degree-2 Walsh has q=3 content (from j=4 layer) and q=5 (from j=6 layer).
  degree-4 Walsh of H(T) involves products of degree-2 Walsh coefficients.

Since H = sum 2^j * alpha_j (OCF), the degree-4 Walsh of H involves:
  - 2^1 * h_4[c_k] for simple cycles (from alpha_1)
  - 2^2 * h_4[alpha_2] for pairs (from alpha_2 = pairs of disjoint cycles)

For alpha_2 (pairs of disjoint 3-cycles for the dominant contribution):
  h_4[alpha_2, {a,b,c,d}] involves the CONVOLUTION of degree-2 Walsh:
  h_2[c_3, {a,b}] * h_2[c_3, {c,d}]

Since c_3 has only q=3 content, this product is eps_3(a,b)*eps_3(c,d).
But c_k for k>=6 has q=5 content, so higher alpha_j contributions mix
q=3 and q=5, generating the "mixed product" terms.

This script:
1. Computes degree-4 Walsh of H for circulant tournaments
2. Tests whether it factorizes as product of degree-2 Walsh coefficients
3. Identifies which alpha_j level controls the dominant degree-4 term
4. Connects the j=4/j=6 layer structure to the product law

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def resonance_sign(ga, gb, p, q):
    if (q * ga - gb) % p == 0:
        return -1
    elif (q * ga + gb) % p == 0:
        return +1
    elif (ga - q * gb) % p == 0:
        return -1
    elif (ga + q * gb) % p == 0:
        return +1
    return 0


def held_karp(A, verts):
    k = len(verts)
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


def compute_H_and_alphas(A, p):
    """Compute H(T) = I(Omega(T), 2) and overlap weights alpha_j."""
    # Count simple cycles by vertex set
    cycle_sets = defaultdict(int)  # frozenset -> directed cycle count

    for k in range(3, p + 1, 2):  # odd cycles only (tournaments)
        if k > p:
            break
        for verts in combinations(range(p), k):
            cnt = held_karp(A, list(verts))
            if cnt > 0:
                cycle_sets[frozenset(verts)] = cnt

    # Independence polynomial of conflict graph
    # Vertices of Omega = cycle vertex-sets, edges = shared vertices
    omega_verts = list(cycle_sets.keys())
    n = len(omega_verts)

    # For small n, enumerate independent sets
    # I(Omega, 2) = sum over independent sets S: 2^|S|
    # Using inclusion-exclusion through alpha_j
    alpha = defaultdict(int)
    alpha[0] = 1

    # Build adjacency for conflict graph
    conflicts = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if omega_verts[i] & omega_verts[j]:  # shared vertex
                conflicts[i][j] = conflicts[j][i] = True

    # Enumerate independent sets up to size limit
    max_j = min(n, 6)  # practical limit

    def count_independent_sets(remaining, selected, last_idx):
        j = len(selected)
        alpha[j] += 1

        for i in range(last_idx + 1, n):
            if i not in remaining:
                continue
            # Check if i conflicts with any selected
            ok = True
            for s in selected:
                if conflicts[i][s]:
                    ok = False
                    break
            if ok and j + 1 <= max_j:
                count_independent_sets(remaining, selected | {i}, i)

    # This is too slow for large n. Use simpler formula for H.
    # H = I(Omega, 2) = sum_{S independent} 2^|S|.
    # For tournaments: H(T) = number of Hamiltonian paths (Redei/OCF).

    # Simpler: use OCF directly. H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # alpha_1 = number of directed odd cycles = sum of cycle_sets values / 2
    # (each cycle counted in both directions? No, held_karp counts directed cycles)
    alpha_1 = sum(cycle_sets.values())

    # alpha_2: pairs of vertex-disjoint cycle sets
    alpha_2 = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not conflicts[i][j]:
                alpha_2 += cycle_sets[omega_verts[i]] * cycle_sets[omega_verts[j]]

    H = 1 + 2 * alpha_1 + 4 * alpha_2
    # Higher alphas omitted for now (they exist but are smaller)
    return H, alpha_1, alpha_2, cycle_sets


def main():
    print("=" * 70)
    print("PRODUCT LAW FROM GAUSS LAYERING")
    print("=" * 70)

    p = 7
    m = (p - 1) // 2
    n_orient = 1 << m
    pairs = [(j, p - j) for j in range(1, m + 1)]
    omega = cmath.exp(2j * cmath.pi / p)

    print(f"\np={p}, m={m}")

    data = []
    for bits in range(n_orient):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))
        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1

        # Compute c_k for k=3,5,7
        c3 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 3))
        c5 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 5))
        c7 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 7))

        # alpha_1 = c3 + c5 + c7
        alpha1 = c3 + c5 + c7

        # alpha_2: count pairs of disjoint cycles
        cycles_by_set = defaultdict(int)
        for k in [3, 5, 7]:
            for verts in combinations(range(p), k):
                cnt = held_karp(A, list(verts))
                if cnt > 0:
                    cycles_by_set[frozenset(verts)] = cnt

        alpha2 = 0
        cycle_list = list(cycles_by_set.items())
        for i in range(len(cycle_list)):
            for j in range(i + 1, len(cycle_list)):
                if not (cycle_list[i][0] & cycle_list[j][0]):
                    alpha2 += cycle_list[i][1] * cycle_list[j][1]

        H = 1 + 2 * alpha1 + 4 * alpha2
        # At p=7, no need for alpha_3 since max 7/3 = 2 disjoint 3-cycles

        data.append({
            'bits': bits, 'c3': c3, 'c5': c5, 'c7': c7,
            'alpha1': alpha1, 'alpha2': alpha2, 'H': H
        })

    sigma = {}
    for bits in range(n_orient):
        sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))

    # Degree-2 Walsh of various quantities
    print(f"\n  Degree-2 Walsh coefficients:")
    chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

    print(f"  {'pair':>8} {'q':>3} {'h_c3':>10} {'h_c5':>10} {'h_c7':>10} "
          f"{'h_a1':>10} {'h_a2':>10} {'h_H':>10}")

    deg2 = {}
    for a, b in chord_pairs:
        ga, gb = a + 1, b + 1
        q = resonance_level(ga, gb, p)
        eps = resonance_sign(ga, gb, p, q) if q <= 5 else 0

        for name in ['c3', 'c5', 'c7', 'alpha1', 'alpha2', 'H']:
            h = sum(data[bits][name] * sigma[bits][a] * sigma[bits][b]
                    for bits in range(n_orient)) / n_orient
            deg2[(name, a, b)] = h

        print(f"  ({a},{b}) q={q:>3}: {deg2[('c3',a,b)]:>10.4f} "
              f"{deg2[('c5',a,b)]:>10.4f} {deg2[('c7',a,b)]:>10.4f} "
              f"{deg2[('alpha1',a,b)]:>10.4f} {deg2[('alpha2',a,b)]:>10.4f} "
              f"{deg2[('H',a,b)]:>10.4f}")

    # Degree-4 Walsh of H
    print(f"\n  Degree-4 Walsh of H:")
    quads = [(a, b, c, d) for a in range(m) for b in range(a+1, m)
             for c in range(b+1, m) for d in range(c+1, m)]

    if not quads:
        # For m=3, there's exactly one degree-4 Walsh (full set)
        # Actually for m=3, degree-4 Walsh needs 4 indices from 0,1,2 = impossible
        # Degree-3 Walsh: one triple (0,1,2)
        print("  (m=3: no degree-4 Walsh, only degree-3)")
        a, b, c = 0, 1, 2
        h3_H = sum(data[bits]['H'] * sigma[bits][a] * sigma[bits][b] * sigma[bits][c]
                   for bits in range(n_orient)) / n_orient
        h3_a1 = sum(data[bits]['alpha1'] * sigma[bits][a] * sigma[bits][b] * sigma[bits][c]
                    for bits in range(n_orient)) / n_orient
        h3_a2 = sum(data[bits]['alpha2'] * sigma[bits][a] * sigma[bits][b] * sigma[bits][c]
                    for bits in range(n_orient)) / n_orient
        h3_c3 = sum(data[bits]['c3'] * sigma[bits][a] * sigma[bits][b] * sigma[bits][c]
                    for bits in range(n_orient)) / n_orient
        h3_c5 = sum(data[bits]['c5'] * sigma[bits][a] * sigma[bits][b] * sigma[bits][c]
                    for bits in range(n_orient)) / n_orient
        h3_c7 = sum(data[bits]['c7'] * sigma[bits][a] * sigma[bits][b] * sigma[bits][c]
                    for bits in range(n_orient)) / n_orient

        print(f"  (0,1,2): h3_H={h3_H:.4f}, h3_a1={h3_a1:.4f}, h3_a2={h3_a2:.4f}")
        print(f"           h3_c3={h3_c3:.4f}, h3_c5={h3_c5:.4f}, h3_c7={h3_c7:.4f}")

    # Now do p=11 for degree-4
    print(f"\n{'='*60}")
    p = 11
    m = (p - 1) // 2
    n_orient = 1 << m
    pairs = [(j, p - j) for j in range(1, m + 1)]
    print(f"p={p}, m={m}")

    data = []
    for bits in range(n_orient):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))
        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1

        c3 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 3))
        c5 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 5))

        # alpha_1 for c3 and c5 only (c7+ too expensive)
        alpha1_35 = c3 + c5

        # For alpha_2: pairs of disjoint 3-cycles only (dominant contribution)
        cycle3_sets = {}
        for verts in combinations(range(p), 3):
            cnt = held_karp(A, list(verts))
            if cnt > 0:
                cycle3_sets[frozenset(verts)] = cnt

        alpha2_33 = 0
        c3_list = list(cycle3_sets.items())
        for i in range(len(c3_list)):
            for j in range(i + 1, len(c3_list)):
                if not (c3_list[i][0] & c3_list[j][0]):
                    alpha2_33 += c3_list[i][1] * c3_list[j][1]

        data.append({
            'bits': bits, 'c3': c3, 'c5': c5,
            'alpha1_35': alpha1_35, 'alpha2_33': alpha2_33
        })

    sigma = {}
    for bits in range(n_orient):
        sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))

    chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

    # Degree-4 Walsh
    quads = [(a, b, c, d) for a in range(m) for b in range(a+1, m)
             for c in range(b+1, m) for d in range(c+1, m)]

    print(f"\n  Degree-4 Walsh of alpha_2 (3-3 pairs) and c_5:")
    print(f"  {'quad':>15} {'q_ab':>5} {'q_cd':>5} {'h4_a2':>12} {'h4_c5':>12} "
          f"{'h2_c3_ab*cd':>15} {'ratio_a2':>12}")

    for a, b, c, d in quads:
        ga, gb, gc, gd = a + 1, b + 1, c + 1, d + 1
        q_ab = resonance_level(ga, gb, p)
        q_cd = resonance_level(gc, gd, p)
        q_ac = resonance_level(ga, gc, p)
        q_bd = resonance_level(gb, gd, p)

        h4_a2 = sum(data[bits]['alpha2_33']
                     * sigma[bits][a] * sigma[bits][b]
                     * sigma[bits][c] * sigma[bits][d]
                     for bits in range(n_orient)) / n_orient

        h4_c5 = sum(data[bits]['c5']
                     * sigma[bits][a] * sigma[bits][b]
                     * sigma[bits][c] * sigma[bits][d]
                     for bits in range(n_orient)) / n_orient

        # Product of degree-2 Walsh of c_3
        h2_c3_ab = sum(data[bits]['c3'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
        h2_c3_cd = sum(data[bits]['c3'] * sigma[bits][c] * sigma[bits][d]
                       for bits in range(n_orient)) / n_orient

        prod_ab_cd = h2_c3_ab * h2_c3_cd

        # Also try other pairings
        h2_c3_ac = sum(data[bits]['c3'] * sigma[bits][a] * sigma[bits][c]
                       for bits in range(n_orient)) / n_orient
        h2_c3_bd = sum(data[bits]['c3'] * sigma[bits][b] * sigma[bits][d]
                       for bits in range(n_orient)) / n_orient
        prod_ac_bd = h2_c3_ac * h2_c3_bd

        h2_c3_ad = sum(data[bits]['c3'] * sigma[bits][a] * sigma[bits][d]
                       for bits in range(n_orient)) / n_orient
        h2_c3_bc = sum(data[bits]['c3'] * sigma[bits][b] * sigma[bits][c]
                       for bits in range(n_orient)) / n_orient
        prod_ad_bc = h2_c3_ad * h2_c3_bc

        if abs(h4_a2) > 0.01 or abs(h4_c5) > 0.01:
            ratio = h4_a2 / prod_ab_cd if abs(prod_ab_cd) > 0.01 else 0
            print(f"  ({a},{b},{c},{d}) q=({q_ab},{q_cd}): "
                  f"{h4_a2:>12.4f} {h4_c5:>12.4f} "
                  f"{prod_ab_cd:>15.4f} {ratio:>12.4f}")

            # Check if other pairings work better
            if abs(prod_ac_bd) > 0.01:
                r2 = h4_a2 / prod_ac_bd
                if abs(r2 - round(r2)) < 0.01:
                    print(f"               ac*bd={prod_ac_bd:.4f}, ratio={r2:.4f}")
            if abs(prod_ad_bc) > 0.01:
                r3 = h4_a2 / prod_ad_bc
                if abs(r3 - round(r3)) < 0.01:
                    print(f"               ad*bc={prod_ad_bc:.4f}, ratio={r3:.4f}")

    # Check: does alpha_2 Walsh factorize as product of c_3 Walsh?
    # For alpha_2 = sum of pairs of disjoint 3-cycles:
    # If cycles at vertex sets U, V are independent (U cap V = empty),
    # then alpha_2 = (1/2) * [alpha_1^2 - sum conflicting pairs]
    # But this isn't quite the convolution...

    # Actually, alpha_2 = sum_{U,V disjoint} c(U)*c(V) where c(U) = directed
    # cycle count on U. The Walsh of alpha_2 is related to the CONVOLUTION
    # of the Walsh of c_3 with itself, modulo the conflict structure.

    # The key question: is h4[alpha_2] approximately proportional to
    # product of h2[c_3] terms?

    print(f"\n  ANALYSIS: Alpha_2 degree-4 Walsh vs c_3 degree-2 products")
    print(f"  Since c_3 is CONSTANT (={data[0]['c3']}), its degree-2 Walsh is 0!")
    print(f"  So the product h2[c3]*h2[c3] = 0 for ALL pairings.")
    print(f"  The nonzero alpha_2 Walsh must come from the STRUCTURE of")
    print(f"  disjoint cycles, not just the total count.")

    # The key insight: alpha_2 depends on WHICH 3-cycles exist (their orientation),
    # not just how many. Even though c_3 is constant, the PAIR structure
    # (which pairs are disjoint) varies by orientation.

    # For a circulant tournament with connection set S:
    # A 3-cycle on vertices {u,v,w} exists iff v-u, w-v, u-w are all in S.
    # The number of directed 3-cycles through a specific vertex set depends on
    # the specific elements of S, hence on the orientation.

    # Wait, c_3 IS constant for circulant tournaments at prime p.
    # But alpha_2 varies because the LOCATION of cycles changes.

    print(f"\n  Alpha_2 values: {sorted(set(d['alpha2_33'] for d in data))}")
    print(f"  c_5 values: {sorted(set(d['c5'] for d in data))}")

    # Key question: is alpha_2 a function of c_5?
    c5_to_a2 = defaultdict(set)
    for d in data:
        c5_to_a2[d['c5']].add(d['alpha2_33'])

    print(f"\n  c_5 -> alpha_2 mapping:")
    for c5_val in sorted(c5_to_a2):
        a2_vals = sorted(c5_to_a2[c5_val])
        print(f"    c_5={c5_val}: alpha_2 = {a2_vals}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
