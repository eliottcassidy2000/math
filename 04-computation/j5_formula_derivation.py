#!/usr/bin/env python3
"""
j5_formula_derivation.py -- Derive the exact formula for J_5(0,d) for Interval

THM-141: J_3(0,d) = min(d, p-d) for Interval S={1,...,m}.
Question: What is J_5(0,d)?

The spectral formula (walks = cycles for k=5 in tournaments) gives:
  J_5(0,d) = sum_{j=1}^{4} [A^j]_{0,d} * [A^{5-j}]_{d,0}

For Interval with Fejer kernel eigenvalues, this should yield an exact formula.

This script:
1. Derives J_5(0,d) combinatorially (counting 5-cycles through 0 and d)
2. Verifies the spectral formula
3. Searches for closed-form expressions
4. Extends the disjointness excess formula (THM-142) to 5-5 pairs
5. Computes the TOTAL alpha_2 excess (all cycle lengths)

Author: kind-pasteur-2026-03-12-S59b
"""

import time
import cmath
from itertools import combinations
from collections import defaultdict
import math


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_cycles(A, verts, start_idx=0):
    """Count directed Hamiltonian cycles on vertex set."""
    k = len(verts)
    dp = {}
    dp[(1 << start_idx, start_idx)] = 1

    for mask in range(1, 1 << k):
        if not (mask & (1 << start_idx)):
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
        if v == start_idx:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start_idx]]:
                total += dp[key]
    return total


def J5_exact(A, p, d):
    """Compute J_5(0,d) = #{directed 5-cycles through both 0 and d}."""
    count = 0
    others = [x for x in range(1, p) if x != d]
    for trio in combinations(others, 3):
        verts = [0, d] + list(trio)
        count += count_ham_cycles(A, verts, 0)
    return count


def J5_combinatorial(p, S, d):
    """Derive J_5(0,d) for Interval by counting vertex placements.

    A directed 5-cycle on {0, d, a, b, c} visits all 5 vertices exactly once
    in some cyclic order. There are (5-1)! = 24 directed cycles on 5 labeled
    vertices, but in a tournament exactly some subset exist.

    For Interval S={1,...,m}: arc i->j iff (j-i) mod p in {1,...,m}.
    """
    S_set = set(S)
    m = (p - 1) // 2

    # Enumerate all 5-cycles through {0, d} by choosing 3 more vertices
    count = 0
    others = [x for x in range(1, p) if x != d]

    for trio in combinations(others, 3):
        verts = [0, d] + list(trio)
        # Count directed Hamiltonian cycles on these 5 vertices
        # A directed Ham cycle is a cyclic permutation (v0, v1, ..., v4)
        # with arcs v0->v1, v1->v2, ..., v4->v0.
        for perm in _cyclic_perms_5(verts):
            valid = True
            for i in range(5):
                diff = (perm[(i+1) % 5] - perm[i]) % p
                if diff not in S_set:
                    valid = False
                    break
            if valid:
                count += 1

    return count


def _cyclic_perms_5(verts):
    """Generate all directed cyclic permutations of 5 elements, fixing first."""
    from itertools import permutations
    # Fix verts[0] = 0, permute the rest
    # Actually, for counting directed Ham cycles, fix start vertex
    start = verts[0]
    rest = verts[1:]
    for perm in permutations(rest):
        yield [start] + list(perm)


def J5_spectral(p, S, d):
    """Compute J_5(0,d) via spectral formula (walks = cycles for k=5).

    J_5(0,d) = sum_{j=1}^{4} [A^j]_{0,d} * [A^{5-j}]_{d,0}
    """
    omega = cmath.exp(2j * cmath.pi / p)
    lam = [sum(omega ** (r * s) for s in S) for r in range(p)]

    def Aj(j, src, dst):
        """[A^j]_{src,dst} for circulant."""
        return sum(lam[r] ** j * omega ** (-r * (dst - src)) for r in range(p)) / p

    total = sum(Aj(j, 0, d) * Aj(5-j, d, 0) for j in range(1, 5))
    return total.real


def J5_interval_closed_form(p, m, d):
    """Attempt to compute J_5(0,d) for Interval using the "feasibility counting" approach.

    A 5-cycle through 0 and d has form (0, a1, a2, a3, a4) -> (a1, a2, a3, a4, 0)
    where d appears among {a1, a2, a3, a4}.

    For Interval: arc i->j iff (j-i) mod p in {1,...,m}.
    An arc goes from i to j if j is "close ahead" of i (within m steps clockwise).

    Vertex d can be at position 1, 2, 3, or 4 in the cycle. By symmetry of the
    cycle (it's a cycle, not a path), we can fix one position and multiply.
    But the cycle has 5 directed edges, so there are 5 rotations, each pinning 0
    at a different position. Since we fix 0, there are 4 positions for d.
    """
    # Approach: enumerate all possible orderings of {0, d, a, b, c} on Z_p
    # that form a directed 5-cycle with Interval arcs.
    #
    # For a 5-cycle (v0, v1, v2, v3, v4) with arcs v0->v1->...->v4->v0,
    # each arc gap g_i = (v_{i+1} - v_i) mod p must be in {1,...,m}.
    # Sum of gaps: g0 + g1 + g2 + g3 + g4 = 0 mod p.
    # Each g_i in {1,...,m}.
    #
    # So we need 5 elements of {1,...,m} summing to 0 mod p = summing to p.
    # (Since each g_i >= 1, sum >= 5. Since each g_i <= m, sum <= 5m = 5(p-1)/2.)
    # For the sum to be 0 mod p: sum = p (for small p) or sum = 2p.
    # 5 <= sum <= 5m. For p=7, m=3: 5 <= sum <= 15. So sum = 7 or 14.
    # For p=11, m=5: 5 <= sum <= 25. So sum = 11 or 22.
    # For p=13, m=6: 5 <= sum <= 30. So sum = 13 or 26.

    # The constraint: vertex d is at position j for some j.
    # If d is at position j: the partial sums g0 + g1 + ... + g_{j-1} = d (mod p).
    # (Since v0 = 0, v_j = g0 + ... + g_{j-1} mod p.)

    # Let me count directly using gap compositions.
    S_set = set(range(1, m + 1))
    count = 0

    # Enumerate 5-tuples (g0,...,g4) in {1,...,m}^5 with sum = 0 mod p
    for g0 in range(1, m + 1):
        for g1 in range(1, m + 1):
            for g2 in range(1, m + 1):
                for g3 in range(1, m + 1):
                    g4 = (-g0 - g1 - g2 - g3) % p
                    if 1 <= g4 <= m:
                        # Check if vertex d appears in the cycle
                        # Partial sums: v1 = g0, v2 = g0+g1, v3 = g0+g1+g2,
                        #               v4 = g0+g1+g2+g3 (all mod p)
                        partial_sums = [g0 % p]
                        partial_sums.append((g0 + g1) % p)
                        partial_sums.append((g0 + g1 + g2) % p)
                        partial_sums.append((g0 + g1 + g2 + g3) % p)

                        # Check that all 5 vertices are distinct
                        verts = [0] + partial_sums
                        if len(set(verts)) == 5:
                            if d in partial_sums:
                                count += 1

    return count


def analyze_J5_structure(p, S, A, name):
    """Deep analysis of J_5(0,d) structure for Interval."""
    m = (p - 1) // 2

    print(f"\n  --- J_5 Deep Structure ({name}, p={p}) ---")

    # Compute exact J_5(0,d) for all d
    j5_profile = [0] * p
    for d in range(p):
        j5_profile[d] = J5_exact(A, p, d)

    print(f"    J_5(0,v) = {j5_profile}")

    # Spectral check
    if True:
        j5_spec = [J5_spectral(p, S, d) for d in range(p)]
        max_err = max(abs(j5_profile[d] - j5_spec[d]) for d in range(1, p))
        print(f"    Spectral match (v>0): max error = {max_err:.6f}")

    # Decompose J_5 by position of d in the cycle
    # d can be at position 1, 2, 3, or 4 in cycle (0, v1, v2, v3, v4)
    S_set = set(S)
    position_counts = {j: [0]*p for j in range(1, 5)}

    for g0 in range(1, m + 1):
        for g1 in range(1, m + 1):
            for g2 in range(1, m + 1):
                for g3 in range(1, m + 1):
                    g4 = (-g0 - g1 - g2 - g3) % p
                    if 1 <= g4 <= m:
                        v1 = g0 % p
                        v2 = (g0 + g1) % p
                        v3 = (g0 + g1 + g2) % p
                        v4 = (g0 + g1 + g2 + g3) % p
                        verts = [0, v1, v2, v3, v4]
                        if len(set(verts)) == 5:
                            for j in range(1, 5):
                                position_counts[j][verts[j]] += 1

    print(f"\n    Position decomposition of J_5(0,d):")
    for d in range(1, (p + 1) // 2 + 1):
        by_pos = [position_counts[j][d] for j in range(1, 5)]
        total_check = sum(by_pos)
        print(f"      d={d}: pos1={by_pos[0]}, pos2={by_pos[1]}, "
              f"pos3={by_pos[2]}, pos4={by_pos[3]}, total={total_check} "
              f"(exact={j5_profile[d]})")

    # Check convolution structure: J_5 = J_3 * something?
    j3_profile = [min(d, p - d) if d > 0 else 0 for d in range(p)]

    # J_5 might decompose as sum of convolutions of J_3
    print(f"\n    Convolution analysis:")
    print(f"    J_3(0,d): {j3_profile}")
    print(f"    J_5(0,d): {j5_profile}")

    # Test: J_5(0,d) = sum_{v} J_3(0,v) * J_3(v,d) - corrections?
    # J_3(0,v) * J_3(v,d) = co-occurrence via intermediate 3-cycle
    # This counts "3-cycle chains" through an intermediate vertex v
    conv_j3 = [0] * p
    for d in range(p):
        for v in range(1, p):
            if v == d:
                continue
            # J_3(v, d) for circulant = J_3(0, d-v mod p)
            gap = (d - v) % p
            conv_j3[d] += j3_profile[v] * j3_profile[gap]

    print(f"    Convolution J_3*J_3: {conv_j3}")
    # Ratio J_5/conv
    for d in range(1, (p + 1) // 2 + 1):
        if conv_j3[d] > 0:
            ratio = j5_profile[d] / conv_j3[d]
            diff = j5_profile[d] - conv_j3[d]
            print(f"      d={d}: J_5={j5_profile[d]}, conv={conv_j3[d]}, "
                  f"ratio={ratio:.4f}, diff={diff}")

    return j5_profile


def compute_total_alpha2(A, p, max_k=None):
    """Compute total alpha_2 (sum over all (k1,k2) disjoint pairs)."""
    if max_k is None:
        max_k = p

    cycles_by_k = defaultdict(list)
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            if count_ham_cycles(A, verts) > 0:
                cycles_by_k[k].append(frozenset(subset))

    # All cycle vertex sets
    all_vs = []
    for k in sorted(cycles_by_k.keys()):
        all_vs.extend(cycles_by_k[k])

    # Also count DIRECTED cycles (for the Omega graph)
    all_directed = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_dir = count_ham_cycles(A, verts)
            for _ in range(n_dir):
                all_directed.append(frozenset(subset))

    # Count disjoint pairs of DIRECTED cycles
    n = len(all_directed)
    alpha2 = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not (all_directed[i] & all_directed[j]):
                alpha2 += 1

    return alpha2, len(all_directed), cycles_by_k


def analyze_55_disjointness(A, p, S, name):
    """Detailed analysis of 5-5 disjointness using co-occurrence."""
    m = (p - 1) // 2
    S_set = set(S)

    # Get 5-cycle vertex sets
    c5_sets = []
    for subset in combinations(range(p), 5):
        verts = list(subset)
        if count_ham_cycles(A, verts) > 0:
            c5_sets.append(frozenset(subset))

    n5 = len(c5_sets)
    print(f"\n  --- 5-5 Disjointness Analysis ({name}, p={p}) ---")
    print(f"    c5 vertex sets: {n5}")

    if n5 == 0 or p < 10:
        print(f"    [Skipping: need k1+k2 <= p for disjoint pairs, need p >= 10]")
        return

    # Co-occurrence for 5-cycles
    co_occ_5 = [0] * p
    for fs in c5_sets:
        if 0 in fs:
            for v in fs:
                if v != 0:
                    co_occ_5[v] += 1

    n0 = sum(1 for fs in c5_sets if 0 in fs)  # sets containing vertex 0
    print(f"    5-cycle sets containing 0: {n0}")
    print(f"    co_occ_5 profile: {co_occ_5}")

    # Co-occurrence by gap
    for d in range(1, (p + 1) // 2 + 1):
        print(f"      gap {d}: {co_occ_5[d]} shared 5-cycle sets")

    # Disjoint pair count
    disjoint_55 = 0
    overlap_dist = defaultdict(int)
    total_pairs = 0
    for i in range(n5):
        for j in range(i + 1, n5):
            ov = len(c5_sets[i] & c5_sets[j])
            overlap_dist[ov] += 1
            total_pairs += 1
            if ov == 0:
                disjoint_55 += 1

    print(f"    Disjoint 5-5 pairs: {disjoint_55}/{total_pairs} "
          f"({100*disjoint_55/total_pairs:.2f}%)")
    print(f"    Overlap dist: {dict(sorted(overlap_dist.items()))}")

    # Apply the same I-E formula as THM-142
    nv = n0  # by circulant symmetry, same for all v
    sum_cnv2 = p * nv * (nv - 1) // 2
    ov2_pairs = sum(overlap_dist.get(ov, 0) for ov in range(2, 6))
    # Actually: sum_v C(n_v,2) = #(ov>=1 pairs counted by vertex) = ov1 + 2*ov2 + 3*ov3 + ...
    weighted_ov = sum(ov * cnt for ov, cnt in overlap_dist.items() if ov > 0)
    check = sum_cnv2  # should equal sum of ov * overlap_dist[ov]

    print(f"    n_v (sets per vertex) = {nv}")
    print(f"    p * C(n_v, 2) = {sum_cnv2}")
    print(f"    sum ov*count = {weighted_ov}")

    # Compute overlap>=2 count
    ov_ge2 = sum(cnt for ov, cnt in overlap_dist.items() if ov >= 2)
    print(f"    |ov>=2| = {ov_ge2}")
    print(f"    I-E check: C(n5,2) - sum_C(nv,2) + |ov>=2| = "
          f"{n5*(n5-1)//2} - {sum_cnv2} + {ov_ge2} = "
          f"{n5*(n5-1)//2 - sum_cnv2 + ov_ge2} (should = {disjoint_55})")

    return disjoint_55


def main():
    print("=" * 70)
    print("J_5 FORMULA DERIVATION + FULL ALPHA_2 ANALYSIS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        tournaments = []
        if p % 4 == 3:
            tournaments = [("Paley", S_qr), ("Interval", S_int)]
        else:
            tournaments = [("Interval", S_int)]

        for name, S in tournaments:
            A = build_adj(p, S)

            print(f"\n{'='*70}")
            print(f"p={p}, {name}, S={S}")
            print(f"{'='*70}")

            # J_5 deep structure
            j5_profile = analyze_J5_structure(p, S, A, name)

            # 5-5 disjointness
            analyze_55_disjointness(A, p, S, name)

    # ====== FULL ALPHA_2 COMPARISON ======
    print(f"\n{'='*70}")
    print("FULL ALPHA_2 COMPARISON (all cycle lengths)")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            A = build_adj(p, S)
            t0 = time.time()
            alpha2, n_dir, cycles_by_k = compute_total_alpha2(A, p)
            t1 = time.time()

            print(f"\n  p={p}, {name}: alpha_1={n_dir}, alpha_2={alpha2} "
                  f"({t1-t0:.1f}s)")
            print(f"    Cycle counts by length: "
                  f"{dict((k, len(v)) for k, v in sorted(cycles_by_k.items()))}")

    # ====== GAP COMPOSITION ANALYSIS ======
    print(f"\n{'='*70}")
    print("GAP COMPOSITION ANALYSIS")
    print("(Exactly which gap tuples contribute to J_5)")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))

        print(f"\n  p={p}, Interval S={S_int}:")

        # Enumerate all valid gap 5-tuples
        valid_gaps = []
        for g0 in range(1, m + 1):
            for g1 in range(1, m + 1):
                for g2 in range(1, m + 1):
                    for g3 in range(1, m + 1):
                        g4 = (-g0 - g1 - g2 - g3) % p
                        if 1 <= g4 <= m:
                            valid_gaps.append((g0, g1, g2, g3, g4))

        print(f"    Valid gap 5-tuples: {len(valid_gaps)}")
        print(f"    (These are DIRECTED cycles; some may have repeated vertices)")

        # Check which have 5 distinct vertices
        simple_count = 0
        non_simple_count = 0
        for g in valid_gaps:
            partial = [0]
            s = 0
            for gi in g[:-1]:
                s = (s + gi) % p
                partial.append(s)
            if len(set(partial)) == 5:
                simple_count += 1
            else:
                non_simple_count += 1

        print(f"    Simple (5 distinct vertices): {simple_count}")
        print(f"    Non-simple (repeated vertex): {non_simple_count}")

        # For simple cycles, tabulate which d values they hit
        d_contributions = defaultdict(int)
        for g in valid_gaps:
            partial = [0]
            s = 0
            for gi in g[:-1]:
                s = (s + gi) % p
                partial.append(s)
            if len(set(partial)) == 5:
                for d in partial[1:]:  # exclude 0
                    d_contributions[d] += 1

        print(f"    d contributions: {dict(sorted(d_contributions.items()))}")
        # Each d should appear J_5(0,d) times
        for d in range(1, (p+1)//2 + 1):
            print(f"      d={d}: count={d_contributions.get(d, 0)}")


if __name__ == '__main__':
    main()
