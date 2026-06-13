#!/usr/bin/env python3
"""
co_occurrence_all_cycles.py -- Joint cycle count J_k(0,v) for ALL cycle lengths

The joint cycle count J_k(0,v) = #{directed k-cycles through both vertex 0 and vertex v}
is the key quantity linking spectral concentration to cycle disjointness.

THM-141 proved: J_3(0,v) = min(v, p-v) for Interval.
This script extends to J_5, J_7, ..., and checks:
  1. Does the gradient persist at higher k?
  2. What is the exact formula for J_k(0,v)?
  3. How does the spectral formula relate to the combinatorial formula?

The spectral representation: for circulant A with eigenvalues lambda_r,
  J_k(0,v) = sum_{j=1}^{k-1} [A^j]_{0v} * [A^{k-j}]_{v0}
  [A^j]_{0v} = (1/p) * sum_r lambda_r^j * omega^{-rv}

Author: kind-pasteur-2026-03-12-S59b
"""

import time
import cmath
from itertools import combinations, permutations
from collections import defaultdict
import math


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def eigenvalues_circulant(p, S):
    """Compute eigenvalues of circulant adjacency matrix.
    lambda_r = sum_{s in S} omega^{rs}, omega = exp(2pi i/p)."""
    omega = cmath.exp(2j * cmath.pi / p)
    lam = []
    for r in range(p):
        val = sum(omega ** (r * s) for s in S)
        lam.append(val)
    return lam


def joint_cycle_count_exact(A, p, k, v):
    """Compute J_k(0,v) = #{directed k-cycles through both 0 and v} exactly.

    Enumerate all directed k-cycles through vertex 0, count those also through v.
    A k-cycle through 0: choose k-1 other vertices, find directed Ham cycles through 0.
    """
    if v == 0:
        # Self: count all k-cycles through 0
        # Each k-cycle on a set containing 0 passes through 0 exactly once
        count = 0
        for other in combinations(range(1, p), k - 1):
            verts = [0] + list(other)
            count += count_directed_ham_cycles_through(A, verts, 0)
        return count

    # Count k-cycles through both 0 and v
    count = 0
    for other in combinations([x for x in range(1, p) if x != v], k - 2):
        verts = [0, v] + list(other)
        count += count_directed_ham_cycles_through_pair(A, verts, 0, v)
    return count


def count_directed_ham_cycles_through(A, verts, start):
    """Count directed Hamiltonian cycles on vertex subset passing through 'start'."""
    k = len(verts)
    idx = {v: i for i, v in enumerate(verts)}
    si = idx[start]

    # DP: paths from start through mask ending at v
    dp = {}
    dp[(1 << si, si)] = 1

    for mask in range(1, 1 << k):
        if not (mask & (1 << si)):
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
        if v == si:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[si]]:
                total += dp[key]
    return total


def count_directed_ham_cycles_through_pair(A, verts, v0, v1):
    """Count directed Ham cycles passing through BOTH v0 and v1."""
    # Same as count_directed_ham_cycles_through, just check that the cycle visits both
    # Since both v0 and v1 are in verts, any Ham cycle on verts visits both.
    return count_directed_ham_cycles_through(A, verts, v0)


def joint_cycle_spectral(p, S, k, v):
    """Compute J_k(0,v) via spectral formula.

    J_k(0,v) ≈ sum_{j=1}^{k-1} [A^j]_{0v} * [A^{k-j}]_{v,0}

    But this counts k-WALKS, not k-cycles. For k=3, walks = cycles (all 3-vertex walks
    returning to start are simple). For k>=5, corrections needed.

    [A^j]_{0v} = (1/p) * sum_r lambda_r^j * omega^{-rv}
    """
    omega = cmath.exp(2j * cmath.pi / p)
    lam = eigenvalues_circulant(p, S)

    # Compute [A^j]_{0v} for each j
    def Aj_0v(j, target):
        return sum(lam[r] ** j * omega ** (-r * target) for r in range(p)) / p

    total = 0
    for j in range(1, k):
        total += Aj_0v(j, v) * Aj_0v(k - j, (p - v) % p)

    return total


def co_occurrence_profile(A, p, S, k, use_spectral=False):
    """Compute J_k(0,v) for v = 0, 1, ..., p-1."""
    profile = []
    for v in range(p):
        if use_spectral:
            val = joint_cycle_spectral(p, S, k, v).real
        else:
            val = joint_cycle_count_exact(A, p, k, v)
        profile.append(val)
    return profile


def co_occurrence_variance(profile, p):
    """Compute variance of co-occurrence profile (excluding v=0 self-count)."""
    vals = [profile[d] for d in range(1, p)]
    mean = sum(vals) / len(vals)
    var = sum((v - mean)**2 for v in vals) / len(vals)
    return mean, var


def disjoint_pairs_from_cycles(A, p, k1, k2):
    """Count pairs of vertex-disjoint k1-cycles and k2-cycles.

    For k1=k2=3, this should match THM-142.
    For mixed lengths (k1=3, k2=5), this is new data.
    """
    # Enumerate k1-cycles and k2-cycles as vertex sets
    cycles_k1 = set()
    cycles_k2 = set()

    for subset in combinations(range(p), k1):
        verts = list(subset)
        if count_directed_ham_cycles_through(A, verts, verts[0]) > 0:
            cycles_k1.add(frozenset(subset))

    if k1 == k2:
        cycles_k2 = cycles_k1
    else:
        for subset in combinations(range(p), k2):
            verts = list(subset)
            if count_directed_ham_cycles_through(A, verts, verts[0]) > 0:
                cycles_k2.add(frozenset(subset))

    # Count disjoint pairs
    list_k1 = sorted(cycles_k1)
    list_k2 = sorted(cycles_k2)

    disjoint = 0
    total = 0
    for i, c1 in enumerate(list_k1):
        start_j = i + 1 if k1 == k2 else 0
        for j in range(start_j, len(list_k2)):
            c2 = list_k2[j]
            if c1 == c2:
                continue
            total += 1
            if not (c1 & c2):
                disjoint += 1

    return len(list_k1), len(list_k2), disjoint, total


def analyze_overlap_weight_by_length(A, p, max_k=None):
    """For each pair of cycle lengths (k1, k2), compute overlap statistics."""
    if max_k is None:
        max_k = p

    # Enumerate all cycle vertex sets by length
    cycles_by_k = defaultdict(list)
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            if count_directed_ham_cycles_through(A, verts, verts[0]) > 0:
                cycles_by_k[k].append(frozenset(subset))

    results = {}
    for k1 in sorted(cycles_by_k.keys()):
        for k2 in sorted(cycles_by_k.keys()):
            if k2 < k1:
                continue
            c1_list = cycles_by_k[k1]
            c2_list = cycles_by_k[k2]

            overlap_dist = defaultdict(int)
            disjoint = 0
            total = 0

            for i, c1 in enumerate(c1_list):
                start_j = i + 1 if k1 == k2 else 0
                for j in range(start_j, len(c2_list)):
                    c2 = c2_list[j]
                    if k1 == k2 and i == j:
                        continue
                    ov = len(c1 & c2)
                    overlap_dist[ov] += 1
                    total += 1
                    if ov == 0:
                        disjoint += 1

            results[(k1, k2)] = {
                'n1': len(c1_list), 'n2': len(c2_list),
                'total_pairs': total, 'disjoint': disjoint,
                'overlap_dist': dict(overlap_dist)
            }

    return results


def fejer_prediction(p, m, k, v):
    """Predict J_k(0,v) using the Fejér kernel model.

    If eigenvalues are lambda_r = F_m(r) (Fejér kernel), then
    [A^j]_{0v} = (1/p) sum_r F_m(r)^j omega^{-rv}

    This should give a peaked profile for Interval.
    """
    omega = cmath.exp(2j * cmath.pi / p)

    # Fejér kernel eigenvalues
    lam = []
    for r in range(p):
        if r == 0:
            lam.append(complex(m))
        else:
            # sin(pi*m*r/p) / sin(pi*r/p)
            val = cmath.sin(cmath.pi * m * r / p) / cmath.sin(cmath.pi * r / p)
            lam.append(val)

    def Aj_0v(j, target):
        return sum(lam[r] ** j * omega ** (-r * target) for r in range(p)) / p

    total = 0
    for j in range(1, k):
        total += Aj_0v(j, v) * Aj_0v(k - j, (p - v) % p)

    return total.real


def main():
    print("=" * 70)
    print("JOINT CYCLE COUNT J_k(0,v) — ALL CYCLE LENGTHS")
    print("Connecting THM-141 (k=3) to general k")
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

            # Eigenvalues
            lam = eigenvalues_circulant(p, S)
            print(f"\n  Eigenvalues |lambda_r|^2:")
            for r in range(p):
                print(f"    r={r}: |lambda|^2 = {abs(lam[r])**2:.4f}")

            # J_k(0,v) for k=3,5,7
            for k in range(3, p + 1, 2):
                t0 = time.time()
                profile = co_occurrence_profile(A, p, S, k)
                t1 = time.time()

                # Also compute spectral walk count for comparison
                walk_profile = []
                for v in range(p):
                    walk_profile.append(joint_cycle_spectral(p, S, k, v).real)

                mean, var = co_occurrence_variance(profile, p)

                print(f"\n  --- J_{k}(0,v) (exact cycle count) ---")
                print(f"    Profile: {profile}")
                print(f"    Mean (v>0) = {mean:.4f}, Var = {var:.4f}")
                print(f"    Time: {t1-t0:.2f}s")

                # Check if gradient pattern: J_k(0,v) = f(min(v, p-v))
                symmetric = all(
                    profile[v] == profile[p - v]
                    for v in range(1, (p + 1) // 2)
                )
                print(f"    Symmetric J(v) = J(p-v): {symmetric}")

                # Check linearity for Interval
                if name == "Interval" and k == 3:
                    linear = all(
                        profile[v] == min(v, p - v)
                        for v in range(1, p)
                    )
                    print(f"    THM-141 check (J_3 = min(v,p-v)): {linear}")

                # Ratio to J_3 (if J_3 > 0)
                if k > 3:
                    ratios = []
                    for v in range(1, (p + 1) // 2 + 1):
                        j3 = min(v, p - v)
                        if j3 > 0:
                            ratios.append(profile[v] / j3)
                    print(f"    J_{k}/J_3 ratios: {[f'{r:.4f}' for r in ratios]}")

                # Walk-cycle discrepancy
                disc = [abs(walk_profile[v] - profile[v]) for v in range(p)]
                max_disc = max(disc)
                if max_disc < 0.01:
                    print(f"    Walk = Cycle (max discrepancy {max_disc:.6f})")
                else:
                    print(f"    Walk != Cycle! Max discrepancy = {max_disc:.2f}")
                    print(f"    Walk profile: {[f'{w:.1f}' for w in walk_profile]}")

            # Cross-length overlap analysis
            print(f"\n  --- Cross-length overlap analysis ---")
            t0 = time.time()
            results = analyze_overlap_weight_by_length(A, p, max_k=p)
            t1 = time.time()
            print(f"  (computed in {t1-t0:.1f}s)")

            for (k1, k2), data in sorted(results.items()):
                disj = data['disjoint']
                total = data['total_pairs']
                pct = 100 * disj / total if total > 0 else 0
                print(f"    ({k1},{k2}): {data['n1']}×{data['n2']} cycles, "
                      f"{disj}/{total} disjoint ({pct:.2f}%)")
                if data['overlap_dist']:
                    print(f"      overlap dist: {dict(sorted(data['overlap_dist'].items()))}")

    # ====== COMPARATIVE: ALPHA_2 DECOMPOSITION BY LENGTH ======
    print(f"\n{'='*70}")
    print("ALPHA_2 DECOMPOSITION: disjoint pairs by (k1,k2) type")
    print("=" * 70)

    for p in [7]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            A = build_adj(p, S)
            results = analyze_overlap_weight_by_length(A, p)

            print(f"\n  p={p}, {name}:")
            total_alpha2 = 0
            for (k1, k2), data in sorted(results.items()):
                disj = data['disjoint']
                total_alpha2 += disj
                if disj > 0:
                    print(f"    alpha_2({k1},{k2}) = {disj}")
            print(f"    TOTAL alpha_2 = {total_alpha2}")

    # ====== J_5 FORMULA INVESTIGATION ======
    print(f"\n{'='*70}")
    print("J_5 FORMULA INVESTIGATION")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)

        print(f"\n  p={p}, Interval:")
        profile_5 = co_occurrence_profile(A, p, S_int, 5)
        print(f"    J_5(0,v) = {profile_5}")

        # Try to find a formula for J_5
        # For k=3: J_3(0,v) = min(v, p-v)
        # For k=5: check if J_5(0,v) = polynomial in min(v, p-v)?
        for v in range(1, (p + 1) // 2 + 1):
            d = min(v, p - v)
            j5 = profile_5[v]
            j3 = d
            print(f"    d={d}: J_5={j5}, J_3={j3}, J_5/J_3={j5/j3:.4f}, "
                  f"J_5-C(d,2)={j5-d*(d-1)//2}")


if __name__ == '__main__':
    main()
