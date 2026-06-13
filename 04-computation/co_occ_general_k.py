#!/usr/bin/env python3
"""
co_occ_general_k.py -- General co-occurrence formula for all cycle lengths

THM-141: co_occ_3(d) = d for Interval (slope 1, intercept 0)
Finding: co_occ_5(d) = 65 + 3d at p=11 (slope 3, intercept 65)

Questions:
1. What is co_occ_k(d) for general k?
2. Is it always LINEAR in d for Interval?
3. What determines the slope and intercept?
4. Connection to Fejer kernel eigenvalues?
5. How does disjointness excess scale with k?

This script computes co_occ_k(d) for vertex sets (not directed cycles)
at multiple primes and cycle lengths.

Author: kind-pasteur-2026-03-12-S59b
"""

import time
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_cycles(A, verts):
    """Count directed Hamiltonian cycles on vertex set."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])

    dp = {}
    dp[(1 << 0, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
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
    for v in range(1, k):
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[0]]:
                total += dp[key]
    return total


def has_ham_cycle(A, verts):
    """Check if vertex set supports at least one directed Ham cycle."""
    return count_ham_cycles(A, verts) > 0


def compute_co_occ_k(A, p, k):
    """Compute co_occ_k(d) = #{k-cycle vertex sets containing both 0 and d}.

    Returns: (profile, n_sets, n_containing_0)
    """
    # Enumerate k-cycle vertex sets
    cycle_sets = []
    for subset in combinations(range(p), k):
        if has_ham_cycle(A, list(subset)):
            cycle_sets.append(frozenset(subset))

    # Count those containing vertex 0
    containing_0 = [fs for fs in cycle_sets if 0 in fs]

    # Co-occurrence profile
    profile = [0] * p
    for fs in containing_0:
        for v in fs:
            if v != 0:
                profile[v] += 1

    return profile, len(cycle_sets), len(containing_0)


def compute_directed_J_k(A, p, k):
    """Compute J_k(0,d) = #{directed k-cycles through both 0 and d}.

    Returns profile indexed by d.
    """
    profile = [0] * p
    for subset in combinations(range(p), k):
        if 0 not in subset:
            continue
        verts = list(subset)
        n_cyc = count_ham_cycles(A, verts)
        if n_cyc > 0:
            for v in verts:
                if v != 0:
                    profile[v] += n_cyc

    return profile


def compute_disjoint_pairs(cycle_sets):
    """Count vertex-disjoint pairs of cycle vertex sets."""
    n = len(cycle_sets)
    disjoint = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not (cycle_sets[i] & cycle_sets[j]):
                disjoint += 1
    return disjoint


def fit_linear(profile, p):
    """Fit co_occ(d) = a + b*d for d=1,...,m where m=(p-1)/2.

    Returns (intercept, slope, max_residual).
    """
    m = (p - 1) // 2
    vals = [profile[d] for d in range(1, m + 1)]
    n = len(vals)

    # Least squares: minimize sum (y_i - a - b*x_i)^2
    sum_x = sum(range(1, n + 1))
    sum_y = sum(vals)
    sum_xx = sum(d*d for d in range(1, n + 1))
    sum_xy = sum(d * vals[d-1] for d in range(1, n + 1))

    b = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x**2)
    a = (sum_y - b * sum_x) / n

    # Residuals
    residuals = [abs(vals[d-1] - (a + b*d)) for d in range(1, n + 1)]
    max_res = max(residuals)

    return a, b, max_res


def main():
    print("=" * 70)
    print("GENERAL co_occ_k(d) ANALYSIS")
    print("=" * 70)

    results_table = []

    for p in [7, 11, 13]:
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

            for k in range(3, min(p, 10), 2):
                t0 = time.time()
                profile, n_sets, n_0 = compute_co_occ_k(A, p, k)
                t1 = time.time()

                # Also get directed cycle J_k profile
                j_profile = compute_directed_J_k(A, p, k)

                # Fit linear
                a, b, max_res = fit_linear(profile, p)

                # Check symmetry
                symmetric = all(
                    profile[d] == profile[p - d]
                    for d in range(1, (p + 1) // 2)
                )

                # Variance
                vals = [profile[d] for d in range(1, p)]
                mean_v = sum(vals) / len(vals)
                var_v = sum((v - mean_v)**2 for v in vals) / len(vals)

                print(f"\n  k={k}: {n_sets} vertex sets, {n_0} through 0 ({t1-t0:.1f}s)")
                print(f"    co_occ_k profile: {profile[1:(p+1)//2+1]} "
                      f"(symmetric: {symmetric})")

                if name == "Interval":
                    print(f"    Linear fit: {a:.2f} + {b:.2f}*d "
                          f"(max residual: {max_res:.4f})")

                    # Check exact linearity
                    is_linear = max_res < 0.01
                    print(f"    EXACTLY LINEAR: {is_linear}")

                    if is_linear:
                        # Report exact slope and intercept
                        # Slope should be integer or simple fraction
                        print(f"    Exact intercept = {round(a)}, slope = {round(b)}")

                else:  # Paley
                    # Check if constant
                    is_const = (var_v < 0.01)
                    print(f"    CONSTANT: {is_const} (value = {mean_v:.1f})")

                print(f"    Variance: {var_v:.4f}")

                # Directed cycles: J_k profile
                j_vals = [j_profile[d] for d in range(1, (p+1)//2 + 1)]
                print(f"    J_k (directed) profile: {j_vals}")

                # Multiplicity: average directed cycles per vertex set
                if n_0 > 0:
                    total_directed_0 = sum(j_profile[d] for d in range(1, p))
                    avg_mult = total_directed_0 / (n_0 * (k - 1))
                    # Each vertex set containing 0 has k-1 other vertices,
                    # total_directed_0 = sum J_k(0,d) = n_0 * (mult * (k-1))
                    print(f"    Avg directed cycles per vertex set: "
                          f"{total_directed_0 / (n_0 * (k-1)):.4f}")

                # Disjoint pairs (if feasible)
                if 2 * k <= p and n_sets <= 5000:
                    t0 = time.time()
                    cycle_list = []
                    for subset in combinations(range(p), k):
                        if has_ham_cycle(A, list(subset)):
                            cycle_list.append(frozenset(subset))
                    disj = compute_disjoint_pairs(cycle_list)
                    t1 = time.time()
                    total_pairs = len(cycle_list) * (len(cycle_list) - 1) // 2
                    pct = 100 * disj / total_pairs if total_pairs > 0 else 0
                    print(f"    Disjoint k-k pairs: {disj}/{total_pairs} "
                          f"({pct:.2f}%) [{t1-t0:.1f}s]")

                results_table.append({
                    'p': p, 'name': name, 'k': k,
                    'n_sets': n_sets, 'slope': round(b, 4),
                    'intercept': round(a, 4), 'max_res': max_res,
                    'variance': var_v
                })

    # ====== SLOPE PATTERN ANALYSIS ======
    print(f"\n{'='*70}")
    print("SLOPE PATTERN: co_occ_k(d) = a_k + b_k * d for Interval")
    print("=" * 70)

    for p in [7, 11, 13]:
        print(f"\n  p={p}:")
        int_results = [r for r in results_table
                       if r['p'] == p and r['name'] == 'Interval']
        for r in int_results:
            print(f"    k={r['k']}: intercept={r['intercept']:.2f}, "
                  f"slope={r['slope']:.2f}, max_res={r['max_res']:.4f}")

        # Look for pattern in slopes
        slopes = [r['slope'] for r in int_results]
        ks = [r['k'] for r in int_results]
        if len(slopes) >= 2:
            print(f"    Slope sequence: {[f'{s:.2f}' for s in slopes]}")
            # Check if slopes follow a pattern
            # For k=3: slope = 1
            # For k=5: slope = 3
            # For k=7: slope = ?
            # Check: slopes = C(k-1, 2)? C(2,2)=1, C(4,2)=6, no.
            # Check: slopes = (k-1)(k-2)/2? k=3: 1, k=5: 6, no.
            # The slope for k=3 is 1, for k=5 is 3.
            # Maybe slope = C(k-2, (k-3)/2)?
            # Or slope relates to the number of gap compositions?

    # ====== CONNECTION TO FEJER KERNEL ======
    print(f"\n{'='*70}")
    print("CONNECTION TO FEJER KERNEL")
    print("=" * 70)

    import cmath

    for p in [11]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))

        # Fejer kernel eigenvalues
        omega = cmath.exp(2j * cmath.pi / p)
        lam_int = [sum(omega ** (r * s) for s in S_int) for r in range(p)]

        # Compute [A^j]_{0,d} from eigenvalues
        def Aj(j, d):
            return sum(lam_int[r] ** j * omega ** (-r * d) for r in range(p)) / p

        print(f"\n  p={p}, Interval eigenvalues |lambda_r|^2:")
        for r in range(p):
            print(f"    r={r}: {abs(lam_int[r])**2:.4f}")

        # Since walks = cycles for k <= 5, the spectral formula gives:
        # J_k(0,d) = sum_{j=1}^{k-1} A^j(0,d) * A^{k-j}(d,0)
        # = sum_{j=1}^{k-1} A^j(0,d) * A^{k-j}(0, p-d)  [by circulant symmetry]

        # Derivative of J_k w.r.t. d (discrete):
        print(f"\n  Spectral derivative analysis (p={p}):")
        for k in [3, 5]:
            for d in range(1, m + 2):
                val = sum(
                    Aj(j, d).real * Aj(k-j, (p-d) % p).real -
                    Aj(j, d).imag * Aj(k-j, (p-d) % p).imag
                    for j in range(1, k)
                )
                val_next = sum(
                    Aj(j, d+1).real * Aj(k-j, (p-d-1) % p).real -
                    Aj(j, d+1).imag * Aj(k-j, (p-d-1) % p).imag
                    for j in range(1, k)
                ) if d < m else 0
                diff = val_next - val if d < m else 0
                print(f"    k={k}, d={d}: J_k={val:.2f}, diff={diff:.2f}")

    # ====== DISJOINTNESS EXCESS FOR HIGHER k ======
    print(f"\n{'='*70}")
    print("DISJOINTNESS EXCESS BY CYCLE LENGTH")
    print("=" * 70)

    for p in [11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for k in [3, 5]:
            if 2 * k > p:
                continue

            disj_data = {}
            for name, S in [("Paley", S_qr), ("Interval", S_int)]:
                A = build_adj(p, S)
                cycle_list = []
                for subset in combinations(range(p), k):
                    if has_ham_cycle(A, list(subset)):
                        cycle_list.append(frozenset(subset))
                disj = compute_disjoint_pairs(cycle_list)
                disj_data[name] = (len(cycle_list), disj)

            n_p, d_p = disj_data['Paley']
            n_i, d_i = disj_data['Interval']
            excess = d_i - d_p

            print(f"\n  p={p}, k={k}:")
            print(f"    Paley: {n_p} vertex sets, {d_p} disjoint pairs")
            print(f"    Interval: {n_i} vertex sets, {d_i} disjoint pairs")
            print(f"    Excess: {excess}")

            # For k=3, verify THM-142
            if k == 3:
                thm142 = p * (p-1) * (p+1) * (p-3) // 192
                print(f"    THM-142 prediction: {thm142} (match: {excess == thm142})")


if __name__ == '__main__':
    main()
