#!/usr/bin/env python3
"""
wrap_structure_analysis.py -- Decompose co_occ_k by number of wraps

KEY INSIGHT: In the Interval tournament S={1,...,m}, a directed k-cycle has
all gaps in {1,...,m}, summing to w*p for some positive integer w (the "wrap number").

For 1-wrap cycles (sum=p), vertices are visited in simple cyclic order on Z_p.
For 2-wrap cycles (sum=2p), the cycle goes around Z_p twice, creating a non-trivial permutation.

Q1: Is co_occ_k(d) linear within each wrap class?
Q2: Does the binomial slope arise from a specific wrap class?
Q3: Connection to composition counting?

Also explores: gap composition structure and the "third vertex sweep" argument.

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations, permutations
from math import comb
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def get_directed_cycles(A, verts):
    """Return all directed Hamiltonian cycles on vertex set as gap tuples."""
    k = len(verts)
    cycles = []

    if k == 3:
        a, b, c = verts
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
        return cycles

    # DP to find all Ham paths starting at verts[0]
    dp = {}
    dp[(1 << 0, 0)] = [[0]]
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    if nkey not in dp:
                        dp[nkey] = []
                    for path in dp[key]:
                        dp[nkey].append(path + [w])

    full = (1 << k) - 1
    for v in range(1, k):
        key = (full, v)
        if key in dp and A[verts[v]][verts[0]]:
            for path in dp[key]:
                cycles.append(tuple(verts[i] for i in path))

    return cycles


def cycle_gaps(cycle, p):
    """Return gap tuple and wrap number for a directed cycle."""
    k = len(cycle)
    gaps = []
    for i in range(k):
        g = (cycle[(i+1) % k] - cycle[i]) % p
        gaps.append(g)
    wrap = sum(gaps) // p
    return tuple(gaps), wrap


def main():
    print("=" * 70)
    print("WRAP STRUCTURE ANALYSIS")
    print("=" * 70)

    for p in [11, 13]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}, Interval")
        print(f"{'='*70}")

        for k in [3, 5, 7]:
            if k > p:
                continue

            print(f"\n--- k={k} ---")

            # Collect all directed cycles through vertex 0
            # and classify by wrap number
            wrap_profile = {}  # wrap -> co_occ profile
            wrap_J_profile = {}  # wrap -> J_k profile

            # Track vertex sets by wrap
            vset_by_wrap = defaultdict(set)  # wrap -> set of frozensets

            t0 = time.time()
            for subset in combinations(range(p), k):
                if 0 not in subset:
                    continue
                verts = list(subset)
                dir_cycles = get_directed_cycles(A, verts)
                if not dir_cycles:
                    continue

                # Classify each directed cycle by wrap
                fs = frozenset(subset)
                wraps_seen = set()
                for cyc in dir_cycles:
                    _, w = cycle_gaps(cyc, p)
                    wraps_seen.add(w)
                    vset_by_wrap[w].add(fs)

                    # J_k contribution
                    if w not in wrap_J_profile:
                        wrap_J_profile[w] = [0] * p
                    for v in subset:
                        if v != 0:
                            wrap_J_profile[w][v] += 1

            t1 = time.time()

            # Compute co_occ profiles by wrap
            for w in sorted(vset_by_wrap.keys()):
                profile = [0] * p
                for fs in vset_by_wrap[w]:
                    for v in fs:
                        if v != 0:
                            profile[v] += 1
                wrap_profile[w] = profile

            # Total profile
            all_vsets = set()
            for vsets in vset_by_wrap.values():
                all_vsets |= vsets
            total_profile = [0] * p
            for fs in all_vsets:
                for v in fs:
                    if v != 0:
                        total_profile[v] += 1

            # Report
            print(f"  {len(all_vsets)} vertex sets through 0 ({t1-t0:.1f}s)")
            print(f"  Total co_occ: {[total_profile[d] for d in range(1, m+1)]}")

            for w in sorted(wrap_profile.keys()):
                prof = wrap_profile[w]
                vals = [prof[d] for d in range(1, m + 1)]
                n_vsets = len(vset_by_wrap[w])

                # Check linearity
                n = len(vals)
                if n >= 2 and max(vals) > 0:
                    sum_x = sum(range(1, n + 1))
                    sum_y = sum(vals)
                    sum_xx = sum(d*d for d in range(1, n + 1))
                    sum_xy = sum(d * vals[d-1] for d in range(1, n + 1))
                    denom = n * sum_xx - sum_x**2
                    if denom != 0:
                        b = (n * sum_xy - sum_x * sum_y) / denom
                        a = (sum_y - b * sum_x) / n
                        residuals = [abs(vals[d-1] - (a + b*d)) for d in range(1, n + 1)]
                        max_res = max(residuals)
                        is_lin = max_res < 0.01
                    else:
                        a, b, max_res, is_lin = 0, 0, 0, True
                else:
                    a, b, max_res, is_lin = 0, 0, 0, True

                print(f"\n  Wrap {w}: {n_vsets} vertex sets")
                print(f"    co_occ: {vals}")
                if is_lin:
                    print(f"    LINEAR: a={a:.1f}, b={b:.1f}")
                else:
                    print(f"    NOT LINEAR (max_res={max_res:.4f})")
                    print(f"    fit: a={a:.2f}, b={b:.2f}")

                # J_k (directed cycles) for this wrap
                if w in wrap_J_profile:
                    j_vals = [wrap_J_profile[w][d] for d in range(1, m + 1)]
                    print(f"    J_k:    {j_vals}")

            # Wrap overlap: vertex sets appearing in multiple wraps
            all_wraps = sorted(vset_by_wrap.keys())
            if len(all_wraps) > 1:
                for i in range(len(all_wraps)):
                    for j in range(i+1, len(all_wraps)):
                        w1, w2 = all_wraps[i], all_wraps[j]
                        overlap = vset_by_wrap[w1] & vset_by_wrap[w2]
                        if overlap:
                            print(f"\n  Vertex sets in BOTH wrap {w1} and {w2}: {len(overlap)}")
                            # These have directed cycles with different wrap numbers
                            for fs in list(overlap)[:3]:
                                verts = sorted(fs)
                                gaps_natural = [(verts[(i+1) % k] - verts[i]) % p for i in range(k)]
                                print(f"    {verts} natural_gaps={gaps_natural}")

    # ====== COMPOSITION COUNTING ANALYSIS ======
    print(f"\n{'='*70}")
    print("COMPOSITION COUNTING vs CO_OCC")
    print("=" * 70)

    for p in [11, 13]:
        m = (p - 1) // 2
        print(f"\np={p}, m={m}")

        for k in [3, 5]:
            # Count compositions of w*p into k parts in {1,...,m} with d as partial sum
            # For each wrap w, N_k(w*p, d) = #{compositions with d among partial sums}
            comp_profile = defaultdict(lambda: [0] * p)

            for target_sum in range(k, k*m + 1):
                if target_sum % p != 0:
                    continue
                w = target_sum // p

                # Generate compositions of target_sum into k parts in {1,...,m}
                def count_comps_with_d(target, k, m, p):
                    """Count compositions by which d values appear as partial sums."""
                    profile = [0] * p
                    n_comps = [0]

                    def gen(remaining_parts, remaining_sum, partial_sums):
                        if remaining_parts == 0:
                            if remaining_sum == 0:
                                n_comps[0] += 1
                                for s in partial_sums:
                                    d = s % p
                                    if d != 0:
                                        profile[d] += 1
                            return
                        for g in range(1, min(m, remaining_sum - remaining_parts + 1) + 1):
                            new_sum = (partial_sums[-1] + g) if partial_sums else g
                            gen(remaining_parts - 1, remaining_sum - g,
                                partial_sums + [new_sum])

                    gen(k, target, [])
                    return profile, n_comps[0]

                prof, nc = count_comps_with_d(target_sum, k, m, p)
                for d in range(p):
                    comp_profile[w][d] += prof[d]

                print(f"\n  k={k}, wrap={w}: {nc} compositions of {target_sum}")
                vals = [comp_profile[w][d] for d in range(1, m + 1)]
                print(f"    J_k (by position) profile: {vals}")

            # Total J_k
            total_J = [0] * p
            for w in comp_profile:
                for d in range(p):
                    total_J[d] += comp_profile[w][d]
            print(f"\n  k={k}, total J_k: {[total_J[d] for d in range(1, m+1)]}")

    # ====== BINOMIAL COEFFICIENT ORIGIN ======
    print(f"\n{'='*70}")
    print("BINOMIAL ORIGIN: Why C(m-j, j)?")
    print("=" * 70)

    print("\nHypothesis: b_k = C(m-(k-1)/2, (k-1)/2) counts the number of")
    print("'boundary compositions' that appear/disappear as d changes by 1.")
    print()

    # For k=5 at p=11: b_5 = C(3,2) = 3
    # When d increases by 1, 3 new vertex sets become feasible
    # and 0 old ones become infeasible (net +3).
    # Let's verify by looking at which vertex sets appear/disappear.

    for p in [11, 13]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)
        k = 5

        print(f"\np={p}, m={m}, k={k}")
        print(f"  Predicted slope = C({m-2},2) = {comb(m-2, 2)}")

        # For each d, find vertex sets containing {0, d}
        vsets_by_d = {}
        for d in range(1, m + 1):
            vsets = set()
            for subset in combinations(range(p), k):
                if 0 in subset and d in subset:
                    verts = list(subset)
                    # Check Ham cycle
                    if has_ham_cycle_fast(A, verts, p, m):
                        vsets.add(frozenset(subset))
            vsets_by_d[d] = vsets

        for d in range(1, m + 1):
            n = len(vsets_by_d[d])
            if d > 1:
                gained = vsets_by_d[d] - vsets_by_d[d-1]
                lost = vsets_by_d[d-1] - vsets_by_d[d]
                # gained and lost don't make sense directly since d changes
                # But we can compute the SYMMETRIC DIFFERENCE
                # Actually, sets for d and d-1 have different "d" vertex,
                # so they can't be directly compared.
                # Instead, just report the count difference.
                delta = n - len(vsets_by_d[d-1])
                print(f"  d={d}: {n} vertex sets (delta = {delta})")
            else:
                print(f"  d={d}: {n} vertex sets")


def has_ham_cycle_fast(A, verts, p, m):
    """Fast Ham cycle check for Interval tournament using gap constraint."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a] + A[a][c] * A[c][b] * A[b][a]) > 0

    # DP
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


if __name__ == '__main__':
    main()
