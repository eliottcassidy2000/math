#!/usr/bin/env python3
"""
conflict_ratio_formula.py -- kind-pasteur-2026-03-13-S60

Analyze the conflict ratio f(k,k')/c_k' pattern.

KEY OBSERVATION from H_ck_theory.py:
- f(k,k')/c_k' = 1 when k+k' > p (pigeonhole: must share a vertex)
- f(k,k')/c_k' < 1 only when k+k' <= p (disjoint placement possible)
- The ratio is SYMMETRIC in k,k' and depends on p,k,k'

When k+k' <= p, the ratio equals the probability that a random k-cycle
and random k'-cycle share a vertex. For UNIFORM random subsets:
  P(share vertex) = 1 - C(p-k, k')/C(p, k')

Let's test: is f(k,k')/c_k' = 1 - C(p-k, k')/C(p, k')?

This would mean the "combinatorial overlap fraction" -- the fraction of
pairs that conflict -- follows the hypergeometric model (as if cycles
were uniformly distributed over vertex sets).
"""

from math import comb
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
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


def compute_all_ratios(p, S, name):
    """Compute conflict ratios and compare to hypergeometric model."""
    A = build_adj(p, S)

    # Enumerate cycles by vertex set
    # A "cycle" here is (frozenset of vertices, count of directed cycles on those vertices)
    vertex_sets_by_k = {}
    directed_counts_by_k = {}

    for k in range(3, p + 1, 2):
        vs_list = []
        dc_list = []
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                for _ in range(nc):
                    vs_list.append(fs)
                    dc_list.append(1)
        vertex_sets_by_k[k] = vs_list
        directed_counts_by_k[k] = dc_list

    c_k = {k: len(v) for k, v in vertex_sets_by_k.items()}

    print(f"\n{'='*60}")
    print(f"  {name}, p={p}, S={S}")
    print(f"{'='*60}")
    for k in sorted(c_k):
        print(f"  c_{k} = {c_k[k]}")

    print(f"\n  Conflict ratio comparison: actual vs hypergeometric model")
    print(f"  {'k':>3} {'k2':>3} {'actual':>10} {'hypergeo':>10} {'delta':>10} {'k+k2':>5} {'<=p?':>5}")

    for k1 in sorted(vertex_sets_by_k):
        for k2 in sorted(vertex_sets_by_k):
            if k1 > k2:
                continue

            # Actual conflict ratio
            cyc1 = vertex_sets_by_k[k1]
            cyc2 = vertex_sets_by_k[k2]

            if not cyc1 or not cyc2:
                continue

            total_conf = 0
            if k1 == k2:
                for i in range(len(cyc1)):
                    for j in range(i + 1, len(cyc1)):
                        if cyc1[i] & cyc1[j]:
                            total_conf += 1
                n_pairs = len(cyc1) * (len(cyc1) - 1) // 2
            else:
                for c1 in cyc1:
                    for c2 in cyc2:
                        if c1 & c2:
                            total_conf += 1
                n_pairs = len(cyc1) * len(cyc2)

            if n_pairs == 0:
                continue

            actual_ratio = total_conf / n_pairs

            # Hypergeometric model: P(share vertex) = 1 - C(p-k1, k2)/C(p, k2)
            if k1 + k2 <= p:
                hyper = 1 - comb(p - k1, k2) / comb(p, k2)
            else:
                hyper = 1.0

            delta = actual_ratio - hyper
            leq = "yes" if k1 + k2 <= p else "no"

            print(f"  {k1:>3} {k2:>3} {actual_ratio:>10.6f} {hyper:>10.6f} {delta:>10.6f} {k1+k2:>5} {leq:>5}")

    # NEW TEST: Is the conflict ratio the SAME for ALL orientations?
    # If so, it's a property of the vertex structure (circulant), not the specific tournament
    return c_k, vertex_sets_by_k


def test_ratio_across_orientations(p):
    """Check if conflict ratios are constant across all orientations."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*60}")
    print(f"  RATIO STABILITY TEST at p={p}")
    print(f"{'='*60}")

    all_ratios = {}

    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)

        vertex_sets_by_k = {}
        for k in range(3, p + 1, 2):
            vs_list = []
            for subset in combinations(range(p), k):
                fs = frozenset(subset)
                nc = count_ham_cycles(A, list(subset))
                for _ in range(nc):
                    vs_list.append(fs)
            vertex_sets_by_k[k] = vs_list

        c_k = {k: len(v) for k, v in vertex_sets_by_k.items()}

        ratios = {}
        for k1 in sorted(vertex_sets_by_k):
            for k2 in sorted(vertex_sets_by_k):
                if k1 > k2:
                    continue
                cyc1 = vertex_sets_by_k[k1]
                cyc2 = vertex_sets_by_k[k2]
                if not cyc1 or not cyc2:
                    continue

                total_conf = 0
                if k1 == k2:
                    for i in range(len(cyc1)):
                        for j in range(i + 1, len(cyc1)):
                            if cyc1[i] & cyc1[j]:
                                total_conf += 1
                    n_pairs = len(cyc1) * (len(cyc1) - 1) // 2
                else:
                    for c1 in cyc1:
                        for c2 in cyc2:
                            if c1 & c2:
                                total_conf += 1
                    n_pairs = len(cyc1) * len(cyc2)

                if n_pairs > 0:
                    ratios[(k1, k2)] = total_conf / n_pairs

        all_ratios[bits] = ratios

    # Check which ratios vary across orientations
    ratio_keys = set()
    for r in all_ratios.values():
        ratio_keys.update(r.keys())

    print(f"\n  {'(k1,k2)':>10} {'min':>10} {'max':>10} {'range':>10} {'constant?':>10}")
    for key in sorted(ratio_keys):
        vals = [all_ratios[bits].get(key, None) for bits in range(N)]
        vals = [v for v in vals if v is not None]
        if vals:
            mn, mx = min(vals), max(vals)
            rng = mx - mn
            const = "YES" if rng < 1e-10 else "NO"
            print(f"  {str(key):>10} {mn:>10.6f} {mx:>10.6f} {rng:>10.6f} {const:>10}")

    return all_ratios


# Main analysis
for p_val in [7, 11]:
    m = (p_val - 1) // 2
    S_qr = sorted(j for j in range(1, p_val) if pow(j, (p_val - 1) // 2, p_val) == 1)
    S_int = list(range(1, m + 1))

    compute_all_ratios(p_val, S_qr, "Paley")
    compute_all_ratios(p_val, S_int, "Interval")

# Stability test
for p_val in [7, 11]:
    test_ratio_across_orientations(p_val)

print("\nDONE.")
