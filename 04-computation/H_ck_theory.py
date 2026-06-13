#!/usr/bin/env python3
"""
H_ck_theory.py -- kind-pasteur-2026-03-13-S60

Theoretical analysis of H = linear(c_k) for circulant tournaments.

KEY INSIGHT: For circulant tournaments on Z_p, all vertices are equivalent.
So the number of k-cycles THROUGH vertex 0 is c_k(0) = k*c_k/p (since each
k-cycle contains k vertices, and by symmetry each is the root equally often).

Wait, that's not quite right. c_k is the number of DIRECTED k-cycles, and
each k-cycle visits exactly k vertices. So the average number of k-cycles
through any vertex is k*c_k/p. By circulant symmetry, this is exact.

Now, two cycles overlap iff they share a vertex. For the conflict graph:
- degree of a k-cycle in the conflict graph = number of other cycles
  sharing at least one vertex with it

For a k-cycle C through vertices v1,...,vk:
- Number of k'-cycles sharing at least vertex vi with C
  depends on the local structure at vi.

For circulant symmetry: the "local structure" at each vertex is identical.

The key quantity: for a cycle C of length k, what is its expected
number of conflicts with k'-cycles?

Let f(k,k') = number of k'-cycles conflicting with a FIXED k-cycle.
Then the total overlap between k-cycles and k'-cycles is:
  overlap(k,k') = c_k * f(k,k') = c_k' * f(k',k)

And disj(k,k') = c_k * c_k' - overlap(k,k') (for k != k')
  disj(k,k) = c_k*(c_k-1)/2 - overlap(k,k)/2

If f(k,k') is a LINEAR function of c_k', then disj(k,k') is linear in c_k'.
And alpha_2 = sum disj(k,k') would be linear in c_k.

Can we derive f(k,k') for circulant tournaments?

This script computes f(k,k') empirically at p=7 and p=11.
"""

from itertools import combinations
from collections import defaultdict
import time


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


def analyze_conflict_structure(p, S, name):
    """Analyze the conflict graph structure."""
    A = build_adj(p, S)

    # Enumerate cycles by length
    cycles_by_k = {}
    all_cycles = []
    for k in range(3, p + 1, 2):
        cyc_k = []
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            for _ in range(n_cyc):
                cyc_k.append(frozenset(subset))
                all_cycles.append((frozenset(subset), k))
        if cyc_k:
            cycles_by_k[k] = cyc_k

    c_k = {k: len(v) for k, v in cycles_by_k.items()}
    n_total = sum(c_k.values())

    print(f"\n  {name}, S={S}")
    print(f"  Total cycles: {n_total}")
    for k in sorted(c_k):
        # k-cycles through vertex 0
        c_k_0 = sum(1 for fs in cycles_by_k[k] if 0 in fs)
        expected = k * c_k[k] / p
        print(f"    c_{k} = {c_k[k]}, through vertex 0: {c_k_0} "
              f"(expected k*c_k/p = {expected:.1f})")

    # For each pair of cycle lengths, compute conflict counts
    print(f"\n  Conflict structure f(k,k'):")
    for k1 in sorted(cycles_by_k):
        for k2 in sorted(cycles_by_k):
            if k1 > k2:
                continue
            # Count total conflicts between k1-cycles and k2-cycles
            total_conf = 0
            if k1 == k2:
                cyc = cycles_by_k[k1]
                for i in range(len(cyc)):
                    for j in range(i + 1, len(cyc)):
                        if cyc[i] & cyc[j]:
                            total_conf += 1
                # Total conflicting pairs (not total conflicts per cycle)
                n_pairs = len(cyc) * (len(cyc) - 1) // 2
                n_disj = n_pairs - total_conf
                # Average conflicts per k1-cycle with other k1-cycles
                avg_conf = 2 * total_conf / len(cyc) if len(cyc) > 0 else 0
            else:
                cyc1 = cycles_by_k[k1]
                cyc2 = cycles_by_k[k2]
                for c1 in cyc1:
                    for c2 in cyc2:
                        if c1 & c2:
                            total_conf += 1
                n_pairs = len(cyc1) * len(cyc2)
                n_disj = n_pairs - total_conf
                avg_conf_per_k1 = total_conf / len(cyc1) if cyc1 else 0
                avg_conf_per_k2 = total_conf / len(cyc2) if cyc2 else 0

            if k1 == k2:
                print(f"    ({k1},{k2}): {total_conf} conflicts / {n_pairs} pairs, "
                      f"disjoint={n_disj}, avg_conf/cycle={avg_conf:.2f}")
            else:
                print(f"    ({k1},{k2}): {total_conf} conflicts / {n_pairs} pairs, "
                      f"disjoint={n_disj}, f({k1},{k2})={avg_conf_per_k1:.2f}, "
                      f"f({k2},{k1})={avg_conf_per_k2:.2f}")

    # Key: f(k,k')/c_k' ratio — is this constant per k?
    print(f"\n  Conflict ratio f(k,k')/c_k':")
    for k1 in sorted(cycles_by_k):
        ratios = []
        for k2 in sorted(cycles_by_k):
            if k1 == k2:
                cyc = cycles_by_k[k1]
                total_conf = 0
                for i in range(len(cyc)):
                    for j in range(i + 1, len(cyc)):
                        if cyc[i] & cyc[j]:
                            total_conf += 1
                avg_conf = 2 * total_conf / len(cyc) if cyc else 0
                ratio = avg_conf / (len(cyc) - 1) if len(cyc) > 1 else 0
            else:
                cyc1 = cycles_by_k[k1]
                cyc2 = cycles_by_k[k2]
                total_conf = sum(1 for c1 in cyc1 for c2 in cyc2 if c1 & c2)
                avg_conf = total_conf / len(cyc1) if cyc1 else 0
                ratio = avg_conf / len(cyc2) if cyc2 else 0
            ratios.append((k2, ratio))
        ratio_str = ', '.join(f'f/c{k2}={r:.6f}' for k2, r in ratios)
        print(f"    k={k1}: {ratio_str}")

    return c_k, cycles_by_k


# Main
print("=" * 70)
print("CONFLICT STRUCTURE ANALYSIS")
print("=" * 70)

for p in [7, 11]:
    m = (p - 1) // 2
    S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
    S_int = list(range(1, m + 1))

    print(f"\n{'='*70}")
    print(f"p = {p}")
    print(f"{'='*70}")

    for name, S in [("Paley", S_qr), ("Interval", S_int)]:
        analyze_conflict_structure(p, S, name)

print("\nDONE.")
