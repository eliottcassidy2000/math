#!/usr/bin/env python3
"""
alpha_directed_p11.py -- Correct OCF decomposition using DIRECTED cycles

The OCF H(T) = I(Omega(T), 2) uses DIRECTED odd cycles as vertices
of the conflict graph Omega. Each vertex set supporting k-Ham-cycles
contributes k distinct vertices to Omega (all mutually adjacent since
they share ALL vertices).

Key correction from alpha_full_p11.py: using vertex sets instead of
directed cycles gave WRONG I(Omega, 2) values.

At p=7:
  Interval: 59 directed cycles → I = 175 ✓
  Paley: 80 directed cycles → I = 189 ✓

This script computes the full decomposition at p=11 using directed cycles.

Author: kind-pasteur-2026-03-12-S59c
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


def count_directed_ham_cycles(A, verts):
    """Count number of directed Hamiltonian cycles on vertex subset."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
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
        if key in dp and dp[key] > 0 and A[verts[v]][verts[start]]:
            total += dp[key]
    return total


def enumerate_directed_cycles(A, p, max_k=None):
    """Enumerate all directed odd cycles. Each vertex set with d directed
    Ham cycles produces d entries in the list.

    Returns: list of (frozenset(verts), length, multiplicity_index)
    Also returns by_k dict: k -> [(frozenset, directed_count)]
    """
    if max_k is None:
        max_k = p

    all_directed = []  # list of frozenset
    by_k_detail = {}   # k -> list of (frozenset, num_directed_on_this_set)

    for k in range(3, max_k + 1, 2):
        sets_k = []
        for subset in combinations(range(p), k):
            d = count_directed_ham_cycles(A, list(subset))
            if d > 0:
                sets_k.append((frozenset(subset), d))
                for _ in range(d):
                    all_directed.append(frozenset(subset))
        by_k_detail[k] = sets_k

    return all_directed, by_k_detail


def main():
    print("=" * 70)
    print("CORRECT OCF DECOMPOSITION USING DIRECTED CYCLES")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)

            t0 = time.time()
            all_directed, by_k = enumerate_directed_cycles(A, p)
            t1 = time.time()

            n = len(all_directed)
            print(f"\n  {name} (S={S}): {n} directed cycles ({t1-t0:.1f}s)")
            for k in sorted(by_k):
                vertex_sets = len(by_k[k])
                directed = sum(d for _, d in by_k[k])
                # Distribution of directed cycles per vertex set
                mult_dist = defaultdict(int)
                for _, d in by_k[k]:
                    mult_dist[d] += 1
                print(f"    k={k}: {vertex_sets} vertex sets, {directed} directed cycles")
                print(f"           multiplicity dist: {dict(sorted(mult_dist.items()))}")

            # Compute I(Omega, 2)
            # Two directed cycles conflict iff they share a vertex.
            # Note: two different directed cycles on the SAME vertex set
            # always conflict (they share ALL vertices).

            # For efficiency, avoid O(n^2) for large n.
            # Instead, use the fact that:
            # I(Omega, 2) = product over cliques of (1 + d_j * x)
            #   where d_j is the multiplicity on vertex set j...
            # NO! That's only if the cliques are disjoint. In general,
            # different vertex sets CAN overlap, creating non-trivial structure.

            # For p=7: n=59/80, feasible for direct computation.
            # For p=11: n could be large. Let me check.

            if n <= 100:
                # Direct computation
                # Build neighbor bitmask
                nbr = [0] * n
                for i in range(n):
                    for j in range(i + 1, n):
                        if all_directed[i] & all_directed[j]:
                            nbr[i] |= (1 << j)
                            nbr[j] |= (1 << i)

                alpha = defaultdict(int)

                def backtrack(v, excluded, size):
                    alpha[size] += 1
                    for w in range(v, n):
                        if not (excluded & (1 << w)):
                            backtrack(w + 1, excluded | nbr[w], size + 1)

                backtrack(-1, 0, 0)

                H = sum(alpha[j] * (2**j) for j in alpha)
                max_j = max(alpha.keys())

                print(f"\n    I(Omega, 2) = H = {H}")
                print(f"    Alpha decomposition:")
                cumulative = 0
                for j in range(max_j + 1):
                    a_j = alpha.get(j, 0)
                    contribution = a_j * (2**j)
                    cumulative += contribution
                    pct = 100 * contribution / H
                    print(f"      alpha_{j} = {a_j:>8}, 2^{j}*alpha_{j} = {contribution:>10} ({pct:>6.2f}%)")

            else:
                # For larger n, compute alpha_1, alpha_2, alpha_3
                print(f"\n    n={n} too large for full backtracking, computing alpha_1,2,3")

                alpha_1 = n

                # alpha_2: count non-adjacent pairs
                t2 = time.time()
                alpha_2 = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        if not (all_directed[i] & all_directed[j]):
                            alpha_2 += 1
                t3 = time.time()

                print(f"    alpha_1 = {alpha_1}")
                print(f"    alpha_2 = {alpha_2} ({t3-t2:.1f}s)")

                # alpha_3
                t4 = time.time()
                alpha_3 = 0
                # Pre-compute non-adjacent pairs for each i
                non_adj = [[] for _ in range(n)]
                for i in range(n):
                    for j in range(i + 1, n):
                        if not (all_directed[i] & all_directed[j]):
                            non_adj[i].append(j)

                for i in range(n):
                    for j in non_adj[i]:
                        union_ij = all_directed[i] | all_directed[j]
                        for l_idx in range(len(non_adj[j])):
                            l = non_adj[j][l_idx]
                            if l <= j:
                                continue
                            if not (all_directed[i] & all_directed[l]):
                                alpha_3 += 1
                t5 = time.time()

                print(f"    alpha_3 = {alpha_3} ({t5-t4:.1f}s)")

                # H estimates
                H_known = {
                    (7, 'Interval'): 175,
                    (7, 'Paley'): 189,
                    (11, 'Interval'): 93027,
                    (11, 'Paley'): 95095,
                }
                H = H_known.get((p, name), 0)

                H_3 = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
                remainder = H - H_3

                print(f"\n    H decomposition:")
                print(f"      Level 0: 1")
                print(f"      Level 1: 2*{alpha_1} = {2*alpha_1}")
                print(f"      Level 2: 4*{alpha_2} = {4*alpha_2}")
                print(f"      Level 3: 8*{alpha_3} = {8*alpha_3}")
                print(f"      H_3 = {H_3}")
                print(f"      H = {H}")
                print(f"      Remainder (Level 4+) = {remainder}")

                # Percentages
                if H > 0:
                    print(f"\n    Percentages:")
                    print(f"      Level 0: {100*1/H:.4f}%")
                    print(f"      Level 1: {100*2*alpha_1/H:.4f}%")
                    print(f"      Level 2: {100*4*alpha_2/H:.4f}%")
                    print(f"      Level 3: {100*8*alpha_3/H:.4f}%")
                    print(f"      Level 4+: {100*remainder/H:.4f}%")

        # Comparison
        if p == 11:
            print(f"\n  {'='*60}")
            print(f"  INTERVAL vs PALEY COMPARISON AT p=11")
            print(f"  {'='*60}")
            print(f"  (see output above for details)")


if __name__ == '__main__':
    main()
