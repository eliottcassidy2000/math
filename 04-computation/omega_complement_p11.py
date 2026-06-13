#!/usr/bin/env python3
"""
omega_complement_p11.py -- Conflict graph complement structure at p=11

At p=7, the Omega complement has exactly p=7 edges for Paley (a perfect matching
on 7 of the 80 vertices). Does this extend to p=11?

The complement edges are pairs of DISJOINT odd cycles -- cycles sharing no vertex.

Key questions:
1. How many complement edges (= disjoint cycle pairs = alpha_2) for each orbit?
2. What is the structure of the complement graph?
3. Is |E(Omega_complement)| a simple function of p for Paley?

At p=7:  Paley alpha_2 = 7 = p (complement has 7 edges)
         Interval alpha_2 = 14 = 2p (complement has 14 edges)

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations, permutations
from collections import defaultdict
import sys


def held_karp_cycles(A, verts):
    """Count directed Hamiltonian cycles on vertex subset using Held-Karp DP."""
    k = len(verts)
    if k < 3:
        return 0
    if k == 3:
        a, b, c = verts
        cnt = 0
        if A[a][b] and A[b][c] and A[c][a]:
            cnt += 1
        if A[a][c] and A[c][b] and A[b][a]:
            cnt += 1
        return cnt

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


def get_all_directed_cycles(A, p):
    """Get all directed odd cycles as (frozenset(vertices), multiplicity) pairs."""
    cycles = []
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_dir = held_karp_cycles(A, verts)
            if n_dir > 0:
                cycles.append((frozenset(subset), n_dir, k))
    return cycles


def main():
    for p in [7, 11]:
        m = (p - 1) // 2
        half = 1 << m
        pairs = [(s, p - s) for s in range(1, m + 1)]

        # Select representative orientations for each orbit
        if p == 7:
            test_bits = [0b000, 0b011]  # Interval, Paley
            labels = ["Interval", "Paley"]
        elif p == 11:
            # 4 orbits at p=11: bits 00000(Interval), 00010, 01010, 01011(Paley)
            test_bits = [0b00000, 0b00010, 0b01010, 0b01011]
            labels = ["Interval", "Orbit-B", "Orbit-C", "Paley"]

        print("=" * 70)
        print(f"OMEGA COMPLEMENT STRUCTURE at p={p}")
        print("=" * 70)

        for idx, (bits, label) in enumerate(zip(test_bits, labels)):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))

            A = [[0]*p for _ in range(p)]
            for v in range(p):
                for s in S:
                    A[v][(v + s) % p] = 1

            # Get all directed odd cycles
            cycles = get_all_directed_cycles(A, p)

            # Expand: each (fset, mult, k) becomes mult copies
            all_cycles = []
            for fset, mult, k in cycles:
                for _ in range(mult):
                    all_cycles.append((fset, k))

            n_total = len(all_cycles)

            # Count alpha_2 = disjoint pairs
            alpha_2 = 0
            disjoint_pairs = []
            for i in range(n_total):
                for j in range(i + 1, n_total):
                    if not (all_cycles[i][0] & all_cycles[j][0]):
                        alpha_2 += 1
                        disjoint_pairs.append((i, j))

            # Analyze disjoint pair structure
            # What lengths participate in disjoint pairs?
            len_pairs = defaultdict(int)
            for i, j in disjoint_pairs:
                k1, k2 = all_cycles[i][1], all_cycles[j][1]
                key = (min(k1, k2), max(k1, k2))
                len_pairs[key] += 1

            # Complement graph degree sequence
            comp_deg = [0] * n_total
            for i, j in disjoint_pairs:
                comp_deg[i] += 1
                comp_deg[j] += 1
            comp_deg_dist = defaultdict(int)
            for d in comp_deg:
                comp_deg_dist[d] += 1

            # Group cycles by length
            by_len = defaultdict(int)
            for fset, k in all_cycles:
                by_len[k] += 1

            # Check if complement is a matching (all degrees <= 1)
            max_comp_deg = max(comp_deg) if comp_deg else 0
            is_matching = max_comp_deg <= 1

            # Alpha_3: count independent triples in Omega
            # (= cliques of size 3 in complement)
            alpha_3 = 0
            if alpha_2 > 0 and n_total <= 200:
                # Build complement adjacency for efficiency
                comp_adj = [set() for _ in range(n_total)]
                for i, j in disjoint_pairs:
                    comp_adj[i].add(j)
                    comp_adj[j].add(i)

                for i in range(n_total):
                    for j in comp_adj[i]:
                        if j > i:
                            for k in comp_adj[i]:
                                if k > j and k in comp_adj[j]:
                                    alpha_3 += 1

            print(f"\n  {label} (bits={bits:0{m}b}, S={S}):")
            print(f"    Total directed cycles (alpha_1): {n_total}")
            print(f"    By length: {dict(sorted(by_len.items()))}")
            print(f"    Disjoint pairs (alpha_2): {alpha_2}")
            print(f"    alpha_2/p = {alpha_2/p:.4f}")
            print(f"    alpha_2 mod p = {alpha_2 % p}")
            print(f"    Disjoint pairs by length: {dict(sorted(len_pairs.items()))}")
            print(f"    Complement degree dist: {dict(sorted(comp_deg_dist.items()))}")
            print(f"    Max complement degree: {max_comp_deg}")
            print(f"    Is matching? {is_matching}")
            if n_total <= 200:
                print(f"    Alpha_3 (disjoint triples): {alpha_3}")

            # H via OCF
            H = 1 + 2 * n_total + 4 * alpha_2 + 8 * alpha_3
            print(f"    H = 1 + 2*{n_total} + 4*{alpha_2} + 8*{alpha_3} = {H}")

            # Check: is alpha_2 = c * p for some integer c?
            if alpha_2 % p == 0:
                print(f"    *** alpha_2 = {alpha_2//p} * p  (p-DIVISIBLE!)")
            else:
                print(f"    alpha_2 NOT p-divisible")

            # Vertex participation in disjoint pairs
            # Which vertex sets participate?
            if alpha_2 > 0 and alpha_2 <= 100:
                fsets_in_disjoint = set()
                for i, j in disjoint_pairs:
                    fsets_in_disjoint.add(all_cycles[i][0])
                    fsets_in_disjoint.add(all_cycles[j][0])
                print(f"    Distinct vertex sets in disjoint pairs: {len(fsets_in_disjoint)}")

                # Check if all disjoint pairs are (3-cycle, k-cycle) or (k-cycle, k-cycle)
                # At p=7: disjoint means they share 0 of 7 vertices, so sizes must sum to <= 7
                # For 3-cycle pairs: 3+3=6 <= 7, so (3,3) and (3,4) and (4,3) possible
                # Wait: at p=7, vertex set is {0,...,6}, so disjoint subsets must sum to <= 7
                # (3,3): 6 <= 7 OK. (3,5): 8 > 7 NOT OK. So only (3,3) and (3,4) if 4-cycles existed
                # But 4-cycles are EVEN — only odd cycles. So only (3,3) at p=7.
                print(f"    Size constraint: vertex sets sum <= {p}")
                for (k1, k2), cnt in sorted(len_pairs.items()):
                    feasible = k1 + k2 <= p
                    print(f"      ({k1},{k2}): {cnt} pairs, sum={k1+k2} {'<= p OK' if feasible else '> p IMPOSSIBLE'}")

        # Summary comparison
        print(f"\n  {'='*60}")
        print(f"  SUMMARY at p={p}:")
        print(f"  Only odd-length pairs with sizes summing to <= p can be disjoint")
        print(f"  At p={p}: possible disjoint length pairs: ", end="")
        possible = []
        for k1 in range(3, p + 1, 2):
            for k2 in range(k1, p + 1, 2):
                if k1 + k2 <= p:
                    possible.append((k1, k2))
        print(possible)

    print("\nDONE.")


if __name__ == '__main__':
    main()
