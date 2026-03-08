#!/usr/bin/env python3
"""
Permanent and Cycle Cover Analysis of Tournaments (kind-pasteur-S34).

For a tournament T with adjacency matrix A:
- per(A) = number of cycle covers (disjoint directed cycles covering all vertices)
- det(A) = signed count of cycle covers
- H(T) = number of Hamiltonian paths

Irving-Omar formula: H(D) = sum_S det(A_bar[S]) * per(A[S^c])
where A_bar[i][j] = 1 - A[i][j] and the sum is over subsets S of vertices.

Tests:
1. per(A) vs H: is there a simple relation?
2. det(A) vs H: correlation?
3. per(A) mod 2: always odd? (Redei says H always odd)
4. Irving-Omar decomposition verification
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, permutations
from collections import defaultdict
import sys


def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def count_ham_paths(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])


def permanent(M, n):
    """Compute permanent of n x n matrix M using Ryser's formula."""
    # Ryser: per(A) = (-1)^n * sum_{S subset [n]} (-1)^{|S|} prod_i sum_{j in S} A[i][j]
    result = 0
    for mask in range(1 << n):
        bits_set = bin(mask).count('1')
        sign = (-1) ** (n - bits_set)
        prod = 1
        for i in range(n):
            row_sum = 0
            for j in range(n):
                if mask & (1 << j):
                    row_sum += M[i][j]
            prod *= row_sum
        result += sign * prod
    return result


def determinant(M, n):
    """Compute determinant via Leibniz formula (brute force for small n)."""
    det = 0
    for perm in permutations(range(n)):
        sign = 1
        # Compute sign of permutation
        visited = [False] * n
        for i in range(n):
            if not visited[i]:
                cycle_len = 0
                j = i
                while not visited[j]:
                    visited[j] = True
                    j = perm[j]
                    cycle_len += 1
                if cycle_len % 2 == 0:
                    sign *= -1
        prod = 1
        for i in range(n):
            prod *= M[i][perm[i]]
        det += sign * prod
    return det


def main():
    print("=== Permanent & Cycle Cover Analysis ===\n")

    for n in [3, 4, 5]:
        print(f"--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        h_to_per = defaultdict(list)
        h_to_det = defaultdict(list)
        per_mod2 = defaultdict(int)

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            per_val = permanent(adj, n)
            det_val = determinant(adj, n)

            h_to_per[H].append(per_val)
            h_to_det[H].append(det_val)
            per_mod2[per_val % 2] += 1

        print(f"  per(A) mod 2: {dict(per_mod2)}")

        print(f"\n  H -> per(A) values:")
        for H in sorted(h_to_per):
            vals = sorted(set(h_to_per[H]))
            unique = "UNIQUE" if len(vals) == 1 else f"{len(vals)} distinct"
            print(f"    H={H:3d}: per(A) = {vals[:5]} ({unique})")

        print(f"\n  H -> det(A) values:")
        for H in sorted(h_to_det):
            vals = sorted(set(h_to_det[H]))
            unique = "UNIQUE" if len(vals) == 1 else f"{len(vals)} distinct"
            print(f"    H={H:3d}: det(A) = {vals[:5]} ({unique})")

        # Check: is per(A) always odd for tournaments?
        all_odd = all(p % 2 == 1 for pp in h_to_per.values() for p in pp)
        print(f"\n  per(A) always odd? {all_odd}")

        # Check: per(A) - H relationship
        print(f"\n  per(A) - H values:")
        diff_set = set()
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            per_val = permanent(adj, n)
            diff_set.add(per_val - H)
        print(f"    per(A) - H in: {sorted(diff_set)}")

        # Check: per(A) / H ratio
        ratio_set = set()
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            per_val = permanent(adj, n)
            if H > 0:
                ratio_set.add(round(per_val / H, 4))
        print(f"    per(A) / H in: {sorted(ratio_set)[:10]}")

        sys.stdout.flush()

    # Part 2: Irving-Omar decomposition at n=4
    print(f"\n\n=== Part 2: Irving-Omar Decomposition (n=4) ===")
    n = 4
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    for bits in range(min(8, max_bits)):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)

        # A_bar[i][j] = 1 - A[i][j] for i != j; A_bar[i][i] = 0
        adj_bar = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    adj_bar[i][j] = 1 - adj[i][j]

        # Irving-Omar: H(T) = sum_S det(A_bar[S]) * per(A[S^c])
        # where S is a subset of vertices
        io_sum = 0
        for mask in range(1 << n):
            S = [i for i in range(n) if mask & (1 << i)]
            Sc = [i for i in range(n) if not (mask & (1 << i))]
            k = len(S)
            kc = len(Sc)

            if k == 0:
                det_S = 1
            else:
                M_S = [[adj_bar[S[i]][S[j]] for j in range(k)] for i in range(k)]
                det_S = determinant(M_S, k)

            if kc == 0:
                per_Sc = 1
            else:
                M_Sc = [[adj[Sc[i]][Sc[j]] for j in range(kc)] for i in range(kc)]
                per_Sc = permanent(M_Sc, kc)

            io_sum += det_S * per_Sc

        match = "OK" if io_sum == H else f"FAIL ({io_sum})"
        print(f"  bits={bits:06b}: H={H}, IO_sum={io_sum}, {match}")


if __name__ == "__main__":
    main()
