#!/usr/bin/env python3
"""
p-adic structure of H(T) beyond p=2 (kind-pasteur-S34).

We know H(T) is always odd (p=2 structure from Redei/OCF).
What about H(T) mod 3, mod 5, mod 7?

Tests:
1. Distribution of H(T) mod p for p=3,5,7
2. Does H mod 3 depend on c3? (H = 1 + 2*alpha_1 + 4*alpha_2 + ... mod 3)
3. Does the 2-adic tower interact with 3-adic structure?
4. Are there impossible residues?
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, permutations
from collections import defaultdict


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


def count_3cycles(adj, n):
    """Count directed 3-cycles (each counted once)."""
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1
    return c3


def main():
    print("=== p-adic Structure of H(T) Beyond p=2 ===\n")

    import sys
    for n in [4, 5, 6]:
        print(f"--- n={n} ---"); sys.stdout.flush()
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        # Distribution of H mod p
        for p in [3, 5, 7]:
            h_mod_p = defaultdict(int)
            c3_to_hmod = defaultdict(lambda: defaultdict(int))

            for bits in range(max_bits):
                adj = tournament_adj(n, bits)
                H = count_ham_paths(adj, n)
                c3 = count_3cycles(adj, n)
                r = H % p
                h_mod_p[r] += 1
                c3_to_hmod[c3][r] += 1

            print(f"\n  H mod {p}:")
            for r in sorted(h_mod_p):
                pct = 100 * h_mod_p[r] / max_bits
                print(f"    H = {r} (mod {p}): {h_mod_p[r]:6d} ({pct:.1f}%)")

            # Check if H mod p depends on c3
            if p == 3:
                print(f"\n  H mod 3 by c3:")
                for c3 in sorted(c3_to_hmod):
                    residues = c3_to_hmod[c3]
                    total = sum(residues.values())
                    parts = ", ".join(f"{r}:{residues[r]}" for r in sorted(residues))
                    unique = len(residues) == 1
                    print(f"    c3={c3:2d}: {parts} {'<-- UNIQUE' if unique else ''}")

        # Check: does H mod 3 have a formula involving c3?
        print(f"\n  Formula check: H mod 3 vs (1 + 2*c3) mod 3:")
        match_count = 0
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            c3 = count_3cycles(adj, n)
            if H % 3 == (1 + 2*c3) % 3:
                match_count += 1
        pct = 100 * match_count / max_bits
        print(f"    Match: {match_count}/{max_bits} ({pct:.1f}%)")

        # Check: H mod 3 vs alpha_1 mod 3
        print(f"\n  Formula check: H mod 3 vs alpha_1-based predictions:")
        # At n=5, Omega is complete, so alpha_1 = sum of all odd cycle counts
        # But alpha_1 is harder to compute. Instead check H mod 3 = (1 + 2*alpha_1) mod 3
        # which = (H mod 4 - 1)/2 mod 3 if we know H = 1 + 2*alpha_1 + 4*alpha_2 + ...
        # Actually H mod 3 = (1 + 2*alpha_1 + alpha_2 + 2*alpha_3) mod 3
        # because 4 = 1, 8 = 2, 16 = 1 mod 3
        # So H mod 3 = (1 + 2*alpha_1 + alpha_2 + 2*alpha_3 + alpha_4 + ...) mod 3
        # = (1 + 2*(alpha_1 + alpha_3 + ...) + (alpha_2 + alpha_4 + ...)) mod 3
        # At n=5: alpha_2 = 0, so H mod 3 = (1 + 2*alpha_1) mod 3
        # At n=6: alpha_3 = 0, so H mod 3 = (1 + 2*alpha_1 + alpha_2) mod 3
        print(f"    (Analysis requires alpha decomposition — see ehrhart_fast.py)")

        print()


if __name__ == "__main__":
    main()
