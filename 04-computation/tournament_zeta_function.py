#!/usr/bin/env python3
"""
Tournament Zeta Function and Cycle Structure (kind-pasteur-S34).

Novel idea: Define a zeta-like function for a tournament T:

Z_T(s) = prod_{k odd} (1 - c_k * x^k)^{-1}  (cycle zeta)

or the Ihara-type zeta function:

zeta_T(u) = prod_{[C] prime cycle} (1 - u^|C|)^{-1}

For tournaments, directed cycles are the "primes". The connection to H(T)
via OCF (H = I(Omega, 2)) suggests Z_T evaluated at specific points
should encode H.

Also tests the "cycle Euler product":
prod_{C odd cycle} (1 + 2^|C|) = ?? (relates to I(Omega, 2) if cycles were independent)

Key question: Is H(T) determined by the cycle count vector (c_3, c_5, c_7, ...)?
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, permutations
from math import comb


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


def count_directed_cycles_per_vset(adj, n, k):
    """Count directed k-cycles for each k-vertex set."""
    results = {}
    for vs in combinations(range(n), k):
        vl = list(vs)
        first = vl[0]
        count = 0
        for perm in permutations(vl[1:]):
            cycle = [first] + list(perm)
            valid = True
            for i in range(k):
                if not adj[cycle[i]][cycle[(i+1) % k]]:
                    valid = False
                    break
            if valid:
                count += 1
        if count > 0:
            results[frozenset(vs)] = count
    return results


def get_cycle_counts(adj, n):
    """Get (c_3, c_5, c_7, ...) — total directed cycles of each odd length."""
    counts = {}
    for k in range(3, n+1, 2):
        total = 0
        for vs in combinations(range(n), k):
            vl = list(vs)
            first = vl[0]
            for perm in permutations(vl[1:]):
                cycle = [first] + list(perm)
                valid = True
                for i in range(k):
                    if not adj[cycle[i]][cycle[(i+1) % k]]:
                        valid = False
                        break
                if valid:
                    total += 1
        counts[k] = total
    return counts


def main():
    print("=== Tournament Zeta Function Analysis ===\n")

    # n=5: Is H determined by (c_3, c_5)?
    print("--- n=5: H vs cycle counts ---")
    n = 5
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    cycle_to_h = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        cc = get_cycle_counts(adj, n)
        key = (cc.get(3, 0), cc.get(5, 0))

        if key not in cycle_to_h:
            cycle_to_h[key] = set()
        cycle_to_h[key].add(H)

    print(f"  (c3, c5) -> H values:")
    for key in sorted(cycle_to_h):
        h_vals = sorted(cycle_to_h[key])
        det = "UNIQUE" if len(h_vals) == 1 else "MULTIPLE"
        print(f"    (c3={key[0]}, c5={key[1]}): H = {h_vals} ({det})")

    # n=6: Is H determined by (c_3, c_5)?
    print("\n\n--- n=6: H vs cycle counts ---")
    n = 6
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    cycle_to_h = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        # Only count c3 and c5 (c7 not possible at n=6... wait, it is not,
        # since we need 7 vertices for a 7-cycle but only have 6)
        # Actually c5 needs only 5 vertices, so it exists at n=6.
        cc3 = 0
        for a, b, c in combinations(range(n), 3):
            if (adj[a][b] and adj[b][c] and adj[c][a]) or \
               (adj[a][c] and adj[c][b] and adj[b][a]):
                cc3 += 1
        cc5 = 0
        for vs in combinations(range(n), 5):
            vl = list(vs)
            first = vl[0]
            for perm in permutations(vl[1:]):
                cycle = [first] + list(perm)
                valid = True
                for i in range(5):
                    if not adj[cycle[i]][cycle[(i+1) % 5]]:
                        valid = False
                        break
                if valid:
                    cc5 += 1

        key = (cc3, cc5)
        if key not in cycle_to_h:
            cycle_to_h[key] = set()
        cycle_to_h[key].add(H)

    non_unique = 0
    for key in sorted(cycle_to_h):
        h_vals = sorted(cycle_to_h[key])
        if len(h_vals) > 1:
            non_unique += 1
            print(f"    (c3={key[0]:2d}, c5={key[1]:2d}): H = {h_vals} (MULTIPLE)")

    print(f"\n  Total (c3,c5) classes: {len(cycle_to_h)}")
    print(f"  Non-unique H: {non_unique} classes")
    if non_unique == 0:
        print(f"  *** H IS DETERMINED by (c3, c5) at n=6! ***")
    else:
        print(f"  H is NOT determined by (c3, c5) at n=6")

    # Key question: does H = f(c3, c5) have a clean formula?
    print("\n\n--- n=5: H formula from cycle counts ---")
    n = 5
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    for bits in range(0, max_bits, max_bits // 20 + 1):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        cc = get_cycle_counts(adj, n)
        c3, c5 = cc.get(3, 0), cc.get(5, 0)
        # Test: H = 1 + 2*(c3 + c5)?
        # At n=5, alpha_1 = c3 + c5 (each vertex set contributes 1 cycle)
        test = 1 + 2*(c3 + c5)
        print(f"  bits={bits:0{num_edges}b}: c3={c3}, c5={c5}, "
              f"H={H}, 1+2*(c3+c5)={test}, match={H==test}")

    print("\n\n--- n=6: Test if H = 1 + 2*alpha1 + 4*alpha2 ---")
    print("  (alpha1 = sum of directed cycle counts, alpha2 = sum over disjoint pairs)")


if __name__ == "__main__":
    main()
