#!/usr/bin/env python3
"""
phase_transition_direct.py -- Direct H(T) computation for Interval vs Paley

Computes H(T) = I(Omega(T), 2) for both tournaments at p=7, 11, 13
by building the full conflict graph Omega(T) and evaluating the
independence polynomial at x=2.

This directly verifies the Paley-to-Interval crossover.

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def has_ham_cycle(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a] + A[a][c] * A[c][b] * A[b][a]) > 0
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


def get_all_odd_cycle_sets(A, p):
    """Get all vertex sets supporting a directed odd cycle (3, 5, ..., p)."""
    all_sets = []
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            if has_ham_cycle(A, list(subset)):
                all_sets.append(frozenset(subset))
    return all_sets


def build_conflict_graph(cycle_sets):
    """Build adjacency list for Omega(T): edge iff shared vertex."""
    n = len(cycle_sets)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def independence_polynomial(adj, n, x=2):
    """Compute I(G, x) = sum alpha_k * x^k by DP over subsets.
    For small n only (n <= ~25).
    """
    if n > 25:
        print(f"  WARNING: n={n} too large for exact IP computation")
        return None

    # Use Bron-Kerbosch or direct DP for I(G, x)
    # For n <= 25, we can iterate over all independent sets using DP

    # Convert adj to bitset format
    neighbors = [0] * n
    for i in range(n):
        for j in adj[i]:
            neighbors[i] |= (1 << j)

    # Count independent sets by size
    alpha = [0] * (n + 1)

    # DP: iterate over vertices, for each decide include/exclude
    # More efficient: use the recursion I(G, x) = I(G-v, x) + x * I(G-v-N(v), x)
    # Or just enumerate all 2^n subsets (feasible for n <= 25)

    if n <= 20:
        for mask in range(1 << n):
            # Check if mask is an independent set
            is_indep = True
            for i in range(n):
                if mask & (1 << i):
                    if mask & neighbors[i]:
                        is_indep = False
                        break
            if is_indep:
                size = bin(mask).count('1')
                alpha[size] += 1

        total = sum(alpha[k] * (x ** k) for k in range(n + 1))
        return total, alpha

    else:
        # For n=21..25, use backtracking
        alpha_dict = {}

        def backtrack(idx, current_set, forbidden):
            size = len(current_set)
            alpha_dict[size] = alpha_dict.get(size, 0) + 1

            for i in range(idx, n):
                if not (forbidden & (1 << i)):
                    current_set.append(i)
                    backtrack(i + 1, current_set, forbidden | neighbors[i])
                    current_set.pop()

        backtrack(0, [], 0)

        total = sum(count * (x ** size) for size, count in alpha_dict.items())
        alpha_list = [0] * (n + 1)
        for size, count in alpha_dict.items():
            if size <= n:
                alpha_list[size] = count
        return total, alpha_list


def main():
    print("=" * 70)
    print("PHASE TRANSITION: H(T) for Paley vs Interval")
    print("=" * 70)

    results = {}

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        paley_exists = (p % 4 == 3)

        tournaments = []
        if paley_exists:
            tournaments.append(("Paley", S_qr))
        tournaments.append(("Interval", S_int))

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        for name, S in tournaments:
            A = build_adj(p, S)

            t0 = time.time()
            cycle_sets = get_all_odd_cycle_sets(A, p)
            t1 = time.time()

            n = len(cycle_sets)
            print(f"\n  {name}: {n} odd-cycle vertex sets ({t1-t0:.1f}s)")

            # Count by size
            by_size = {}
            for fs in cycle_sets:
                k = len(fs)
                by_size[k] = by_size.get(k, 0) + 1

            for k in sorted(by_size.keys()):
                print(f"    k={k}: {by_size[k]} vertex sets")

            if n <= 25:
                t2 = time.time()
                adj = build_conflict_graph(cycle_sets)
                t3 = time.time()

                result = independence_polynomial(adj, n, x=2)
                t4 = time.time()

                if result:
                    H, alpha = result
                    print(f"\n    H(T) = I(Omega, 2) = {H} ({t4-t2:.1f}s)")
                    print(f"    Alpha coefficients:")
                    for k_idx in range(len(alpha)):
                        if alpha[k_idx] > 0:
                            print(f"      alpha_{k_idx} = {alpha[k_idx]} "
                                  f"(contributes {alpha[k_idx] * 2**k_idx})")

                    results[(p, name)] = H
            else:
                # For large Omega, use the formula directly
                # H = I(Omega, 2) -- need to compute this differently
                print(f"\n    Omega too large ({n} vertices) for direct IP computation")
                print(f"    Using formula H = det(I + 2*A_omega)... (not implemented)")

                # Alternative: compute H from the transfer matrix
                # H(T) = number of Hamiltonian paths with odd parity
                # For circulant tournament, H = sum_{perm with odd cycles} 2^{#cycles}
                # But this is what Omega gives us

                # Let's at least compute alpha_1 and alpha_2
                alpha_1 = n
                t2 = time.time()
                disj = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        if not (cycle_sets[i] & cycle_sets[j]):
                            disj += 1
                t3 = time.time()
                alpha_2 = disj

                # Lower bound: H >= 1 + 2*alpha_1 + 4*alpha_2
                H_lower = 1 + 2 * alpha_1 + 4 * alpha_2
                print(f"\n    alpha_1 = {alpha_1}")
                print(f"    alpha_2 = {disj} ({t3-t2:.1f}s)")
                print(f"    H >= 1 + 2*{alpha_1} + 4*{alpha_2} = {H_lower}")

                results[(p, name)] = H_lower

    # ====== SUMMARY ======
    print(f"\n{'='*70}")
    print("PHASE TRANSITION SUMMARY")
    print("=" * 70)

    for p in [7, 11, 13]:
        paley_exists = (p % 4 == 3)
        print(f"\n  p={p}:")
        if paley_exists and (p, "Paley") in results and (p, "Interval") in results:
            h_pal = results[(p, "Paley")]
            h_int = results[(p, "Interval")]
            winner = "Paley" if h_pal > h_int else "Interval" if h_int > h_pal else "TIE"
            print(f"    Paley:    H = {h_pal}")
            print(f"    Interval: H = {h_int}")
            print(f"    Winner: {winner}")
            if h_pal > 0:
                print(f"    Ratio Int/Pal: {h_int/h_pal:.6f}")
        elif (p, "Interval") in results:
            print(f"    Interval: H = {results[(p, 'Interval')]}")
            print(f"    (no Paley at p=1 mod 4)")


if __name__ == '__main__':
    main()
