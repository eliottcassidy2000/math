#!/usr/bin/env python3
"""
max_H_circulant_p13.py -- Find the H-maximizing circulant tournament at p=13

At p=13 (1 mod 4), no Paley tournament exists. There are 2^6 = 64 circulant
tournaments on Z_13 (choose one element from each pair {s, 13-s}).

Which one maximizes H(T)?

Also computes for p=7, 11 to verify Paley maximality.

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


def count_ham_paths(A, n):
    """Held-Karp DP for Hamiltonian path count."""
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def all_circulant_tournaments(p):
    """Generate all circulant tournament connection sets on Z_p."""
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    # Choose one from each pair
    result = []
    for bits in range(1 << m):
        S = []
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.append(a)
            else:
                S.append(b)
        result.append(sorted(S))
    return result


def main():
    print("=" * 70)
    print("H-MAXIMIZING CIRCULANT TOURNAMENT")
    print("=" * 70)

    for p in [7, 11, 13, 17]:
        m = (p - 1) // 2
        tournaments = all_circulant_tournaments(p)
        n_tours = len(tournaments)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}: {n_tours} circulant tournaments")
        print(f"{'='*70}")

        # Known tournaments
        S_int = list(range(1, m + 1))
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        paley_exists = (p % 4 == 3) and len(S_qr) == m

        if p <= 13:
            # Compute H for all circulant tournaments
            results = []
            t0 = time.time()
            for S in tournaments:
                A = build_adj(p, S)
                H = count_ham_paths(A, p)
                results.append((H, S))
            t1 = time.time()

            results.sort(reverse=True)

            print(f"\nAll H values computed in {t1-t0:.1f}s")
            print(f"\nTop 10:")
            for i, (H, S) in enumerate(results[:10]):
                labels = []
                if S == S_int:
                    labels.append("INTERVAL")
                if paley_exists and S == S_qr:
                    labels.append("PALEY")
                label = " <-- " + ", ".join(labels) if labels else ""
                print(f"  {i+1}. H={H:>10}, S={S}{label}")

            print(f"\nBottom 5:")
            for i, (H, S) in enumerate(results[-5:]):
                rank = n_tours - 4 + i
                labels = []
                if S == S_int:
                    labels.append("INTERVAL")
                if paley_exists and S == S_qr:
                    labels.append("PALEY")
                label = " <-- " + ", ".join(labels) if labels else ""
                print(f"  {rank}. H={H:>10}, S={S}{label}")

            # Find Interval and Paley ranks
            for i, (H, S) in enumerate(results):
                if S == S_int:
                    print(f"\nInterval: rank {i+1}/{n_tours}, H={H}")
                if paley_exists and S == S_qr:
                    print(f"Paley:    rank {i+1}/{n_tours}, H={H}")

            # H distribution statistics
            H_values = [H for H, _ in results]
            print(f"\nH statistics:")
            print(f"  Max: {max(H_values)}")
            print(f"  Min: {min(H_values)}")
            print(f"  Mean: {sum(H_values)/len(H_values):.1f}")
            print(f"  Distinct values: {len(set(H_values))}")

            # Check if maximizer is unique
            max_H = results[0][0]
            maximizers = [(H, S) for H, S in results if H == max_H]
            print(f"\n  Number of maximizers: {len(maximizers)}")
            for H, S in maximizers:
                labels = []
                if S == S_int:
                    labels.append("INTERVAL")
                if paley_exists and S == S_qr:
                    labels.append("PALEY")
                # Check if S is QR for some other structure
                label = " <-- " + ", ".join(labels) if labels else ""
                print(f"    S={S}{label}")

        else:
            # p=17: too many tournaments (2^8=256), but each is fast
            results = []
            t0 = time.time()
            for S in tournaments:
                A = build_adj(p, S)
                H = count_ham_paths(A, p)
                results.append((H, S))
            t1 = time.time()
            results.sort(reverse=True)

            print(f"\nAll H values computed in {t1-t0:.1f}s")
            print(f"\nTop 5:")
            for i, (H, S) in enumerate(results[:5]):
                labels = []
                if S == S_int:
                    labels.append("INTERVAL")
                if paley_exists and S == S_qr:
                    labels.append("PALEY")
                label = " <-- " + ", ".join(labels) if labels else ""
                print(f"  {i+1}. H={H:>15}, S={S}{label}")

            for i, (H, S) in enumerate(results):
                if S == S_int:
                    print(f"\nInterval: rank {i+1}/{n_tours}, H={H}")
                if paley_exists and S == S_qr:
                    print(f"Paley:    rank {i+1}/{n_tours}, H={H}")

            max_H = results[0][0]
            maximizers = [(H, S) for H, S in results if H == max_H]
            print(f"\nNumber of maximizers: {len(maximizers)}")
            for H, S in maximizers:
                labels = []
                if S == S_int:
                    labels.append("INTERVAL")
                if paley_exists and S == S_qr:
                    labels.append("PALEY")
                label = " <-- " + ", ".join(labels) if labels else ""
                print(f"  S={S}{label}")


if __name__ == '__main__':
    main()
