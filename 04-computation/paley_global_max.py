"""
paley_global_max.py — Check if Paley maximizes H among ALL tournaments (not just circulants)

At n=7, Paley has H=189 among 2^((n-1)/2) = 8 circulants.
But there are n!/|Aut| ≈ 5040/21 = 240 distinct tournaments up to relabeling
(total labeled tournaments = 2^C(n,2) = 2^21 = 2097152).

We enumerate all labeled tournaments at n=7 and check:
1. Is Paley globally H-maximal?
2. What is the distribution of H values?
3. Do non-circulant tournaments ever beat Paley?

Also: at n=5, H=15 for all circulants. Is 15 the global max?

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    """Count Hamiltonian paths by bitmask DP."""
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1

    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]

    return sum(dp[full][v] for v in range(n))


def score_sequence(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))


def is_circulant(A, n):
    """Check if tournament is circulant (invariant under rotation)."""
    # Check if A[i][j] = A[(i+1)%n][(j+1)%n] for all i,j
    for i in range(n):
        for j in range(n):
            if A[i][j] != A[(i+1)%n][(j+1)%n]:
                return False
    return True


def enumerate_all_tournaments(n):
    """Enumerate ALL labeled tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def sample_tournaments(n, count, rng):
    """Sample random tournaments."""
    for _ in range(count):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        yield A


def main():
    print("PALEY GLOBAL MAXIMALITY CHECK")
    print("=" * 70)

    # n=5: exhaustive
    n = 5
    print(f"\nn={n}: Exhaustive enumeration of all 2^C({n},2) = {2**((n*(n-1))//2)} tournaments")
    t0 = time.time()
    H_dist = Counter()
    max_H = 0
    max_A = None
    total = 0
    for A in enumerate_all_tournaments(n):
        H = hamiltonian_paths_dp(A, n)
        H_dist[H] += 1
        if H > max_H:
            max_H = H
            max_A = A
        total += 1
    elapsed = time.time() - t0
    print(f"  {total} tournaments, {elapsed:.1f}s")
    print(f"  H distribution: {dict(sorted(H_dist.items()))}")
    print(f"  Max H = {max_H} (expected 15 from circulant analysis)")
    print(f"  Maximizer count: {H_dist[max_H]}")
    if max_A:
        print(f"  Maximizer score seq: {score_sequence(max_A, n)}")
        print(f"  Maximizer circulant: {is_circulant(max_A, n)}")

    # n=7: exhaustive (2^21 = 2097152 tournaments)
    n = 7
    print(f"\nn={n}: Exhaustive enumeration of all 2^C({n},2) = {2**((n*(n-1))//2)} tournaments")
    t0 = time.time()
    H_dist = Counter()
    max_H = 0
    max_count = 0
    max_circ = 0
    max_score_seqs = Counter()
    total = 0
    for A in enumerate_all_tournaments(n):
        H = hamiltonian_paths_dp(A, n)
        H_dist[H] += 1
        if H > max_H:
            max_H = H
            max_count = 1
            max_circ = 1 if is_circulant(A, n) else 0
            max_score_seqs = Counter()
            max_score_seqs[score_sequence(A, n)] = 1
        elif H == max_H:
            max_count += 1
            if is_circulant(A, n):
                max_circ += 1
            max_score_seqs[score_sequence(A, n)] += 1
        total += 1
        if total % 500000 == 0:
            print(f"  ... {total}/{2**21} done, current max H={max_H}")
    elapsed = time.time() - t0
    print(f"  {total} tournaments, {elapsed:.1f}s")
    print(f"  Distinct H values: {len(H_dist)}")
    print(f"  Max H = {max_H} (expected 189 from Paley circulant)")
    print(f"  Maximizers: {max_count} labeled, {max_circ} circulant")
    print(f"  Maximizer score seqs: {dict(max_score_seqs)}")

    # Top 10 H values
    top = sorted(H_dist.items(), key=lambda x: -x[0])[:10]
    print(f"  Top H values: {top}")

    # For the global maximum, check if all maximizers are relabelings of Paley
    # Paley at n=7: QR = {1,2,4}
    print(f"\n  Checking Paley structure of maximizers...")
    paley_S = frozenset([1, 2, 4])
    paley_A = [[0]*7 for _ in range(7)]
    for i in range(7):
        for s in paley_S:
            paley_A[i][(i+s)%7] = 1
    paley_H = hamiltonian_paths_dp(paley_A, 7)
    print(f"  Paley (S={set(paley_S)}) H = {paley_H}")

    # n=9: sampling (too large for exhaustive)
    n = 9
    print(f"\nn={n}: Sampling 500000 random tournaments")
    rng = np.random.RandomState(42)
    t0 = time.time()
    H_dist_9 = Counter()
    max_H_9 = 0
    total_9 = 0
    for A in sample_tournaments(n, 500000, rng):
        H = hamiltonian_paths_dp(A, n)
        H_dist_9[H] += 1
        if H > max_H_9:
            max_H_9 = H
        total_9 += 1
        if total_9 % 100000 == 0:
            print(f"  ... {total_9} done, current max H={max_H_9}")
    elapsed = time.time() - t0
    print(f"  Max H from sampling = {max_H_9} (circulant max was 3357)")
    print(f"  Top sampled: {sorted(H_dist_9.items(), key=lambda x: -x[0])[:5]}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
