"""
circulant_reduction.py — Test Step A: Does Z_p-averaging reduce to circulants?

Question: For tournament T on Z_p, define the "Z_p-averaged" tournament T̄ by:
  T̄(i→j) = 1 iff |{a ∈ Z_p : T(i+a → j+a)}| > p/2

This is the "majority vote" circulant. If H(T̄) ≥ H(T) for all T,
then max_T H(T) = max_{circulant} H(T).

But H might not be convex under this operation!

Alternative: The question is whether the global H-maximizer is always circulant
(or at least isomorphic to a circulant). This was verified exhaustively at n=7
(all 240 maximizers are Paley relabelings). Test at n=5 too.

For n=9 (composite), we can sample.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import defaultdict, Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
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


def zp_average(A, n):
    """Compute the Z_n-averaged circulant tournament from A.

    For each difference d ∈ {1,...,n-1}, count how many (i,j) pairs with j-i≡d have A[i][j]=1.
    If count > n/2, include d in the connection set.
    """
    S = set()
    for d in range(1, n):
        count = sum(A[i][(i+d)%n] for i in range(n))
        if count > n / 2:
            S.add(d)
        elif count == n / 2:
            # Tie: need to break somehow. Include d iff d < n-d (consistent orientation)
            if d < n - d:
                S.add(d)
            else:
                S.add(n - d)  # Actually this means n-d is included; skip d
    # Verify |S| = (n-1)/2
    if len(S) != (n-1)//2:
        # Fix: tournament requires exactly one of {d, n-d} in S
        S_fixed = set()
        for d in range(1, (n+1)//2):
            count_d = sum(A[i][(i+d)%n] for i in range(n))
            count_nd = sum(A[i][(i+n-d)%n] for i in range(n))
            if count_d >= count_nd:
                S_fixed.add(d)
            else:
                S_fixed.add(n-d)
        S = S_fixed

    return frozenset(S)


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def random_tournament(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def is_circulant(A, n):
    for i in range(n):
        for j in range(n):
            if A[i][j] != A[(i+1)%n][(j+1)%n]:
                return False
    return True


def score_sequence(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))


def main():
    print("CIRCULANT REDUCTION: DOES Z_p-AVERAGING INCREASE H?")
    print("=" * 75)

    # n=5: exhaustive
    n = 5
    print(f"\nn={n}: Exhaustive check (2^{n*(n-1)//2} = {2**(n*(n-1)//2)} tournaments)")
    t0 = time.time()
    increased = 0
    decreased = 0
    equal = 0
    total = 0
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    max_H_orig = 0
    max_H_avg = 0
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        H_orig = hamiltonian_paths_dp(A, n)
        S = zp_average(A, n)
        A_avg = circulant_adj(n, S)
        H_avg = hamiltonian_paths_dp(A_avg, n)
        if H_avg > H_orig:
            increased += 1
        elif H_avg < H_orig:
            decreased += 1
        else:
            equal += 1
        max_H_orig = max(max_H_orig, H_orig)
        max_H_avg = max(max_H_avg, H_avg)
        total += 1

    elapsed = time.time() - t0
    print(f"  {total} tournaments, {elapsed:.1f}s")
    print(f"  H_avg > H_orig: {increased} ({100*increased/total:.1f}%)")
    print(f"  H_avg < H_orig: {decreased} ({100*decreased/total:.1f}%)")
    print(f"  H_avg = H_orig: {equal} ({100*equal/total:.1f}%)")
    print(f"  Max H_orig = {max_H_orig}, Max H_avg = {max_H_avg}")

    # n=7: sampling
    for n in [7, 9, 11]:
        print(f"\nn={n}: Sampling 50000 random tournaments")
        rng = np.random.RandomState(42)
        increased = 0
        decreased = 0
        equal = 0
        max_H_orig = 0
        max_H_avg = 0
        t0 = time.time()

        for trial in range(50000):
            A = random_tournament(n, rng)
            H_orig = hamiltonian_paths_dp(A, n)
            S = zp_average(A, n)
            A_avg = circulant_adj(n, S)
            H_avg = hamiltonian_paths_dp(A_avg, n)

            if H_avg > H_orig:
                increased += 1
            elif H_avg < H_orig:
                decreased += 1
            else:
                equal += 1
            max_H_orig = max(max_H_orig, H_orig)
            max_H_avg = max(max_H_avg, H_avg)

            if trial % 10000 == 9999:
                elapsed = time.time() - t0
                print(f"  ... {trial+1} done ({elapsed:.1f}s): "
                      f"+{increased} -{decreased} ={equal}")

        elapsed = time.time() - t0
        total = 50000
        print(f"  {total} tournaments, {elapsed:.1f}s")
        print(f"  H_avg > H_orig: {increased} ({100*increased/total:.1f}%)")
        print(f"  H_avg < H_orig: {decreased} ({100*decreased/total:.1f}%)")
        print(f"  H_avg = H_orig: {equal} ({100*equal/total:.1f}%)")
        print(f"  Max H_orig = {max_H_orig}, Max H_avg = {max_H_avg}")

        if decreased > 0:
            print(f"  *** Z_n-averaging can DECREASE H! Not a valid reduction. ***")


if __name__ == '__main__':
    main()
    print("\nDONE.")
