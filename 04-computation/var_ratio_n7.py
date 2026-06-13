"""
var_ratio_n7.py -- kind-pasteur-2026-03-14-S105f
Compute EXACT Var(H)/Mean(H)^2 at n=7 using Held-Karp DP.

n=7 has m=21 arcs, 2^21 = 2,097,152 tournaments.
Each Held-Karp computation is O(2^n * n^2) = O(128 * 49) = O(6272).
Total: ~13 billion operations. Might take a while.

Actually: let's use numpy arrays and batch processing.
"""

import sys, math
import time
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def count_ham_paths(adj, n):
    """Count Hamiltonian paths using Held-Karp DP."""
    dp = [0] * (n * (1 << n))
    for v in range(n):
        dp[v * (1 << n) + (1 << v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[v * (1 << n) + mask]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[u * (1 << n) + (mask | (1 << u))] += val
    full = (1 << n) - 1
    total = 0
    for v in range(n):
        total += dp[v * (1 << n) + full]
    return total

def main():
    n = 7
    m = n*(n-1)//2  # 21
    total_tournaments = 1 << m  # 2^21 = 2,097,152

    print(f"Computing EXACT Var(H)/Mean(H)^2 at n={n}")
    print(f"m = {m}, total tournaments = {total_tournaments}")
    print(f"This will take a while...")

    sum_h = 0
    sum_h2 = 0
    start = time.time()
    checkpoint = total_tournaments // 20

    for bits in range(total_tournaments):
        adj = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
                idx += 1

        h = count_ham_paths(adj, n)
        sum_h += h
        sum_h2 += h * h

        if (bits + 1) % checkpoint == 0:
            elapsed = time.time() - start
            pct = 100 * (bits + 1) / total_tournaments
            est_total = elapsed / (bits + 1) * total_tournaments
            print(f"  {pct:.0f}% done, elapsed {elapsed:.1f}s, "
                  f"est total {est_total:.1f}s, "
                  f"sum_h so far = {sum_h}")

    elapsed = time.time() - start
    print(f"\nDone in {elapsed:.1f}s")

    N = total_tournaments
    mean = Fraction(sum_h, N)
    mean_h2 = Fraction(sum_h2, N)
    var = mean_h2 - mean * mean
    ratio = var / (mean * mean)

    level2 = Fraction(2*(n-2), n*(n-1))
    e4_plus = ratio - level2

    print(f"\n  sum(H) = {sum_h}")
    print(f"  sum(H^2) = {sum_h2}")
    print(f"  Mean = {mean} = {float(mean):.6f}")
    print(f"  Mean(H^2) = {mean_h2} = {float(mean_h2):.6f}")
    print(f"  Var = {var} = {float(var):.6f}")
    print(f"  Var/Mean^2 = {ratio} = {float(ratio):.10f}")
    print(f"  Level-2 = {level2} = {float(level2):.10f}")
    print(f"  E_4+/E_0 = {e4_plus} = {float(e4_plus):.10f}")

    diff = Fraction(1, 3) - ratio
    print(f"\n  Deviation from 1/3 = {diff} = {float(diff):.10f}")
    print(f"  In 180ths: {diff * 180}")
    print(f"  Predicted (k=3): 15/180 = {15/180:.10f}")
    print(f"  Actual: {float(diff * 180):.6f}/180")

if __name__ == '__main__':
    main()
