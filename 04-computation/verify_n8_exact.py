"""
verify_n8_exact.py -- kind-pasteur-2026-03-14-S107g
Compute EXACT Var(H)/Mean(H)^2 at n=8 to confirm prediction 59/252.

n=8: m=28 arcs, 2^28 = 268,435,456 tournaments.
TOO MANY for exhaustive. Use Monte Carlo with large sample.
"""

import sys, math, random, time
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def count_ham_paths_fast(adj, n):
    """Held-Karp DP using flat array."""
    size = n * (1 << n)
    dp = [0] * size
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
                if adj[v * n + u]:
                    dp[u * (1 << n) + (mask | (1 << u))] += val
    full = (1 << n) - 1
    return sum(dp[v * (1 << n) + full] for v in range(n))

def random_tournament_flat(n):
    adj = [0] * (n * n)
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i * n + j] = 1
            else:
                adj[j * n + i] = 1
    return adj

def main():
    n = 8
    N_samples = 50000
    random.seed(12345)

    print(f"Computing Var/Mean^2 at n={n} via Monte Carlo ({N_samples} samples)")
    print(f"Prediction from grand formula: 59/252 = {float(Fraction(59,252)):.10f}")

    start = time.time()
    sum_h = 0
    sum_h2 = 0

    for i in range(N_samples):
        adj = random_tournament_flat(n)
        h = count_ham_paths_fast(adj, n)
        sum_h += h
        sum_h2 += h * h
        if (i+1) % 10000 == 0:
            elapsed = time.time() - start
            mean = sum_h / (i+1)
            mean_h2 = sum_h2 / (i+1)
            var = mean_h2 - mean**2
            ratio = var / mean**2 if mean > 0 else 0
            print(f"  {i+1}/{N_samples} ({elapsed:.1f}s): "
                  f"mean={mean:.2f}, Var/Mean^2={ratio:.6f}")

    elapsed = time.time() - start
    mean = sum_h / N_samples
    mean_h2 = sum_h2 / N_samples
    var = mean_h2 - mean**2
    ratio = var / mean**2

    predicted = float(Fraction(59, 252))

    print(f"\nFinal ({elapsed:.1f}s):")
    print(f"  Mean(H) = {mean:.4f} (expected: {math.factorial(8)/2**7:.4f})")
    print(f"  Var(H)  = {var:.4f}")
    print(f"  Var/Mean^2 = {ratio:.8f}")
    print(f"  Predicted  = {predicted:.8f}")
    print(f"  Difference = {abs(ratio - predicted):.8f}")
    print(f"  Within 1%: {abs(ratio - predicted) / predicted < 0.01}")

if __name__ == '__main__':
    main()
