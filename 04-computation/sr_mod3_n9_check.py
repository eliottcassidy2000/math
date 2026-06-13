"""
sr_mod3_n9_check.py
kind-pasteur-2026-03-07-S36

Check: does S_r = 0 mod 3 for all r still hold at n=9?
If YES: the Eulerian-number approach is insufficient (all A(9,k) = 1 mod 3).
If NO: mod 9 universality may fail at n=9, despite v_3(9!) = 4.

Also check n=10 (A(10,k) also has no zeros mod 3).

KEY ALGEBRAIC FACT (derived in S36):
  If S_r = 0 mod 3 for all r AND v_3(n!) >= 2, then F(T,omega) = 0 mod 9.
  The proof uses palindrome + sum = n! to show the "quotient" is always = 0 mod 3.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import random
from collections import defaultdict


def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def compute_F_dp(adj, n):
    """Compute F(T,x) coefficients using DP."""
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v][0] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            for fwd in range(n):
                if dp[mask][last][fwd] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    new_mask = mask | (1 << nxt)
                    if adj[last][nxt]:
                        dp[new_mask][nxt][fwd + 1] += dp[mask][last][fwd]
                    else:
                        dp[new_mask][nxt][fwd] += dp[mask][last][fwd]

    full = (1 << n) - 1
    F = [0] * n
    for last in range(n):
        for fwd in range(n):
            F[fwd] += dp[full][last][fwd]
    return F


def check_sr_mod3(n, num_samples=5000):
    m = n * (n - 1) // 2
    num_T = 1 << m

    random.seed(42)

    sr_mod3_dist = defaultdict(int)
    fk_mod3_always0 = [True] * n
    fk_mod3_counts = [[0,0,0] for _ in range(n)]
    failures = 0
    total = 0

    for _ in range(num_samples):
        bits = random.randint(0, num_T - 1)
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)

        S = [0, 0, 0]
        for k in range(n):
            S[k % 3] += F[k]
            r = F[k] % 3
            fk_mod3_counts[k][r] += 1
            if r != 0:
                fk_mod3_always0[k] = False

        key = tuple(s % 3 for s in S)
        sr_mod3_dist[key] += 1
        if any(s % 3 != 0 for s in S):
            failures += 1
        total += 1

    print(f"\nn={n} ({num_samples} samples):")
    print(f"  (S_0, S_1, S_2) mod 3 distribution:")
    for key in sorted(sr_mod3_dist.keys()):
        pct = 100 * sr_mod3_dist[key] / total
        print(f"    {key}: {sr_mod3_dist[key]} ({pct:.1f}%)")
    print(f"  S_r = 0 mod 3 for all r: {total - failures}/{total} "
          f"({'UNIVERSAL' if failures == 0 else f'{failures} FAILURES'})")

    print(f"  Individual F_k always = 0 mod 3:")
    for k in range(n):
        div3_pct = 100 * fk_mod3_counts[k][0] / total
        status = "YES" if fk_mod3_always0[k] else "no"
        print(f"    F_{k}: {status} ({div3_pct:.1f}%)")


# Check n=8 (should work, Eulerian zeros at k=2,5)
check_sr_mod3(8, 5000)

# Check n=9 (CRITICAL: no Eulerian zeros)
check_sr_mod3(9, 3000)

# Check n=10 (no Eulerian zeros)
# n=10 is expensive (36 bits, huge DP). Let's try with fewer samples.
check_sr_mod3(10, 500)

print("\nDONE")
