"""
taylor_cj_mod3_analysis.py
kind-pasteur-2026-03-07-S37

Compute c_j(T) = sum_k C(k,j) F_k(T) for tournaments T at small n.
Check which c_j are always = 0 mod 3.

Key question: at n=7, is c_4(T) = 0 mod 3 for ALL tournaments?
(Palindrome forces c_3=0, and c_5 = c_4 mod 3 from palindrome.)
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import random
from math import comb


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


def compute_cj(F, n):
    """Compute Taylor coefficients c_j = sum_k C(k,j) F_k."""
    c = [0] * n
    for j in range(n):
        for k in range(j, n):
            c[j] += comb(k, j) * F[k]
    return c


def v3(n):
    if n == 0:
        return float('inf')
    v = 0
    n = abs(n)
    while n % 3 == 0:
        v += 1
        n //= 3
    return v


# Exhaustive check for n <= 6
print("=" * 70)
print("Taylor coefficients c_j(T) mod 3: exhaustive check")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    num_T = 1 << m

    always_0_mod3 = [True] * n
    cj_mod3_dist = [[0, 0, 0] for _ in range(n)]
    min_v3_cj = [float('inf')] * n

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c = compute_cj(F, n)

        for j in range(n):
            r = c[j] % 3
            cj_mod3_dist[j][r] += 1
            if r != 0:
                always_0_mod3[j] = False
            vv = v3(c[j])
            if vv < min_v3_cj[j]:
                min_v3_cj[j] = vv

    print(f"\n  n={n} (exhaustive, {num_T} tournaments):")
    for j in range(n):
        pct0 = 100 * cj_mod3_dist[j][0] / num_T
        v3_str = f"min_v3={min_v3_cj[j]}" if min_v3_cj[j] < float('inf') else "all=0"
        status = "ALWAYS 0 mod 3" if always_0_mod3[j] else "NOT always 0"
        print(f"    c_{j}: {pct0:5.1f}% div by 3  {v3_str:12s}  {status}")


# Sampled check at n=7, 8, 9
print("\n" + "=" * 70)
print("Sampled check at n=7, 8, 9")
print("=" * 70)

for n, num_samples in [(7, 10000), (8, 5000), (9, 2000)]:
    m = n * (n - 1) // 2
    num_T = 1 << m
    random.seed(42)

    always_0_mod3 = [True] * n
    cj_mod3_dist = [[0, 0, 0] for _ in range(n)]

    for _ in range(num_samples):
        bits = random.randint(0, num_T - 1)
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c = compute_cj(F, n)

        for j in range(n):
            r = c[j] % 3
            cj_mod3_dist[j][r] += 1
            if r != 0:
                always_0_mod3[j] = False

    print(f"\n  n={n} ({num_samples} samples):")
    for j in range(n):
        pct0 = 100 * cj_mod3_dist[j][0] / num_samples
        status = "ALWAYS 0 mod 3" if always_0_mod3[j] else "NOT always 0"
        print(f"    c_{j}: {pct0:5.1f}% div by 3  {status}")


# Summary: pattern of c_j = 0 mod 3 universality
print("\n" + "=" * 70)
print("Pattern: c_j always = 0 mod 3 for ALL T?")
print("=" * 70)
print("(Y = yes for all tested tournaments, N = no)")
print(f"{'n':>3s}: " + " ".join(f"c_{j}" for j in range(10)))
print("-" * 60)

# From exhaustive data
for n in range(3, 7):
    m = n * (n - 1) // 2
    num_T = 1 << m
    always_0 = [True] * n
    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c = compute_cj(F, n)
        for j in range(n):
            if c[j] % 3 != 0:
                always_0[j] = False

    row = [("Y" if always_0[j] else "N") for j in range(n)]
    row += [" "] * (10 - n)
    print(f"{n:3d}: " + "  ".join(f"{r:>2s}" for r in row))

# From sampled data
for n, num_samples in [(7, 10000), (8, 5000), (9, 2000)]:
    m = n * (n - 1) // 2
    num_T = 1 << m
    random.seed(42)
    always_0 = [True] * n
    for _ in range(num_samples):
        bits = random.randint(0, num_T - 1)
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c = compute_cj(F, n)
        for j in range(n):
            if c[j] % 3 != 0:
                always_0[j] = False

    row = [("Y" if always_0[j] else "N") for j in range(n)]
    row += [" "] * (10 - n)
    print(f"{n:3d}: " + "  ".join(f"{r:>2s}" for r in row))


# Desarmenien factorization: A_n(x) mod 3 = product of A_{d_i+1}(x^{3^i}) mod 3
# where n-1 = sum d_i * 3^i. Check what this predicts for the c_j structure.
print("\n" + "=" * 70)
print("Desarm. factorization: (x-1)-adic valuation of A_n(x) mod 3")
print("=" * 70)

def eulerian_numbers(n):
    if n == 1:
        return [1]
    prev = eulerian_numbers(n - 1)
    A = [0] * n
    for k in range(n):
        A[k] = (k + 1) * prev[k] if k < len(prev) else 0
        if k > 0:
            A[k] += (n - k) * prev[k - 1]
    return A

def base3_digits(n):
    """Return base-3 digits of n, least significant first."""
    if n == 0:
        return [0]
    digits = []
    while n > 0:
        digits.append(n % 3)
        n //= 3
    return digits

for n in range(3, 25):
    A = eulerian_numbers(n)
    # Compute (x-1)-adic valuation of A_n(x) mod 3
    # = smallest j such that c_j(A_n) != 0 mod 3
    c_trans = [0] * n
    for j in range(n):
        for k in range(j, n):
            c_trans[j] += comb(k, j) * A[k]

    val = n  # if all c_j = 0 mod 3
    for j in range(n):
        if c_trans[j] % 3 != 0:
            val = j
            break

    digits = base3_digits(n - 1)
    min_digit = min(d for d in digits if d > 0) if any(d > 0 for d in digits) else 0

    eulerian_zeros = len([k for k in range(n) if A[k] % 3 == 0])
    nonzero_count = n - eulerian_zeros

    print(f"  n={n:2d}: n-1={n-1:3d} (base3={''.join(str(d) for d in reversed(digits))}), "
          f"(x-1)-val={val}, #nonzero={nonzero_count}/{n}, "
          f"digits={digits}")


print("\nDONE")
