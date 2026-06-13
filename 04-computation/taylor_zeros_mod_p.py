"""
taylor_zeros_mod_p.py
kind-pasteur-2026-03-07-S38

Investigate THM-086 generalization: for prime p, what is the (x-1)-adic
valuation of F(T,x) mod p? Does it match the Eulerian valuation?

For mod 3: val(n) = 2*floor((n-1)/2), matching A_n(x) mod 3.
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


def eulerian_val(n, p):
    """(x-1)-adic valuation of A_n(x) mod p."""
    A = eulerian_numbers(n)
    for j in range(n):
        cj = sum(comb(k, j) * A[k] for k in range(n))
        if cj % p != 0:
            return j
    return n


# ============================================================
# Part 1: Eulerian (x-1)-adic valuations
# ============================================================
print("=" * 70)
print("(x-1)-adic valuation of A_n(x) mod p")
print("=" * 70)

for p in [2, 3, 5, 7]:
    vals = []
    print(f"\n  p={p}:")
    for n in range(2, 20):
        v = eulerian_val(n, p)
        vals.append(v)
        print(f"    n={n:2d}: val = {v}")


# ============================================================
# Part 2: Universal Taylor zeros of F(T,x) mod p (n<=7 only)
# ============================================================
print("\n" + "=" * 70)
print("Universal Taylor zeros of F(T,x) mod p")
print("=" * 70)

for p in [2, 5, 7]:
    print(f"\n  === p = {p} ===")

    for n in range(3, 8):
        m_bits = n * (n - 1) // 2
        num_T = 1 << m_bits
        random.seed(42)

        if n <= 6:
            samples = list(range(num_T))
            method = "exhaustive"
        else:
            samples = [random.randint(0, num_T - 1) for _ in range(3000)]
            method = "3000 samples"

        num_samples = len(samples)
        failures = [0] * n

        for bits in samples:
            adj = tournament_from_bits(n, bits)
            F = compute_F_dp(adj, n)
            c = [sum(comb(k, j) * F[k] for k in range(n)) for j in range(n)]
            for j in range(n):
                if c[j] % p != 0:
                    failures[j] += 1

        threshold = 0
        for j in range(n):
            if failures[j] == 0:
                threshold = j + 1
            else:
                break

        eul_val = eulerian_val(n, p)
        match = "MATCH" if threshold == eul_val else f"MISMATCH (eul={eul_val})"

        print(f"    n={n} ({method}): universal zeros = {threshold}, "
              f"Eulerian val = {eul_val}  {match}")

        # Show detail for mismatches
        if threshold != eul_val:
            for j in range(min(n, max(threshold, eul_val) + 2)):
                status = "ALWAYS 0" if failures[j] == 0 else f"{failures[j]} fail ({100*failures[j]/num_samples:.1f}%)"
                print(f"      c_{j}: {status}")


# ============================================================
# Part 3: Mod 2 detail (F(T,x) mod 2 patterns)
# ============================================================
print("\n" + "=" * 70)
print("Mod 2: F(T,x) mod 2 patterns")
print("=" * 70)

for n in range(3, 7):
    m_bits = n * (n - 1) // 2
    num_T = 1 << m_bits
    patterns = {}
    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        F_mod2 = tuple(f % 2 for f in F)
        patterns[F_mod2] = patterns.get(F_mod2, 0) + 1

    print(f"\n  n={n}: {len(patterns)} distinct F mod 2 out of {num_T}")
    for pat, count in sorted(patterns.items(), key=lambda x: -x[1])[:6]:
        print(f"    {pat}: {count} ({100*count/num_T:.1f}%)")


# ============================================================
# Part 4: Palindrome and universal zeros interaction
# ============================================================
print("\n" + "=" * 70)
print("For p=2: what forces c_j = 0 mod 2?")
print("=" * 70)

# c_0 = n!, c_1 = n!(n-1)/2
# For p=2: v_2(n!) grows fast. c_0 = n! always even for n>=2.
# c_1 = n!(n-1)/2: for n>=2, this is always an integer. Is it always even?
# n=3: c_1 = 6*1 = 6, even. n=4: c_1 = 24*3/2 = 36, even.
# Actually c_1 = n!(n-1)/2. v_2(c_1) = v_2(n!) + v_2(n-1) - 1.

from math import factorial

for n in range(3, 12):
    c0 = factorial(n)
    c1 = factorial(n) * (n-1) // 2
    v2_c0 = 0
    t = c0
    while t % 2 == 0:
        v2_c0 += 1; t //= 2
    v2_c1 = 0
    t = c1
    while t % 2 == 0:
        v2_c1 += 1; t //= 2
    print(f"  n={n}: v_2(c_0)={v2_c0}, v_2(c_1)={v2_c1}")


print("\nDONE")
