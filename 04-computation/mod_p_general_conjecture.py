"""
mod_p_general_conjecture.py
kind-pasteur-2026-03-07-S38

General conjecture: For prime p and n >= p+2,
  p | c_j(T) for all j < val_p(n) and all tournaments T.

Verify at:
- p=5, n=7,8: should match val_5(7)=4, val_5(8)=4
- p=7, n=9: should match val_7(9)=6
- p=2, n=8: should match val_2(8)=7
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
    A = eulerian_numbers(n)
    for j in range(n):
        cj = sum(comb(k, j) * A[k] for k in range(n))
        if cj % p != 0:
            return j
    return n


def test_universal_zeros(n, p, num_samples=2000):
    m_bits = n * (n - 1) // 2
    num_T = 1 << m_bits
    random.seed(42)

    samples = [random.randint(0, num_T - 1) for _ in range(num_samples)]
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

    return threshold, failures


# ============================================================
# Test cases
# ============================================================
print("=" * 70)
print("General mod-p conjecture verification")
print("=" * 70)

tests = [
    (2, 8, 1000, "THM-094 extension"),
    (5, 7, 2000, "p=5, n=p+2"),
    (5, 8, 1000, "p=5, n=p+3"),
]

for p, n, ns, desc in tests:
    eul_val = eulerian_val(n, p)
    threshold, failures = test_universal_zeros(n, p, ns)
    match = "MATCH" if threshold == eul_val else f"MISMATCH (threshold={threshold})"

    print(f"\n  p={p}, n={n} ({desc}, {ns} samples):")
    print(f"    Eulerian val = {eul_val}, universal zeros = {threshold}  {match}")
    for j in range(min(n, max(threshold, eul_val) + 2)):
        status = "ALWAYS 0" if failures[j] == 0 else f"{failures[j]} fail ({100*failures[j]/ns:.1f}%)"
        print(f"    c_{j}: {status}")


# ============================================================
# p=2 at n=8: verify F mod 2 = A mod 2 (individual F_k)
# ============================================================
print("\n" + "=" * 70)
print("p=2, n=8: individual F_k mod 2")
print("=" * 70)

n = 8
A = eulerian_numbers(n)
A_mod2 = [a % 2 for a in A]
print(f"  A(8,k) mod 2 = {A_mod2}")
print(f"  C(7,k) mod 2 = {[comb(7,k) % 2 for k in range(8)]}")

m = n * (n - 1) // 2
random.seed(42)
num_samples = 1000
failures = 0

for _ in range(num_samples):
    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_from_bits(n, bits)
    F = compute_F_dp(adj, n)
    F_mod2 = [f % 2 for f in F]
    if F_mod2 != A_mod2:
        failures += 1

print(f"  {num_samples} samples: {failures} failures")
if failures == 0:
    print(f"  F_k(T) = A(8,k) mod 2 for all T: CONFIRMED")


# ============================================================
# For p=5 at n=7: check individual F_k mod 5 vs A(7,k) mod 5
# ============================================================
print("\n" + "=" * 70)
print("p=5, n=7: individual F_k mod 5 vs Eulerian")
print("=" * 70)

n = 7
p = 5
A = eulerian_numbers(n)
A_mod5 = [a % 5 for a in A]
print(f"  A(7,k) mod 5 = {A_mod5}")

# Find positions where A(7,k) = 0 mod 5
euler_zeros = [k for k in range(n) if A[k] % 5 == 0]
print(f"  Eulerian zeros mod 5 at k = {euler_zeros}")

m = n * (n - 1) // 2
random.seed(42)
num_samples = 3000
fk_failures = [0] * n

for _ in range(num_samples):
    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_from_bits(n, bits)
    F = compute_F_dp(adj, n)
    for k in euler_zeros:
        if F[k] % 5 != 0:
            fk_failures[k] += 1

for k in euler_zeros:
    status = "ALWAYS 0 mod 5" if fk_failures[k] == 0 else f"{fk_failures[k]} fail"
    print(f"  F_{k}: {status}")

if euler_zeros and all(fk_failures[k] == 0 for k in euler_zeros):
    print("  Eulerian conjecture mod 5 at n=7: HOLDS")


print("\nDONE")
