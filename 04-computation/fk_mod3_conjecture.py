"""
fk_mod3_conjecture.py
kind-pasteur-2026-03-07-S36

CONJECTURE: 3 | F_k(T) for all tournaments T iff 3 | A(n,k) (Eulerian number).

Equivalently: F_k(T) = A(n,k) mod 3 for all T? (stronger)
Or just: 3 | A(n,k) => 3 | F_k(T) for all T? (weaker, but sufficient)

Eulerian numbers mod 3:
n=3: 1,1,1 (none divisible by 3)
n=4: 1,2,2,1 (none)
n=5: 1,2,0,2,1 (k=2)
n=6: 1,0,2,2,0,1 (k=1,4)
n=7: 1,0,0,1,0,0,1 (k=1,2,4,5)
n=8: 1,?,?,?,?,?,?,1

This script:
1. Compute Eulerian numbers mod 3 for n up to 12
2. Verify the conjecture at n=8 (sampled)
3. Study the structure of which k have A(n,k) = 0 mod 3
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import random
from collections import defaultdict


def eulerian_numbers(n):
    """Compute Eulerian numbers A(n,k) for k=0,...,n-1."""
    if n == 1:
        return [1]
    # Recurrence: A(n,k) = (k+1)*A(n-1,k) + (n-k)*A(n-1,k-1)
    prev = eulerian_numbers(n - 1)
    A = [0] * n
    for k in range(n):
        A[k] = (k + 1) * prev[k] if k < len(prev) else 0
        if k > 0:
            A[k] += (n - k) * prev[k - 1]
    return A


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


# Step 1: Eulerian numbers mod 3
print("=" * 60)
print("Eulerian numbers A(n,k) mod 3")
print("=" * 60)
for n in range(2, 13):
    A = eulerian_numbers(n)
    mod3 = [a % 3 for a in A]
    zero_positions = [k for k in range(n) if mod3[k] == 0]
    print(f"  n={n:2d}: {mod3}  zeros at k={zero_positions}")


# Step 2: Verify at n=5,6 (exhaustive) and n=7,8 (sampled)
print("\n" + "=" * 60)
print("Verification of conjecture: 3|A(n,k) => 3|F_k(T) for all T")
print("=" * 60)

for n in [5, 6]:
    A = eulerian_numbers(n)
    zero_k = [k for k in range(n) if A[k] % 3 == 0]
    m = n * (n - 1) // 2
    num_T = 1 << m

    failures = {k: 0 for k in zero_k}
    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        for k in zero_k:
            if F[k] % 3 != 0:
                failures[k] += 1

    print(f"\n  n={n} (exhaustive, {num_T} tournaments):")
    print(f"    A(n,k) mod 3 zeros at k={zero_k}")
    for k in zero_k:
        print(f"    F_{k} always = 0 mod 3: {'YES' if failures[k] == 0 else f'NO ({failures[k]} failures)'}")

for n, num_samples in [(7, 10000), (8, 5000)]:
    A = eulerian_numbers(n)
    zero_k = [k for k in range(n) if A[k] % 3 == 0]
    m = n * (n - 1) // 2
    num_T = 1 << m

    random.seed(42)
    failures = {k: 0 for k in zero_k}
    for _ in range(num_samples):
        bits = random.randint(0, num_T - 1)
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        for k in zero_k:
            if F[k] % 3 != 0:
                failures[k] += 1

    print(f"\n  n={n} ({num_samples} samples):")
    print(f"    A(n,k) mod 3 zeros at k={zero_k}")
    for k in zero_k:
        print(f"    F_{k} always = 0 mod 3: {'YES (0 failures)' if failures[k] == 0 else f'NO ({failures[k]} failures)'}")


# Step 3: Stronger conjecture? F_k(T) = A(n,k) mod 3 for ALL T and ALL k?
print("\n" + "=" * 60)
print("Stronger test: F_k(T) = A(n,k) mod 3 for ALL T?")
print("=" * 60)

for n in [5, 6]:
    A = eulerian_numbers(n)
    m = n * (n - 1) // 2
    num_T = 1 << m

    # For each k, check if F_k(T) mod 3 is always = A(n,k) mod 3
    all_match = [True] * n
    Fk_residues = [set() for _ in range(n)]
    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        for k in range(n):
            Fk_residues[k].add(F[k] % 3)
            if F[k] % 3 != A[k] % 3:
                all_match[k] = False

    print(f"\n  n={n}:")
    for k in range(n):
        print(f"    F_{k}: A(n,k) mod 3 = {A[k]%3}, "
              f"observed residues = {sorted(Fk_residues[k])}, "
              f"always = A(n,k) mod 3: {all_match[k]}")


# Step 4: Pattern in Eulerian number zeros mod 3
print("\n" + "=" * 60)
print("Pattern: which k have A(n,k) = 0 mod 3?")
print("=" * 60)

for n in range(2, 16):
    A = eulerian_numbers(n)
    mod3 = [a % 3 for a in A]
    zero_k = [k for k in range(n) if mod3[k] == 0]
    nonzero_k = [k for k in range(n) if mod3[k] != 0]
    num_zeros = len(zero_k)
    print(f"  n={n:2d}: {num_zeros}/{n} zeros. "
          f"Nonzero at k={nonzero_k} (mod 3 residues: {[mod3[k] for k in nonzero_k]})")


# Step 5: Check if the mod-3 Eulerian pattern matches F_k for random n=9
print("\n" + "=" * 60)
print("n=9 verification (sampled)")
print("=" * 60)

n = 9
A = eulerian_numbers(n)
zero_k = [k for k in range(n) if A[k] % 3 == 0]
m = n * (n - 1) // 2
num_T = 1 << m

print(f"  A(9,k) mod 3: {[a%3 for a in A]}")
print(f"  Zeros at k={zero_k}")

random.seed(42)
num_samples = 2000
failures = {k: 0 for k in zero_k}
for _ in range(num_samples):
    bits = random.randint(0, num_T - 1)
    adj = tournament_from_bits(n, bits)
    F = compute_F_dp(adj, n)
    for k in zero_k:
        if F[k] % 3 != 0:
            failures[k] += 1

for k in zero_k:
    print(f"  F_{k}: {'PASS' if failures[k] == 0 else f'FAIL ({failures[k]})'} ({num_samples} samples)")


print("\nDONE")
