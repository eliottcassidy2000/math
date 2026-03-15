#!/usr/bin/env python3
"""
pi_spectral_H_89c.py — Can we express H(T) spectrally?
opus-2026-03-14-S89c

For Paley tournaments, all eigenvalues are ±i√p.
H(P_p) has been computed for p=3,7,11,19,23.
Can H be expressed as a function of the spectrum?

Also investigate: the PERMANENT of the adjacency matrix A (where A[i][j]=1 if i→j, else 0).
The permanent is ≥ H since it counts all not-necessarily-simple directed walks.
Actually perm(A) counts perfect matchings in the bipartite graph, not H.
We know perm correlates with H (0.936 at n=6).
"""

from itertools import permutations
from math import factorial, sqrt, pi, e as euler_e
from fractions import Fraction

def legendre(a, p):
    if a % p == 0: return 0
    r = pow(a, (p-1)//2, p)
    return -1 if r == p-1 else r

def paley_adj(p):
    """0-1 adjacency matrix: A[i][j]=1 if i→j"""
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and legendre(j-i, p) == 1:
                A[i][j] = 1
    return A

def count_hp_dp(adj, n):
    """Count Hamiltonian paths via DP"""
    dp = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in dp[mask]:
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u) == 0 and adj[v][u] == 1:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full].values())

def permanent(M, n):
    """Permanent via inclusion-exclusion (Ryser formula)"""
    total = 0
    for S in range(1 << n):
        # S encodes subset of columns
        bits = bin(S).count('1')
        sign = (-1) ** (n - bits)
        # Product of row sums restricted to S
        prod = 1
        for i in range(n):
            row_sum = 0
            for j in range(n):
                if S & (1 << j):
                    row_sum += M[i][j]
            prod *= row_sum
        total += sign * prod
    return total

print("=" * 70)
print("SPECTRAL AND ALGEBRAIC ANALYSIS OF H(P_p)")
print("=" * 70)

# Known H values
H_data = {3: 4, 7: 189, 11: 95095, 19: 119409657765, 23: 15760206976379349}

for p in [3, 7, 11, 19, 23]:
    if p > 19:  # skip expensive computation
        H = H_data.get(p, None)
        if H is None:
            continue
        print(f"\n--- p = {p} (using known H) ---")
    else:
        A = paley_adj(p)
        H = count_hp_dp(A, p)
        print(f"\n--- p = {p} ---")
        print(f"  H(P_{p}) = {H}")

    # H = p × (p-1)/2 × (a+b) decomposition
    # Known: H/p is always an integer (THM-212)
    print(f"  H/p = {H // p}")
    half = (p-1) // 2
    print(f"  H/(p×(p-1)/2) = {H // (p * half)}")

    # E[H] = n!/2^{n-1}
    EH = factorial(p) / 2**(p-1)
    ratio = H / EH
    print(f"  E[H] = {EH:.2f}")
    print(f"  H/E[H] = {ratio:.6f}")

    # Spectral data: eigenvalues are 0, +i√p (×(p-1)/2), -i√p (×(p-1)/2)
    # The characteristic polynomial of S is:
    # x × (x² + p)^{(p-1)/2}
    #
    # Can H be expressed in terms of p alone?
    # H(P_3) = 4, H(P_7) = 189, H(P_11) = 95095, H(P_19) = 119409657765

    # Let's look at H/p!
    # 4/6 = 2/3
    # 189/5040 = 189/5040 = 27/720 = 3/80
    # 95095/39916800 = ...
    print(f"  H/p! = {Fraction(H, factorial(p))}")

    # H × 2^{p-1} / p!
    val = Fraction(H * 2**(p-1), factorial(p))
    print(f"  H × 2^(p-1) / p! = {val} = {float(val):.8f}")

    # What about H / (p × Subfactorial(p))?
    # Subfactorial: D(n) = n! × Σ (-1)^k/k!
    D_p = sum((-1)**k * factorial(p) // factorial(k) for k in range(p+1))
    print(f"  D({p}) = {D_p}")
    if D_p > 0:
        print(f"  H / D({p}) = {H / D_p:.8f}")

print()
print("=" * 70)
print("PERMANENT vs H FOR PALEY TOURNAMENTS")
print("=" * 70)

for p in [3, 7, 11]:
    A = paley_adj(p)
    H = count_hp_dp(A, p)
    perm_A = permanent(A, p)
    print(f"\n  p={p}: H = {H}, perm(A) = {perm_A}")
    print(f"    perm/H = {perm_A / H:.6f}")
    print(f"    perm/p! = {Fraction(perm_A, factorial(p))}")

    # The permanent counts: number of ways to assign each row to a distinct column
    # where A[i][σ(i)]=1 for all i. This is the number of permutations σ such that
    # i → σ(i) in the tournament. These are the permutation matrices contained in A.
    # This is related to the number of directed perfect matchings.

    # H counts Hamiltonian paths = orderings π with A[π(i), π(i+1)]=1 for all i.
    # Permanent counts permutation matrices in A.
    # They're related but different.

print()
print("=" * 70)
print("H(P_p) FACTORED AND ANALYZED")
print("=" * 70)

# Factor H values
def factor(n):
    if n <= 1: return [(n, 1)]
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors

for p in [3, 7, 11, 19, 23]:
    H = H_data[p]
    f = factor(H)
    print(f"  H(P_{p}) = {H} = {'×'.join(map(str, f))}")

    # Check divisibility by p, (p-1)/2, p!
    print(f"    div by p={p}: {H % p == 0}")
    print(f"    div by (p-1)/2={(p-1)//2}: {H % ((p-1)//2) == 0}")

print()
print("=" * 70)
print("THE H/E[H] SEQUENCE AND ITS LIMIT")
print("=" * 70)

# H(P_p)/E[H_p] for Paley:
print("\n  Ratios H(P_p)/E[H]:")
for p in [3, 7, 11, 19, 23]:
    H = H_data[p]
    EH = factorial(p) / 2**(p-1)
    r = H / EH
    # Compare with various constants
    log_r = __import__('math').log(r)
    print(f"  p={p:2d}: H/E[H] = {r:.8f}, ln(H/E[H]) = {log_r:.6f}")

# The ratios: 2.667, 2.4, 2.441, 2.525, 2.559
# These are growing slowly. What's the limit?
# Previous analysis showed gap ≈ 0.645/p^0.756
# So H/E[H] → e? → √(2π)? → something else?

import math

ratios = []
for p in [3, 7, 11, 19, 23]:
    H = H_data[p]
    EH = factorial(p) / 2**(p-1)
    ratios.append(H / EH)

print(f"\n  Candidate limits:")
print(f"    e = {euler_e:.6f}")
print(f"    √(2π) = {sqrt(2*pi):.6f}")
print(f"    π = {pi:.6f}")
print(f"    Current max ratio (p=23): {ratios[-1]:.6f}")

# By Stirling: E[H] = n!/2^{n-1} ≈ √(2πn) × (n/2e)^n
# If H(P_p) ≈ C × E[H], what is C?
# Actually for large p, the Paley tournament is "almost regular"
# and H should be close to the mean, so C → 1??
# But the ratio is > 2 and slowly growing...

# Wait, let me reconsider. E[H] is the average over ALL tournaments.
# The Paley tournament is special. Maybe H(P_p) grows faster than E[H].

# Let me check: H(P_p) / p × (permanent of Paley) ?
# The permanent might be the right normalization.

print()
print("=" * 70)
print("STIRLING ANALYSIS OF H(P_p)")
print("=" * 70)

for p in [3, 7, 11, 19, 23]:
    H = H_data[p]
    # Stirling: n! ≈ √(2πn) (n/e)^n
    stirling = sqrt(2*pi*p) * (p/euler_e)**p
    # E[H] = p!/2^{p-1}
    EH_exact = factorial(p) / 2**(p-1)
    # H/stirling
    print(f"  p={p}: log₂(H) = {math.log2(H):.2f}, log₂(p!) = {math.log2(factorial(p)):.2f}")
    print(f"    H × 2^(p-1) / p! = {H * 2**(p-1) / factorial(p):.8f}")
    # This ratio is H/E[H] which we already know

# The key question: what is the combinatorial significance of H(P_p)?
# H = p × (p-1)/2 × (a+b) where a,b are endpoint-class counts.
# For p=11: a+b = 1729 = 7×13×19 (taxicab!)
# For p=7: a+b = 27 = 3³

ab_vals = {3: H_data[3] // (3 * 1), 7: H_data[7] // (7 * 3),
           11: H_data[11] // (11 * 5), 19: H_data[19] // (19 * 9),
           23: H_data[23] // (23 * 11)}

print(f"\n  a+b values (H = p × (p-1)/2 × (a+b)):")
for p in [3, 7, 11, 19, 23]:
    ab = ab_vals[p]
    f = factor(ab)
    print(f"    p={p}: a+b = {ab} = {'×'.join(map(str, f))}")
    # Divide by (p-2)!
    if p <= 11:
        r = Fraction(ab, factorial(p-2))
        print(f"      (a+b)/(p-2)! = {r} = {float(r):.8f}")

print("\nDone!")
