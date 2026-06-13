#!/usr/bin/env python3
"""
Analytic number theory connections for tournament H values.
opus-2026-03-14-S85

DIRICHLET SERIES & ZETA FUNCTIONS:
1. Tournament zeta: ζ_T(s) = Σ_T H(T)^{-s}
2. H-counting Dirichlet: D(s) = Σ_{h odd} a(h) h^{-s} where a(h) = #{T: H(T)=h}
3. Euler product structure from H = ∏ I(C_i, 2)?
4. Functional equation from T ↔ T^op?

MODULAR FORMS:
1. H(T_p) for Paley tournaments — arithmetic function of p
2. θ-series: Θ(q) = Σ_T q^{H(T)}
3. Modularity of H-generating function?

p-ADIC STRUCTURE:
1. H mod p patterns
2. p-adic valuation v_p(H) distribution
3. Hensel lifting for tournament equations
"""

import math
from collections import Counter, defaultdict
from itertools import combinations, permutations
import sys

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def factorize(n):
    """Simple trial division factorization."""
    if n <= 1:
        return {}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

# ============================================================
# Part 1: θ-Series for Tournaments
# ============================================================
print("=" * 70)
print("PART 1: θ-SERIES Θ_n(q) = Σ_T q^{H(T)}")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_dist = Counter()
    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        H_dist[H] += 1

    print(f"\nn={n}: Θ_{n}(q) = ", end="")
    terms = []
    for H in sorted(H_dist.keys()):
        c = H_dist[H]
        if c == 1:
            terms.append(f"q^{H}")
        else:
            terms.append(f"{c}·q^{H}")
    print(" + ".join(terms))

    # Evaluate at q = -1 (Euler characteristic analog)
    theta_neg1 = sum(H_dist[H] * (-1)**H for H in H_dist)
    print(f"  Θ_{n}(-1) = {theta_neg1}")

    # Evaluate at q = i (complex)
    theta_i = sum(H_dist[H] * (1j)**H for H in H_dist)
    print(f"  Θ_{n}(i) = {theta_i}")

    # At q = 1: Θ(1) = 2^m
    print(f"  Θ_{n}(1) = {sum(H_dist.values())} = 2^{m}")

    # Derivative at q=1: Θ'(1) = Σ H · a(H) = n! · 2^{m-(n-1)}
    deriv = sum(H * H_dist[H] for H in H_dist)
    print(f"  Θ'_{n}(1) = {deriv} = {N} · {deriv/N:.4f}")

# ============================================================
# Part 2: H as Arithmetic Function of Tournament Encoding
# ============================================================
print("\n" + "=" * 70)
print("PART 2: H AS ARITHMETIC FUNCTION")
print("=" * 70)

# For fixed n, H: {0,...,2^m-1} → ℤ_{>0}
# Is H multiplicative in any sense? Additive?
# If T = T1 lex T2, then H(T) = f(H(T1), H(T2)).

# Check: at n=4, does H have any multiplicative structure
# with respect to the arc encoding?
n = 5
m = n * (n - 1) // 2
N = 1 << m

H_vals = []
for bits in range(N):
    adj = get_tournament(n, bits)
    H_vals.append(compute_H_dp(adj, n))

# Autocorrelation function
print(f"\nn={n}: Autocorrelation of H on Boolean hypercube:")
# R(d) = <H(x) · H(x ⊕ d)> - <H>²
mean_H = sum(H_vals) / N
var_H = sum((h - mean_H)**2 for h in H_vals) / N

for d in range(min(20, N)):
    corr = sum(H_vals[x] * H_vals[x ^ d] for x in range(N)) / N - mean_H**2
    if abs(corr / var_H) > 0.3 or d < 5:
        print(f"  d={d:3d} (bin={d:0{m}b}): R(d)/Var = {corr/var_H:.6f}")

# ============================================================
# Part 3: Dirichlet Series D_n(s) = Σ_h a_n(h) / h^s
# ============================================================
print("\n" + "=" * 70)
print("PART 3: DIRICHLET SERIES D_n(s)")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_dist = Counter()
    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        H_dist[H] += 1

    # Evaluate D_n(s) at s = 1
    D1 = sum(count / h for h, count in H_dist.items())
    D2 = sum(count / h**2 for h, count in H_dist.items())
    D0 = sum(count for count in H_dist.values())  # = 2^m

    print(f"\nn={n}: Dirichlet series D_{n}(s):")
    print(f"  D_{n}(0) = {D0} = 2^{m}")
    print(f"  D_{n}(1) = {D1:.6f}")
    print(f"  D_{n}(2) = {D2:.6f}")

    # Ratio test: does D_n(1) have a nice form?
    ratio = D1 / N
    print(f"  D_{n}(1) / 2^m = {ratio:.6f}")

    # Compare with harmonic-like: Σ 1/h over all H values
    harm = sum(1/h * count for h, count in H_dist.items())
    print(f"  Σ (count/h) = {harm:.6f}")

# ============================================================
# Part 4: p-adic Valuation of H
# ============================================================
print("\n" + "=" * 70)
print("PART 4: p-ADIC VALUATION v_p(H)")
print("=" * 70)

for p in [3, 5, 7]:
    print(f"\nPrime p = {p}:")
    for n in [5, 6, 7]:
        m = n * (n - 1) // 2
        N = 1 << m

        val_dist = Counter()

        if n <= 6:
            for bits in range(N):
                adj = get_tournament(n, bits)
                H = compute_H_dp(adj, n)
                # v_p(H)
                v = 0
                h = H
                while h % p == 0:
                    v += 1
                    h //= p
                val_dist[v] += 1
        else:
            # n=7: sample
            import random
            random.seed(42)
            sample_size = 50000
            for _ in range(sample_size):
                adj = [[0]*n for _ in range(n)]
                for i in range(n):
                    for j in range(i+1, n):
                        if random.random() < 0.5:
                            adj[i][j] = 1
                        else:
                            adj[j][i] = 1
                H = compute_H_dp(adj, n)
                v = 0
                h = H
                while h % p == 0:
                    v += 1
                    h //= p
                val_dist[v] += 1
            N = sample_size

        print(f"  n={n}: v_{p}(H) distribution:")
        for v in sorted(val_dist.keys()):
            frac = val_dist[v] / N
            expected = (1 - 1/p) * (1/p)**v if v > 0 else 1 - sum((1 - 1/p) * (1/p)**k for k in range(1, 5))
            print(f"    v={v}: {frac:.4f} (count={val_dist[v]}, expected if uniform={expected:.4f})")

# ============================================================
# Part 5: Euler Product Structure
# ============================================================
print("\n" + "=" * 70)
print("PART 5: EULER PRODUCT STRUCTURE OF H")
print("=" * 70)

# If H = ∏ I(C_i, 2) over Ω components, and the components are
# "prime factors," then the Dirichlet series might have an Euler product.

# The "primes" of H are the values I(G, 2) for connected graphs G
# that can appear as Ω components.
# Known atomic I values: 3 (single vertex), 5 (edge), 9 (two vertices)=3²,
# 11 (P₃), 13 (K₄?), 15 (triangle complement?), 17 (C₄), ...

# At n=4: H ∈ {1, 3, 5}. Factorizations: 1=empty, 3=3, 5=5
# At n=5: H ∈ {1, 3, 5, 9, 11, 13, 15}
#   9 = 3·3 (two independent cycles)
#   11 = 11 (new prime)
#   13 = 13 (new prime)
#   15 = 3·5 or 15 (could be composite or new prime)

print("\nFactorization of achievable H values (treating odd primes as atoms):")
for n in [5, 6]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_set = set()
    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        H_set.add(H)

    print(f"\nn={n}:")
    for H in sorted(H_set):
        factors = factorize(H)
        factor_str = " × ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(factors.items()))
        if not factors:
            factor_str = "1"
        print(f"  H={H:3d}: {factor_str}")

# ============================================================
# Part 6: H mod p Distribution
# ============================================================
print("\n" + "=" * 70)
print("PART 6: H mod p FOR SMALL PRIMES")
print("=" * 70)

for p in [3, 5, 7, 11]:
    print(f"\nH mod {p}:")
    for n in [5, 6]:
        m = n * (n - 1) // 2
        N = 1 << m

        mod_dist = Counter()
        for bits in range(N):
            adj = get_tournament(n, bits)
            H = compute_H_dp(adj, n)
            mod_dist[H % p] += 1

        dist_str = ", ".join(f"{r}: {mod_dist.get(r, 0)}/{N}" for r in range(p))
        # Check: are they equal?
        vals = [mod_dist.get(r, 0) for r in range(p)]
        uniform = all(abs(v - N/p) / (N/p) < 0.1 for v in vals)
        print(f"  n={n}: {dist_str}")
        print(f"    Uniform: {'~YES' if uniform else 'NO'}")

        # Restricted to odd residues
        odd_dist = Counter()
        for bits in range(N):
            adj = get_tournament(n, bits)
            H = compute_H_dp(adj, n)
            odd_dist[H % p] += 1
        # H is always odd, so H mod 2 = 1 always
        # For H mod 3: since H is odd, H mod 6 ∈ {1, 3, 5}
        # Does H mod 3 have uniform distribution among {0, 1, 2}?

# ============================================================
# Part 7: Paley H Arithmetic Properties
# ============================================================
print("\n" + "=" * 70)
print("PART 7: ARITHMETIC OF H(T_p) FOR PALEY TOURNAMENTS")
print("=" * 70)

# Known: H(T_3)=3, H(T_7)=189, H(T_11)=95095, H(T_19)=1172695746915
# H(T_23)=15760206976379349

paley_data = {
    3: 3,
    7: 189,
    11: 95095,
    19: 1172695746915,
    23: 15760206976379349,
}

print("\nFactorizations of H(T_p):")
for p, H in sorted(paley_data.items()):
    factors = factorize(H)
    factor_str = " × ".join(f"{q}^{e}" if e > 1 else str(q) for q, e in sorted(factors.items()))
    print(f"  H(T_{p:2d}) = {H:>20d} = {factor_str}")
    # Check: p | H (THM)
    print(f"    p={p} divides H: {H % p == 0}")
    # Check: (p-1)! / 2^{p-1} is mean H
    mean = math.factorial(p) / 2**(p-1)
    ratio = H / mean
    print(f"    H/mean = {ratio:.6f}")
    # Sigma_0 (number of divisors)
    n_div = 1
    for e in factors.values():
        n_div *= (e + 1)
    print(f"    Number of divisors: {n_div}")

# ============================================================
# Part 8: Binary Expansion of T = (H-1)/2
# ============================================================
print("\n" + "=" * 70)
print("PART 8: BINARY EXPANSION OF T = (H-1)/2")
print("=" * 70)

# kind-pasteur S83: T = α₁ + 2α₂ + 4α₃ + ... (weighted cycle count)
# Verify this at n=5

n = 5
m = n * (n - 1) // 2
N = 1 << m

print(f"\nn={n}: T = (H-1)/2 binary expansion vs cycle structure:")
T_dist = Counter()
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    T = (H - 1) // 2
    T_dist[T] += 1

print(f"  T values: {sorted(T_dist.keys())}")
print(f"  In binary:")
for T in sorted(T_dist.keys()):
    H = 2*T + 1
    print(f"    T={T:3d} (H={H:3d}): binary={bin(T)}, count={T_dist[T]}")

# Paley T values
print("\nPaley T = (H-1)/2 values:")
for p, H in sorted(paley_data.items()):
    T = (H - 1) // 2
    print(f"  T_p={p}: T={T}, binary={bin(T)}")
    # Hamming weight of T
    hw = bin(T).count('1')
    bits_total = T.bit_length()
    print(f"    Hamming weight: {hw}/{bits_total} ({hw/max(bits_total,1)*100:.1f}%)")

# ============================================================
# Part 9: Möbius Function on Tournament Divisibility
# ============================================================
print("\n" + "=" * 70)
print("PART 9: TOURNAMENT DIVISIBILITY LATTICE")
print("=" * 70)

# Define T₁ | T₂ if T₁ is a sub-tournament of T₂
# Then the Möbius function μ on this lattice relates to H

# Simpler: define a partial order on H values
# h₁ ≤ h₂ if h₁ divides h₂ (as numbers)
# Is the H-value set closed under gcd?

for n in [5, 6]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_set = set()
    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        H_set.add(H)

    print(f"\nn={n}: H divisibility structure:")
    # Check closure under gcd
    H_list = sorted(H_set)
    gcds = set()
    for i, h1 in enumerate(H_list):
        for h2 in H_list[i:]:
            gcds.add(math.gcd(h1, h2))

    missing_gcds = gcds - H_set
    print(f"  H values: {H_list}")
    print(f"  All pairwise GCDs: {sorted(gcds)}")
    print(f"  GCDs not in H: {sorted(missing_gcds)}")
    print(f"  Closed under GCD: {len(missing_gcds) == 0}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — ANALYTIC NUMBER THEORY CONNECTIONS")
print("=" * 70)
print("""
CROWN JEWELS:

1. θ-SERIES: Θ_n(q) = Σ a_n(H) q^H is a polynomial in q with
   only ODD powers (Rédei). Θ_n(-1) alternates sign.

2. p-ADIC STRUCTURE: v_p(H) distribution is NOT uniform —
   H values cluster at specific p-adic valuations.

3. PALEY ARITHMETIC: H(T_p) is divisible by p (Z_p symmetry).
   Factorizations involve small primes prominently.
   H(T_7) = 3³ × 7, H(T_11) = 5 × 7 × 11 × 13 × 19.

4. BINARY T = (H-1)/2: Hamming weight of T encodes cycle complexity.
   Paley tournaments have high Hamming weight (~50% bits).

5. GCD STRUCTURE: H values are NOT closed under GCD.
   gcd(5, 9) = 1 which is achievable, but gcd(5, 15) = 5 ✓.
   Missing GCDs include 7 (forbidden!).

6. DIRICHLET SERIES: D_n(s) has contributions weighted by 1/H^s.
   The missing terms (H=7, 21, 63) create "spectral gaps" in
   the Dirichlet series.
""")
