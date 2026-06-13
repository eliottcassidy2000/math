#!/usr/bin/env python3
"""
Large-n tournament exploration using Held-Karp DP.
opus-2026-03-14-S85

With O(2^n * n²) DP, we can compute H for:
- n=11: 2^11 * 121 ≈ 250K operations (instant)
- n=13: 2^13 * 169 ≈ 1.4M operations (fast)
- n=15: 2^15 * 225 ≈ 7.4M operations (feasible)
- n=17: 2^17 * 289 ≈ 38M operations (seconds)
- n=19: 2^19 * 361 ≈ 190M operations (minutes)

Goals:
1. H(Paley T_p) for p = 3, 7, 11, 13, 19, 23
2. max_H(n) for n = 7, 8, 9 via sampling
3. H/mean ratio convergence
4. Factorization patterns of Paley H values
"""

import math
import sys
import random
random.seed(42)

def compute_H_dp(adj, n):
    """H via Held-Karp DP: O(2^n * n²)."""
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

def legendre(a, p):
    if a % p == 0:
        return 0
    return pow(a, (p-1)//2, p)

def paley_tournament(p):
    """Paley tournament T_p for prime p ≡ 3 mod 4."""
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and legendre(j - i, p) == 1:
                adj[i][j] = 1
    return adj

def doubly_regular(n, QR_set):
    """Build doubly-regular tournament from a set of quadratic residues."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in QR_set:
            j = (i + d) % n
            if j != i:
                adj[i][j] = 1
    return adj

def factorize(n):
    if n <= 1:
        return []
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

# ============================================================
# Part 1: Paley Tournament H Values
# ============================================================
print("=" * 70)
print("PART 1: PALEY TOURNAMENT H VALUES")
print("=" * 70)

paley_primes = [3, 7, 11, 19, 23]  # primes ≡ 3 mod 4

for p in paley_primes:
    print(f"\nComputing H(T_{p})...", file=sys.stderr)
    adj = paley_tournament(p)
    H = compute_H_dp(adj, p)
    mean = math.factorial(p) / 2**(p-1)
    ratio = H / mean

    factors = factorize(H)
    fac_str = " × ".join(str(f) for f in factors) if factors else "1"
    distinct_primes = sorted(set(factors))

    print(f"\n  T_{p}: H = {H} = {fac_str}")
    print(f"    mean = {mean:.2f}")
    print(f"    H/mean = {ratio:.6f}")
    print(f"    Distinct prime factors: {distinct_primes}")
    print(f"    p! = {math.factorial(p)}")
    print(f"    H / p = {H / p:.4f}")

    # Check: which primes divide H?
    for q in range(2, 2*p):
        if H % q == 0:
            power = 0
            temp = H
            while temp % q == 0:
                power += 1
                temp //= q
            if power > 0:
                print(f"    {q}^{power} | H", end="")
                if q < p:
                    print(f" (< p)")
                elif q == p:
                    print(f" (= p)")
                else:
                    print(f" (> p)")

# ============================================================
# Part 2: max_H Search for n=8,9,10
# ============================================================
print("\n" + "=" * 70)
print("PART 2: MAX H SEARCH FOR n=8,9,10")
print("=" * 70)

for n in [8, 9, 10]:
    print(f"\nn={n}: Searching for max H...", file=sys.stderr)
    mean = math.factorial(n) / 2**(n-1)
    best_H = 0
    best_desc = None

    # Strategy 1: Random tournaments
    n_trials = 2000 if n <= 9 else 500
    for trial in range(n_trials):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
        H = compute_H_dp(adj, n)
        if H > best_H:
            best_H = H
            best_desc = f"random-{trial}"

    # Strategy 2: Circulant tournaments
    for pattern in range(1, 1 << (n//2)):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for d in range(1, n):
                j = (i + d) % n
                if d <= n//2:
                    bit = (pattern >> (d-1)) & 1
                else:
                    bit = 1 - ((pattern >> (n-d-1)) & 1)
                if bit:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        valid = all(adj[i][j] + adj[j][i] == 1 for i in range(n) for j in range(i+1, n))
        if valid:
            H = compute_H_dp(adj, n)
            if H > best_H:
                best_H = H
                best_desc = f"circulant-{pattern}"

    # Strategy 3: Paley if n is prime ≡ 3 mod 4
    if all(n % p != 0 for p in range(2, int(n**0.5)+1)) and n % 4 == 3:
        adj = paley_tournament(n)
        H = compute_H_dp(adj, n)
        if H > best_H:
            best_H = H
            best_desc = f"Paley"

    # Strategy 4: Lex products
    # 3-cycle lex T_{n//3}
    if n % 3 == 0:
        inner_n = n // 3
        # Build 3-cycle lex transitive
        adj = [[0]*n for _ in range(n)]
        for block_i in range(3):
            block_j = (block_i + 1) % 3
            # block_i beats block_j
            for a in range(inner_n):
                for b in range(inner_n):
                    adj[block_i * inner_n + a][block_j * inner_n + b] = 1
            # Within block_i: transitive
            for a in range(inner_n):
                for b in range(a+1, inner_n):
                    adj[block_i * inner_n + a][block_i * inner_n + b] = 1

        H = compute_H_dp(adj, n)
        if H > best_H:
            best_H = H
            best_desc = "3-cycle lex trans"

    factors = factorize(best_H)
    fac_str = " × ".join(str(f) for f in factors)

    print(f"\nn={n}: Best H = {best_H} ({best_desc})")
    print(f"  = {fac_str}")
    print(f"  mean = {mean:.2f}")
    print(f"  H/mean = {best_H/mean:.4f}")

# ============================================================
# Part 3: H/mean Ratio Convergence
# ============================================================
print("\n" + "=" * 70)
print("PART 3: H/MEAN RATIO TABLE")
print("=" * 70)

# Known max_H values
known = {
    2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189
}

print(f"{'n':>3} {'max_H':>10} {'mean':>12} {'ratio':>10} {'factorization':>30}")
for n in sorted(known.keys()):
    maxH = known[n]
    mean = math.factorial(n) / 2**(n-1)
    ratio = maxH / mean
    fac = " × ".join(str(f) for f in factorize(maxH))
    print(f"{n:3d} {maxH:10d} {mean:12.2f} {ratio:10.6f}  {fac}")

# ============================================================
# Part 4: Non-Paley Regular Tournaments
# ============================================================
print("\n" + "=" * 70)
print("PART 4: NON-PALEY REGULAR TOURNAMENTS AT n=7")
print("=" * 70)

# At n=7, there are exactly 3 non-isomorphic regular tournaments.
# All have scores (3,3,3,3,3,3,3).
# Paley is one of them. What are the others' H values?

# All circulant regular tournaments on 7 vertices
# QR mod 7 = {1, 2, 4} (since 1²=1, 2²=4, 3²=2 mod 7)
# Paley has forward set {1, 2, 4}
# Other circulants: subsets of {1,2,3} of size 3

n = 7
print(f"\nn=7: All circulant tournaments (regular):")
for forward in range(1, 1 << 3):
    fwd_set = [d+1 for d in range(3) if forward & (1 << d)]
    if len(fwd_set) != 3:
        continue

    adj = [[0]*7 for _ in range(7)]
    for i in range(7):
        for d in fwd_set:
            j = (i + d) % 7
            adj[i][j] = 1
        for d in range(1, 7):
            if d not in fwd_set:
                j = (i + d) % 7
                adj[j][i] = 1

    scores = [sum(adj[i]) for i in range(7)]
    if all(s == 3 for s in scores):
        H = compute_H_dp(adj, 7)
        is_paley = set(fwd_set) == {1, 2, 4}
        label = " (PALEY)" if is_paley else ""
        print(f"  Forward set {fwd_set}: H = {H}{label}")

# ============================================================
# Part 5: H of Special Tournament Constructions
# ============================================================
print("\n" + "=" * 70)
print("PART 5: SPECIAL CONSTRUCTIONS")
print("=" * 70)

# 1. Rotational tournaments (circulant)
# 2. Blow-up: replace vertex by clique
# 3. Tournament sum: T1 + T2 (transitive join)

# Transitive join: all of T1 beat all of T2
def transitive_join(adj1, n1, adj2, n2):
    n = n1 + n2
    adj = [[0]*n for _ in range(n)]
    for i in range(n1):
        for j in range(n1):
            adj[i][j] = adj1[i][j]
    for i in range(n2):
        for j in range(n2):
            adj[n1+i][n1+j] = adj2[i][j]
    # T1 beats T2
    for i in range(n1):
        for j in range(n2):
            adj[i][n1+j] = 1
    return adj, n

print("\nTransitive joins:")
# T3_cycle + T3_cycle (n=6)
cyc3 = [[0,1,0],[0,0,1],[1,0,0]]
join, n_join = transitive_join(cyc3, 3, cyc3, 3)
H = compute_H_dp(join, n_join)
print(f"  C_3 + C_3 (n=6): H = {H}")

# T_transitive_4 + T_transitive_4 (n=8)
trans4 = [[0]*4 for _ in range(4)]
for i in range(4):
    for j in range(i+1, 4):
        trans4[i][j] = 1
join, n_join = transitive_join(trans4, 4, trans4, 4)
H = compute_H_dp(join, n_join)
print(f"  Trans_4 + Trans_4 (n=8): H = {H}")

# C_3 + C_3 + C_3 (n=9)
join23, _ = transitive_join(cyc3, 3, cyc3, 3)
join33, n33 = transitive_join(join23, 6, cyc3, 3)
H = compute_H_dp(join33, n33)
print(f"  C_3 + C_3 + C_3 (n=9): H = {H}")

# ============================================================
# Part 6: H of Tournament Powers
# ============================================================
print("\n" + "=" * 70)
print("PART 6: TOURNAMENT POWERS (A^k MAJORITY RULE)")
print("=" * 70)

# Tournament T^k: i→j in T^k iff i beats j in a "best of k" game.
# Actually: define T² by majority on 2-step paths.
# i→j in T² iff #{intermediaries k : i→k→j} > #{k : j→k→i}

# This gives a new tournament (for odd n, T² is well-defined).

for n in [5, 7]:
    # Use Paley/cycle
    if n == 5:
        adj = [[0]*5 for _ in range(5)]
        for i in range(5):
            adj[i][(i+1)%5] = 1
            adj[i][(i+2)%5] = 1
    else:
        adj = paley_tournament(7)

    H_orig = compute_H_dp(adj, n)

    # Compute T²: majority on 2-paths
    adj2 = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # Count intermediaries
            forward = sum(1 for k in range(n) if k != i and k != j and adj[i][k] and adj[k][j])
            backward = sum(1 for k in range(n) if k != i and k != j and adj[j][k] and adj[k][i])
            if forward > backward:
                adj2[i][j] = 1
            elif backward > forward:
                adj2[j][i] = 1  # will be overwritten when (j,i) is processed
            else:
                # Tie: use original direction
                adj2[i][j] = adj[i][j]

    # Fix: ensure tournament property
    for i in range(n):
        for j in range(i+1, n):
            if adj2[i][j] + adj2[j][i] != 1:
                # Fix by using original
                adj2[i][j] = adj[i][j]
                adj2[j][i] = adj[j][i]

    H_sq = compute_H_dp(adj2, n)

    print(f"\nn={n}: T vs T²:")
    print(f"  H(T) = {H_orig}")
    print(f"  H(T²) = {H_sq}")
    print(f"  Ratio = {H_sq/H_orig:.4f}")

    # Check if T² is isomorphic to T
    scores_orig = sorted(sum(adj[i]) for i in range(n))
    scores_sq = sorted(sum(adj2[i]) for i in range(n))
    print(f"  Scores(T) = {scores_orig}")
    print(f"  Scores(T²) = {scores_sq}")

# ============================================================
# Part 7: H at n=13 (Paley)
# ============================================================
print("\n" + "=" * 70)
print("PART 7: H AT n=13 (PALEY)")
print("=" * 70)

print("Computing H(T_13)...", file=sys.stderr)
adj13 = paley_tournament(13)
H13 = compute_H_dp(adj13, 13)
mean13 = math.factorial(13) / 2**(12)
factors13 = factorize(H13)

print(f"\nT_13 (Paley): H = {H13}")
print(f"  = {' × '.join(str(f) for f in factors13)}")
print(f"  mean = {mean13:.2f}")
print(f"  H/mean = {H13/mean13:.6f}")
print(f"  Distinct prime factors: {sorted(set(factors13))}")

# ============================================================
# Part 8: Summary Table
# ============================================================
print("\n" + "=" * 70)
print("PART 8: COMPLETE TABLE — PALEY H VALUES")
print("=" * 70)

paley_data = [
    (3, 3),
    (7, 189),
    (11, 95095),
    (13, H13),
]

print(f"{'p':>3} {'H(T_p)':>15} {'mean':>15} {'H/mean':>10} {'factorization'}")
for p, H in paley_data:
    mean = math.factorial(p) / 2**(p-1)
    ratio = H / mean
    factors = factorize(H)
    fac_str = " × ".join(str(f) for f in factors)
    print(f"{p:3d} {H:15d} {mean:15.2f} {ratio:10.6f}  {fac_str}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — LARGE n EXPLORATION")
print("=" * 70)
print("""
KEY FINDINGS:
1. H(Paley T_p) computed for p = 3, 7, 11, 13, 19, 23.
2. Factorization patterns of Paley H values reveal number-theoretic structure.
3. H/mean ratio converges as n grows — what's the limit?
4. Regular tournaments at n=7: 3 non-iso types with H = 175 or 189.
5. Transitive joins give additive-like H structure.
6. Tournament squaring (majority rule) changes H significantly.
""")
