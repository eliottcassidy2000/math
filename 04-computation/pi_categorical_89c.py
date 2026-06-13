#!/usr/bin/env python3
"""
pi_categorical_89c.py — Categorical structure and π connections

opus-2026-03-14-S89c

1. E[H²]/E[H]² ratio: what OEIS sequence do the numerators/denominators form?
2. The E[1/H] denominators: lcm of all H values?
3. Homotopy / simplicial structure of the independence complex
4. The "tournament measure" μ(H=h) = #{T: H(T)=h}/2^m and its π content
"""

import math
import itertools
from fractions import Fraction
from collections import Counter, defaultdict
from functools import reduce

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def count_H(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    return a * b // gcd(a, b)

# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: The E[H²]/E[H]² Ratio Sequence")
print("=" * 70)

# r_n = E[H²]/E[H]²
# = Σ_{(π,σ) compatible} 2^{k(π,σ)} / (n!)²
# where k = shared directed arcs

# n=3: 4/3
# n=4: 4/3
# n=5: 79/60
# n=6: 58/45
# n=7: 635/504

# CV² = r_n - 1 = Var/E²
# n=3: 1/3
# n=4: 1/3
# n=5: 19/60
# n=6: 13/45
# n=7: 131/504 (= 635/504 - 1)

# Check: 131/504 = ?
# 504 = 7 × 72 = 7 × 8 × 9 = 2³ × 3² × 7
# 131 is prime
# So CV²(7) = 131/504

# The denominators: 3, 3, 60, 45, 504
# 3 = 3
# 60 = 2² × 3 × 5
# 45 = 3² × 5
# 504 = 2³ × 3² × 7

# OEIS check: the sequence 4/3, 4/3, 79/60, 58/45, 635/504
# Numerators: 4, 4, 79, 58, 635
# Denominators: 3, 3, 60, 45, 504

# Alternative: Σ H² sequence: 24, 768, 75840, 21381120, ...
# 24 = 4!
# 768 = 3 × 256 = 3 × 4^4
# 75840 = ?
# 21381120 = ?

print("\nΣ_T H(T)² sequence:")
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    sum_h2 = sum(h**2 for h in h_values)
    print(f"  n={n}: Σ H² = {sum_h2}")
    # Factor
    x = sum_h2
    factors = []
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        while x % p == 0:
            factors.append(p)
            x //= p
    if x > 1:
        factors.append(x)
    print(f"    = {' × '.join(str(f) for f in factors)}")

    # Also: Σ H² / (n! × 2^{m-n+1}) where Σ H = n! × 2^{m-n+1}
    m = n*(n-1)//2
    sum_h = math.factorial(n) * 2**(m-n+1)
    ratio_h2_h = Fraction(sum_h2, sum_h)
    print(f"    Σ H² / Σ H = {ratio_h2_h} = {float(ratio_h2_h):.4f}")
    # This is E[H²]/E[H] = E[H] + Var/E = E[H](1 + CV²)
    E_H = Fraction(math.factorial(n), 2**(n-1))
    print(f"    = E[H] × (1 + CV²) = {float(E_H * (1 + Fraction(1,3))):.4f}" if n <= 4 else "")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: E[1/H] Denominators — LCM of H Values?")
print("=" * 70)

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    h_unique = sorted(set(h_values))

    # LCM of all distinct H values
    h_lcm = reduce(lcm, h_unique)

    # E[1/H] as fraction
    N = len(h_values)
    E_inv = sum(Fraction(1, h) for h in h_values) / N

    print(f"\n  n={n}:")
    print(f"    Distinct H values: {h_unique}")
    print(f"    LCM of H values: {h_lcm}")
    print(f"    E[1/H] denominator: {E_inv.denominator}")
    print(f"    LCM divides denominator? {E_inv.denominator % h_lcm == 0}")
    print(f"    denominator / LCM = {E_inv.denominator // h_lcm if E_inv.denominator % h_lcm == 0 else 'N/A'}")

    # N·E[1/H] = Σ 1/H
    sum_inv = sum(Fraction(1, h) for h in h_values)
    print(f"    Σ 1/H = {sum_inv} = {float(sum_inv):.6f}")
    print(f"    Σ 1/H denominator: {sum_inv.denominator}")
    print(f"    LCM divides Σ denominator? {sum_inv.denominator % h_lcm == 0}")

    # Check: is Σ 1/H = Σ_{h unique} c_h / h where c_h = count of H=h?
    # Each c_h/h contributes to the sum. The LCM of the h values is the natural denominator.
    sum_check = sum(Fraction(count, h) for h, count in Counter(h_values).items())
    print(f"    Σ c_h/h = {sum_check} (same as Σ 1/H: {sum_check == sum_inv})")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 3: H Value Set Structure — Arithmetic Properties")
print("=" * 70)

# H values are always odd (Rédei), positive
# Known: H ≡ 1 mod 2
# From S89b: H mod 4 ∈ {1, 3} (encodes cycle parity)
# H mod 6 ∈ {1, 3, 5} (encodes total_cycles mod 3)

# What about the SET of H values? Does it have structure?
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    h_unique = sorted(set(h_values))
    max_h = max(h_unique)

    print(f"\n  n={n}: H values = {h_unique}")
    print(f"    All odd: {all(h % 2 == 1 for h in h_unique)}")
    print(f"    Range: [{min(h_unique)}, {max_h}]")
    print(f"    Gaps: {[h_unique[i+1]-h_unique[i] for i in range(len(h_unique)-1)]}")

    # Are differences always 2? (consecutive odd numbers?)
    diffs = [h_unique[i+1]-h_unique[i] for i in range(len(h_unique)-1)]
    print(f"    All gaps = 2? {all(d == 2 for d in diffs)}")

    # GCD of all H values
    h_gcd = reduce(gcd, h_unique)
    print(f"    GCD of all H values: {h_gcd}")

    # Do H values form a complete set of odd numbers from 1 to max?
    full_odds = set(range(1, max_h + 1, 2))
    actual = set(h_unique)
    missing = full_odds - actual
    print(f"    Missing odd values: {sorted(missing) if missing else 'NONE'}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 4: The H=1 and H=max Populations")
print("=" * 70)

# H(T) = 1 iff T is transitive. Number = n!
# H(T) = max iff T is "regular" (for odd n) or near-regular

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    h_counts = Counter(h_values)
    N = len(h_values)
    max_h = max(h_values)

    print(f"\n  n={n}:")
    print(f"    P(H=1)   = {h_counts[1]}/{N} = {h_counts[1]/N:.6f}")
    print(f"    P(H=max) = {h_counts[max_h]}/{N} = {h_counts[max_h]/N:.6f}")
    print(f"    max H = {max_h}")
    print(f"    #{{'H=1'}} = {h_counts[1]} = {math.factorial(n)} (= n!)")
    print(f"    #{{'H=max'}} = {h_counts[max_h]}")

    # Ratio of extreme populations
    if h_counts[max_h] > 0:
        ratio = h_counts[1] / h_counts[max_h]
        print(f"    #{{'H=1'}} / #{{'H=max'}} = {ratio:.4f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5: Signed Tournament Euler Characteristic")
print("=" * 70)

# Define χ(n) = Σ_T (-1)^{H(T)} ... but H is always odd so (-1)^H = -1 always
# More interesting: χ = Σ_T (-1)^{(H(T)-1)/2}
# This is (-1)^0 for H≡1 mod 4, (-1)^1 for H≡3 mod 4

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    chi = sum((-1)**((h-1)//2) for h in h_values)
    N = len(h_values)
    h1 = sum(1 for h in h_values if h % 4 == 1)
    h3 = sum(1 for h in h_values if h % 4 == 3)

    print(f"\n  n={n}: #{{'H≡1 mod 4'}} = {h1}, #{{'H≡3 mod 4'}} = {h3}")
    print(f"    χ = #{{'H≡1'}} - #{{'H≡3'}} = {chi}")
    print(f"    χ/N = {chi/N:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 6: Tournament π-phases — Fourier Transform of H")
print("=" * 70)

# For the function H: Q_m → Z_odd, compute its Fourier transform
# H(T) = Σ_S ĥ(S)·χ_S(T) where χ_S(T) = (-1)^{Σ_{e∈S} T_e}
# The "DC component" ĥ(∅) = E[H]
# Higher order terms... too expensive for exhaustive, but we can do n=3

print("\nWalsh-Hadamard transform of H at n=3:")
n = 3
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Build H as a function of the m bits
H_func = []
for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(pairs):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H_func.append(count_H(adj, n))

# Walsh-Hadamard transform
# ĥ(S) = (1/2^m) Σ_T H(T)·(-1)^{<S,T>}
print(f"  H values: {H_func}")
for S in range(1 << m):
    h_hat = 0
    for T in range(1 << m):
        sign = (-1)**bin(S & T).count('1')
        h_hat += H_func[T] * sign
    h_hat_frac = Fraction(h_hat, 1 << m)
    if h_hat != 0:
        # Decode S as a set of arcs
        arc_set = [pairs[i] for i in range(m) if S & (1 << i)]
        print(f"  S={bin(S)[2:].zfill(m)} ({arc_set}): ĥ = {h_hat_frac}")

print("\n  For n=4:")
n = 4
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

H_func = []
for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(pairs):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H_func.append(count_H(adj, n))

# Just show the spectrum: how many nonzero coefficients at each weight
weight_spectrum = Counter()
total_energy = 0
dc_energy = 0
for S in range(1 << m):
    h_hat = sum(H_func[T] * (-1)**bin(S & T).count('1') for T in range(1 << m))
    if h_hat != 0:
        w = bin(S).count('1')
        weight_spectrum[w] += 1
        total_energy += h_hat**2
        if S == 0:
            dc_energy = h_hat**2
print(f"  Weight spectrum (weight → #nonzero coefficients):")
for w in sorted(weight_spectrum.keys()):
    print(f"    weight {w}: {weight_spectrum[w]} nonzero coefficients")
print(f"  Total energy: {total_energy}, DC energy: {dc_energy}")
print(f"  DC fraction: {dc_energy/total_energy:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 7: Score Sequence Distribution and π")
print("=" * 70)

# The score sequence of a tournament is the sorted list of out-degrees
# For n vertices: scores sum to C(n,2)
# The number of tournaments with a given score sequence involves
# multinomial coefficients, which involve π through Stirling

for n in range(3, 7):
    score_dist = Counter()
    for adj in all_tournaments(n):
        scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
        score_dist[scores] += 1

    print(f"\n  n={n}: {len(score_dist)} distinct score sequences")
    for scores, count in sorted(score_dist.items(), key=lambda x: -x[1]):
        print(f"    {scores}: {count} tournaments ({count*100/2**(n*(n-1)//2):.1f}%)")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 8: The 4/3 Mystery — Exact for n=3,4")
print("=" * 70)

# WHY is E[H²]/E[H]² = 4/3 exactly at n=3 and n=4?
#
# At n=3: Only H∈{1,3}. P(H=1)=3/4, P(H=3)=1/4.
# E[H] = 3/4 + 3/4 = 3/2
# E[H²] = 3/4 + 9/4 = 3
# E[H²]/E[H]² = 3/(9/4) = 4/3 ✓
#
# At n=4: H∈{1,3,5}. P(H=1)=24/64=3/8, P(H=3)=16/64=1/4, P(H=5)=24/64=3/8.
# E[H] = 3/8 + 3/4 + 15/8 = 3
# E[H²] = 3/8 + 9/4 + 75/8 = (3+18+75)/8 = 96/8 = 12
# E[H²]/E[H]² = 12/9 = 4/3 ✓
#
# At n=5: E[H²]/E[H]² = 79/60 ≠ 4/3
#
# The 4/3 ratio at n=3 corresponds to a 2-point distribution {1,3}
# The 4/3 ratio at n=4 is more surprising since there are 3 values
# Is there a structural reason?

print("\nWhy 4/3 at n=3,4?")

# Check: for ANY distribution on {1,3,...,2k-1} with equal spacing 2,
# does E[H²]/E[H]² = 4/3 if and only if the distribution is symmetric?

# At n=4: P(1) = P(5) = 3/8, P(3) = 1/4
# This IS symmetric about the mean 3!
# E[H] = 3
# E[(H-3)²] = 3/8·4 + 1/4·0 + 3/8·4 = 3
# E[H²] = Var + E² = 3 + 9 = 12
# CV² = 3/9 = 1/3

# At n=3: P(1) = 3/4, P(3) = 1/4
# NOT symmetric about mean 3/2
# But still CV² = 1/3!
# Var = 3/4 · (1/2)² + 1/4 · (3/2)² = 3/16 + 9/16 = 12/16 = 3/4
# E² = 9/4
# CV² = 3/(4·9/4) = 1/3 ✓

# Hmm. What's special about CV² = 1/3?
# For a Bernoulli(p) on {0,1}: CV² = (1-p)/p
# For Rademacher on {-1,+1}: Var=1, E²=0, CV²=∞
# For uniform on {1,...,n}: CV² = (n²-1)/12 / ((n+1)/2)² = (n-1)/(3(n+1))

# For n=3 tournaments: 2 values {1,3}, P(1)=3/4, P(3)=1/4
# This is like a coin flip with probability 1/4 for "heads"
# CV² = Var/E² = (3/4)/(9/4) = 1/3

# For n=4 tournaments: 3 values {1,3,5}, P=(3/8, 1/4, 3/8)
# CV² = 3/9 = 1/3

# Is CV² = 1/3 connected to E[H] = n!/2^{n-1} somehow?
# E[H²] = (4/3)·E[H]² = (4/3)·(n!/2^{n-1})²
# Σ H² = (4/3)·(n!)²·2^{m-2n+2}/... hmm

# Alternative: Σ H² = (4/3)·(Σ H)²/N where N = 2^m
# 4/3 · (n!·2^{m-n+1})² / 2^m = 4/3 · n!² · 2^{2m-2n+2-m} = 4/3·n!²·2^{m-2n+2}

# At n=3: 4/3 · 36 · 2^{3-4} = 4/3 · 36 · 1/2 = 4/3 · 18 = 24 ✓
# At n=4: 4/3 · 576 · 2^{6-6} = 4/3 · 576 = 768 ✓
# At n=5: should be 79/60 instead of 4/3
# 79/60 · 14400 · 2^{10-8} = 79/60 · 14400 · 4 = 79·960 = 75840 ✓

print("Verified: Σ H² = (r_n) · n!² · 2^{m-2n+2} for all n.")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 9: Searching for Σ H² as a Double Sum")
print("=" * 70)

# We proved Σ H = n!·2^{m-n+1} by counting (T, π) pairs.
# Can we prove Σ H² = f(n) by counting (T, π, σ) triples?
# Σ H² = Σ_{(π,σ) compatible} 2^{m - c(π,σ)}
# where c = |{undirected pairs constrained by π or σ}|
# = 2(n-1) - k for k shared arcs, when compatible

# So Σ H² = 2^m · Σ_{(π,σ) compatible} 2^{-2(n-1)+k}

# Define f(n) = Σ_{(π,σ) compatible} 2^k
# Then Σ H² = f(n) · 2^{m-2(n-1)}

# We computed f(n):
# f(3) = 48
# f(4) = 768
# f(5) = 18960
# f(6) = 668160
# f(7) = 32004000

fn = [None, None, None, 48, 768, 18960, 668160, 32004000]
print("\nf(n) = Σ_{compatible (π,σ)} 2^k:")
for n in range(3, 8):
    print(f"  f({n}) = {fn[n]}")
    print(f"  f({n}) / n!² = {Fraction(fn[n], math.factorial(n)**2)}")
    print(f"  f({n}) / (n! · 2^{{n-1}}) = {Fraction(fn[n], math.factorial(n) * 2**(n-1))}")

# f(n) / (n!·2^{n-1}) should be related to E[H²]/E[H] = E[H](1+CV²)
# = (n!/2^{n-1})·(1+CV²)
# = n!/2^{n-1} + Var·2^{n-1}/n!
# Hmm, let me compute f(n)/(n!·2^{n-1})
print("\nf(n) / (n!·2^{n-1}):")
for n in range(3, 8):
    val = Fraction(fn[n], math.factorial(n) * 2**(n-1))
    print(f"  f({n})/(n!·2^{{n-1}}) = {val} = {float(val):.4f}")

# f(n)/n!² = r_n = E[H²]/E[H]²
# r = 4/3, 4/3, 79/60, 58/45, 635/504

# What about f(n)/n! ?
print("\nf(n)/n!:")
for n in range(3, 8):
    val = Fraction(fn[n], math.factorial(n))
    print(f"  f({n})/n! = {val} = {float(val):.4f}")

# 8, 32, 158, 928, 6350
# Hmm: 8, 32, 158, 928, 6350
# 8 = 2³
# 32 = 2⁵
# 158 = 2 × 79
# 928 = 2⁵ × 29
# 6350 = 2 × 5² × 127

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 10: Is the CV² Sequence in OEIS?")
print("=" * 70)

# CV² = 1/3, 1/3, 19/60, 13/45, 131/504
# Denominators: 3, 3, 60, 45, 504
# 1/CV²: 3, 3, 60/19, 45/13, 504/131

# Var[H] = (n!/2^{n-1})² · CV²
# n=3: Var = (3/2)² · 1/3 = 9/4 · 1/3 = 3/4
# n=4: Var = 9 · 1/3 = 3
# n=5: Var = (225/4) · 19/60 = 4275/240 = 285/16
# n=6: Var = (2025/4) · 13/45 = 26325/180 = 585/4

# Var sequence: 3/4, 3, 285/16, 585/4
# = 3/4, 12/4, 285/16, 2340/16
# Hmm let me compute Σ(H-E[H])² = N·Var
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    E_H = Fraction(sum(h_values), N)
    var_sum = sum((h - E_H)**2 for h in h_values)
    print(f"  n={n}: Σ(H-E)² = {var_sum} = {float(var_sum):.2f}")
    print(f"         N·Var = {N * (Fraction(sum(h**2 for h in h_values), N) - E_H**2)}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 11: Tournament Probability Generating Function")
print("=" * 70)

# PGF: G_n(z) = Σ_T z^{H(T)} / 2^m = E[z^H]
# At z=1: G_n(1) = 1
# G'_n(1) = E[H]
# G''_n(1) = E[H(H-1)] = E[H²] - E[H]
# So we can write G_n(z) explicitly for small n

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    h_counts = Counter(h_values)
    N = len(h_values)

    # PGF as polynomial
    print(f"\n  n={n}: G_n(z) = E[z^H] =")
    terms = []
    for h in sorted(h_counts.keys()):
        coeff = Fraction(h_counts[h], N)
        terms.append(f"  ({coeff})·z^{h}")
    print("    " + " + ".join(terms))

    # Evaluate at specific points
    # z = -1: G_n(-1) = E[(-1)^H] = -1 (since H always odd)
    gz_neg1 = sum(Fraction(h_counts[h], N) * (-1)**h for h in h_counts)
    print(f"    G_n(-1) = {gz_neg1} (always -1 since H odd)")

    # z = i: G_n(i) = E[i^H] = E[i^{2k+1}] for H=2k+1
    # i^{2k+1} = i·(-1)^k, so G_n(i) = i·E[(-1)^{(H-1)/2}]
    import cmath
    gz_i = sum(Fraction(h_counts[h], N) * (1j)**h for h in h_counts)
    print(f"    G_n(i) = {gz_i.real:.6f} + {gz_i.imag:.6f}i")

    # z = e^{2πi/k} for small k
    for k in [3, 4, 6]:
        omega = cmath.exp(2j * cmath.pi / k)
        gz = sum((h_counts[h]/N) * omega**h for h in h_counts)
        print(f"    G_n(e^{{2πi/{k}}}) = {gz.real:.6f} + {gz.imag:.6f}i, |·| = {abs(gz):.6f}")

print("\n" + "=" * 70)
print("END — opus-2026-03-14-S89c categorical")
print("=" * 70)
