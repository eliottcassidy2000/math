#!/usr/bin/env python3
"""
five_recurrence_deep.py — opus-2026-03-14-S73
Deep exploration of 5 in the recurrence lattice.

Key threads:
1. The Q(√5) bridge: Fibonacci ↔ x=11 ↔ x=31 ↔ ...
2. 5 as the second-order recurrence discriminant
3. The pentanacci (k=5) as optimal approximation
4. 5-adic structure of tournament numbers
5. The role of 5 in modular H arithmetic
"""

import numpy as np
from itertools import combinations, permutations
from fractions import Fraction
from math import gcd, comb, sqrt, log2
import random

def banner(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")

# ─────────────────────────────────────────────────────────────────────
# PART 1: The Q(√5) Tower — all x values with √5 discriminant
# ─────────────────────────────────────────────────────────────────────
banner("PART 1: THE Q(√5) TOWER")

print("Second-order recurrence f(n) = f(n-1) + x·f(n-2)")
print("Discriminant Δ = 1 + 4x")
print("Roots = (1 ± √Δ) / 2")
print()
print("x values where Δ involves √5 (i.e., Δ = 5m² for integer m):")
print("  x = (5m² - 1)/4, requiring m ≡ 1 or 3 (mod 4) for integer x")
print()

# Find all integer x with Δ = 5m²
q5_tower = []
for m in range(1, 200):
    val = 5*m*m - 1
    if val % 4 == 0:
        x = val // 4
        delta = 1 + 4*x
        sqrt_delta = m * sqrt(5)
        r1 = (1 + sqrt_delta) / 2
        r2 = (1 - sqrt_delta) / 2
        q5_tower.append((m, x, delta, r1, r2))

print(f"{'m':>4} {'x':>8} {'Δ':>10} {'√Δ':>12} {'r₊':>12} {'r₋':>12}  {'r₊ in terms of φ':>20}")
print("-" * 90)
phi = (1 + sqrt(5)) / 2
for m, x, delta, r1, r2 in q5_tower[:20]:
    # Express r1 in terms of φ: r1 = (1 + m√5)/2 = a + bφ
    # φ = (1+√5)/2, so √5 = 2φ-1
    # r1 = (1 + m(2φ-1))/2 = (1 - m + 2mφ)/2 = (1-m)/2 + mφ
    a_coeff = (1 - m) / 2  # rational part
    b_coeff = m  # φ coefficient
    if m % 2 == 1:  # integer coefficients
        phi_expr = f"{int((1-m)//2)} + {m}φ" if (1-m)//2 != 0 else f"{m}φ"
    else:
        phi_expr = f"({1-m}/2 + {m}φ)"
    print(f"{m:4d} {x:8d} {delta:10d} {m:3d}√5={sqrt(delta):7.3f} {r1:12.6f} {r2:12.6f}  {phi_expr:>20}")

print()
print("PATTERN: The Q(√5) tower has x = 1, 11, 31, 61, 101, 151, ...")
print("Differences: ", [q5_tower[i+1][1] - q5_tower[i][1] for i in range(min(10, len(q5_tower)-1))])
print("Second differences: ", [q5_tower[i+2][1] - 2*q5_tower[i+1][1] + q5_tower[i][1] for i in range(min(9, len(q5_tower)-2))])

# Check if x values follow a pattern
xs = [t[1] for t in q5_tower[:15]]
print(f"\nx values: {xs}")
print(f"These are x = (5m²-1)/4 for m = {[t[0] for t in q5_tower[:15]]}")

# Check: are these related to centered pentagonal numbers?
print("\nComparison with known sequences:")
for m, x, _, _, _ in q5_tower[:10]:
    print(f"  m={m}: x={x}, 5·T(m-1)+m-1 = {5*m*(m-1)//2 + m - 1}, Pent = {m*(3*m-1)//2}")

# ─────────────────────────────────────────────────────────────────────
# PART 2: Fibonacci ↔ x=11 explicit bridge
# ─────────────────────────────────────────────────────────────────────
banner("PART 2: THE FIBONACCI ↔ x=11 BRIDGE")

print("Fibonacci (x=1): f(n) = f(n-1) + f(n-2)")
print("  Roots: φ = (1+√5)/2 ≈ 1.618, ψ = (1-√5)/2 ≈ -0.618")
print("  F(n) = (φⁿ - ψⁿ)/√5")
print()
print("x=11 recurrence: f(n) = f(n-1) + 11·f(n-2)")
print("  Roots: (1+3√5)/2 = 3φ-1 ≈ 3.854, (1-3√5)/2 = -3/φ ≈ -2.854")
print("  G(n) = ((3φ-1)ⁿ - (-(3/φ))ⁿ) / (3√5)")
print()

# Verify the φ-connection
r11_plus = (1 + 3*sqrt(5))/2
r11_minus = (1 - 3*sqrt(5))/2
print(f"  3φ - 1 = {3*phi - 1:.10f}")
print(f"  (1+3√5)/2 = {r11_plus:.10f}")
print(f"  Match: {abs(3*phi - 1 - r11_plus) < 1e-12}")
print()

# Generate both sequences and look for relationships
fib = [0, 1]
g11 = [0, 1]
for i in range(2, 20):
    fib.append(fib[-1] + fib[-2])
    g11.append(g11[-1] + 11*g11[-2])

print("  n   Fib(n)     G₁₁(n)    G₁₁/Fib     G₁₁ mod 5")
for n in range(1, 16):
    ratio = g11[n] / fib[n] if fib[n] != 0 else float('inf')
    print(f"  {n:2d}  {fib[n]:8d}  {g11[n]:10d}  {ratio:10.4f}  {g11[n] % 5:5d}")

print()
# Check: G₁₁(n) mod 5
print("G₁₁(n) mod 5:", [g11[n] % 5 for n in range(1, 16)])
print("Fib(n) mod 5:", [fib[n] % 5 for n in range(1, 16)])
print()

# The key: both live in Q(√5), so there should be algebraic relations
# G₁₁(n) = ((3φ-1)ⁿ - (2-3φ)ⁿ) / (3√5)
# F(n) = (φⁿ - (-1/φ)ⁿ) / √5
# So G₁₁(n) / F(n) → (3φ-1)ⁿ / (3φⁿ) = ((3φ-1)/(φ))ⁿ / 3
# (3φ-1)/φ = 3 - 1/φ = 3 - (√5-1)/2 = (7-√5)/2

ratio_limit = (3*phi - 1) / phi
print(f"Asymptotic ratio G₁₁(n)/F(n) → ((3φ-1)/φ)ⁿ / 3")
print(f"  (3φ-1)/φ = 3 - 1/φ = {ratio_limit:.10f}")
print(f"  = (7-√5)/2 = {(7-sqrt(5))/2:.10f}")
print()

# ─────────────────────────────────────────────────────────────────────
# PART 3: 5-adic valuation of tournament numbers
# ─────────────────────────────────────────────────────────────────────
banner("PART 3: 5-ADIC STRUCTURE OF TOURNAMENT NUMBERS")

# Jacobsthal numbers: T(n) = (2^n + (-1)^n) / 3 for n ≥ 0... 
# Actually tournament counts: t(n) = 2^{n(n-1)/2} / n! * ... 
# Let's look at Jacobsthal mod 5

def jacobsthal(n):
    """J(n) = (2^n - (-1)^n) / 3"""
    return (2**n - (-1)**n) // 3

def v5(n):
    """5-adic valuation"""
    if n == 0: return float('inf')
    v = 0
    while n % 5 == 0:
        v += 1
        n //= 5
    return v

print("Jacobsthal numbers and their 5-adic valuations:")
print(f"{'n':>4} {'J(n)':>15} {'v₅(J)':>6} {'J mod 5':>8} {'J mod 25':>8}")
for n in range(1, 30):
    j = jacobsthal(n)
    print(f"{n:4d} {j:15d} {v5(j):6d} {j%5:8d} {j%25:8d}")

print()
# When is J(n) ≡ 0 (mod 5)?
# J(n) = (2^n - (-1)^n) / 3
# 5 | J(n) iff 2^n ≡ (-1)^n (mod 15)
# ord(2 mod 15) = 4: 2,4,8,16≡1, so period 4
# (-1)^n mod 15: 14,1,14,1,...
# 2^n mod 15: 2,4,8,1,2,4,8,1,...
# Need 2^n ≡ (-1)^n (mod 15):
# n=1: 2 ≡ 14? No
# n=2: 4 ≡ 1? No
# n=3: 8 ≡ 14? No
# n=4: 1 ≡ 1? Yes!
# So J(n) ≡ 0 (mod 5) iff n ≡ 0 (mod 4)?

j_mod5 = [(n, jacobsthal(n) % 5) for n in range(1, 25)]
zeros = [n for n, jm in j_mod5 if jm == 0]
print(f"J(n) ≡ 0 (mod 5) at n = {zeros}")
print(f"Pattern: n ≡ {[z % 4 for z in zeros]} (mod 4)")

# Now look at 2^(n(n-1)/2) mod 5
print("\n2^{n(n-1)/2} mod 5 (tournament-related):")
for n in range(1, 20):
    edges = n*(n-1)//2
    val = pow(2, edges, 5)
    print(f"  n={n:2d}: edges={edges:3d}, 2^edges mod 5 = {val}, v₅(2^edges-1) = ", end="")
    # Don't compute 2^edges - 1 for large edges, just check mod
    v = 0
    temp = pow(2, edges) - 1 if edges < 60 else None
    if temp is not None:
        print(f"{v5(temp)}")
    else:
        print("(too large)")

# ─────────────────────────────────────────────────────────────────────
# PART 4: The pentanacci as optimal 5-step lookback
# ─────────────────────────────────────────────────────────────────────
banner("PART 4: PENTANACCI — OPTIMAL 5-STEP LOOKBACK")

print("The k-nacci root φ_k satisfies z^{k+1} - 2z^k + 1 = 0 (z≠1).")
print("At k=5: z⁶ - 2z⁵ + 1 = 0")
print()

# Factor z^6 - 2z^5 + 1
# We know (z-1) divides z^{k+1} - 2z^k + 1 if we include it
# Actually z^{k+1} - 2z^k + 1 = (z-1)(z^k - z^{k-1} - ... - 1)
# For k=5: z^6-2z^5+1 = (z-1)(z^5-z^4-z^3-z^2-z-1)

coeffs_5nacci = [1, -1, -1, -1, -1, -1]  # z^5 - z^4 - z^3 - z^2 - z - 1
roots_5 = np.roots(coeffs_5nacci)
print(f"Pentanacci characteristic roots:")
for i, r in enumerate(sorted(roots_5, key=lambda x: -abs(x))):
    print(f"  r_{i+1} = {r:.10f}  |r| = {abs(r):.10f}")

dom = max(roots_5, key=abs).real
print(f"\nDominant root φ₅ = {dom:.15f}")
print(f"Gap to 2: {2-dom:.15e}")
print(f"1/2⁵ = {1/32:.15e}")
print(f"Ratio (gap)/(1/32) = {(2-dom)*32:.15f}")
print()

# The pentanacci has a special property: k=5 is where 
# the gap first drops below the Fibonacci gap at k=2
print("Gap comparison across k:")
print(f"{'k':>4} {'φ_k':>15} {'gap=2-φ_k':>15} {'1/2^k':>15} {'gap/Fib_gap':>15}")
fib_gap = 2 - phi  # gap at k=2
for k in range(2, 12):
    coeffs = [1] + [-1]*k
    roots = np.roots(coeffs)
    dom_root = max(r.real for r in roots if abs(r.imag) < 1e-10)
    gap = 2 - dom_root
    print(f"{k:4d} {dom_root:15.12f} {gap:15.12e} {1/2**k:15.12e} {gap/fib_gap:15.12f}")

# ─────────────────────────────────────────────────────────────────────
# PART 5: The golden ratio in tournament H-values
# ─────────────────────────────────────────────────────────────────────
banner("PART 5: GOLDEN RATIO φ IN H-VALUES")

def adj_matrix(n, perm_idx):
    """Generate tournament adjacency matrix from index."""
    m = n*(n-1)//2
    A = [[0]*n for _ in range(n)]
    bits = perm_idx
    for i in range(n):
        for j in range(i+1, n):
            if bits & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            bits >>= 1
    return A

def count_ham_paths(A, n):
    """Count Hamiltonian paths via inclusion-exclusion on permutations."""
    count = 0
    for p in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

# At n=5, check H values mod 5 and relationship to φ
print("H(T) distribution at n=5 (all 1024 tournaments):")
h_counts = {}
h_vals = []
for idx in range(1024):
    A = adj_matrix(5, idx)
    h = count_ham_paths(A, 5)
    h_vals.append(h)
    h_counts[h] = h_counts.get(h, 0) + 1

for h in sorted(h_counts.keys()):
    mod5 = h % 5
    mod_phi = ""  # conceptual
    print(f"  H={h:3d}: count={h_counts[h]:4d}, H mod 5 = {mod5}")

print(f"\nH mod 5 distribution:")
mod5_dist = {}
for h in h_vals:
    m = h % 5
    mod5_dist[m] = mod5_dist.get(m, 0) + 1
for m in sorted(mod5_dist.keys()):
    print(f"  H ≡ {m} (mod 5): {mod5_dist[m]} ({mod5_dist[m]/1024:.3f})")

# Check: is H always ≡ 1 (mod 2)? (Rédei)
print(f"\nH mod 2: all odd? {all(h % 2 == 1 for h in h_vals)}")
print(f"H mod 3: distribution = {sorted(set(h % 3 for h in h_vals))}")
print(f"H mod 5: distribution = {sorted(set(h % 5 for h in h_vals))}")
print(f"H mod 7: distribution = {sorted(set(h % 7 for h in h_vals))}")

# ─────────────────────────────────────────────────────────────────────
# PART 6: The (3,5)-packing at n=8 — detailed analysis
# ─────────────────────────────────────────────────────────────────────
banner("PART 6: (3,5)-PACKING AT n=8")

print("At n=8 = 2³, the first cross-level (3,5) packing appears.")
print("A disjoint (3-cycle, 5-cycle) pair uses 3+5=8 = ALL vertices.")
print("So at n=8, a (3,5) pair is a PARTITION of the vertex set.\n")

# For n=8 we can't enumerate all 2^28 tournaments, sample instead
random.seed(42)
n = 8
num_samples = 500

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    count += 1
    return count

def has_5cycle_on_vertices(A, verts):
    """Check if there's a directed 5-cycle on exactly these 5 vertices."""
    for perm in permutations(verts):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            return True
    return False

def count_35_pairs(A, n):
    """Count vertex-disjoint (3-cycle, 5-cycle) pairs."""
    count = 0
    # Choose 3 vertices for the 3-cycle
    for triple in combinations(range(n), 3):
        i, j, k = triple
        has3 = False
        if A[i][j] and A[j][k] and A[k][i]:
            has3 = True
        if A[i][k] and A[k][j] and A[j][i]:
            has3 = True
        if has3:
            remaining = [v for v in range(n) if v not in triple]
            if len(remaining) == 5:
                if has_5cycle_on_vertices(A, remaining):
                    count += 1
    return count

print(f"Sampling {num_samples} random tournaments at n={n}...")
import time
t0 = time.time()

has_35_pair = 0
pair_counts = []
for trial in range(num_samples):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    
    c35 = count_35_pairs(A, n)
    pair_counts.append(c35)
    if c35 > 0:
        has_35_pair += 1
    
    if trial % 100 == 0:
        print(f"  trial {trial}: {time.time()-t0:.1f}s")

print(f"Done: {time.time()-t0:.1f}s")
print(f"\n(3,5)-pair statistics at n=8:")
print(f"  Tournaments with ≥1 (3,5)-pair: {has_35_pair}/{num_samples} ({has_35_pair/num_samples:.1%})")
print(f"  Mean (3,5)-pair count: {np.mean(pair_counts):.2f}")
print(f"  Max (3,5)-pair count: {max(pair_counts)}")
print(f"  Distribution of pair counts:")
from collections import Counter
pc_dist = Counter(pair_counts)
for c in sorted(pc_dist.keys()):
    print(f"    {c} pairs: {pc_dist[c]} tournaments ({pc_dist[c]/num_samples:.1%})")

# ─────────────────────────────────────────────────────────────────────
# PART 7: 5 in the master polynomial at x=2
# ─────────────────────────────────────────────────────────────────────
banner("PART 7: MASTER POLYNOMIAL SPECIALIZATIONS AT x=2")

print("Master polynomial: z^{k+1} - (1+x)z^k + x^k = 0")
print("At x=2: z^{k+1} - 3z^k + 2^k = 0")
print()

# For each k, find the Q_k polynomial (after removing z=x=2 root)
# Q_k(z) = (z^{k+1} - 3z^k + 2^k) / (z - 2)
for k in range(2, 8):
    # Polynomial division: (z^{k+1} - 3z^k + 2^k) / (z-2)
    # Coefficients of z^{k+1} - 3z^k + 0·z^{k-1} + ... + 2^k
    poly = [0] * (k+2)
    poly[0] = 1  # z^{k+1}
    poly[1] = -3  # z^k
    poly[k+1] = 2**k  # constant
    
    # Synthetic division by (z-2)
    quotient = [0] * (k+1)
    quotient[0] = poly[0]
    for i in range(1, k+1):
        quotient[i] = poly[i] + 2 * quotient[i-1]
    remainder = poly[k+1] + 2 * quotient[k]
    
    # Find roots of Q_k
    q_roots = np.roots(quotient)
    dom_root = max(r.real for r in q_roots if abs(r.imag) < 1e-8)
    
    q_coeffs_str = " + ".join(f"{int(c)}z^{k-i}" if c != 0 else "" for i, c in enumerate(quotient))
    print(f"k={k}: Q_{k} coefficients = {[int(round(c)) for c in quotient]}")
    print(f"     remainder = {int(round(remainder))} (should be 0)")
    print(f"     dom root of Q_{k} = {dom_root:.10f} (→ 3 as k→∞)")
    
    # Check: does 5 appear in Q_k?
    has_5 = any(abs(abs(c) - 5) < 0.5 for c in quotient)
    if has_5:
        print(f"     *** Q_{k} has coefficient ±5! ***")
    print()

# ─────────────────────────────────────────────────────────────────────
# PART 8: 5-cycle contribution to H — the exact formula
# ─────────────────────────────────────────────────────────────────────
banner("PART 8: 5-CYCLE CONTRIBUTION TO H AT n=5")

print("At n=5: H = 1 + 2α₁ where α₁ = dc3 + dc5")
print("(No α₂ possible since 3+3=6 > 5)")
print()

# Compute all n=5 tournaments with full cycle data
results = []
for idx in range(1024):
    A = adj_matrix(5, idx)
    h = count_ham_paths(A, 5)
    
    # Count directed 3-cycles
    dc3 = 0
    for triple in combinations(range(5), 3):
        i, j, k = triple
        if A[i][j] and A[j][k] and A[k][i]: dc3 += 1
        if A[i][k] and A[k][j] and A[j][i]: dc3 += 1
    
    # Count directed 5-cycles (Hamiltonian cycles)
    dc5 = 0
    for perm in permutations(range(5)):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            dc5 += 1
    dc5 //= 5  # each 5-cycle counted 5 times (cyclic)
    
    alpha1 = dc3 + dc5
    h_pred = 1 + 2*alpha1
    results.append((idx, h, dc3, dc5, alpha1, h_pred))

# Check the formula
mismatches = sum(1 for r in results if r[1] != r[5])
print(f"H = 1 + 2(dc3 + dc5) verified: {1024 - mismatches}/1024 correct")
if mismatches > 0:
    print(f"  (Mismatches at: {[(r[0], r[1], r[5]) for r in results if r[1] != r[5]][:5]})")
print()

# Contribution analysis
print("Relative contribution of 5-cycles vs 3-cycles to α₁:")
total_dc3 = sum(r[2] for r in results)
total_dc5 = sum(r[3] for r in results)
total_alpha1 = sum(r[4] for r in results)
print(f"  Total dc3 across all tournaments: {total_dc3}")
print(f"  Total dc5 across all tournaments: {total_dc5}")
print(f"  Total α₁ = dc3+dc5: {total_alpha1}")
print(f"  dc5/α₁ = {total_dc5/total_alpha1:.4f}")
print(f"  dc3/α₁ = {total_dc3/total_alpha1:.4f}")
print(f"  dc5/dc3 = {total_dc5/total_dc3:.4f}")
print()

# What fraction of H comes from 5-cycles?
print("H contribution analysis:")
print(f"  Mean H = {np.mean([r[1] for r in results]):.4f}")
print(f"  Mean contribution from dc3: 2·mean(dc3) = {2*np.mean([r[2] for r in results]):.4f}")
print(f"  Mean contribution from dc5: 2·mean(dc5) = {2*np.mean([r[3] for r in results]):.4f}")
print(f"  Fraction of (H-1) from 5-cycles: {2*total_dc5 / sum(r[1]-1 for r in results if r[1]>1):.4f}")

# ─────────────────────────────────────────────────────────────────────
# PART 9: The 5-fold symmetry group
# ─────────────────────────────────────────────────────────────────────
banner("PART 9: Z₅ AND CYCLIC TOURNAMENT STRUCTURE")

print("The unique tournament on Z₅ with arc i→j iff j-i ∈ {1,2} (mod 5):")
A_Z5 = [[0]*5 for _ in range(5)]
QR = {1, 2}  # quadratic residues mod 5: 1²=1, 2²=4≡-1, but for a tournament we just pick
# Actually QR mod 5 = {1, 4}. Let's use {1, 2} as the "forward" set
for i in range(5):
    for j in range(5):
        if i != j and (j - i) % 5 in QR:
            A_Z5[i][j] = 1

print("Adjacency matrix:")
for row in A_Z5:
    print(f"  {row}")

h_z5 = count_ham_paths(A_Z5, 5)
print(f"\nH(Z₅ circulant) = {h_z5}")

# Count cycles
dc3_z5 = 0
for triple in combinations(range(5), 3):
    i, j, k = triple
    if A_Z5[i][j] and A_Z5[j][k] and A_Z5[k][i]: dc3_z5 += 1
    if A_Z5[i][k] and A_Z5[k][j] and A_Z5[j][i]: dc3_z5 += 1

dc5_z5 = 0
for perm in permutations(range(5)):
    if all(A_Z5[perm[i]][perm[(i+1)%5]] for i in range(5)):
        dc5_z5 += 1
dc5_z5 //= 5

print(f"dc3 = {dc3_z5}, dc5 = {dc5_z5}")
print(f"α₁ = {dc3_z5 + dc5_z5}")
print(f"H = 1 + 2·{dc3_z5 + dc5_z5} = {1 + 2*(dc3_z5+dc5_z5)}")
print()

# The QR tournament on Z₅
print("The QUADRATIC RESIDUE tournament on Z₅:")
print("QR(5) = {1, 4} = {1, -1} mod 5")
A_QR5 = [[0]*5 for _ in range(5)]
QR5 = {1, 4}
for i in range(5):
    for j in range(5):
        if i != j and (j - i) % 5 in QR5:
            A_QR5[i][j] = 1

print("Adjacency matrix:")
for row in A_QR5:
    print(f"  {row}")

h_qr5 = count_ham_paths(A_QR5, 5)
dc3_qr5 = 0
for triple in combinations(range(5), 3):
    i, j, k = triple
    if A_QR5[i][j] and A_QR5[j][k] and A_QR5[k][i]: dc3_qr5 += 1
    if A_QR5[i][k] and A_QR5[k][j] and A_QR5[j][i]: dc3_qr5 += 1

dc5_qr5 = 0
for perm in permutations(range(5)):
    if all(A_QR5[perm[i]][perm[(i+1)%5]] for i in range(5)):
        dc5_qr5 += 1
dc5_qr5 //= 5

print(f"\nH(QR₅) = {h_qr5}")
print(f"dc3 = {dc3_qr5}, dc5 = {dc5_qr5}")
print(f"α₁ = {dc3_qr5 + dc5_qr5}")

# ─────────────────────────────────────────────────────────────────────
# PART 10: Synthesis — 5 as the discriminant modulus
# ─────────────────────────────────────────────────────────────────────
banner("PART 10: SYNTHESIS — 5 AS DISCRIMINANT MODULUS")

print("THEOREM (informal): 5 is the discriminant modulus of the")
print("tournament recurrence tower.\n")

print("Evidence:")
print("1. The second-order recurrence f(n) = f(n-1) + x·f(n-2)")
print("   has discriminant Δ = 1+4x.")
print("   At x=1: Δ = 5 (the Fibonacci case)")
print("   At x=2: Δ = 9 = 3² (the Jacobsthal/tournament case)")
print()
print("2. The transition from irrational (√5) to rational (3)")
print("   discriminant happens at x = 1 → x = 2.")
print("   This is the SAME transition as Fibonacci → Jacobsthal.")
print()
print("3. The integer-discriminant points are x = k(k-1):")
print("   x=0: Δ=1, x=2: Δ=9, x=6: Δ=25, x=12: Δ=49, ...")
print("   The discriminants are: 1², 3², 5², 7², 9², 11², ...")
print("   AND 5 APPEARS IN THIS SEQUENCE (at x=6, Δ=25=5²).")
print()
print("4. The Q(√5) field connects:")
print("   x=1 (Fibonacci, Δ=5) ↔ x=11 (Δ=45=9·5) ↔ x=31 (Δ=125=25·5)")
print("   These x-values are (5m²-1)/4 for m=1,3,5,...")
print()
print("5. At x=11: roots = (1±3√5)/2 = {3φ-1, 2-3φ}")
print("   So x=11 roots are EXACTLY 3× Fibonacci roots, shifted by -1.")
print("   The 'decimal pair' (10,11) connects to the golden ratio!")
print()
print("6. The J₆ sequence (x=6 level, denominator 5) generates 7:")
print("   J₆(3) = (3³-(-2)³)/5 = (27+8)/5 = 35/5 = 7")
print("   So 5 GENERATES the (7,8) transition threshold.")
print()
print("CONCLUSION: 5 is the algebraic number that makes the")
print("Fibonacci→Jacobsthal transition WORK. It is the discriminant")
print("of the base case (x=1), and through Q(√5) it threads through")
print("the entire hierarchy: 1→5→11→31→... All levels where the")
print("golden ratio echoes in the tournament structure.")

print("\nDone.")
