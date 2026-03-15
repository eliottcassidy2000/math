#!/usr/bin/env python3
"""
grand_reframe_S72b.py — opus-2026-03-15-S72b

GETTING A HANDLE ON EVERYTHING: THE UNIFIED PICTURE
=====================================================

This project has five interlocking programs. This script reframes them
as facets of ONE underlying structure, then explores the connections.

THE FIVE PROGRAMS:
  1. OCF: H(T) = I(Omega(T), 2)  — counting ↔ topology of cycles
  2. Walsh-Fourier: H_hat[S] decomposition — spectral structure
  3. Signed permanent S(T) — universal congruences, 2-adic tower
  4. GLMY path homology: β₀=1, β₁∈{0,1}, β₂=0, β₃ onset — holes in path complex
  5. Paley/circulant structure: QR tournaments, eigenspace decomposition — number theory

THE UNIFYING IDEA:
  A tournament encodes a BINARY CHOICE on each edge pair.
  These choices create LOCAL structure (3-cycles, scores) and GLOBAL structure
  (Hamiltonian paths, homological holes, orientability).

  The five programs measure different aspects:
  - OCF: how local cycles ADD UP to global path count
  - Walsh: how each binary choice INDEPENDENTLY contributes
  - Signed perm: what's UNIVERSAL vs tournament-specific
  - Homology: what TOPOLOGY the directed paths create
  - Paley: what happens when choices have NUMBER-THEORETIC structure

THE KEY DISCOVERIES:
  - β₂=0 (PROVED): tournament completeness kills 2-holes
  - β₁ = orientability obstruction (NEW, this session): the twist
  - H=7 and H=21 are permanent gaps (PROVED): poisoning graph
  - Walsh proof of OCF (THM-077): elementary, bypasses P-partitions
  - φ·τ ≈ 3: Fibonacci × tribonacci ≈ ternary weight limit
  - F₇ = T₇ = 13: Fibonacci = tribonacci at n=7 (first β₃ onset)

WHAT'S MISSING — THE GAPS IN UNDERSTANDING:
  1. No unified theory connecting β₁ to Walsh spectrum
  2. No algebraic proof of β₁·β₃=0
  3. No formula for β₁ or β₃ in terms of known invariants
  4. The "twist" (β₁=1 with full coverage) has no clean characterization
  5. H=7 gap and H=21 gap: is there a TOPOLOGICAL explanation?
  6. Paley maximizes H — but WHY? (only verified, not proved)

This script explores the connections between these gaps.
"""

import numpy as np
from math import factorial, comb, sqrt, gcd
import sys
import time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 32749

def adj_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def rank_mod(M, p):
    M = M.copy() % p
    nrows, ncols = M.shape
    r = 0
    for col in range(ncols):
        nz = np.where(M[r:, col] != 0)[0]
        if len(nz) == 0:
            continue
        pivot = nz[0] + r
        if pivot != r:
            M[[r, pivot]] = M[[pivot, r]]
        inv = pow(int(M[r, col]), p - 2, p)
        M[r] = M[r] * inv % p
        factors = M[:, col].copy()
        factors[r] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
        r += 1
    return r

def compute_H(A, n):
    """DP Hamiltonian path count."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            count = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def compute_beta1(A, n):
    """Fast β₁ via rank formula."""
    edges = []
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    ne = len(edges)

    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[i][k]:
                    triples.append((i, j, k))
    nt = len(triples)

    d1 = np.zeros((n, ne), dtype=np.int64)
    for idx, (i, j) in enumerate(edges):
        d1[j, idx] += 1
        d1[i, idx] -= 1

    d2 = np.zeros((ne, nt), dtype=np.int64)
    for t_idx, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx: d2[edge_idx[(j, k)], t_idx] += 1
        if (i, k) in edge_idx: d2[edge_idx[(i, k)], t_idx] -= 1
        if (i, j) in edge_idx: d2[edge_idx[(i, j)], t_idx] += 1

    rk1 = rank_mod(d1 % PRIME, PRIME)
    rk2 = rank_mod(d2 % PRIME, PRIME)
    return ne - rk1 - rk2

def t3_count(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += A[i][j]*A[j][k]*A[k][i] + A[j][i]*A[i][k]*A[k][j]
    return count

def score_seq(A, n):
    return tuple(sorted(int(sum(A[i])) for i in range(n)))

print("=" * 72)
print("GRAND REFRAME: THE UNIFIED PICTURE")
print("opus-2026-03-15-S72b")
print("=" * 72)

# =====================================================================
# CONNECTION 1: β₁ and H — what's the relationship?
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 1: β₁ AND H")
print("=" * 60)

n = 5
total = 1 << (n*(n-1)//2)

h_by_beta1 = {0: [], 1: []}
t3_by_beta1 = {0: [], 1: []}
score_by_beta1 = {0: {}, 1: {}}

for bits in range(total):
    A = adj_from_bits(n, bits)
    H = compute_H(A, n)
    b1 = compute_beta1(A, n)
    t3 = t3_count(A, n)
    s = score_seq(A, n)

    h_by_beta1[b1].append(H)
    t3_by_beta1[b1].append(t3)
    score_by_beta1[b1][s] = score_by_beta1[b1].get(s, 0) + 1

print(f"\nn=5 ({total} tournaments):")
for b1 in [0, 1]:
    vals = h_by_beta1[b1]
    t3s = t3_by_beta1[b1]
    print(f"  β₁={b1}: {len(vals)} tournaments")
    print(f"    H: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")
    print(f"    t3: min={min(t3s)}, max={max(t3s)}, mean={np.mean(t3s):.1f}")
    print(f"    scores: {score_by_beta1[b1]}")

# Is there a formula? H = f(t3, beta1)?
print(f"\n  H vs t3 by β₁:")
for b1 in [0, 1]:
    # Group by t3
    ht3 = {}
    for bits in range(total):
        A = adj_from_bits(n, bits)
        H = compute_H(A, n)
        b = compute_beta1(A, n)
        t = t3_count(A, n)
        if b == b1:
            if t not in ht3:
                ht3[t] = set()
            ht3[t].add(H)
    print(f"  β₁={b1}:")
    for t in sorted(ht3.keys()):
        print(f"    t3={t}: H ∈ {sorted(ht3[t])}")

# =====================================================================
# CONNECTION 2: β₁ and the signed permanent S(T)
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 2: β₁ AND SIGNED PERMANENT")
print("=" * 60)

n = 5
# S(T) = sum_P prod B[P_i, P_{i+1}] where B = 2A - J
s_by_beta1 = {0: [], 1: []}

for bits in range(total):
    A = adj_from_bits(n, bits)
    B = 2 * A - np.ones((n, n), dtype=int)
    np.fill_diagonal(B, 0)

    # Compute S(T)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + dp[(mask, v)] * B[v][u]

    full = (1 << n) - 1
    S = sum(dp.get((full, v), 0) for v in range(n))
    b1 = compute_beta1(A, n)
    s_by_beta1[b1].append(S)

for b1 in [0, 1]:
    vals = s_by_beta1[b1]
    from collections import Counter
    dist = Counter(vals)
    print(f"  β₁={b1}: S(T) distribution = {dict(sorted(dist.items()))}")

# c0 = S / 2^{n-1} = S / 16
print(f"\n  c₀ = S/16 distribution:")
for b1 in [0, 1]:
    c0s = [s / 16 for s in s_by_beta1[b1]]
    from collections import Counter
    dist = Counter(c0s)
    print(f"  β₁={b1}: c₀ ∈ {sorted(set(c0s))}")

# =====================================================================
# CONNECTION 3: β₁ and Walsh spectrum
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 3: β₁ AND WALSH SPECTRUM")
print("=" * 60)

# At n=5, compute the degree-2 Walsh energy for each tournament
# Degree-2 Walsh = sum of H_hat[S]^2 over |S|=2

# First, compute Walsh transform of H for ALL tournaments
n = 5
m = n*(n-1)//2  # 10

# Compute H for all tournaments
H_all = np.zeros(1 << m, dtype=np.int64)
for bits in range(1 << m):
    A = adj_from_bits(n, bits)
    H_all[bits] = compute_H(A, n)

# Walsh-Hadamard transform
H_hat = H_all.copy().astype(float)
for i in range(m):
    half = 1 << i
    for j in range(0, 1 << m, 2 * half):
        for k in range(half):
            a, b = H_hat[j + k], H_hat[j + k + half]
            H_hat[j + k] = a + b
            H_hat[j + k + half] = a - b
H_hat /= (1 << m)

# Compute degree-2 Walsh energy per tournament
# Actually we need the Walsh coefficients for each tournament, not the global average
# H_hat[S] = (1/2^m) sum_T (-1)^{T·S} H(T) is a single number
# What we want is: for each tournament T, its Walsh "signature"

# Better: the Walsh spectrum tells us about the AVERAGE tournament
# For individual tournaments, we should look at the local Walsh contributions

# Actually let me think about this differently.
# At n=5, m=10 bits. For each tournament T (bit string t),
# H(t) = sum_S H_hat[S] * (-1)^{t·chi_S} (Walsh reconstruction)
# The "Walsh energy at degree d" for tournament t is:
# E_d(t) = sum_{|S|=d} H_hat[S] * (-1)^{t·chi_S}

# Compute this for each tournament
from itertools import combinations

# Get subsets of each size
edges = []
for i in range(n):
    for j in range(i+1, n):
        edges.append((i,j))
# edges[k] is the k-th edge pair, bit k of tournament encoding

walsh_deg2_energy = np.zeros(1 << m)
walsh_deg4_energy = np.zeros(1 << m)

for S_mask in range(1 << m):
    popcount = bin(S_mask).count('1')
    if popcount == 2:
        coeff = H_hat[S_mask]
        for bits in range(1 << m):
            sign = 1 - 2 * (bin(bits & S_mask).count('1') % 2)
            walsh_deg2_energy[bits] += coeff * sign
    elif popcount == 4:
        coeff = H_hat[S_mask]
        for bits in range(1 << m):
            sign = 1 - 2 * (bin(bits & S_mask).count('1') % 2)
            walsh_deg4_energy[bits] += coeff * sign

# Now compare β₁ with Walsh energies
wd2_by_b1 = {0: [], 1: []}
wd4_by_b1 = {0: [], 1: []}

for bits in range(total):
    A = adj_from_bits(n, bits)
    b1 = compute_beta1(A, n)
    wd2_by_b1[b1].append(walsh_deg2_energy[bits])
    wd4_by_b1[b1].append(walsh_deg4_energy[bits])

print(f"\n  Walsh degree-2 energy by β₁:")
for b1 in [0, 1]:
    vals = wd2_by_b1[b1]
    print(f"    β₁={b1}: mean={np.mean(vals):.4f}, std={np.std(vals):.4f}, "
          f"range=[{min(vals):.2f}, {max(vals):.2f}]")

print(f"\n  Walsh degree-4 energy by β₁:")
for b1 in [0, 1]:
    vals = wd4_by_b1[b1]
    print(f"    β₁={b1}: mean={np.mean(vals):.4f}, std={np.std(vals):.4f}, "
          f"range=[{min(vals):.2f}, {max(vals):.2f}]")

# =====================================================================
# CONNECTION 4: β₁ and OCF decomposition α₁, α₂
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 4: β₁ AND OCF DECOMPOSITION")
print("=" * 60)

# H = 1 + 2*α₁ + 4*α₂ + ...
# At n=5: α₁ = t3 (3-cycles), α₂ = i2 (#independent pairs of 3-cycles)
# Is β₁ determined by (α₁, α₂)?

alpha_by_b1 = {0: {}, 1: {}}
for bits in range(total):
    A = adj_from_bits(n, bits)
    b1 = compute_beta1(A, n)
    t3 = t3_count(A, n)

    # Count independent pairs of 3-cycles (alpha_2)
    cycles3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i] + A[j][i]*A[i][k]*A[k][j] > 0:
                    cycles3.append((i, j, k))

    alpha2 = 0
    for i in range(len(cycles3)):
        for j in range(i+1, len(cycles3)):
            s1 = set(cycles3[i])
            s2 = set(cycles3[j])
            if not s1 & s2:
                alpha2 += 1

    key = (t3, alpha2)
    if key not in alpha_by_b1[b1]:
        alpha_by_b1[b1][key] = 0
    alpha_by_b1[b1][key] += 1

print(f"\n  (t3, α₂) distribution by β₁:")
all_keys = sorted(set(list(alpha_by_b1[0].keys()) + list(alpha_by_b1[1].keys())))
for key in all_keys:
    c0 = alpha_by_b1[0].get(key, 0)
    c1 = alpha_by_b1[1].get(key, 0)
    H_val = 1 + 2*key[0] + 4*key[1]
    print(f"    (t3={key[0]}, α₂={key[1]}): H={H_val}, "
          f"β₁=0: {c0}, β₁=1: {c1}" +
          (" ← AMBIGUOUS" if c0 > 0 and c1 > 0 else ""))

# =====================================================================
# CONNECTION 5: H=7 gap and β₁
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 5: H=7 GAP AND β₁")
print("=" * 60)

# H=7 requires α₁=3, α₂=0 → t3=3, no disjoint pair of 3-cycles
# At n=5: t3=3 gives H=7? Let's check
print(f"\n  n=5: tournaments with t3=3:")
count_t3_3 = 0
H_vals_t3_3 = {}
for bits in range(total):
    A = adj_from_bits(n, bits)
    t3 = t3_count(A, n)
    if t3 == 3:
        H = compute_H(A, n)
        b1 = compute_beta1(A, n)
        key = (H, b1)
        H_vals_t3_3[key] = H_vals_t3_3.get(key, 0) + 1
        count_t3_3 += 1

print(f"    Total: {count_t3_3}")
for key in sorted(H_vals_t3_3.keys()):
    print(f"    H={key[0]}, β₁={key[1]}: {H_vals_t3_3[key]}")

# H=7 = 1 + 2*3 + 4*0. But at n=5, t3=3 gives H=9 or H=11, not 7!
# Because α₂ ≥ 1 when t3=3 at n=5 (3 disjoint cycles on 5 vertices
# can't avoid sharing... wait, 3 cycles × 3 vertices = 9 > 5)
# So any 3 3-cycles on 5 vertices must share vertices, but they might
# still be disjoint in pairs (Omega graph).
# Let me check: 2 disjoint 3-cycles = 6 vertices. n=5 < 6. So α₂ = 0!

# Then H = 1 + 2*3 = 7. But we see H=9 for some t3=3 tournaments.
# That means OCF gives H = 1 + 2*t3 + 4*α₂ where at n=5, α₁ = t3 + t5.
# At n=5, 5-cycles exist! So α₁ = t3 + t5, not just t3.

# Recompute with full α₁
print(f"\n  Correcting: α₁ = t3 + t5 (all odd cycles)")
for bits in range(total):
    A = adj_from_bits(n, bits)
    t3 = t3_count(A, n)
    H = compute_H(A, n)
    if t3 == 3:
        # Count 5-cycles
        t5 = 0
        from itertools import permutations
        for perm in permutations(range(n)):
            is_cycle = True
            for i in range(5):
                if not A[perm[i]][perm[(i+1) % 5]]:
                    is_cycle = False
                    break
            if is_cycle:
                t5 += 1
        t5 //= 5  # each 5-cycle counted 5 times (cyclic shifts)

        alpha1 = t3 + t5
        # H = 1 + 2*alpha1 + 4*alpha2 + ...
        # So alpha2 = (H - 1 - 2*alpha1) / 4 (if no higher terms)
        remaining = H - 1 - 2*alpha1
        print(f"    bits={bits}: t3={t3}, t5={t5}, α₁={alpha1}, H={H}, "
              f"H-1-2α₁={remaining}, α₂={remaining//4 if remaining%4==0 else '?'}")
        if count_t3_3 > 10:
            break

# =====================================================================
# CONNECTION 6: The tau-phi clock and β₁ onset
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 6: TAU-PHI CLOCK AND BETTI ONSET")
print("=" * 60)

phi = (1 + sqrt(5)) / 2
tau = 1.8392867552141612

print(f"""
  KEY OBSERVATION: F₇ = T₇ = 13.

  At n=7, Fibonacci and tribonacci are EQUAL.
  This is also the first n where β₃ > 0 appears.

  Is this a coincidence? Let's check the Betti onset sizes:
    β₁ first appears at n=3  (3-cycles)
    β₃ first appears at n=6  (but rare; common at n=7)
    β₄ first appears at n=7  (Paley T₇ has β₄=6)
    β₅ first appears at n=8

  Fibonacci numbers:  F₃=2, F₆=8, F₇=13, F₈=21
  Tribonacci numbers: T₃=1, T₆=7, T₇=13, T₈=24

  The CROSSOVER point F_n = T_n occurs at n=7.
  Before n=7: F_n > T_n (Fibonacci dominates → paths dominate cycles)
  After n=7: T_n > F_n (tribonacci dominates → cycles dominate paths)

  This matches the homological picture:
    n < 7: β₁ (1-holes from cycles) is the main phenomenon
    n = 7: transition point — β₃, β₄ emerge
    n > 7: higher Betti numbers proliferate

  The tau-phi clock at n=7:
    phi^7 = {phi**7:.2f}  ≈ F₇ * phi = 13 * 1.618 = {13*phi:.2f}
    tau^7 = {tau**7:.2f}  ≈ T₇ * K = 13 * K
    Ratio: tau^7/phi^7 = (tau/phi)^7 = {(tau/phi)**7:.6f}

  The tau-phi clock crosses over at n=7: after this, tribonacci > Fibonacci.
""")

def fib_mod_period(m):
    a, b = 0, 1
    for i in range(1, 6*m*m):
        a, b = b, (a+b) % m
        if a == 0 and b == 1:
            return i
    return -1

def trib_mod_period(m):
    a, b, c = 0, 0, 1
    for i in range(1, 6*m*m*m):
        a, b, c = b, c, (a+b+c) % m
        if a == 0 and b == 0 and c == 1:
            return i
    return -1

# =====================================================================
# CONNECTION 7: What determines β₁? The complete picture at n=5
# =====================================================================
print("\n" + "=" * 60)
print("CONNECTION 7: COMPLETE DETERMINATION OF β₁ AT n=5")
print("=" * 60)

# From our earlier work:
# β₁=0: 720 tournaments, β₁=1: 304 tournaments
# β₁=1 has two subtypes: twisted (144) and untwisted (160)
# β₁=0 implies full edge coverage (THM-222)
# Untwisted: has uncovered edge (source with out-deg 1, target with out-deg n-2)

# Is β₁ determined by (t3, t5)?
t3t5_by_b1 = {0: {}, 1: {}}
for bits in range(total):
    A = adj_from_bits(n, bits)
    b1 = compute_beta1(A, n)
    t3 = t3_count(A, n)

    # Count 5-cycles
    t5 = 0
    for perm in permutations(range(n)):
        is_cycle = True
        for i in range(5):
            if not A[perm[i]][perm[(i+1) % 5]]:
                is_cycle = False
                break
        if is_cycle:
            t5 += 1
    t5 //= 5

    key = (t3, t5)
    t3t5_by_b1[b1][key] = t3t5_by_b1[b1].get(key, 0) + 1

print(f"\n  (t3, t5) → β₁ at n=5:")
all_keys = sorted(set(list(t3t5_by_b1[0].keys()) + list(t3t5_by_b1[1].keys())))
ambiguous = 0
for key in all_keys:
    c0 = t3t5_by_b1[0].get(key, 0)
    c1 = t3t5_by_b1[1].get(key, 0)
    tag = ""
    if c0 > 0 and c1 > 0:
        tag = " ← AMBIGUOUS"
        ambiguous += 1
    print(f"    (t3={key[0]}, t5={key[1]}): β₁=0: {c0:4d}, β₁=1: {c1:4d}{tag}")

print(f"\n  Ambiguous (t3,t5) classes: {ambiguous}")
if ambiguous == 0:
    print(f"  ★ (t3, t5) DETERMINES β₁ at n=5! ★")
else:
    print(f"  (t3, t5) does NOT determine β₁ at n=5.")

# =====================================================================
# SYNTHESIS
# =====================================================================
print("\n" + "=" * 60)
print("GRAND SYNTHESIS")
print("=" * 60)

print(f"""
  THE UNIFIED PICTURE OF TOURNAMENT THEORY:

  ┌──────────────────────────────────────────────────────────┐
  │                    TOURNAMENT T                          │
  │              (n vertices, C(n,2) binary choices)         │
  └──────────────┬──────────────┬──────────────┬────────────┘
                 │              │              │
          ┌──────▼──────┐ ┌────▼────┐  ┌──────▼──────┐
          │ LOCAL       │ │ GLOBAL  │  │ TOPOLOGICAL │
          │ (subgraphs) │ │ (paths) │  │ (homology)  │
          └──────┬──────┘ └────┬────┘  └──────┬──────┘
                 │              │              │
  t3, t5, ...   │     H(T)    │     β₀,β₁,...│
  score seq     │     M[a,b]  │     Ω dims   │
  α₁, α₂       │     S(T)    │     twist    │
                 │              │              │
          ┌──────▼──────────────▼──────────────▼──────┐
          │              OCF BRIDGE                   │
          │    H(T) = I(Ω(T), 2) = 1 + 2α₁ + 4α₂   │
          │    Links local (cycles) to global (paths) │
          └──────────────────┬────────────────────────┘
                             │
          ┌──────────────────▼────────────────────────┐
          │          WALSH-FOURIER LENS               │
          │  Decomposes everything by binary choices   │
          │  H_hat[S] = ε·2^r·(n-2k)!/2^(n-1)       │
          │  Each S is an even-length path union       │
          └──────────────────┬────────────────────────┘
                             │
          ┌──────────────────▼────────────────────────┐
          │           NUMBER THEORY LENS              │
          │  Paley: QR structure → eigenspace decomp  │
          │  THM-125: n× speedup from circulant       │
          │  φ·τ ≈ 3: Fibonacci × tribonacci ≈ ternary│
          └───────────────────────────────────────────┘

  THE KEY CONNECTIONS DISCOVERED:

  1. β₁ = ORIENTABILITY OBSTRUCTION
     The path complex can be locally but not globally orientable.
     This is invisible to subgraph counts (THM-222, HYP-1578).

  2. β₂ = 0 ALWAYS (THM-108/109)
     Tournament completeness kills all 2-dimensional holes.
     The twin vertex mechanism: every counterexample needs twins.

  3. β₁ · β₃ = 0 (verified through n=7)
     SLACK-DENSITY DUALITY: twisted ↔ can't wrap higher,
     dense ↔ can't twist. The τ-φ clock desynchronization.

  4. H=7 and H=21 are PERMANENT GAPS (PROVED)
     The poisoning graph argument. OCF algebraic constraint.

  5. PALEY MAXIMIZES H at primes p ≡ 3 mod 4 (verified p=3,7,11)
     QR structure creates maximal cycle coverage.

  6. τ/φ has minimal polynomial of degree 6
     The two growth rates are algebraically entangled
     but not simply related.

  WHAT'S NEXT — the deepest open questions:

  A. PROVE β₁·β₃=0 algebraically.
     The slack-density argument is conceptual but not rigorous.
     Needs: translate "twist prevents higher wrapping" into
     a statement about chain complex dimensions.

  B. FIND a formula for β₁ in terms of known invariants.
     At n=5: (t3,t5) determines β₁ (if confirmed above).
     At n=6: (c3v,sc4) gives 99.93% accuracy (HYP-1568).
     But no EXACT formula is known for any n.

  C. CONNECT β₁ to the Walsh spectrum.
     The twist is about SIGNS of triple boundaries.
     Walsh coefficients are about signs of tournament encodings.
     There should be a deep connection.

  D. EXPLAIN why Paley maximizes H.
     The QR structure gives maximal SYMMETRY.
     The eigenspace decomposition (THM-125) compresses computation.
     But why does symmetry → max H?
""")

print("DONE — grand_reframe_S72b.py complete")
