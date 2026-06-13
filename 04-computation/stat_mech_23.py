#!/usr/bin/env python3
"""
Statistical mechanics of tournaments — partition functions and phase transitions.
opus-2026-03-14-S85

ISING-TYPE MODEL ON TOURNAMENTS:
- Arc variables σ_{ij} ∈ {+1, -1} (i<j), σ=+1 means i→j, σ=-1 means j→i.
- Energy E(T) = -H(T) (minimize energy = maximize H).
- Partition function Z(β) = Σ_T exp(-β E(T)) = Σ_T exp(β H(T)).
- Free energy F(β) = -log Z(β) / β.
- Phase transitions: does F(β) have singularities?

CONNECTION TO POTTS MODEL:
- Tournament = orientation of complete graph K_n.
- H(T) = # Hamiltonian paths = sum over permutations of arc products.
- This is like a Potts model with n! states, but restricted to tournament structure.

INFORMATION GEOMETRY:
- The exponential family {P_β(T) = exp(β H(T))/Z(β)} defines a
  1-dimensional statistical manifold.
- Fisher information I(β) = Var_β(H) = d²F/dβ².
- Critical β: where I(β) is maximized (phase transition).
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
import sys
import numpy as np

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

# ============================================================
# Part 1: Partition Function Z(β) and Free Energy
# ============================================================
print("=" * 70)
print("PART 1: PARTITION FUNCTION Z(β)")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    # Compute all H values
    H_all = []
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_all.append(compute_H_dp(adj, n))

    H_dist = Counter(H_all)
    print(f"\nn={n}: H distribution: {dict(sorted(H_dist.items()))}")

    # Compute Z(β) for various β
    betas = np.linspace(-2, 2, 41)
    print(f"\n  {'β':>6s} {'Z(β)':>15s} {'<H>_β':>10s} {'Var(H)_β':>10s} {'F(β)':>12s}")
    for beta in betas:
        # Z(β) = Σ_h #(H=h) * exp(β*h)
        Z = sum(count * math.exp(beta * h) for h, count in H_dist.items())
        # <H>_β = Σ h * P_β(h)
        mean_H = sum(h * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        # <H²>_β
        mean_H2 = sum(h**2 * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        var_H = mean_H2 - mean_H**2
        F = -math.log(Z) / beta if beta != 0 else 0

        if abs(beta) < 0.01 or abs(beta - 1) < 0.01 or abs(beta + 1) < 0.01 or abs(beta - 0.5) < 0.01:
            print(f"  {beta:6.2f} {Z:15.2f} {mean_H:10.4f} {var_H:10.4f} {F:12.4f}")

    # Detailed: β=1 (Boltzmann temperature)
    beta = 1.0
    Z = sum(count * math.exp(beta * h) for h, count in H_dist.items())
    probs = {h: count * math.exp(beta * h) / Z for h, count in H_dist.items()}
    print(f"\n  β=1 Boltzmann distribution:")
    for h in sorted(probs.keys()):
        print(f"    P(H={h}) = {probs[h]:.6f}")

# ============================================================
# Part 2: Specific Heat (β-derivative of Energy)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: SPECIFIC HEAT C(β) = β² Var_β(H)")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_all = []
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_all.append(compute_H_dp(adj, n))

    H_dist = Counter(H_all)

    betas = np.linspace(0.01, 3, 60)
    print(f"\nn={n}: Specific heat:")
    max_C = 0
    max_beta = 0
    for beta in betas:
        Z = sum(count * math.exp(beta * h) for h, count in H_dist.items())
        mean_H = sum(h * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        mean_H2 = sum(h**2 * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        var_H = mean_H2 - mean_H**2
        C = beta**2 * var_H

        if C > max_C:
            max_C = C
            max_beta = beta

    print(f"  Max specific heat: C = {max_C:.4f} at β = {max_beta:.3f}")
    print(f"  Critical temperature T_c = 1/β_c = {1/max_beta:.4f}")

    # Fine scan near critical point
    betas_fine = np.linspace(max(0.01, max_beta - 0.5), max_beta + 0.5, 100)
    print(f"\n  Fine scan near β_c:")
    for beta in betas_fine:
        Z = sum(count * math.exp(beta * h) for h, count in H_dist.items())
        mean_H = sum(h * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        mean_H2 = sum(h**2 * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        C = beta**2 * (mean_H2 - mean_H**2)
        if abs(beta - max_beta) < 0.02:
            print(f"    β={beta:.4f}: C={C:.4f}, <H>={mean_H:.4f}")

# ============================================================
# Part 3: Entropy of Tournament Distribution
# ============================================================
print("\n" + "=" * 70)
print("PART 3: ENTROPY S(β) = log Z + β<H>")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_all = []
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_all.append(compute_H_dp(adj, n))

    H_dist = Counter(H_all)

    print(f"\nn={n}: Entropy at various β:")
    for beta in [0, 0.1, 0.5, 1.0, 2.0, 5.0]:
        if beta == 0:
            S = math.log(N)
            print(f"  β=0: S = log(2^m) = {S:.4f} = {m} log 2 = {m * math.log(2):.4f}")
            continue

        Z = sum(count * math.exp(beta * h) for h, count in H_dist.items())
        mean_H = sum(h * count * math.exp(beta * h) for h, count in H_dist.items()) / Z
        S = math.log(Z) + beta * mean_H  # Wait, S = log Z - β<E> = log Z + β<H>
        # Actually S = -Σ P log P
        S_exact = 0
        for h, count in H_dist.items():
            p = count * math.exp(beta * h) / Z
            if p > 0:
                S_exact -= p * math.log(p)

        print(f"  β={beta:.1f}: S = {S_exact:.4f} (max = {math.log(N):.4f})")

# ============================================================
# Part 4: Lee-Yang Zeros of Partition Function
# ============================================================
print("\n" + "=" * 70)
print("PART 4: LEE-YANG ZEROS")
print("=" * 70)

# Z(β) = Σ_h N_h * exp(βh)
# As a function of z = exp(β), this is a polynomial:
# Z = Σ_h N_h * z^h
# The Lee-Yang theorem asks: where are the zeros of Z(z)?

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_all = []
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_all.append(compute_H_dp(adj, n))

    H_dist = Counter(H_all)
    max_H = max(H_dist.keys())

    # Build polynomial coefficients: Z(z) = Σ_h N_h * z^h
    coeffs = [0] * (max_H + 1)
    for h, count in H_dist.items():
        coeffs[h] = count

    # Find zeros
    # Reverse coefficients for np.roots (highest degree first)
    np_coeffs = list(reversed(coeffs))
    while np_coeffs and np_coeffs[0] == 0:
        np_coeffs.pop(0)

    zeros = np.roots(np_coeffs)

    # Analyze zeros
    print(f"\nn={n}: Lee-Yang zeros of Z(z):")
    print(f"  Degree of Z: {max_H}")
    print(f"  Number of zeros: {len(zeros)}")

    # How many are on/near the unit circle?
    on_circle = sum(1 for z in zeros if abs(abs(z) - 1) < 0.01)
    near_circle = sum(1 for z in zeros if abs(abs(z) - 1) < 0.1)
    print(f"  Zeros on unit circle (|z|≈1): {on_circle}")
    print(f"  Zeros near unit circle (||z|-1|<0.1): {near_circle}")

    # Negative real zeros
    neg_real = sum(1 for z in zeros if z.real < 0 and abs(z.imag) < 0.01)
    print(f"  Negative real zeros: {neg_real}")

    # Print zeros
    for z in sorted(zeros, key=lambda z: (abs(z), z.real)):
        if abs(z.imag) < 0.001:
            print(f"    z = {z.real:.6f} (|z|={abs(z):.4f})")
        else:
            print(f"    z = {z.real:.6f} + {z.imag:.6f}i (|z|={abs(z):.4f})")

# ============================================================
# Part 5: Boltzmann Machine — H-Optimal Arc Selection
# ============================================================
print("\n" + "=" * 70)
print("PART 5: SIMULATED ANNEALING FOR MAX-H")
print("=" * 70)

import random
random.seed(42)

# Use simulated annealing to find max-H tournaments at larger n
for n in [7, 8, 9]:
    m = n * (n - 1) // 2

    # Start from random tournament
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

    H = compute_H_dp(adj, n)
    best_H = H
    best_adj = [row[:] for row in adj]

    # Annealing schedule
    T0 = 10.0
    T_final = 0.1
    n_steps = 5000
    for step in range(n_steps):
        T = T0 * (T_final / T0) ** (step / n_steps)

        # Flip a random arc
        i = random.randint(0, n-2)
        j = random.randint(i+1, n-1)
        # Flip
        adj[i][j], adj[j][i] = adj[j][i], adj[i][j]
        H_new = compute_H_dp(adj, n)

        # Accept?
        dH = H_new - H
        if dH > 0 or random.random() < math.exp(dH / T):
            H = H_new
            if H > best_H:
                best_H = H
                best_adj = [row[:] for row in adj]
        else:
            # Reject — flip back
            adj[i][j], adj[j][i] = adj[j][i], adj[i][j]

    mean_H = math.factorial(n) / 2**(n-1)
    print(f"\nn={n}: Simulated annealing max H = {best_H}")
    print(f"  mean H = {mean_H:.2f}")
    print(f"  H/mean = {best_H/mean_H:.4f}")
    scores = sorted(sum(best_adj[i][j] for j in range(n) if j != i) for i in range(n))
    print(f"  Scores: {scores}")

# ============================================================
# Part 6: Tournament Ising Model — Correlation Functions
# ============================================================
print("\n" + "=" * 70)
print("PART 6: CORRELATION FUNCTIONS <σ_ij σ_kl>")
print("=" * 70)

# In the Boltzmann distribution P_β(T) = exp(βH(T))/Z(β),
# what are the arc-arc correlations?
# <σ_{ij}> = expected direction of arc (i,j)
# <σ_{ij} σ_{kl}> = correlation between two arcs

n = 4
m = n * (n - 1) // 2
N = 1 << m
arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

H_all = []
for bits in range(N):
    adj = get_tournament(n, bits)
    H_all.append(compute_H_dp(adj, n))

for beta in [0, 0.5, 1.0, 2.0]:
    # Compute Boltzmann weights
    Z = sum(math.exp(beta * h) for h in H_all)
    probs = [math.exp(beta * h) / Z for h in H_all]

    # Arc expectations
    arc_means = [0.0] * m
    for bits in range(N):
        for k in range(m):
            sigma = 2 * ((bits >> k) & 1) - 1  # +1 or -1
            arc_means[k] += probs[bits] * sigma

    # Arc-arc correlations
    corr_matrix = [[0.0]*m for _ in range(m)]
    for bits in range(N):
        for k1 in range(m):
            for k2 in range(m):
                s1 = 2 * ((bits >> k1) & 1) - 1
                s2 = 2 * ((bits >> k2) & 1) - 1
                corr_matrix[k1][k2] += probs[bits] * s1 * s2

    print(f"\nn=4, β={beta}:")
    print(f"  Arc means: {[f'{a:.4f}' for a in arc_means]}")
    print(f"  All means zero: {all(abs(a) < 1e-10 for a in arc_means)}")

    # Off-diagonal correlations
    off_diag = []
    for k1 in range(m):
        for k2 in range(k1+1, m):
            connected = corr_matrix[k1][k2] - arc_means[k1] * arc_means[k2]
            off_diag.append(connected)

    mean_corr = sum(off_diag) / len(off_diag)
    max_corr = max(abs(c) for c in off_diag)
    print(f"  Mean connected correlation: {mean_corr:.6f}")
    print(f"  Max |connected correlation|: {max_corr:.6f}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — STATISTICAL MECHANICS")
print("=" * 70)
print("""
KEY FINDINGS:
1. PARTITION FUNCTION: Z(β) = Σ N_h exp(βh) is a polynomial in z=exp(β).
2. SPECIFIC HEAT: Has a maximum at critical β_c, indicating a crossover
   (finite-size analog of phase transition).
3. ENTROPY: Goes from log(2^m) at β=0 (uniform) to 0 as β→∞ (max-H only).
4. LEE-YANG ZEROS: Zeros of Z(z) have specific patterns.
   Not all on unit circle (unlike true Ising).
5. SIMULATED ANNEALING: Finds near-optimal H for larger n.
6. CORRELATION FUNCTIONS: At β=0, all arc correlations are zero.
   At β>0, arcs become correlated via the H landscape.
""")
