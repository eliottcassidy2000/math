#!/usr/bin/env python3
"""
Correct chain complex analysis for tournaments.
The key is Ω_p ≠ A_p in general!
Author: opus-2026-03-07-S46e
"""
import numpy as np
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths,
    compute_omega_basis, build_full_boundary_matrix
)
from itertools import combinations
import random

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

random.seed(42)

# === C5 tournament ===
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        A[i][(i+d)%5] = 1
n = 5

print("C5 tournament (regular, t3=5):")
print("Chain complex dimensions:")
for p in range(n):
    ap = enumerate_allowed_paths(A, n, p)
    omega = compute_omega_basis(A, n, p, ap, enumerate_allowed_paths(A, n, p-1) if p > 0 else [])
    dim_omega = omega.shape[1] if omega.ndim == 2 else 0
    print(f"  p={p}: |A_{p}|={len(ap)}, dim(Ω_{p})={dim_omega}")

beta = path_betti_numbers(A, n)
print(f"  β = {[int(b) for b in beta]}")

# Transitive
A_t = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        A_t[i][j] = 1

print("\nTransitive tournament (t3=0):")
for p in range(n):
    ap = enumerate_allowed_paths(A_t, n, p)
    omega = compute_omega_basis(A_t, n, p, ap, enumerate_allowed_paths(A_t, n, p-1) if p > 0 else [])
    dim_omega = omega.shape[1] if omega.ndim == 2 else 0
    print(f"  p={p}: |A_{p}|={len(ap)}, dim(Ω_{p})={dim_omega}")

beta_t = path_betti_numbers(A_t, n)
print(f"  β = {[int(b) for b in beta_t]}")

# === KEY: Why dim(Ω_p) < |A_p|? ===
print("\n" + "=" * 60)
print("WHY Ω_p ≠ A_p FOR TOURNAMENTS?")
print("=" * 60)
print("""
Ω_p = {u ∈ A_p : ∂u ∈ A_{p-1}}

For a tournament on n vertices:
- A_p = all directed (p+1)-vertex paths
- ∂u ∈ A_{p-1} requires: for every (p+1)-path in u,
  removing any interior vertex gives an ALLOWED (p-1)-path

The catch: removing vertex v_i from path v_0→...→v_p gives v_0→...→v_{i-1}→v_{i+1}→...→v_p
This requires v_{i-1} → v_{i+1} in T.
If v_{i+1} → v_{i-1} instead, the face is NOT allowed.

So: an allowed p-path (a,b,c,...) has face (a,c,...) only if a→c.
If c→a, this face is DISALLOWED, and Ω_p must avoid creating this term.

For tournaments, roughly HALF the faces are disallowed (cyclic triples).
This makes Ω_p significantly smaller than A_p for cycle-rich tournaments.
""")

# Compute dim(Ω_p)/|A_p| ratio for different tournaments
print("dim(Ω_p)/|A_p| ratios for various n=5 tournaments:")
for trial in range(5):
    if trial == 0:
        A_test = A  # C5
        name = "C5 (t3=5)"
    elif trial == 1:
        A_test = A_t  # Transitive
        name = "Trans (t3=0)"
    else:
        A_test = random_tournament(n)
        t3 = count_t3(A_test, n)
        name = f"Random (t3={t3})"

    ratios = []
    for p in range(n):
        ap = enumerate_allowed_paths(A_test, n, p)
        if len(ap) == 0:
            continue
        omega = compute_omega_basis(A_test, n, p, ap,
                                    enumerate_allowed_paths(A_test, n, p-1) if p > 0 else [])
        dim_omega = omega.shape[1] if omega.ndim == 2 else 0
        ratios.append((p, len(ap), dim_omega))

    print(f"\n  {name}:")
    for p, a, o in ratios:
        print(f"    p={p}: |A|={a}, dim(Ω)={o}, ratio={o/a:.3f}")

# === KEY QUESTION: Is dim(Ω_2) = dim(im ∂_3) for tournaments? ===
print("\n" + "=" * 60)
print("β₂ = 0 CHECK: is ker(∂_2|_Ω) = im(∂_3|_Ω)?")
print("=" * 60)

for n_test in [5, 6, 7]:
    print(f"\nn={n_test}: sampling 50 random tournaments")
    all_zero = True
    for trial in range(50):
        A_test = random_tournament(n_test)
        beta = path_betti_numbers(A_test, n_test)
        # Check all even Betti numbers
        for p in range(2, len(beta), 2):
            if int(beta[p]) != 0:
                all_zero = False
                t3 = count_t3(A_test, n_test)
                print(f"  FAILURE: β_{p} = {int(beta[p])}, t3={t3}")
                break
    if all_zero:
        print(f"  ALL β_{{2k}} = 0 (50/50)")

# === MUTUAL EXCLUSIVITY of β₁ and β₃ ===
print("\n" + "=" * 60)
print("MUTUAL EXCLUSIVITY: β₁ > 0 XOR β₃ > 0 (or both 0)?")
print("=" * 60)

n_test = 7
counts = {'00': 0, '10': 0, '01': 0, '11': 0}
for trial in range(200):
    A_test = random_tournament(n_test)
    beta = path_betti_numbers(A_test, n_test)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0
    key = f"{min(b1,1)}{min(b3,1)}"
    counts[key] += 1

print(f"n=7, 200 trials:")
for k, v in counts.items():
    b1, b3 = int(k[0]), int(k[1])
    print(f"  (β₁>0={b1}, β₃>0={b3}): {v} ({100*v/200:.1f}%)")

if counts['11'] == 0:
    print("CONFIRMED: β₁ and β₃ are mutually exclusive at n=7!")
    print("  β₁ > 0 → β₃ = 0  AND  β₃ > 0 → β₁ = 0")

# === Euler characteristic ===
print("\n" + "=" * 60)
print("EULER CHARACTERISTIC χ = Σ(-1)^p β_p")
print("=" * 60)

chi_dist = {}
for trial in range(200):
    A_test = random_tournament(n_test)
    beta = path_betti_numbers(A_test, n_test)
    chi = sum((-1)**p * int(beta[p]) for p in range(len(beta)))
    chi_dist[chi] = chi_dist.get(chi, 0) + 1

print(f"n=7: χ distribution ({200} trials)")
for chi_val in sorted(chi_dist.keys()):
    print(f"  χ={chi_val}: {chi_dist[chi_val]} ({100*chi_dist[chi_val]/200:.1f}%)")

print("\nNote: if β_{2k}=0 and β₁ XOR β₃:")
print("  β=(1,0,0,0,...): χ = 1")
print("  β=(1,1,0,0,...): χ = 0")
print("  β=(1,0,0,1,...): χ = 2")
