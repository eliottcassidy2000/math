#!/usr/bin/env python3
"""
interval_betti_comparison.py — Compare Betti numbers of Interval vs Paley

Paley T_7: β = [1, 0, 0, 0, 6, 0, 0], χ = 7
What is β for Interval T_7?

This connects the spectral/Walsh framework to topology:
- Paley has FLAT spectrum → β_4 = 6 (high-dim homology)
- Interval has PEAKED spectrum → β_? = ?
- Does the spectral structure predict Betti numbers?

Author: opus-2026-03-12-S67
"""

import sys
sys.path.insert(0, '04-computation')

from path_homology_v2 import path_betti_numbers

def make_circulant(n, S):
    """Build adjacency matrix for circulant tournament on Z_n with connection set S."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in S:
            adj[i][(i+d)%n] = 1
    return adj

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

print("=" * 70)
print("BETTI NUMBER COMPARISON: INTERVAL vs PALEY vs OTHER CIRCULANTS")
print("=" * 70)

# p = 7
p = 7
m = (p-1)//2

# All circulant orientations
orientations = {}
for bits in range(1 << m):
    sigma = tuple(1 if (bits >> j) & 1 else -1 for j in range(m))
    S = set()
    for k in range(1, m+1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    orientations[sigma] = sorted(S)

# Known H values
H_values = {
    (1,1,1): 175,    # Interval
    (1,1,-1): 189,   # Paley
    (-1,-1,1): 189,  # anti-Paley
    (-1,-1,-1): 175, # anti-Interval
    (1,-1,1): 175,
    (-1,1,-1): 175,
    (1,-1,-1): 175,
    (-1,1,1): 175,
}

print(f"\np={p}: Computing Betti numbers for all {len(orientations)} circulant tournaments")
print()

# Identify special orientations
sigma_int = (1,1,1)
sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))

for sigma in sorted(orientations, key=lambda s: -H_values.get(s, 0)):
    S = orientations[sigma]
    adj = make_circulant(p, S)

    tag = ""
    if sigma == sigma_int: tag = " [INTERVAL]"
    elif sigma == sigma_pal: tag = " [PALEY]"
    elif sigma == tuple(-x for x in sigma_pal): tag = " [anti-PALEY]"
    elif sigma == tuple(-x for x in sigma_int): tag = " [anti-INTERVAL]"

    betti = path_betti_numbers(adj, p, max_dim=p-1)

    chi = sum((-1)**i * b for i, b in enumerate(betti))
    H = H_values.get(sigma, '?')

    print(f"  σ={sigma}, S={S}, H={H}")
    print(f"    β = {betti}, χ = {chi}{tag}")
    print()

print("=" * 70)
print("KEY QUESTION: Does spectral structure predict Betti numbers?")
print("=" * 70)
print("""
If Paley (flat spectrum) has β_4=6 and Interval (peaked spectrum) has
a DIFFERENT Betti profile, this would show that spectral concentration
directly influences topological structure.

The connection would be:
  Eigenvalue distribution → cycle structure → chain complex → Betti numbers

Each step in this chain is known:
  1. Eigenvalues determine cycle counts (trace formula)
  2. Cycle structure determines Ω(T) (odd-cycle conflict graph)
  3. Ω(T) structure determines the chain complex via GLMY
  4. Chain complex determines Betti numbers

So spectral structure MUST predict Betti numbers — the question is
whether the prediction is sharp (exact Betti) or approximate (bounds).
""")

# Also compare at p=11 if feasible
print("=" * 70)
print("p=11: Interval vs Paley Betti (if feasible)")
print("=" * 70)
p = 11
m = 5
sigma_int = tuple([1]*m)
sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))

for name, sigma in [("Interval", sigma_int), ("Paley", sigma_pal)]:
    S = set()
    for k in range(1, m+1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    adj = make_circulant(p, sorted(S))

    print(f"\n  {name}: S={sorted(S)}")
    print(f"  Computing Betti numbers (this may take a while for p=11)...")
    try:
        betti = path_betti_numbers(adj, p, max_dim=min(p-1, 6))
        chi_partial = sum((-1)**i * b for i, b in enumerate(betti))
        print(f"    β (up to degree 6) = {betti}")
        print(f"    Partial χ = {chi_partial}")
    except Exception as e:
        print(f"    Error: {e}")

print("\nDONE.")
