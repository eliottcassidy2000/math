#!/usr/bin/env python3
"""
s_phase_maximizer_n7.py — Are H-maximizers always S-phase?

At n=6: ALL 240 H-maximizers (H=45, score (2,2,2,3,3,3)) have beta_3>0.
At n=7: Are the 240 H-maximizers (H=189, regular score (3,3,3,3,3,3,3)) S-phase?

This would connect H-maximization to path homology topology.
If H-max => S-phase, then:
  H-max => near-degenerate skew eigenvalues => small spectral gap
  H-max => large |Pf| (at even n)

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def paley_tournament(p):
    """Paley tournament on Z/pZ for p = 3 mod 4."""
    qr = set((a*a) % p for a in range(1, p))
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def skew_adj(A, n):
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def spectral_gap(S):
    eigs = np.linalg.eigvals(S)
    pos = sorted([abs(e.imag) for e in eigs if abs(e.imag) > 0.01])
    if len(pos) < 2:
        return 0.0
    return pos[-1] - pos[0]

# ===== Test 1: Paley T_7 =====
print("=" * 70)
print("PALEY T_7: path homology")
print("=" * 70)

A7 = paley_tournament(7)
H7 = H_tournament(A7, 7)
print(f"H(T_7) = {H7}")

t0 = time.time()
beta7 = path_betti_numbers(A7, 7, max_dim=5)
print(f"beta = {beta7} (computed in {time.time()-t0:.1f}s)")

S7 = skew_adj(A7, 7)
gap7 = spectral_gap(S7)
print(f"spectral gap = {gap7:.6f}")
print(f"Phase: {'S' if beta7[3] > 0 else ('C' if beta7[1] > 0 else 'P')}")

# ===== Test 2: Generate n=7 maximizers =====
print("\n" + "=" * 70)
print("n=7 H-MAXIMIZER SEARCH + BETTI NUMBERS")
print("=" * 70)

# Generate tournaments by trying random regular ones
# Regular n=7: score (3,3,3,3,3,3,3)
# All maximizers have H=189 (known from OEIS)

# Strategy: generate random n=7 tournaments, keep those with H >= 170 (near-max)
n = 7
high_H = []
t0 = time.time()

for trial in range(50000):
    A = random_tournament(n)
    H = H_tournament(A, n)
    if H >= 160:  # Near-maximizers
        S = skew_adj(A, n)
        gap = spectral_gap(S)
        try:
            beta = path_betti_numbers(A, n, max_dim=4)
        except:
            beta = [1, 0, 0, 0, 0]
        b1 = int(beta[1]) if len(beta) > 1 else 0
        b3 = int(beta[3]) if len(beta) > 3 else 0
        phase = 'S' if b3 > 0 else ('C' if b1 > 0 else 'P')
        high_H.append({'H': H, 'gap': gap, 'phase': phase, 'beta': beta[:5]})

elapsed = time.time() - t0
print(f"\nSearched 50000 random tournaments in {elapsed:.1f}s")
print(f"Found {len(high_H)} with H >= 160")

# Results
print("\nH >= 160 tournaments:")
for d in sorted(high_H, key=lambda x: -x['H'])[:30]:
    print(f"  H={d['H']}, gap={d['gap']:.4f}, phase={d['phase']}, beta={d['beta']}")

# Statistics
print("\n--- Phase distribution for H >= 160 ---")
phase_counts = Counter(d['phase'] for d in high_H)
print(f"Phases: {dict(phase_counts)}")

# H = 189 (maximizers)
max_H = [d for d in high_H if d['H'] == 189]
if max_H:
    print(f"\nH=189 (max): {len(max_H)} found")
    phase_189 = Counter(d['phase'] for d in max_H)
    print(f"  Phases: {dict(phase_189)}")
    gaps_189 = [d['gap'] for d in max_H]
    print(f"  Spectral gap: {sorted(set(round(g,4) for g in gaps_189))}")

# H = 175 (second class)
h175 = [d for d in high_H if d['H'] == 175]
if h175:
    print(f"\nH=175: {len(h175)} found")
    phase_175 = Counter(d['phase'] for d in h175)
    print(f"  Phases: {dict(phase_175)}")

# H = 171 (third class)
h171 = [d for d in high_H if d['H'] == 171]
if h171:
    print(f"\nH=171: {len(h171)} found")
    phase_171 = Counter(d['phase'] for d in h171)
    print(f"  Phases: {dict(phase_171)}")

# ===== Test 3: check near-min H tournaments =====
print("\n" + "=" * 70)
print("LOW H TOURNAMENTS: phase distribution")
print("=" * 70)

low_H = []
random.seed(1000)
for trial in range(20000):
    A = random_tournament(n)
    H = H_tournament(A, n)
    if H <= 30:
        try:
            beta = path_betti_numbers(A, n, max_dim=4)
        except:
            beta = [1, 0, 0, 0, 0]
        b1 = int(beta[1]) if len(beta) > 1 else 0
        b3 = int(beta[3]) if len(beta) > 3 else 0
        phase = 'S' if b3 > 0 else ('C' if b1 > 0 else 'P')
        low_H.append({'H': H, 'phase': phase})

print(f"Found {len(low_H)} with H <= 30")
if low_H:
    phase_low = Counter(d['phase'] for d in low_H)
    print(f"Phases: {dict(phase_low)}")
    for H_val in sorted(set(d['H'] for d in low_H)):
        cnt = sum(1 for d in low_H if d['H'] == H_val)
        phases = Counter(d['phase'] for d in low_H if d['H'] == H_val)
        print(f"  H={H_val}: {cnt}, phases={dict(phases)}")
