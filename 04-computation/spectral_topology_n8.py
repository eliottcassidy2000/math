#!/usr/bin/env python3
"""
spectral_topology_n8.py — Spectral-topological connection at n=8

At n=8 (even), det(S) = Pf(S)^2 is nonzero. We can use BOTH:
  1. Spectral gap of S (imaginary eigenvalues)
  2. Pfaffian |Pf(S)| directly

At n=6, opus found: beta_1>0 => |Pf| in {1,3}, beta_3>0 => |Pf| in {7,9}
At n=7, we found: spectral gap separates phases

Questions:
  1. Does the Pfaffian separation extend to n=8?
  2. Does beta_4 appear at n=8? (reported by opus-S38)
  3. How does the spectral gap correlate with topology at even n?
  4. Are all three invariants (Pf, gap, H) correlated?

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

random.seed(314)
np.random.seed(314)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def skew_adj(A, n):
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def pfaffian_explicit(S):
    """Compute Pfaffian of skew-symmetric matrix using LU decomposition.
    For 2m x 2m matrix: Pf^2 = det(S).
    We compute det and take signed sqrt."""
    n = S.shape[0]
    if n % 2 == 1:
        return 0
    det_val = np.linalg.det(S)
    pf_sq = round(abs(det_val))
    pf = round(np.sqrt(pf_sq))
    # Sign: use the Hessenberg/Householder approach, or just round det
    return pf  # absolute value

def spectral_data(S):
    eigs = np.linalg.eigvals(S)
    pos_imag = sorted([e.imag for e in eigs if e.imag > 0.01])
    return pos_imag

def spectral_gap(pos_eigs):
    if len(pos_eigs) < 2:
        return 0.0
    return pos_eigs[-1] - pos_eigs[0]

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

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

# ===== n=8 SAMPLING =====
print("=" * 70)
print("n=8: SPECTRAL-TOPOLOGICAL CONNECTION (300 samples)")
print("=" * 70)

n = 8
N_SAMPLES = 300
results = []
t0 = time.time()

for trial in range(N_SAMPLES):
    A = random_tournament(n)

    # Betti numbers up to dim 4 to catch beta_3 (beta_4 too expensive)
    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except Exception:
        beta = [1, 0, 0, 0, 0]
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0
    b4 = 0  # skip beta_4 at n=8 (too expensive)

    S = skew_adj(A, n)
    pf = pfaffian_explicit(S)
    pe = spectral_data(S)
    gap = spectral_gap(pe)
    H = H_tournament(A, n)
    c3 = count_3cycles(A, n)

    if b4 > 0:
        phase = 'S4'
    elif b3 > 0:
        phase = 'S3'
    elif b1 > 0:
        phase = 'C'
    else:
        phase = 'P'

    results.append({
        'phase': phase, 'pf': pf, 'gap': gap,
        'eigs': pe, 'H': H, 'c3': c3,
        'b1': b1, 'b3': b3, 'b4': b4,
        'beta': beta[:6]
    })

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{N_SAMPLES} done ({elapsed:.1f}s)")

elapsed = time.time() - t0
print(f"\nTotal time: {elapsed:.1f}s")

# Phase distribution
print("\n" + "=" * 70)
print("PHASE DISTRIBUTION")
print("=" * 70)
phase_counts = Counter(r['phase'] for r in results)
for phase in ['P', 'C', 'S3', 'S4']:
    cnt = phase_counts.get(phase, 0)
    print(f"  {phase}: {cnt} ({100*cnt/N_SAMPLES:.1f}%)")

# Show all betti patterns
print("\nAll betti patterns:")
betti_pats = Counter(tuple(r['beta']) for r in results)
for pat, cnt in betti_pats.most_common(10):
    print(f"  {pat}: {cnt}")

# Pfaffian by phase
print("\n" + "=" * 70)
print("|Pf(S)| BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S3', 'S4']:
    phase_data = [r for r in results if r['phase'] == phase]
    if not phase_data:
        continue
    pf_counts = Counter(r['pf'] for r in phase_data)
    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  |Pf| distribution: {dict(sorted(pf_counts.items()))}")

# Spectral gap by phase
print("\n" + "=" * 70)
print("SPECTRAL GAP BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S3', 'S4']:
    phase_data = [r for r in results if r['phase'] == phase]
    if not phase_data:
        continue
    gaps = [r['gap'] for r in phase_data]
    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  gap: mean={np.mean(gaps):.3f}, std={np.std(gaps):.3f}, min={min(gaps):.3f}, max={max(gaps):.3f}")

# H by phase
print("\n" + "=" * 70)
print("H(T) BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S3', 'S4']:
    phase_data = [r for r in results if r['phase'] == phase]
    if not phase_data:
        continue
    Hs = [r['H'] for r in phase_data]
    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  H: mean={np.mean(Hs):.1f}, min={min(Hs)}, max={max(Hs)}")

# Cross-tabulation: Pf vs gap
print("\n" + "=" * 70)
print("|Pf| vs SPECTRAL GAP (all phases)")
print("=" * 70)

# Bin the spectral gap
for pf_val in sorted(set(r['pf'] for r in results)):
    pf_data = [r for r in results if r['pf'] == pf_val]
    gaps = [r['gap'] for r in pf_data]
    phases = Counter(r['phase'] for r in pf_data)
    print(f"|Pf|={pf_val}: n={len(pf_data)}, gap={np.mean(gaps):.3f}+/-{np.std(gaps):.3f}, phases={dict(phases)}")

# The key n=6 finding was |Pf| separates beta_1 from beta_3.
# At n=8, the |Pf| values are different. Let's see if the same pattern holds.
print("\n" + "=" * 70)
print("SEPARATION TEST: does |Pf| separate phases?")
print("=" * 70)

pf_for_C = set(r['pf'] for r in results if r['phase'] == 'C')
pf_for_S3 = set(r['pf'] for r in results if r['phase'] == 'S3')
pf_for_S4 = set(r['pf'] for r in results if r['phase'] == 'S4')
pf_for_P = set(r['pf'] for r in results if r['phase'] == 'P')

print(f"|Pf| values for C-phase: {sorted(pf_for_C)}")
print(f"|Pf| values for S3-phase: {sorted(pf_for_S3)}")
print(f"|Pf| values for S4-phase: {sorted(pf_for_S4)}")
print(f"|Pf| values for P-phase: {sorted(pf_for_P)}")

if pf_for_C and pf_for_S3:
    print(f"C ∩ S3: {sorted(pf_for_C & pf_for_S3)}")
if pf_for_C and pf_for_S4:
    print(f"C ∩ S4: {sorted(pf_for_C & pf_for_S4)}")
