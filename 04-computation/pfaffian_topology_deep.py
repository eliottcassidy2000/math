#!/usr/bin/env python3
"""
pfaffian_topology_deep.py — Deep investigation of Pf-topology connection

KEY FINDINGS SO FAR:
n=6: |Pf| PERFECTLY separates beta_1 from beta_3
  beta_1>0 => |Pf| in {1,3}
  beta_3>0 => |Pf| in {7,9}
n=7 (odd): det(S)=0, but spectral gap separates phases
n=8: |Pf| does NOT perfectly separate, but C-phase has |Pf| in {1,3,5} (small)

HYPOTHESIS: At EVEN n, C-phase (beta_1>0) requires |Pf(S)| <= n-1

This script tests:
1. Does beta_4>0 require specific |Pf| values?
2. Is there a THRESHOLD on |Pf| separating C from S3?
3. What invariant replaces |Pf| at odd n?

Also: investigate the specific beta_4>0 tournaments at n=8.

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

random.seed(999)
np.random.seed(999)

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

def pfaffian_abs(S):
    """Absolute value of Pfaffian (from det)."""
    n = S.shape[0]
    if n % 2 == 1:
        return 0
    det_val = np.linalg.det(S)
    return round(np.sqrt(abs(round(det_val))))

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

# ======= n=6 EXHAUSTIVE (verify opus findings) =======
print("=" * 70)
print("n=6: EXHAUSTIVE Pf-Betti analysis (32768 tournaments)")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
t0 = time.time()

pf_by_phase = {'P': Counter(), 'C': Counter(), 'S': Counter()}
phase_counts = Counter()
gap_by_phase = {'P': [], 'C': [], 'S': []}

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except Exception:
        continue

    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    S = skew_adj(A, n)
    pf = pfaffian_abs(S)
    pe = spectral_data(S)
    gap = spectral_gap(pe)

    phase = 'S' if b3 > 0 else ('C' if b1 > 0 else 'P')
    pf_by_phase[phase][pf] += 1
    phase_counts[phase] += 1
    gap_by_phase[phase].append(gap)

    if (bits + 1) % 10000 == 0:
        print(f"  {bits+1}/{total} ({time.time()-t0:.1f}s)")

elapsed = time.time() - t0
print(f"\nTotal: {elapsed:.1f}s")

print(f"\nPhase counts: P={phase_counts['P']}, C={phase_counts['C']}, S={phase_counts['S']}")

for phase in ['P', 'C', 'S']:
    print(f"\n{phase} ({phase_counts[phase]}):")
    print(f"  |Pf| distribution: {dict(sorted(pf_by_phase[phase].items()))}")
    gaps = gap_by_phase[phase]
    if gaps:
        print(f"  spectral gap: mean={np.mean(gaps):.3f}, min={min(gaps):.3f}, max={max(gaps):.3f}")

# Verify opus finding
if pf_by_phase['C']:
    c_pf = set(pf_by_phase['C'].keys())
    print(f"\nC-phase |Pf| values: {sorted(c_pf)}")
if pf_by_phase['S']:
    s_pf = set(pf_by_phase['S'].keys())
    print(f"S-phase |Pf| values: {sorted(s_pf)}")
if pf_by_phase['C'] and pf_by_phase['S']:
    print(f"C intersect S: {sorted(c_pf & s_pf)}")
    if not (c_pf & s_pf):
        print("CONFIRMED: |Pf| PERFECTLY SEPARATES C from S at n=6!")

# ======= Spectral gap at n=6 =======
print("\n" + "=" * 70)
print("n=6: SPECTRAL GAP BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    gaps = gap_by_phase[phase]
    if not gaps:
        continue
    print(f"\n{phase} ({len(gaps)}):")
    print(f"  gap: mean={np.mean(gaps):.3f}, std={np.std(gaps):.3f}, min={min(gaps):.3f}, max={max(gaps):.3f}")

# Check: is spectral gap also a perfect separator at n=6?
if gap_by_phase['C'] and gap_by_phase['S']:
    c_gap_range = (min(gap_by_phase['C']), max(gap_by_phase['C']))
    s_gap_range = (min(gap_by_phase['S']), max(gap_by_phase['S']))
    print(f"\nC gap range: {c_gap_range}")
    print(f"S gap range: {s_gap_range}")
    if c_gap_range[0] > s_gap_range[1]:
        print("SPECTRAL GAP also PERFECTLY SEPARATES C from S!")
    elif c_gap_range[1] < s_gap_range[0]:
        print("SPECTRAL GAP also PERFECTLY SEPARATES C from S! (reversed)")
    else:
        print("Spectral gap ranges OVERLAP — not a perfect separator")
