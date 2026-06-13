#!/usr/bin/env python3
"""
pfaffian_betti_n7.py — Verify Pfaffian-Betti connection at n=7

At n=6 (exhaustive, 32768 tournaments, opus-S46e):
  β₁>0 ⟹ |Pf(S)| ∈ {1,3}
  β₃>0 ⟹ |Pf(S)| ∈ {7,9}
  Pfaffian COMPLETELY SEPARATES β₁ from β₃

At n=7 (odd), det(S)=0 (always for odd-dimensional skew-symmetric).
So we need different spectral invariants. Options:
  (a) Pfaffian of the (n-1)×(n-1) skew-adj of vertex-deleted T-v
  (b) Principal Pfaffians (sub-Pfaffians)
  (c) Skew eigenvalues of S
  (d) tr(S^k) for various k
  (e) The invariant polynomial det(xI - S) = x * prod(x^2 + lambda_k^2)

Since n=7 is odd, we'll use:
  1. Sub-Pfaffians: for each v, compute Pf(S[V\v]) where S[V\v] is 6×6 sub-matrix
  2. Skew spectrum: eigenvalues of iS (real)
  3. (x-1) Taylor coefficients of F(T,x)

Also check: does the *product* of sub-Pfaffians relate to β?
For even n, we can also check n=8 structure.

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

import io
old_stdout = sys.stdout
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from importlib import import_module
# Suppress path_homology_v2 print output
_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
ph_mod = import_module('path_homology_v2')
sys.stdout = _saved
path_betti_numbers = ph_mod.path_betti_numbers

random.seed(42)
np.random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def skew_adj(A, n):
    """S[i][j] = A[i][j] - A[j][i] = ±1 for i≠j, 0 on diagonal"""
    S = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def pfaffian_6(S6):
    """Compute Pfaffian of a 6×6 skew-symmetric matrix.
    Pf = sum over perfect matchings of K_6, with sign and product."""
    # Perfect matchings of {0,1,2,3,4,5}:
    # Standard formula: fix 0, pair with j, then recurse on remaining 4
    pf = 0
    remaining = [1,2,3,4,5]
    for idx_j, j in enumerate(remaining):
        rem4 = remaining[:idx_j] + remaining[idx_j+1:]
        # Sign for pairing (0,j): (-1)^idx_j (bringing j to position 1)
        sign_oj = (-1)**idx_j
        # Now compute Pf of 4×4 sub-matrix on rem4
        a, b, c, d = rem4
        # Pf of 4×4: S[a,b]*S[c,d] - S[a,c]*S[b,d] + S[a,d]*S[b,c]
        pf4 = S6[a,b]*S6[c,d] - S6[a,c]*S6[b,d] + S6[a,d]*S6[b,c]
        pf += sign_oj * S6[0,j] * pf4
    return pf

def sub_pfaffians(S, n):
    """For each vertex v, compute Pf of S with row/col v deleted."""
    result = []
    for v in range(n):
        idx = [i for i in range(n) if i != v]
        sub = S[np.ix_(idx, idx)]
        pf = pfaffian_6(sub)
        result.append(pf)
    return result

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def H_tournament(A, n):
    """Count Hamiltonian paths using bitmask DP."""
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

# ===== n=7 SAMPLING =====
print("=" * 70)
print("n=7: Pfaffian-Betti Connection (Sampling)")
print("=" * 70)

n = 7
N_SAMPLES = 500

# Collect data
data = []
t0 = time.time()

for trial in range(N_SAMPLES):
    A = random_tournament(n)

    # Betti numbers (max_dim=4 to catch β₃)
    beta = path_betti_numbers(A, n, max_dim=4)
    b1 = beta[1] if len(beta) > 1 else 0
    b3 = beta[3] if len(beta) > 3 else 0

    # Skew adjacency
    S = skew_adj(A, n)

    # Sub-Pfaffians (Pf of each 6×6 minor)
    spf = sub_pfaffians(S, n)

    # Skew eigenvalues
    eigs = np.linalg.eigvals(S * 1.0)
    skew_eigs = sorted([abs(e.imag) for e in eigs if abs(e.imag) > 0.01])

    # Cycle counts
    c3 = count_3cycles(A, n)
    H = H_tournament(A, n)

    # Score sequence
    scores = sorted([sum(A[i]) for i in range(n)])

    phase = 'S' if b3 > 0 else ('C' if b1 > 0 else 'P')

    data.append({
        'phase': phase, 'b1': b1, 'b3': b3,
        'spf': spf, 'spf_abs': sorted([abs(p) for p in spf]),
        'skew_eigs': skew_eigs,
        'c3': c3, 'H': H, 'scores': tuple(scores)
    })

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{N_SAMPLES} done ({elapsed:.1f}s)")

elapsed = time.time() - t0
print(f"\nTotal time: {elapsed:.1f}s")

# ===== ANALYSIS =====
print("\n" + "=" * 70)
print("PHASE DISTRIBUTION")
print("=" * 70)
phase_counts = Counter(d['phase'] for d in data)
for phase in ['P', 'C', 'S']:
    print(f"  {phase}: {phase_counts[phase]} ({100*phase_counts[phase]/N_SAMPLES:.1f}%)")

# Sub-Pfaffian analysis
print("\n" + "=" * 70)
print("SUB-PFAFFIAN (|Pf(S\\v)|) BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    phase_data = [d for d in data if d['phase'] == phase]
    if not phase_data:
        print(f"\n{phase}: no data")
        continue

    # Collect all |Pf| values across all vertices
    all_pf_abs = Counter()
    sorted_pf_patterns = Counter()
    for d in phase_data:
        for p in d['spf']:
            all_pf_abs[abs(p)] += 1
        sorted_pf_patterns[tuple(d['spf_abs'])] += 1

    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  |Pf| value distribution: {dict(sorted(all_pf_abs.items()))}")
    print(f"  Unique sorted |Pf| patterns: {len(sorted_pf_patterns)}")

    # Show most common patterns
    for pat, cnt in sorted_pf_patterns.most_common(5):
        print(f"    {pat}: {cnt}")

# Sum and product of sub-Pfaffians
print("\n" + "=" * 70)
print("SUM / PRODUCT OF SUB-PFAFFIANS BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    phase_data = [d for d in data if d['phase'] == phase]
    if not phase_data:
        continue

    sums = Counter()
    abs_sums = Counter()
    products = Counter()
    for d in phase_data:
        sums[sum(d['spf'])] += 1
        abs_sums[sum(abs(p) for p in d['spf'])] += 1

    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  sum(Pf(S\\v)) values: {dict(sorted(sums.items()))}")
    print(f"  sum(|Pf(S\\v)|) values: {dict(sorted(abs_sums.items()))}")

# Skew eigenvalue analysis
print("\n" + "=" * 70)
print("SKEW EIGENVALUE SPECTRUM BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    phase_data = [d for d in data if d['phase'] == phase]
    if not phase_data:
        continue

    # Eigenvalue profile: sort the 3 positive imaginary eigenvalues
    eigen_profiles = []
    for d in phase_data:
        eig3 = d['skew_eigs'][-3:] if len(d['skew_eigs']) >= 3 else d['skew_eigs']
        eigen_profiles.append(tuple(round(e, 2) for e in eig3))

    # Statistics
    max_eigs = [d['skew_eigs'][-1] if d['skew_eigs'] else 0 for d in phase_data]
    min_eigs = [d['skew_eigs'][0] if d['skew_eigs'] else 0 for d in phase_data]

    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  max|skew_eig|: mean={np.mean(max_eigs):.3f}, std={np.std(max_eigs):.3f}")
    print(f"  min|skew_eig|: mean={np.mean(min_eigs):.3f}, std={np.std(min_eigs):.3f}")
    print(f"  spectral gap (max-min): mean={np.mean([mx-mn for mx,mn in zip(max_eigs,min_eigs)]):.3f}")

# H and c3 by phase
print("\n" + "=" * 70)
print("H(T) AND c3 BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    phase_data = [d for d in data if d['phase'] == phase]
    if not phase_data:
        continue

    Hs = [d['H'] for d in phase_data]
    c3s = [d['c3'] for d in phase_data]

    print(f"\n{phase} ({len(phase_data)} tournaments):")
    print(f"  H: mean={np.mean(Hs):.1f}, min={min(Hs)}, max={max(Hs)}")
    print(f"  c3: mean={np.mean(c3s):.1f}, min={min(c3s)}, max={max(c3s)}")

# Connection: does sum of sub-Pfaffians = ±H?
print("\n" + "=" * 70)
print("COFACTOR EXPANSION: sum Pf vs H")
print("=" * 70)

# At odd n, det(S)=0. But cofactor expansion:
# sum_v (-1)^v * S[0,v] * Pf(S_0v) relates to the expansion.
# Actually, for a row expansion of the Pfaffian (which doesn't apply at odd n).
# Let's check: is sum_v Pf(S\v) related to H?

match_count = 0
for d in data[:50]:
    spf_sum = sum(d['spf'])
    H = d['H']
    print(f"  H={H:5d}, sum(Pf)={spf_sum:6d}, ratio={spf_sum/H:.4f}" if H > 0 else f"  H=0")
    if spf_sum == H or spf_sum == -H:
        match_count += 1

print(f"\n  Exact match sum(Pf) = ±H: {match_count}/50")

# Check: signed cofactor expansion
print("\n  Signed cofactor: sum_v (-1)^v * Pf(S\\v)")
for d in data[:20]:
    signed_sum = sum((-1)**v * d['spf'][v] for v in range(n))
    print(f"  H={d['H']:5d}, signed_sum={signed_sum:6d}, phase={d['phase']}")
