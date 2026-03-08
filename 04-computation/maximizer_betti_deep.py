#!/usr/bin/env python3
"""
maximizer_betti_deep.py — Deep analysis of H-maximizer Betti split at n=6

At n=6: 480 H-maximizers (H=45), split:
  240 with beta = [1,1,0,0,0,0] (C-phase)
  240 with beta = [1,0,0,1,0,0] (S-phase)

Questions:
1. What score sequences does each group have?
2. What is the 3-cycle count c3 for each group?
3. Are these isomorphic to each other? (complement? reversal?)
4. What's the Pfaffian / spectral gap distribution?

At n=7: ALL 240 maximizers have beta_4 = 6.
Why does the split vanish at n=7?

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

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def skew_adj(A, n):
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def pfaffian_abs(S):
    n = S.shape[0]
    if n % 2 == 1: return 0
    det_val = np.linalg.det(S)
    return round(np.sqrt(abs(round(det_val))))

def spectral_gap(S):
    eigs = np.linalg.eigvals(S)
    pos = sorted([abs(e.imag) for e in eigs if abs(e.imag) > 0.01])
    if len(pos) < 2: return 0.0
    return pos[-1] - pos[0]

def complement(A, n):
    """Complement tournament: reverse all arcs."""
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                C[i][j] = 1 - A[i][j]
    return C

# ===== n=6 EXHAUSTIVE =====
print("=" * 70)
print("n=6 EXHAUSTIVE: H=45 maximizer analysis")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
t0 = time.time()

c_max = []  # C-phase maximizers
s_max = []  # S-phase maximizers

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

    H = H_tournament(A, n)
    if H != 45:
        continue

    try:
        beta = path_betti_numbers(A, n, max_dim=5)
    except:
        continue

    beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
    c3 = count_3cycles(A, n)
    sc = score_seq(A, n)
    S = skew_adj(A, n)
    pf = pfaffian_abs(S)
    gap = spectral_gap(S)

    entry = {
        'bits': bits, 'beta': beta_list, 'c3': c3,
        'score': sc, 'pf': pf, 'gap': gap, 'A': A
    }

    if beta_list[1] > 0:
        c_max.append(entry)
    elif beta_list[3] > 0:
        s_max.append(entry)
    else:
        print(f"  UNEXPECTED: H=45 with beta={beta_list}")

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")
print(f"C-phase (beta_1>0): {len(c_max)}")
print(f"S-phase (beta_3>0): {len(s_max)}")

# Score sequences
print("\n--- Score sequences ---")
c_scores = Counter(d['score'] for d in c_max)
s_scores = Counter(d['score'] for d in s_max)
print(f"C-phase scores: {dict(c_scores)}")
print(f"S-phase scores: {dict(s_scores)}")

# c3 distribution
print("\n--- 3-cycle counts ---")
c_c3 = Counter(d['c3'] for d in c_max)
s_c3 = Counter(d['c3'] for d in s_max)
print(f"C-phase c3: {dict(sorted(c_c3.items()))}")
print(f"S-phase c3: {dict(sorted(s_c3.items()))}")

# Pfaffian
print("\n--- |Pf| distribution ---")
c_pf = Counter(d['pf'] for d in c_max)
s_pf = Counter(d['pf'] for d in s_max)
print(f"C-phase |Pf|: {dict(sorted(c_pf.items()))}")
print(f"S-phase |Pf|: {dict(sorted(s_pf.items()))}")

# Spectral gap
print("\n--- Spectral gap ---")
c_gaps = [d['gap'] for d in c_max]
s_gaps = [d['gap'] for d in s_max]
if c_gaps:
    print(f"C-phase gap: mean={np.mean(c_gaps):.4f}, min={min(c_gaps):.4f}, max={max(c_gaps):.4f}")
if s_gaps:
    print(f"S-phase gap: mean={np.mean(s_gaps):.4f}, min={min(s_gaps):.4f}, max={max(s_gaps):.4f}")

# Check complement relationship
print("\n--- Complement analysis ---")
# For each C-phase maximizer, is its complement an S-phase maximizer?
c_bits = set(d['bits'] for d in c_max)
s_bits = set(d['bits'] for d in s_max)

comp_pairs = 0
for d in c_max:
    A = d['A']
    C = complement(A, n)
    # Encode complement
    c_bits_val = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if C[i][j]:
                c_bits_val |= (1 << idx)
            idx += 1
    if c_bits_val in s_bits:
        comp_pairs += 1
    elif c_bits_val in c_bits:
        pass  # complement is also C-phase

print(f"C-phase maximizers whose complement is S-phase: {comp_pairs}/{len(c_max)}")

# Check if complement of C-max is even in the maximizer set
comp_in_max = 0
comp_in_s = 0
comp_in_c = 0
comp_not_max = 0
for d in c_max:
    A = d['A']
    C = complement(A, n)
    HC = H_tournament(C, n)
    if HC == 45:
        comp_in_max += 1
        c_bits_val = 0
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if C[i][j]:
                    c_bits_val |= (1 << idx)
                idx += 1
        if c_bits_val in s_bits:
            comp_in_s += 1
        else:
            comp_in_c += 1
    else:
        comp_not_max += 1

print(f"Complement of C-max: max={comp_in_max} (S={comp_in_s}, C={comp_in_c}), not-max={comp_not_max}")

# Also for S-phase
comp_of_s_in_max = 0
comp_of_s_in_c = 0
for d in s_max:
    A = d['A']
    C = complement(A, n)
    HC = H_tournament(C, n)
    if HC == 45:
        comp_of_s_in_max += 1
    c_bits_val = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if C[i][j]:
                c_bits_val |= (1 << idx)
            idx += 1
    if c_bits_val in c_bits:
        comp_of_s_in_c += 1

print(f"Complement of S-max is C-max: {comp_of_s_in_c}/{len(s_max)}")
print(f"Complement of S-max has H=45: {comp_of_s_in_max}/{len(s_max)}")

# H values for complements
print("\n--- H of complements ---")
comp_H_of_c = []
for d in c_max[:50]:
    C = complement(d['A'], n)
    comp_H_of_c.append(H_tournament(C, n))
print(f"H of complements of C-max (first 50): {Counter(comp_H_of_c)}")

comp_H_of_s = []
for d in s_max[:50]:
    C = complement(d['A'], n)
    comp_H_of_s.append(H_tournament(C, n))
print(f"H of complements of S-max (first 50): {Counter(comp_H_of_s)}")

# ===== Check self-complementary maximizers =====
print("\n--- Self-complementary maximizers? ---")
sc_count = 0
for d in c_max + s_max:
    A = d['A']
    C = complement(A, n)
    c_bits_val = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if C[i][j]:
                c_bits_val |= (1 << idx)
            idx += 1
    if c_bits_val == d['bits']:
        sc_count += 1
        print(f"  Self-comp: bits={d['bits']}, phase={'C' if d['beta'][1]>0 else 'S'}, score={d['score']}")

print(f"Self-complementary maximizers: {sc_count}/{len(c_max)+len(s_max)}")

# ===== Quick look at which H values have which Betti =====
print("\n" + "=" * 70)
print("ALL H VALUES: Betti distribution at n=6")
print("=" * 70)

betti_by_H = {}
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

    H = H_tournament(A, n)
    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except:
        continue

    beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(5)]
    if H not in betti_by_H:
        betti_by_H[H] = []
    betti_by_H[H].append(tuple(beta_list))

print(f"\nH -> Betti distribution:")
for H_val in sorted(betti_by_H.keys(), reverse=True):
    entries = betti_by_H[H_val]
    betti_counts = Counter(entries)
    counts_str = ", ".join(f"{list(b)}:{cnt}" for b, cnt in betti_counts.most_common(3))
    print(f"  H={H_val:3d} ({len(entries):5d}): {counts_str}")
