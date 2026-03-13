#!/usr/bin/env python3
"""
overlap_weight_spectrum.py -- kind-pasteur-2026-03-13-S60

The overlap weight w(C) determines the "penalty" of including cycle C in
an independent set. Cycles with LOWER overlap weight conflict with fewer 
others, making them easier to include in independent sets.

KEY QUESTION: Does the Paley tournament maximize H because its overlap
weight spectrum is most FAVORABLE for independent set formation?

ANALYSIS:
1. Compare overlap weight distributions across orientations at p=7 and p=11
2. Check: does H correlate with overlap weight UNIFORMITY or SKEWNESS?
3. The Paley tournament has pair_qr = pair_nqr (edge-type symmetry).
   Does this symmetry translate to overlap weight uniformity?
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total

def compute_H_heldkarp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


# At p=7: compute overlap weight distribution for all 8 orientations
p = 7
m = (p - 1) // 2
N = 1 << m

print(f"{'='*70}")
print(f"OVERLAP WEIGHT SPECTRUM at p={p}")
print(f"{'='*70}")

results = []
for bits in range(N):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))
    
    A = build_adj(p, S)
    H = compute_H_heldkarp(A, p)
    
    # For each cycle, compute its overlap weight = N_total - 1 - C_odd(comp)
    cycles = []  # (vertex_set, length, directed_cycle_count)
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                cycles.append((frozenset(subset), k, nc))
    
    N_total = sum(nc for _, _, nc in cycles)
    c_k = defaultdict(int)
    for _, k, nc in cycles:
        c_k[k] += nc
    
    # Overlap weight for each active vertex set
    weights = []
    for fs, k, nc in cycles:
        comp = sorted(frozenset(range(p)) - fs)
        c_odd_comp = 0
        for kk in range(3, len(comp)+1, 2):
            for sub in combinations(comp, kk):
                c_odd_comp += count_ham_cycles(A, list(sub))
        w = N_total - 1 - c_odd_comp
        weights.append((w, k, nc))
    
    # Weight statistics
    all_w = [w for w, _, _ in weights]
    w_range = max(all_w) - min(all_w) if all_w else 0
    w_var = sum((w - sum(all_w)/len(all_w))**2 for w in all_w) / len(all_w) if all_w else 0
    
    # Unique weights
    unique_w = len(set(all_w))
    
    S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
    is_paley = sorted(S) == S_qr
    
    results.append({
        'bits': bits, 'S': S, 'H': H, 'N_total': N_total,
        'c_k': dict(c_k), 'n_cycles': len(cycles),
        'w_range': w_range, 'w_var': w_var, 'unique_w': unique_w,
        'weights': all_w, 'is_paley': is_paley,
    })
    
    paley_tag = " <-- PALEY" if is_paley else ""
    print(f"  bits={bits}: S={S}, H={H}, N={N_total}, c_k={dict(c_k)}, "
          f"w_range={w_range}, var={w_var:.1f}, unique_w={unique_w}{paley_tag}")

# Sort by H
results.sort(key=lambda r: -r['H'])
print(f"\nRanked by H:")
for r in results:
    paley_tag = " <-- PALEY" if r['is_paley'] else ""
    print(f"  H={r['H']}: N={r['N_total']}, w_range={r['w_range']}, "
          f"var={r['w_var']:.1f}, unique_w={r['unique_w']}{paley_tag}")

# KEY: Does H correlate with N_total? With w_var?
print(f"\nCorrelation analysis:")
Hs = [r['H'] for r in results]
Ns = [r['N_total'] for r in results]
import numpy as np
if len(set(Hs)) > 1 and len(set(Ns)) > 1:
    corr = np.corrcoef(Hs, Ns)[0,1]
    print(f"  corr(H, N_total) = {corr:.4f}")

wvars = [r['w_var'] for r in results]
if len(set(wvars)) > 1:
    corr2 = np.corrcoef(Hs, wvars)[0,1]
    print(f"  corr(H, w_var) = {corr2:.4f}")

# At p=11: sample a few orientations
print(f"\n\n{'='*70}")
print(f"OVERLAP WEIGHT SPECTRUM at p=11 (selected orientations)")
print(f"{'='*70}")

p = 11
m = (p - 1) // 2
S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)

# Test Paley and Interval, plus a few others
test_cases = [
    ("Paley", S_qr),
    ("Interval", list(range(1, m+1))),
    ("Alt1", [1, 10, 3, 8, 5]),  # alternating
    ("Alt2", [10, 2, 8, 4, 5]),  # another pattern
]

for name, S in test_cases:
    A = build_adj(p, S)
    H = compute_H_heldkarp(A, p)
    
    # Count all cycles with their overlap weights
    cycle_data = []
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                cycle_data.append((frozenset(subset), k, nc))
    
    N_total = sum(nc for _, _, nc in cycle_data)
    c_k = defaultdict(int)
    for _, k, nc in cycle_data:
        c_k[k] += nc
    
    # Overlap weights (only compute for k=3 and k=5 to save time)
    w_by_k = defaultdict(list)
    for fs, k, nc in cycle_data:
        if k > 5:
            continue  # skip large k (complement structure determined by smaller cycles anyway)
        comp = sorted(frozenset(range(p)) - fs)
        c_odd_comp = 0
        for kk in range(3, len(comp)+1, 2):
            for sub in combinations(comp, kk):
                c_odd_comp += count_ham_cycles(A, list(sub))
        w = N_total - 1 - c_odd_comp
        w_by_k[k].append(w)
    
    print(f"\n  {name}: S={S}, H={H}, N_total={N_total}")
    print(f"    c_k: {dict(sorted(c_k.items()))}")
    for k in sorted(w_by_k):
        ws = w_by_k[k]
        w_unique = sorted(set(ws))
        w_range = max(ws) - min(ws) if ws else 0
        print(f"    k={k}: {len(ws)} cycles, w in {w_unique[:5]}{'...' if len(w_unique) > 5 else ''}, "
              f"range={w_range}")

print("\nDONE.")
