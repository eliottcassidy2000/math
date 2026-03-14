#!/usr/bin/env python3
"""
determination_power.py — opus-2026-03-14-S71d

What determines what? Using CORRECT directed cycle counts.

Given H, what can we determine about (dc3, dc5, dc7)?
Given (H, α₁), what more do we learn?
Given (H, α₁, α₂), is everything determined?

The fundamental relations at n≤8 (α₃=0):
  H = 1 + 2α₁ + 4α₂
  α₁ = dc3 + dc5 + dc7
  α₂ = #{disjoint pairs of directed cycles}

Also: the det(I+xA) coefficients give DIFFERENT cycle information:
  c_3^det = #{3-cycle vertex-sets} = dc3
  c_5^det = #{5-cycle covers minus paired terms}

Can we use det(I+xA) + OCF together to extract everything?
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
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

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {}
    dp[(1 << 0, 0)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

def det_coefficients(A, n):
    """Compute coefficients of det(I + xA) via principal minors."""
    M = A.astype(float)
    coeffs = [1]  # c_0 = 1
    for k in range(1, n+1):
        total = 0
        for combo in combinations(range(n), k):
            sub = M[np.ix_(combo, combo)]
            total += np.linalg.det(sub)
        coeffs.append(int(round(total)))
    return coeffs

# ======================================================================
# n=5 exhaustive: what determines what?
# ======================================================================
print("=" * 70)
print("DETERMINATION HIERARCHY AT n=5 (exhaustive)")
print("=" * 70)

n = 5
tb = n*(n-1)//2

data = []
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_ham_cycles(A, n)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    coeffs = det_coefficients(A, n)
    data.append({
        'H': H, 'dc3': dc3, 'dc5': dc5, 'a1': a1, 'a2': a2,
        'c3': coeffs[3], 'c4': coeffs[4], 'c5': coeffs[5],
        'bits': bits
    })

# Group by H alone
groups_H = defaultdict(set)
for d in data:
    groups_H[d['H']].add((d['dc3'], d['dc5']))

print(f"\n  H alone → (dc3, dc5):")
for H in sorted(groups_H.keys()):
    vals = groups_H[H]
    amb = "AMBIGUOUS" if len(vals) > 1 else "unique"
    print(f"    H={H:3d}: {sorted(vals)} ({amb})")

# Group by (H, c3) — since c3 = dc3
groups_Hc3 = defaultdict(set)
for d in data:
    groups_Hc3[(d['H'], d['c3'])].add((d['dc3'], d['dc5']))

amb_Hc3 = sum(1 for g in groups_Hc3.values() if len(g) > 1)
print(f"\n  (H, c3) → (dc3, dc5): {amb_Hc3}/{len(groups_Hc3)} ambiguous")

# Since c3 = dc3, this is really (H, dc3) → dc5
groups_Hdc3 = defaultdict(set)
for d in data:
    groups_Hdc3[(d['H'], d['dc3'])].add(d['dc5'])

amb_Hdc3 = sum(1 for g in groups_Hdc3.values() if len(g) > 1)
print(f"  (H, dc3) → dc5: {amb_Hdc3}/{len(groups_Hdc3)} ambiguous")
for key, vals in sorted(groups_Hdc3.items()):
    if len(vals) > 1:
        print(f"    (H={key[0]}, dc3={key[1]}): dc5 ∈ {sorted(vals)}")

# Group by det(I+xA) coefficients
groups_det = defaultdict(set)
for d in data:
    groups_det[(d['c3'], d['c4'], d['c5'])].add((d['dc3'], d['dc5'], d['H']))

print(f"\n  det(I+xA) coefficients → (dc3, dc5, H):")
for key in sorted(groups_det.keys()):
    vals = groups_det[key]
    if len(vals) > 1:
        print(f"    c=(0,0,{key[0]},{key[1]},{key[2]}): {sorted(vals)} (AMBIGUOUS)")

# Does det(I+xA) determine H?
groups_det_H = defaultdict(set)
for d in data:
    groups_det_H[tuple(d[f'c{k}'] for k in [3,4,5])].add(d['H'])

amb_det_H = sum(1 for g in groups_det_H.values() if len(g) > 1)
print(f"\n  det(I+xA) → H: {amb_det_H}/{len(groups_det_H)} ambiguous")

# ======================================================================
# n=7 sample: what determines what?
# ======================================================================
print("\n" + "=" * 70)
print("DETERMINATION HIERARCHY AT n=7 (200 samples)")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

data7 = []
t0 = time.time()
for trial in range(200):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4

    data7.append({
        'H': H, 'dc3': dc3, 'dc5': dc5, 'dc7': dc7,
        'a1': a1, 'a2': a2
    })
    if trial % 50 == 0:
        print(f"  trial {trial}: {time.time()-t0:.1f}s")

# Group by H
groups_H7 = defaultdict(set)
for d in data7:
    groups_H7[d['H']].add((d['dc3'], d['dc5'], d['dc7']))

amb_H = sum(1 for g in groups_H7.values() if len(g) > 1)
print(f"\n  H → (dc3, dc5, dc7): {amb_H}/{len(groups_H7)} ambiguous")

# Group by (H, dc3)
groups_Hdc3_7 = defaultdict(set)
for d in data7:
    groups_Hdc3_7[(d['H'], d['dc3'])].add((d['dc5'], d['dc7']))

amb = sum(1 for g in groups_Hdc3_7.values() if len(g) > 1)
print(f"  (H, dc3) → (dc5, dc7): {amb}/{len(groups_Hdc3_7)} ambiguous")

# Group by (H, α₁)
groups_Ha1 = defaultdict(set)
for d in data7:
    groups_Ha1[(d['H'], d['a1'])].add((d['dc3'], d['dc5'], d['dc7']))

amb = sum(1 for g in groups_Ha1.values() if len(g) > 1)
print(f"  (H, α₁) → (dc3, dc5, dc7): {amb}/{len(groups_Ha1)} ambiguous")

# Group by (α₁, α₂)
groups_a12 = defaultdict(set)
for d in data7:
    groups_a12[(d['a1'], d['a2'])].add((d['dc3'], d['dc5'], d['dc7']))

amb = sum(1 for g in groups_a12.values() if len(g) > 1)
print(f"  (α₁, α₂) → (dc3, dc5, dc7): {amb}/{len(groups_a12)} ambiguous")

# Group by (dc3, dc5, dc7)
groups_full = defaultdict(set)
for d in data7:
    groups_full[(d['dc3'], d['dc5'], d['dc7'])].add(d['H'])

amb = sum(1 for g in groups_full.values() if len(g) > 1)
print(f"  (dc3, dc5, dc7) → H: {amb}/{len(groups_full)} ambiguous")

# ======================================================================
# KEY QUESTION: Does (dc3, dc5, dc7) determine H?
# ======================================================================
print("\n" + "=" * 70)
print("DOES (dc3, dc5, dc7) DETERMINE H?")
print("=" * 70)

# H = 1 + 2(dc3+dc5+dc7) + 4α₂
# So H is determined by (dc3, dc5, dc7) iff α₂ is determined by them
# α₂ = #{disjoint pairs of directed cycles}

# At n=7: disjoint pairs must use 6 vertices (3+3)
# So α₂ = #{pairs of vertex-disjoint directed 3-cycles}
# This depends on the STRUCTURE of 3-cycles, not just their count!

# Two tournaments with the same dc3 can have different α₂
# (e.g., 3 vertex-disjoint 3-cycles vs 3 overlapping 3-cycles)

print(f"\n  (dc3, dc5, dc7) → α₂ ambiguity shows α₂ is NOT cycle-count-determined:")
for key, vals in sorted(groups_full.items()):
    if len(vals) > 1:
        d_ex = [d for d in data7 if (d['dc3'], d['dc5'], d['dc7']) == key]
        a2_vals = sorted(set(d['a2'] for d in d_ex))
        H_vals = sorted(vals)
        print(f"    (dc3={key[0]}, dc5={key[1]}, dc7={key[2]}): α₂ ∈ {a2_vals}, H ∈ {H_vals}")

print("\nDone.")
