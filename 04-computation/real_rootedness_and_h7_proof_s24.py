#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Three investigations:
1. Real-rootedness of I(Ω(T), x) — is Ω(T) always claw-free?
2. Structural proof that H=7 is forbidden
3. H=21 investigation at n=7

Key finding so far:
- ALL zeros of I(Ω(T),x) are real and negative at n=3..6 (exhaustive)
- At n=7, α₁=3 always has α₂=2 → H=15 not H=7
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import time

# ─── Utilities ───

def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def count_3cycles(T):
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def hamiltonian_path_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def find_all_directed_odd_cycles(T):
    n = len(T)
    cycles = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(length):
                    if T[path[i]][path[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    mi = path.index(min(path))
                    canon = tuple(path[(mi + j) % length] for j in range(length))
                    cycles.add(canon)
                    break
    return list(cycles)

def build_conflict_graph(cycles):
    k = len(cycles)
    adj = [[False]*k for _ in range(k)]
    vsets = [set(c) for c in cycles]
    for i in range(k):
        for j in range(i+1, k):
            if vsets[i] & vsets[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def independence_polynomial_coeffs(adj):
    k = len(adj)
    if k == 0: return [1]
    coeffs = [0] * (k + 1)
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[len(verts)] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def has_claw(adj):
    """Check if graph has an induced claw (K_{1,3})."""
    k = len(adj)
    for center in range(k):
        neighbors = [j for j in range(k) if adj[center][j]]
        if len(neighbors) < 3:
            continue
        for triple in combinations(neighbors, 3):
            a, b, c = triple
            # Check if a,b,c are pairwise non-adjacent (independent)
            if not adj[a][b] and not adj[a][c] and not adj[b][c]:
                return True, center, triple
    return False, None, None

# ═══════════════════════════════════════════
# PART 1: CLAW-FREE CHECK
# ═══════════════════════════════════════════
print("="*60)
print("PART 1: Is Ω(T) always claw-free?")
print("="*60)

# At n≤6, exhaustive
for n in [3, 4, 5, 6]:
    has_claw_count = 0
    total = 0
    for T in all_tournaments(n):
        cycles = find_all_directed_odd_cycles(T)
        if len(cycles) < 4:
            total += 1
            continue
        adj = build_conflict_graph(cycles)
        total += 1
        result = has_claw(adj)
        if result[0]:
            has_claw_count += 1

    print(f"n={n}: {has_claw_count}/{total} tournaments have claw in Ω(T)")

# At n=7, sampling
print("\nn=7: sampling 5000 tournaments...")
has_claw_7 = 0
checked_7 = 0
for trial in range(5000):
    T = random_tournament(7)
    cycles = find_all_directed_odd_cycles(T)
    if len(cycles) >= 4:
        adj = build_conflict_graph(cycles)
        result = has_claw(adj)
        if result[0]:
            has_claw_7 += 1
            if has_claw_7 <= 3:
                print(f"  Claw found! center={result[1]}, leaves={result[2]}, #cycles={len(cycles)}")
                # Analyze the claw structure
                cycle_at = {i: cycles[i] for i in [result[1]] + list(result[2])}
                for idx, cyc in cycle_at.items():
                    print(f"    cycle {idx}: {cyc} (len={len(cyc)})")
    checked_7 += 1

print(f"n=7: {has_claw_7}/{checked_7} have claw ({100*has_claw_7/checked_7:.1f}%)")

# ═══════════════════════════════════════════
# PART 2: REAL-ROOTEDNESS AT n=7
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 2: Real-rootedness of I(Ω(T), x) at n=7")
print("="*60)

n_real = 0
n_not_real = 0
n_neg_real = 0

for trial in range(3000):
    T = random_tournament(7)
    cycles = find_all_directed_odd_cycles(T)
    adj = build_conflict_graph(cycles)
    coeffs = independence_polynomial_coeffs(adj)

    if len(coeffs) > 1:
        np_coeffs = list(reversed(coeffs))
        zeros = np.roots(np_coeffs)
        all_r = all(abs(z.imag) < 1e-8 for z in zeros)
        all_nr = all(abs(z.imag) < 1e-8 and z.real < 1e-8 for z in zeros)
        if all_r:
            n_real += 1
        else:
            n_not_real += 1
            if n_not_real <= 3:
                print(f"  Non-real zeros! coeffs={coeffs}, zeros={zeros}")
        if all_nr:
            n_neg_real += 1
    else:
        n_real += 1
        n_neg_real += 1

total = n_real + n_not_real
print(f"n=7 ({total} samples): all-real={n_real} ({100*n_real/total:.1f}%), "
      f"all-neg-real={n_neg_real} ({100*n_neg_real/total:.1f}%)")

# ═══════════════════════════════════════════
# PART 3: H=7 PROOF ANALYSIS
# ═══���═══════════════════════════════════════
print(f"\n{'='*60}")
print("PART 3: Why H=7 is forbidden — structural analysis")
print("="*60)

# H=7 = 1 + 2·α₁ + 4·α₂ + ...
# Only solution with α_k ≥ 0 and α_k ≤ C(α₁, k):
# α₁=3, α₂=0 (and all higher = 0)
# So we need to prove: α₁=3 → α₂ ≥ 1 for all n.

# At n=7: α₁=3 always has α₂=2. Let's check the STRUCTURE.
print("\nStructure of α₁=3 tournaments at n=7:")
structures = Counter()
n_alpha1_3 = 0
target = 200

random.seed(42)
for trial in range(1000000):
    T = random_tournament(7)
    c3 = count_3cycles(T)
    if c3 != 3:
        continue
    cycles = find_all_directed_odd_cycles(T)
    if len(cycles) != 3:
        continue

    n_alpha1_3 += 1
    vsets = [set(c) for c in cycles]

    # Classify the intersection pattern
    intersections = []
    for i in range(3):
        for j in range(i+1, 3):
            inter = vsets[i] & vsets[j]
            intersections.append(len(inter))
    pattern = tuple(sorted(intersections))
    structures[pattern] += 1

    if n_alpha1_3 <= 5:
        print(f"  cycles: {cycles}")
        print(f"  vertex sets: {vsets}")
        print(f"  intersections: {intersections}")
        union = vsets[0] | vsets[1] | vsets[2]
        others = set(range(7)) - union
        print(f"  union={union}, others={others}")

    if n_alpha1_3 >= target:
        break

print(f"\nIntersection patterns (pair sizes) for α₁=3 at n=7 ({n_alpha1_3} found):")
for pattern, count in sorted(structures.items()):
    pct = 100 * count / n_alpha1_3
    print(f"  {pattern}: {count} ({pct:.1f}%)")

# Now analyze: which patterns give α₂=0 (needed for H=7)?
# α₂=0 means ALL pairs share a vertex.
# Pattern (0, a, b) with a,b ≥ 1: one pair disjoint → α₂ ≥ 1 → H ≥ 11
# Pattern (a, b, c) with all ≥ 1: no disjoint pair → α₂ = 0 → H = 7
print(f"\nPatterns with all intersections ≥ 1 (would give α₂=0, H=7):")
for pattern, count in sorted(structures.items()):
    if all(x >= 1 for x in pattern):
        print(f"  {pattern}: {count}")
if not any(all(x >= 1 for x in p) for p in structures):
    print(f"  NONE! α₂=0 with α���=3 is impossible at n=7!")

# ═══════════════��═══════════════════════════
# PART 4: H=7 AT n=8 SAMPLING
# ═════════════════════════���═════════════════
print(f"\n{'='*60}")
print("PART 4: H=7 at n=8 (sampling)")
print("="*60)

h7_n8 = 0
h21_n8 = 0
h_near_7 = Counter()
for trial in range(100000):
    T = random_tournament(8)
    H = hamiltonian_path_count(T)
    if H == 7:
        h7_n8 += 1
    if H == 21:
        h21_n8 += 1
    if H <= 30:
        h_near_7[H] += 1

print(f"100K samples at n=8: H=7 found {h7_n8}, H=21 found {h21_n8}")
print(f"H ≤ 30 distribution: {dict(sorted(h_near_7.items()))}")
missing = sorted(set(range(1, 50, 2)) - set(h_near_7.keys()) - {h for h in range(0, 50, 2)})
print(f"Missing odd H ≤ 49: {missing}")

# ═══════��═══════════════════════════════════
# PART 5: WHY H=21 IS FORBIDDEN
# ════���══════════════════════════════════════
print(f"\n{'='*60}")
print("PART 5: H=21 analysis")
print("="*60)

# H=21 = 1 + 2α₁ + 4α₂ + 8α₃ + ...
# Solutions: (α₁,α₂,α₃): (10,0,0), (8,1,0), (6,2,0), (4,3,0), (2,4,0), (0,5,0)
# Also: (6,0,1), (4,1,1), (2,2,1), (0,3,1), (2,0,2), (0,1,2)

# At n=6: forbidden. At n=7: also forbidden (from sampling).
# Check the possible α decompositions at n=7.

# Sample tournaments near H=21 at n=7
print("Sampling n=7 tournaments near H=21...")
h_near_21 = Counter()
for trial in range(50000):
    T = random_tournament(7)
    H = hamiltonian_path_count(T)
    if 19 <= H <= 23:
        h_near_21[H] += 1

print(f"H near 21: {dict(sorted(h_near_21.items()))}")
print(f"H=21 found: {h_near_21.get(21, 0)}")

# Check: what (α₁, α₂) give H near 21?
print("\nPolynomial structure for H=19, 23 at n=7 (first few found):")
poly_for_h = defaultdict(list)
for trial in range(200000):
    T = random_tournament(7)
    H = hamiltonian_path_count(T)
    if H in [19, 21, 23] and len(poly_for_h[H]) < 10:
        cycles = find_all_directed_odd_cycles(T)
        adj = build_conflict_graph(cycles)
        coeffs = independence_polynomial_coeffs(adj)
        poly_for_h[H].append(tuple(coeffs))

for H_val in [19, 21, 23]:
    polys = Counter(poly_for_h[H_val])
    print(f"  H={H_val}: {dict(polys)}")

# ════════════���══════════════════════════════
# PART 6: PROVE THE PATTERN
# ════════════════��══════════════════════════
print(f"\n{'='*60}")
print("PART 6: α₁ vs α₂ relationship — is there a formula?")
print("="*60)

# For each (α₁, α₂) pair at n=6, what is the achievable H?
a1_a2_6 = Counter()
for T in all_tournaments(6):
    cycles = find_all_directed_odd_cycles(T)
    adj = build_conflict_graph(cycles)
    coeffs = independence_polynomial_coeffs(adj)
    a1 = coeffs[1] if len(coeffs) > 1 else 0
    a2 = coeffs[2] if len(coeffs) > 2 else 0
    H = hamiltonian_path_count(T)
    a1_a2_6[(a1, a2, H)] += 1

print("(α₁, α₂) → H at n=6:")
by_a1 = defaultdict(list)
for (a1, a2, H), count in sorted(a1_a2_6.items()):
    by_a1[a1].append((a2, H, count))

for a1 in sorted(by_a1):
    entries = by_a1[a1]
    # Group by a2
    by_a2 = defaultdict(list)
    for a2, H, count in entries:
        by_a2[a2].append((H, count))
    for a2 in sorted(by_a2):
        h_vals = [h for h, c in by_a2[a2]]
        print(f"  α₁={a1:2d}, α₂={a2}: H = {sorted(set(h_vals))}")

# Key check: H = 1 + 2α₁ + 4α₂ always?
print("\nVerify H = 1 + 2α₁ + 4α₂ (i.e., α_k=0 for k≥3):")
for (a1, a2, H), count in sorted(a1_a2_6.items()):
    predicted = 1 + 2*a1 + 4*a2
    if predicted != H:
        print(f"  MISMATCH: α��={a1}, α₂={a2}, H={H}, predicted={predicted}")
print("  (no output means all match)")

print("\n\nDONE.")
