#!/usr/bin/env python3
"""
path_homology_ocf_connection.py - Explore OCF-path homology connection.

KEY FINDINGS from previous script:
1. β_{2k} = 0 for all tournaments (even Betti vanishing)
2. β₁ and β₃ are MUTUALLY EXCLUSIVE
3. χ(T) = 1 - β₁ + β₃ - β₅ + ... ∈ {0, 1} at n=5

HYPOTHESIS TO TEST:
- χ(T) = 1 for "acyclic-like" tournaments (few cycles)
- χ(T) = 0 for "cycle-rich" tournaments
- Is χ(T) = H(T) mod 2? (always 1 by Redei)
  No, since χ=0 exists. But maybe χ(T) = 1 - β₁ always?

DEEPER QUESTION:
OCF: H(T) = 1 + 2·α₁ + 4·α₂ + ...
χ(T) = 1 - β₁ + β₃ - β₅ + ...

If β₁ = indicator(some cycle condition) and β₃ = indicator(some higher condition):
  χ = 1    if no cycles "detected" by path homology
  χ = 0    if β₁ = 1 (one "directed hole" from 3-cycle structure)
  χ = 0    if β₃ = 1 (replaces the β₁ = 1 case?)

The mutual exclusivity suggests a PHASE TRANSITION in topology.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import factorial

def count_t3(A, n):
    return sum(1 for i, j, k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def ham_path_count(A, n):
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            count += 1
    return count

def count_alpha1(A, n):
    """Count directed odd cycles (alpha_1 in OCF)."""
    total = 0
    for L in range(3, n+1, 2):  # odd L
        for combo in combinations(range(n), L):
            for perm in permutations(combo):
                if all(A[perm[i]][perm[(i+1)%L]] for i in range(L)):
                    total += 1
                    break  # Count each vertex set once? No...
    # Actually alpha_1 counts independent sets of size 1 in Omega
    # = number of directed odd cycles
    # Let me count each directed odd cycle separately
    total = 0
    for L in range(3, n+1, 2):
        count = 0
        for combo in combinations(range(n), L):
            for perm in permutations(combo):
                if all(A[perm[i]][perm[(i+1)%L]] for i in range(L)):
                    count += 1
        total += count // L  # Each cycle counted L times
    return total

# === Exhaustive n=5 ===
print("=" * 70)
print("n=5: EXHAUSTIVE β, H, t₃, α₁ CORRELATION")
print("=" * 70)

n = 5
m = n*(n-1)//2
results = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, n)
    t3 = count_t3(A, n)
    H = ham_path_count(A, n)
    chi = sum((-1)**p * int(beta[p]) for p in range(len(beta)))
    results.append({
        'beta': tuple(int(b) for b in beta),
        't3': t3, 'H': H, 'chi': chi,
        'b1': int(beta[1]) if len(beta) > 1 else 0,
        'b3': int(beta[3]) if len(beta) > 3 else 0
    })

# Cross-tabulate (β₁, β₃) with H mod 4
print("\n(β₁, β₃) vs H mod 4:")
tab = defaultdict(Counter)
for r in results:
    tab[(r['b1'], r['b3'])][r['H'] % 4] += 1
for key in sorted(tab.keys()):
    print(f"  (β₁={key[0]}, β₃={key[1]}): H mod 4 = {dict(tab[key])}")

# Cross-tabulate β with t₃
print("\nβ₁ vs t₃:")
tab2 = defaultdict(Counter)
for r in results:
    tab2[r['t3']][r['b1']] += 1
for t3 in sorted(tab2.keys()):
    print(f"  t₃={t3}: β₁ distribution = {dict(tab2[t3])}")

# Check: is β₁ = 1 iff t₃ >= threshold?
print("\nIs β₁ a threshold function of t₃?")
for threshold in range(1, 6):
    correct = sum(1 for r in results if (r['b1'] == 1) == (r['t3'] >= threshold))
    print(f"  β₁ = 1 iff t₃ >= {threshold}: {correct}/{len(results)} ({100*correct/len(results):.1f}%)")

# Exhaustive n=4
print("\n" + "=" * 70)
print("n=4: EXHAUSTIVE")
print("=" * 70)

n = 4
m = n*(n-1)//2
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, n)
    t3 = count_t3(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0

# At n=4: β₁=1 iff t₃=2 (the only non-transitive 4-tournament has 2 directed 3-cycles??)
# Actually n=4 has 0,1,2,4 as possible t₃ values. Let me check.
results4 = []
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, n)
    t3 = count_t3(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    results4.append({'t3': t3, 'b1': b1})

tab4 = defaultdict(Counter)
for r in results4:
    tab4[r['t3']][r['b1']] += 1
print("n=4: β₁ vs t₃")
for t3 in sorted(tab4.keys()):
    print(f"  t₃={t3}: β₁ = {dict(tab4[t3])}")

# n=4: max t₃ = C(4,3) - 0 = 4 (all triples cyclic). Actually max is...
# t₃ = C(n,3) - sum C(s_i, 2). For score (1,1,2,2): t₃ = 4 - 2 = 2.
# For score (0,1,2,3): t₃ = 0. Max t₃ at n=4 is 2 (not 4).
# Wait n=4 has C(4,3)=4 triples. Max t₃ = 4 - min(sum C(s_i,2)).

# === n=6 exhaustive (64 tournaments? No, 2^15 = 32768) ===
# That's feasible
print("\n" + "=" * 70)
print("n=6: EXHAUSTIVE (32768 tournaments)")
print("=" * 70)

n = 6
m = n*(n-1)//2
betti_dist6 = Counter()
b_combos = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, n)
    beta_tup = tuple(int(b) for b in beta)
    betti_dist6[beta_tup] += 1
    t3 = count_t3(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0
    b_combos[(b1, b3, t3)] += 1

print(f"Total: {1 << m}")
print("Betti distribution:")
for k, v in sorted(betti_dist6.items(), key=lambda x: -x[1]):
    print(f"  β={list(k)}: {v} ({100*v/(1<<m):.2f}%)")

# β₁ vs β₃ mutual exclusivity
print("\n(β₁, β₃) at n=6:")
b1b3 = Counter()
for (b1, b3, t3), count in b_combos.items():
    b1b3[(b1, b3)] += count
for k, v in sorted(b1b3.items()):
    print(f"  (β₁={k[0]}, β₃={k[1]}): {v}")

# β₃ vs t₃
print("\nβ₃ > 0 by t₃:")
b3_by_t3 = defaultdict(lambda: [0, 0])
for (b1, b3, t3), count in b_combos.items():
    if b3 > 0:
        b3_by_t3[t3][1] += count
    else:
        b3_by_t3[t3][0] += count
for t3 in sorted(b3_by_t3.keys()):
    total = sum(b3_by_t3[t3])
    b3_yes = b3_by_t3[t3][1]
    if b3_yes > 0:
        print(f"  t₃={t3}: β₃>0 = {b3_yes}/{total} ({100*b3_yes/total:.1f}%)")
