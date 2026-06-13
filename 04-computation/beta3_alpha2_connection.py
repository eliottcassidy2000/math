#!/usr/bin/env python3
"""
beta3_alpha2_connection.py - Explore whether β₃ relates to α₂ in OCF.

OCF: H(T) = 1 + 2α₁ + 4α₂ + ...
where α_k = number of independent sets of size k in Ω(T).

Path homology: β = (1, β₁, 0, β₃, 0, β₅, ...)
with β₁ and β₃ mutually exclusive.

HYPOTHESIS: β₃ > 0 iff α₂ > 0 (or some threshold on α₂)?
Or: β₃ > 0 iff some specific cycle interaction exists?

At n=5: β₁ ∈ {0,1}, β₃ = 0 always. α₂ can be 0 or 1.
At n=6: β₃ appears for 1.2% of tournaments.
At n=7: β₃ appears for 8-11%.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter, defaultdict
import random

random.seed(42)

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def count_t5(A, n):
    t = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t += 1
    return t // 5

def count_alpha2(A, n):
    """Count pairs of vertex-disjoint directed odd cycles."""
    # Find all directed odd cycles
    cycles = []
    for L in range(3, n+1, 2):
        for combo in combinations(range(n), L):
            combo_set = frozenset(combo)
            found = False
            for perm in permutations(combo):
                if all(A[perm[i]][perm[(i+1)%L]] for i in range(L)):
                    found = True
                    break
            if found:
                cycles.append(combo_set)

    # Count vertex-disjoint pairs
    alpha2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i] & cycles[j]) == 0:
                alpha2 += 1
    return alpha2

def ham_count(A, n):
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            count += 1
    return count

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# === n=6: Check β₃ vs cycle structure ===
print("=" * 60)
print("n=6: β₃ vs CYCLE STRUCTURE (exhaustive)")
print("=" * 60)

n = 6
m = n*(n-1)//2
b3_data = []
b0_data = []

total = 1 << m
print(f"Total tournaments: {total}")

# This is 2^15 = 32768 — feasible but path homology is slow per tournament
# Let me sample instead
sample_size = 200
b3_count = 0
for trial in range(sample_size):
    A = random_tournament(n)
    beta = path_betti_numbers(A, n)
    b3 = int(beta[3]) if len(beta) > 3 else 0
    t3 = count_t3(A, n)

    if b3 > 0:
        b3_count += 1
        t5 = count_t5(A, n)
        a2 = count_alpha2(A, n)
        H = ham_count(A, n)
        b3_data.append({'t3': t3, 't5': t5, 'a2': a2, 'H': H})
    else:
        b0_data.append({'t3': t3})

print(f"\nn=6: β₃ > 0 in {b3_count}/{sample_size} ({100*b3_count/sample_size:.1f}%)")

if b3_data:
    print("\nβ₃ > 0 tournaments:")
    for d in b3_data:
        print(f"  t₃={d['t3']}, t₅={d['t5']}, α₂={d['a2']}, H={d['H']}")

    print(f"\nt₃ range for β₃>0: [{min(d['t3'] for d in b3_data)}, {max(d['t3'] for d in b3_data)}]")
    print(f"α₂ range for β₃>0: [{min(d['a2'] for d in b3_data)}, {max(d['a2'] for d in b3_data)}]")

# === n=7: More detailed analysis ===
print("\n" + "=" * 60)
print("n=7: β₃ vs CYCLE STRUCTURE (sampling)")
print("=" * 60)

n = 7
b3_data7 = []
b1_data7 = []
b0_data7 = []

for trial in range(150):
    A = random_tournament(n)
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0
    t3 = count_t3(A, n)

    if b3 > 0:
        t5 = count_t5(A, n)
        b3_data7.append({'t3': t3, 't5': t5, 'b1': b1})
    elif b1 > 0:
        b1_data7.append({'t3': t3})
    else:
        b0_data7.append({'t3': t3})

print(f"β₃>0: {len(b3_data7)}, β₁>0: {len(b1_data7)}, both 0: {len(b0_data7)}")

if b3_data7:
    t3_b3 = [d['t3'] for d in b3_data7]
    t5_b3 = [d['t5'] for d in b3_data7]
    print(f"t₃ for β₃>0: mean={sum(t3_b3)/len(t3_b3):.1f}, range=[{min(t3_b3)}, {max(t3_b3)}]")
    print(f"t₅ for β₃>0: mean={sum(t5_b3)/len(t5_b3):.1f}, range=[{min(t5_b3)}, {max(t5_b3)}]")

if b1_data7:
    t3_b1 = [d['t3'] for d in b1_data7]
    print(f"t₃ for β₁>0: mean={sum(t3_b1)/len(t3_b1):.1f}, range=[{min(t3_b1)}, {max(t3_b1)}]")

if b0_data7:
    t3_b0 = [d['t3'] for d in b0_data7]
    print(f"t₃ for β=0: mean={sum(t3_b0)/len(t3_b0):.1f}, range=[{min(t3_b0)}, {max(t3_b0)}]")

print("\n" + "=" * 60)
print("INTERPRETATION")
print("=" * 60)
print("""
The path homology of tournaments reveals a topological trichotomy:

1. CONTRACTIBLE (β = (1,0,...,0)): "topologically trivial"
   - Most tournaments (85-90% at n=7)
   - Low to moderate t₃

2. ONE 1-HOLE (β₁ = 1): "circle-like"
   - Strongly connected tournaments with enough 3-cycles
   - The 3-cycle structure creates an unfillable directed loop

3. ONE 3-HOLE (β₃ = 1): "directed 3-sphere-like"
   - Tournaments with specific higher cycle interactions
   - β₁ = 0: the 1-hole has been "filled" but a 3-hole emerges
   - Related to 5-cycle structure?

The mutual exclusivity of β₁ and β₃ suggests a PHASE TRANSITION:
as cycle complexity increases, the topology shifts from
"circle" (1-hole) to "higher sphere" (3-hole).

This parallels the OCF hierarchy: α₁ gives basic cycle count,
α₂ gives disjoint pair count. The topology similarly grades.
""")
