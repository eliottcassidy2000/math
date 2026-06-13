#!/usr/bin/env python3
"""
H-CONNECTION: Path homology and Hamiltonian path count

Key observation: S-phase (β_3=1) tournaments have mean H=114.8
while P-phase has mean H=75.3 and C-phase has H=73.5.
"""
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import random
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    path_betti_numbers, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# ===== n=5: Exhaustive H vs beta =====
print("=" * 70)
print("n=5: H(T) vs TOPOLOGICAL TYPE (EXHAUSTIVE)")
print("=" * 70)

n = 5
H_by_phase = defaultdict(list)
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

for mask in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    betti = path_betti_numbers(A, n, max_dim=n-1)
    H = ham_path_count(A, n)
    phase = "P" if betti[1]==0 else "C"
    H_by_phase[phase].append(H)

print(f"\nPhase P: {len(H_by_phase['P'])} tournaments")
print(f"  H values: {sorted(Counter(H_by_phase['P']).items())}")
print(f"  Mean H: {np.mean(H_by_phase['P']):.2f}")

print(f"\nPhase C: {len(H_by_phase['C'])} tournaments")
print(f"  H values: {sorted(Counter(H_by_phase['C']).items())}")
print(f"  Mean H: {np.mean(H_by_phase['C']):.2f}")

# ===== n=7: H by phase =====
print("\n\n" + "=" * 70)
print("n=7: H(T) vs TOPOLOGICAL TYPE (500 samples)")
print("=" * 70)

n = 7
H_by_phase = defaultdict(list)
for trial in range(500):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    H = ham_path_count(A, n)
    phase = "P" if betti[1]==0 and betti[3]==0 else ("C" if betti[1]>0 else "S")
    H_by_phase[phase].append(H)
    if trial % 100 == 99:
        print(f"  ... {trial+1} done", flush=True)

for phase in ['P', 'C', 'S']:
    if H_by_phase[phase]:
        vals = H_by_phase[phase]
        print(f"\nPhase {phase}: {len(vals)} tournaments")
        print(f"  H range: [{min(vals)}, {max(vals)}]")
        print(f"  H mean: {np.mean(vals):.1f}, median: {np.median(vals):.0f}")

print(f"\nH thresholds:")
for H_thresh in [40, 80, 100, 120, 140]:
    counts = {p: sum(1 for h in H_by_phase[p] if h >= H_thresh) for p in ['P','C','S']}
    parts = ", ".join(f"{p}:{counts[p]}/{len(H_by_phase[p])}" for p in ['P','C','S'] if H_by_phase[p])
    print(f"  H >= {H_thresh}: {parts}")

# ===== OCF cycle counts vs phase =====
print("\n\n" + "=" * 70)
print("OCF CYCLE COUNTS vs PHASE (n=7)")
print("=" * 70)

def count_kcycles(A, n, k):
    from itertools import permutations
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            is_cycle = True
            for i in range(k):
                if A[perm[i]][perm[(i+1)%k]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count // k

n = 7
cycle_data = defaultdict(list)
for trial in range(200):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    phase = "P" if betti[1]==0 and betti[3]==0 else ("C" if betti[1]>0 else "S")

    t3 = count_3cycles(A, n)
    t5 = count_kcycles(A, n, 5)
    H = ham_path_count(A, n)
    cycle_data[phase].append((t3, t5, H))

    if trial % 50 == 49:
        print(f"  ... {trial+1} done", flush=True)

for phase in ['P', 'C', 'S']:
    if cycle_data[phase]:
        d = cycle_data[phase]
        print(f"\n  Phase {phase} ({len(d)} tournaments):")
        print(f"    t3: mean={np.mean([x[0] for x in d]):.1f}")
        print(f"    t5: mean={np.mean([x[1] for x in d]):.1f}")
        print(f"    H:  mean={np.mean([x[2] for x in d]):.1f}")

print("\nDone.")
