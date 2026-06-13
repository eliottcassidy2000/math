#!/usr/bin/env python3
"""
beta4_classes_n7.py — Compare Betti for all 3 regular tournament classes at n=7

Regular n=7 has 3 rigid classes:
  H=189 (BIBD, 240 tours): beta_4 = 6
  H=175 (720 tours): beta_1 = 1 (C-phase)
  H=171 (1680 tours): contractible

Questions:
1. What is the path complex dimension for each class?
2. Why does BIBD get beta_4 = 6 but the others don't?
3. Is there a simple combinatorial reason?

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time, random
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers, enumerate_allowed_paths
sys.stdout = _saved

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

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

def paley_tournament(p):
    qr = set((a*a) % p for a in range(1, p))
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A

n = 7

# Find representatives of each regular class
classes = {}
for trial in range(200000):
    A = random_tournament(n)
    scores = [sum(A[i]) for i in range(n)]
    if all(s == 3 for s in scores):
        H = H_tournament(A, n)
        if H not in classes:
            classes[H] = [row[:] for row in A]
        if len(classes) >= 3:
            break

# Add Paley if not found
if 189 not in classes:
    A7 = paley_tournament(7)
    classes[189] = A7

print("=" * 70)
print("REGULAR n=7 TOURNAMENT CLASSES: Path complex comparison")
print("=" * 70)

for H_val in sorted(classes.keys(), reverse=True):
    A = classes[H_val]

    # Path complex dimensions
    omega_dims = []
    for dim in range(8):
        paths = enumerate_allowed_paths(A, n, dim)
        omega_dims.append(len(paths))

    # Betti numbers
    try:
        beta = path_betti_numbers(A, n, max_dim=6)
        beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(7)]
    except:
        beta_list = "FAILED"

    print(f"\nH={H_val}:")
    print(f"  Omega dims: {omega_dims}")
    print(f"  Betti: {beta_list}")

    # Per-subset Ham path counts for 5- and 6-vertex subsets
    from itertools import combinations
    from collections import Counter

    # 5-vertex subsets
    h5 = []
    for combo in combinations(range(n), 5):
        sub_A = [[0]*5 for _ in range(5)]
        for i in range(5):
            for j in range(5):
                sub_A[i][j] = A[combo[i]][combo[j]]
        h5.append(H_tournament(sub_A, 5))
    h5_dist = Counter(h5)
    print(f"  5-vertex sub H: {dict(sorted(h5_dist.items()))}")

    # 6-vertex subsets
    h6 = []
    for combo in combinations(range(n), 6):
        sub_A = [[0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                sub_A[i][j] = A[combo[i]][combo[j]]
        h6.append(H_tournament(sub_A, 6))
    h6_dist = Counter(h6)
    print(f"  6-vertex sub H: {dict(sorted(h6_dist.items()))}")

    # Betti of 6-vertex subtournaments
    print(f"  6-vertex sub Betti:")
    for combo in combinations(range(n), 6):
        sub_A = [[0]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                sub_A[i][j] = A[combo[i]][combo[j]]
        try:
            sub_beta = path_betti_numbers(sub_A, 6, max_dim=4)
            sub_beta_list = [int(sub_beta[k]) if k < len(sub_beta) else 0 for k in range(5)]
        except:
            sub_beta_list = "FAILED"
        H_sub = H_tournament(sub_A, 6)
        print(f"    {combo}: H={H_sub}, beta={sub_beta_list}")

print("\n" + "=" * 70)
print("KEY OBSERVATION")
print("=" * 70)
print("""
If H=189 (BIBD) has ALL 6-vertex sub-tournaments as H=45 maximizers,
and 5-vertex sub-tournaments with specific H values, this constrains
the path complex strongly.

beta_4 = 6 means 6 independent 4-cycles in the path complex.
The uniformity of the BIBD structure (vertex-transitive, all subtournaments
identical) creates these topological features.

The non-BIBD regular classes have heterogeneous subtournament structure,
which kills the high-dimensional homology.
""")
