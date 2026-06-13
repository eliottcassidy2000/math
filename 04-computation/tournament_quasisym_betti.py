#!/usr/bin/env python3
"""
tournament_quasisym_betti.py - Does β₃>0 have a quasisymmetric function signature?

Grinberg-Stanley: U_T ∈ N[p₁, 2p₃, 2p₅, ...] where p_k are power-sum QSFs.
The coefficient of p₁^a · (2p₃)^b · (2p₅)^c · ... encodes cycle decomposition.

For β₃>0 tournaments at n=6:
- All have t₃=8 (max), t₅=12 (max), α₂=1, H=45 (max)
- This means specific coefficients in U_T are nonzero

QUESTION: What's special about the quasisymmetric expansion for β₃>0?

Actually, let me think about this differently. The descent set composition of
a tournament T is encoded by F(T,x). But U_T encodes MORE: it's the
generating function over ALL permutations weighted by tournament arcs.

At n=6, let's examine the DESCENT SET distribution for β₃>0 vs β₁>0.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import permutations, combinations
from collections import Counter, defaultdict
import random

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def descent_set_distribution(A, n):
    """For each permutation, compute its descent set w.r.t. tournament T.
    Descent at position i if A[perm[i]][perm[i+1]] = 0 (i.e., perm[i+1]→perm[i])."""
    dist = Counter()
    for perm in permutations(range(n)):
        # Forward edges: positions where perm[i]→perm[i+1]
        fwd_set = frozenset(i for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd_set] += 1
    return dist

def descent_type_distribution(A, n):
    """Compute distribution of # forward edges (= F_k counts)."""
    counts = [0] * n
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        counts[fwd] += 1
    return counts

# === n=5: Descent pattern for β₁=0 vs β₁=1 ===
print("=" * 60)
print("n=5: F-polynomial by β₁ (exhaustive)")
print("=" * 60)

n = 5
m = n*(n-1)//2
fpoly_by_b1 = defaultdict(Counter)

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
    b1 = int(beta[1]) if len(beta) > 1 else 0
    F = tuple(descent_type_distribution(A, n))
    fpoly_by_b1[b1][F] += 1

print(f"\nβ₁=0: {sum(fpoly_by_b1[0].values())} tournaments, {len(fpoly_by_b1[0])} distinct F-polynomials")
print(f"β₁=1: {sum(fpoly_by_b1[1].values())} tournaments, {len(fpoly_by_b1[1])} distinct F-polynomials")

# Show the F-polynomials for β₁=1
print("\nF-polynomials for β₁=1:")
for F, count in sorted(fpoly_by_b1[1].items(), key=lambda x: -x[1]):
    print(f"  F = {list(F)}: {count} tournaments")

# Show the F-polynomials for β₁=0 (top 5)
print("\nTop F-polynomials for β₁=0:")
for F, count in sorted(fpoly_by_b1[0].items(), key=lambda x: -x[1])[:5]:
    print(f"  F = {list(F)}: {count} tournaments")

# === n=6: Sample β₃>0 tournaments and get F-poly ===
print("\n" + "=" * 60)
print("n=6: F-polynomial by β type (sampling)")
print("=" * 60)

n = 6
fpoly_by_type = {'b0': Counter(), 'b1': Counter(), 'b3': Counter()}
count_by_type = Counter()

for trial in range(3000):
    A = random_tournament(n)
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    key = 'b3' if b3 > 0 else ('b1' if b1 > 0 else 'b0')
    count_by_type[key] += 1

    F = tuple(descent_type_distribution(A, n))
    fpoly_by_type[key][F] += 1

    if (trial + 1) % 500 == 0:
        print(f"  ... {trial+1}/3000 (b3: {count_by_type['b3']})")

print(f"\nCounts: b0={count_by_type['b0']}, b1={count_by_type['b1']}, b3={count_by_type['b3']}")

if count_by_type['b3'] > 0:
    print(f"\nβ₃>0 F-polynomials ({len(fpoly_by_type['b3'])} distinct):")
    for F, count in sorted(fpoly_by_type['b3'].items(), key=lambda x: -x[1]):
        print(f"  F = {list(F)}: {count}")

    print(f"\nβ₁>0 F-polynomials (top 5 of {len(fpoly_by_type['b1'])}):")
    for F, count in sorted(fpoly_by_type['b1'].items(), key=lambda x: -x[1])[:5]:
        print(f"  F = {list(F)}: {count}")

# KEY: Do β₃>0 tournaments have a UNIQUE F-polynomial?
# If so, the path homology is determined by the descent statistics!
b3_F_set = set(fpoly_by_type['b3'].keys())
b1_F_set = set(fpoly_by_type['b1'].keys())
overlap = b3_F_set & b1_F_set
print(f"\nF-poly overlap between β₃>0 and β₁>0: {len(overlap)} polynomials")
if overlap:
    print(f"  Overlapping: {[list(F) for F in overlap]}")

b3_b0_overlap = b3_F_set & set(fpoly_by_type['b0'].keys())
print(f"F-poly overlap between β₃>0 and β=0: {len(b3_b0_overlap)}")

if not overlap and not b3_b0_overlap:
    print("\n*** β₃>0 has UNIQUE F-polynomial(s) — path homology detectable from F(T,x)! ***")
