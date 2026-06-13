#!/usr/bin/env python3
"""
VERIFY: D_S mod 2^{n-1-k} is universal for ALL position sets S of size k,
not just non-adjacent ones.

If this holds, the universal congruence proof is COMPLETE and SIMPLE.
"""
from itertools import permutations, combinations
from math import factorial, comb
import random

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_DS(A, n, positions):
    count = 0
    for P in permutations(range(n)):
        ok = all(A[P[pos]][P[pos+1]] for pos in positions)
        if ok:
            count += 1
    return count

# Test at n=5 (fast, exhaustive for each S)
print("=== n=5: All D_S mod 2^{n-1-k} ===")
n = 5
edge_positions = list(range(n-1))

for k in range(1, n):
    mod = 2**(n-1-k)
    if mod == 0: mod = 1
    print(f"\n  k={k}, mod = 2^{n-1-k} = {mod}:")

    for S in combinations(edge_positions, k):
        DS_mod_vals = set()
        for trial in range(50):
            A = random_tournament(n)
            DS = compute_DS(A, n, S)
            DS_mod_vals.add(DS % mod)

        universal = len(DS_mod_vals) == 1
        val = DS_mod_vals.pop() if universal else None
        print(f"    S={S}: {'UNIVERSAL = ' + str(val) if universal else 'NOT UNIVERSAL: ' + str(DS_mod_vals)}")

# Test at n=7 for higher k
print("\n\n=== n=7: D_S mod 2^{n-1-k} for k=4,5 ===")
n = 7
edge_positions = list(range(n-1))

for k in [4, 5]:
    mod = 2**(n-1-k)
    if mod == 0: mod = 1
    print(f"\n  k={k}, mod = 2^{n-1-k} = {mod}:")

    all_universal = True
    for S in combinations(edge_positions, k):
        DS_mod_vals = set()
        for trial in range(20):
            A = random_tournament(n)
            DS = compute_DS(A, n, S)
            DS_mod_vals.add(DS % mod)

        universal = len(DS_mod_vals) == 1
        if not universal:
            all_universal = False
            print(f"    S={S}: NOT UNIVERSAL: {DS_mod_vals}")
        else:
            val = DS_mod_vals.pop()
            print(f"    S={S}: UNIVERSAL = {val}")

    if all_universal:
        print(f"  ALL position sets universal for k={k}!")

# Test k=6 (single position set: all 6 edge positions)
print(f"\n  k=6: D_{{0,1,2,3,4,5}} = H(T)")
S = tuple(range(6))
mod = 2**(n-1-6)  # = 2^0 = 1
print(f"    mod = 1, so trivially universal")
# But H mod 2 = 1 always (H is always odd), which is non-trivial
DS_vals = set()
for trial in range(20):
    A = random_tournament(n)
    DS = compute_DS(A, n, S)
    DS_vals.add(DS % 2)
print(f"    H mod 2 values: {DS_vals}")
