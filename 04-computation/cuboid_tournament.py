#!/usr/bin/env python3
"""
Cuboid Tournament Investigation
opus-2026-03-14-S71f

Simplex tournaments have I(Ω,x) = (1+x)^m, giving H = 3^m.
Question: Do "cuboid tournaments" exist with I(Ω,x) = (1+2x)^m, giving H = 5^m?

For m=1: H=5, I=(1+2x). Need α₁=2, α₂=0.
  → 2 directed 3-cycles, 0 disjoint pairs. So the 2 cycles must share a vertex.
  → This means n≤5 (two triangles sharing a vertex on 5 vertices).

For m=2: H=25, I=(1+2x)² = 1+4x+4x². Need α₁=4, α₂=4.
  → 4 directed 3-cycles, 4 disjoint pairs.

Let's check computationally.
"""
import random
from itertools import permutations, combinations

def count_hp(A, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            count += 1
        if A[i][k] and A[k][j] and A[j][i]:
            count += 1
    return count

def find_all_3cycles(A, n):
    """Find all directed 3-cycles as frozensets of vertices."""
    cycles = []
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.append(frozenset([i,j,k]))
        if A[i][k] and A[k][j] and A[j][i]:
            cycles.append(frozenset([i,j,k]))
    return cycles

def compute_alpha(A, n):
    """Compute α₁, α₂ from 3-cycle structure."""
    cycles = find_all_3cycles(A, n)
    # Remove duplicates (each undirected triangle appears once as a directed 3-cycle)
    unique = list(set(cycles))
    alpha1 = len(unique)
    
    # Count disjoint pairs
    alpha2 = 0
    for i in range(len(unique)):
        for j in range(i+1, len(unique)):
            if unique[i].isdisjoint(unique[j]):
                alpha2 += 1
    return alpha1, alpha2

# Part 1: Check n=5 for H=5 (cuboid m=1)
print("="*60)
print("CUBOID TOURNAMENT SEARCH")
print("="*60)

# m=1: H=5, need α₁=2, α₂=0
print("\n--- m=1: Looking for H=5 with α₁=2, α₂=0 ---")
found_cuboid1 = []
for n in [4, 5]:
    count = 0
    total = 0
    for trial in range(100000):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        h = count_hp(A, n)
        if h == 5:
            a1, a2 = compute_alpha(A, n)
            if a1 == 2 and a2 == 0:
                count += 1
                if count <= 3:
                    found_cuboid1.append((n, A))
            total += 1
    print(f"  n={n}: H=5 found {total}/100k, with (α₁=2,α₂=0): {count}")

# Part 2: Check the I(Ω,x) decomposition for H=5 tournaments  
print("\n--- Analyzing H=5 tournaments ---")
for n in [4, 5]:
    for trial in range(200000):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        h = count_hp(A, n)
        if h == 5:
            a1, a2 = compute_alpha(A, n)
            print(f"  n={n}: H={h}, α₁={a1}, α₂={a2}, H_check={1+2*a1+4*a2}")
            break

# Part 3: m=2: H=25, need α₁=4, α₂=4
print("\n--- m=2: Looking for H=25 with α₁=4, α₂=4 ---")
for n in [6, 7]:
    count = 0
    total_h25 = 0
    alpha_dist = {}
    for trial in range(200000):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        h = count_hp(A, n)
        if h == 25:
            a1, a2 = compute_alpha(A, n)
            key = (a1, a2)
            alpha_dist[key] = alpha_dist.get(key, 0) + 1
            if a1 == 4 and a2 == 4:
                count += 1
            total_h25 += 1
    print(f"  n={n}: H=25 found {total_h25}/200k, with (4,4): {count}")
    if alpha_dist:
        print(f"    α distributions: {sorted(alpha_dist.items())}")

# Part 4: What I(Ω,x) looks like for H=25
print("\n--- I(Ω,x) for H=25 at n=6 ---")
found = 0
for trial in range(500000):
    A = [[0]*6 for _ in range(6)]
    for i in range(6):
        for j in range(i+1, 6):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    h = count_hp(A, 6)
    if h == 25:
        a1, a2 = compute_alpha(A, 6)
        print(f"  Found H=25: α₁={a1}, α₂={a2}, I(Ω,x)={a2}x²+{a1}x+1, I(Ω,2)={1+2*a1+4*a2}")
        found += 1
        if found >= 10:
            break

# Part 5: Check if (1+2x)^m factorization exists for any H
print("\n--- Checking (1+2x)^m pattern ---")
print("(1+2x)^1 = 1+2x → α₁=2, α₂=0 → H=5")
print("(1+2x)^2 = 1+4x+4x² → α₁=4, α₂=4 → H=25")
print("(1+2x)^3 = 1+6x+12x²+8x³ → α₁=6, α₂=12, α₃=8 → H=125")

# Part 6: General (1+cx)^m analysis
print("\n--- General (1+cx)^m pattern ---")
for c in [1, 2, 3]:
    print(f"\n  (1+{c}x)^m:")
    for m in range(1, 5):
        from math import comb
        alphas = [comb(m, k) * c**k for k in range(1, m+1)]
        H = sum(2**k * comb(m, k) * c**k for k in range(m+1))
        print(f"    m={m}: H={H}, αs={alphas}")
        # Check: H should be (1+2c)^m
        assert H == (1+2*c)**m, f"Mismatch: {H} != {(1+2*c)**m}"

print("\n\n" + "="*60)
print("SIMPLEX-CUBOID-TESSERACT HIERARCHY")
print("="*60)
print(f"Simplex  (1+x)^m:  H = 3^m  ({[3**m for m in range(1,6)]})")
print(f"Cuboid   (1+2x)^m: H = 5^m  ({[5**m for m in range(1,6)]})")
print(f"Tesseract(1+3x)^m: H = 7^m  ({[7**m for m in range(1,6)]})")
print(f"  But H=7 is FORBIDDEN! So tesseracts DON'T EXIST for m=1!")
print(f"  What about m≥2? (1+3x)^2 = 1+6x+9x² → H=49")

# Part 7: Check H=49 at n=7
print("\n--- Can H=49 (= 7²) exist? ---")
count_49 = 0
for trial in range(200000):
    A = [[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(i+1, 7):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    h = count_hp(A, 7)
    if h == 49:
        count_49 += 1
        if count_49 <= 3:
            a1, a2 = compute_alpha(A, 7)
            print(f"  H=49: α₁={a1}, α₂={a2}")
print(f"  H=49 at n=7: found {count_49}/200k")

# Part 8: What about mixed products?
print("\n--- Mixed products: (1+x)(1+2x) = 1+3x+2x² → H=15 ---")
# H=15 = 3*5 = simplex × cuboid
print("  (1+x)(1+2x): α₁=3, α₂=2 → H=15")
print("  (1+x)²(1+2x): α₁=4, α₂=5, α₃=2 → H=45 = 9*5")
print("  (1+x)(1+2x)²: α₁=5, α₂=10, α₃=4 → H=75 = 3*25")

