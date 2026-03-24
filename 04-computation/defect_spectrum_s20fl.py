#!/usr/bin/env python3
"""
defect_spectrum_s20fl.py — Full defect spectrum D(d) for n=4,5,6
kind-pasteur-2026-03-24-S20fl

D(d) = #{(labeled T, sigma in S_n) : defect(T, sigma) = d}

Properties:
  D(0) = sum_sigma Fix(sigma) = n! * V(n) (Burnside)
  D(1) = my formula (kind-pasteur-S20el)
  D(M) = ? (anti-automorphisms, related to SC count)
  sum D(d) = n! * 2^{C(n,2)} (total pairs)
  D(d) = D(d) under complement (defect is complement-invariant)
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FULL DEFECT SPECTRUM")
print("  kind-pasteur-2026-03-24-S20fl")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [4, 5, 6]:
    t0 = time.time()
    N = n
    M = comb(n, 2)
    TILES = get_tiles(n)
    m = len(TILES)
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    print(f"\n{'#'*60}")
    print(f"  n = {n}, M = C({n},2) = {M}, 2^m = {2**m}")
    print(f"{'#'*60}")

    # Compute full defect spectrum
    spectrum = Counter()

    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)

        for sigma in all_perms:
            si = [0]*N
            for i in range(N): si[sigma[i]] = i

            d = 0
            for i in range(N):
                for j in range(i+1, N):
                    if A[si[i]][si[j]] != A[i][j]:
                        d += 1
            spectrum[d] += 1

    total_pairs = sum(spectrum.values())
    expected_total = factorial(N) * (2 ** comb(N, 2))

    print(f"\n  Total (T, sigma) pairs: {total_pairs} (expected {expected_total})")
    print(f"  {'Defect':>6} {'Raw count':>12} {'Count/n!':>10} {'Fraction':>10}")

    for d in range(M+1):
        count = spectrum.get(d, 0)
        per_nf = count / factorial(N)
        frac = count / total_pairs
        marker = ""
        if d == 0: marker = " <- Fix(sigma) = n!*V"
        elif d == 1: marker = " <- D(n)"
        elif d == M-1: marker = " <- anti-D"
        elif d == M: marker = " <- anti-aut (SC-related)"
        print(f"  {d:6d} {count:12d} {per_nf:10.1f} {frac:10.6f}{marker}")

    # Symmetry check: is D(d) = D(M-d)?
    print(f"\n  SYMMETRY CHECK: D(d) vs D(M-d)")
    symmetric = True
    for d in range(M+1):
        if spectrum.get(d, 0) != spectrum.get(M-d, 0):
            symmetric = False
            print(f"    D({d}) = {spectrum.get(d,0)} vs D({M-d}) = {spectrum.get(M-d,0)}: ASYMMETRIC")
    if symmetric:
        print(f"    SYMMETRIC! D(d) = D(M-d) for all d.")
    else:
        # What is the asymmetry ratio?
        ratios = []
        for d in range(M//2 + 1):
            a = spectrum.get(d, 0)
            b = spectrum.get(M-d, 0)
            if b > 0:
                ratios.append(a/b)
        if ratios:
            print(f"    Asymmetry ratios D(d)/D(M-d): {[f'{r:.2f}' for r in ratios[:6]]}")

    # Even/odd pattern
    even_sum = sum(spectrum.get(d, 0) for d in range(0, M+1, 2))
    odd_sum = sum(spectrum.get(d, 0) for d in range(1, M+1, 2))
    print(f"\n  Even defect sum: {even_sum} ({even_sum/total_pairs*100:.1f}%)")
    print(f"  Odd defect sum:  {odd_sum} ({odd_sum/total_pairs*100:.1f}%)")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
