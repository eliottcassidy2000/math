#!/usr/bin/env python3
"""
sc_vs_overlap_s20fi.py — Is the odd-n vanishing explained by SC filtering?
kind-pasteur-2026-03-24-S20fi

HYPOTHESIS: At odd n, defect = M-1 occurs ONLY for SC tournaments.
If T is SC, [T^c] = [T], so the "overlap" is just a self-loop.
Only NON-SC tournaments with defect = M-1 create actual overlap edges.

At n=5: 360 cases of defect=9. Are they ALL for SC tournaments?
At n=7: 4 cases of defect=20. Are they ALL for SC tournaments?
At n=6: 305 cases of defect=14. Some must be non-SC (we found 5 overlap pairs).
"""

import sys
import random
from math import factorial, comb
from itertools import permutations
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  SC vs NON-SC IN HIGH-DEFECT CASES")
print("  kind-pasteur-2026-03-24-S20fi")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6, 7]:
    t0 = time.time()
    N = n
    M = comb(n, 2)
    TILES = get_tiles(n)
    m = len(TILES)
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    tile_verts = []
    for i, (x, y) in enumerate(TILES):
        tile_verts.append((VERTS.index(x), VERTS.index(y)))

    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tile_verts[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def adj_complement(A):
        return [[1-A[i][j] if i!=j else 0 for j in range(N)] for i in range(N)]

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def is_SC(A):
        cn = canonicalize(A)
        cn_comp = canonicalize(adj_complement(A))
        return cn == cn_comp

    def compute_defect(A, sigma):
        sigma_inv = [0]*N
        for i in range(N): sigma_inv[sigma[i]] = i
        defect = 0
        for i in range(N):
            for j in range(i+1, N):
                if A[sigma_inv[i]][sigma_inv[j]] != A[i][j]:
                    defect += 1
        return defect

    print(f"\n{'#'*60}")
    print(f"  n = {n}, M = C({n},2) = {M}, target defect = {M-1}")
    print(f"{'#'*60}")

    # Enumerate ALL tilings at n=5,6 (small enough)
    # At n=7: sample
    if n <= 6:
        masks = list(range(1 << m))
    else:
        masks = [random.randint(0, 2**m - 1) for _ in range(10000)]

    sc_at_target = 0
    nsc_at_target = 0
    total_checked = 0

    for mask in masks:
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)

        # Check if palindromic scores (necessary)
        scores = tuple(sorted(sum(A[i]) for i in range(N)))
        palindromic = all(scores[i] + scores[N-1-i] == N-1 for i in range(N//2))
        if not palindromic:
            continue

        total_checked += 1

        # Check defect for all permutations (or sample)
        if n <= 6:
            sigmas = all_perms
        else:
            sigmas = [tuple(random.sample(range(N), N)) for _ in range(200)]

        found_target = False
        for sigma in sigmas:
            d = compute_defect(A, sigma)
            if d == M - 1:
                found_target = True
                break

        if found_target:
            if is_SC(A):
                sc_at_target += 1
            else:
                nsc_at_target += 1

    print(f"  Palindromic tilings checked: {total_checked}")
    print(f"  Tilings with defect = {M-1}:")
    print(f"    Self-complementary (SC): {sc_at_target}")
    print(f"    Non-self-complementary (NSC): {nsc_at_target}")
    print(f"    Total: {sc_at_target + nsc_at_target}")

    if sc_at_target + nsc_at_target > 0:
        print(f"    SC fraction: {sc_at_target / (sc_at_target + nsc_at_target) * 100:.1f}%")
    print(f"    NSC count = actual overlap edges (new!)")

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

print("\n" + "=" * 60)
print("CONCLUSION")
print("=" * 60)
print("""
If NSC count = 0 at odd n (5, 7): The odd-n vanishing is explained
by SC filtering. All high-defect cases at odd n are SC tournaments,
where the complement IS the same class (self-loop, not new edge).

If NSC count > 0 at odd n: The vanishing has a different cause.
Some non-SC tournaments achieve defect = M-1, but these specific
classes happen to not overlap (the complement class is NOT connected
by a wiggly edge to the original class).
""")

print("DONE.")
print("=" * 80)
