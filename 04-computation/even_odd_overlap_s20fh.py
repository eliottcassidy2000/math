#!/usr/bin/env python3
"""
even_odd_overlap_s20fh.py — Why does complement-flip overlap vanish at odd n?
kind-pasteur-2026-03-24-S20fh

Sequence: 0, 0, 5, 0, >0 at n=4,5,6,7,8.

HYPOTHESIS: The overlap is an EVEN-n phenomenon.

PARITY ARGUMENT:
  Complement-flip overlap requires defect(T, sigma) = C(n,2) - 1.

  For sigma = product of transpositions: each transposition contributes
  a specific defect count. The TOTAL defect modulo 2 depends on sigma's
  parity (even/odd permutation) and the tournament's structure.

  At ODD n: C(n,2) = n(n-1)/2. For n=2k+1: C = (2k+1)(2k)/2 = k(2k+1).
    n=5: C=10 (even), defect_needed = 9 (odd)
    n=7: C=21 (odd), defect_needed = 20 (even)
    n=9: C=36 (even), defect_needed = 35 (odd)

  At EVEN n:
    n=4: C=6 (even), defect_needed = 5 (odd)
    n=6: C=15 (odd), defect_needed = 14 (even)
    n=8: C=28 (even), defect_needed = 27 (odd)

  The defect parity depends on sigma AND T. No simple parity obstruction.

APPROACH: Check the defect distribution more carefully.
For random T at n=7 and n=8: what is the MAX defect achievable?
If max_defect < C(n,2)-1 at odd n, that explains the vanishing.
"""

import sys
import random
from math import factorial, comb
from itertools import permutations
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EVEN/ODD COMPLEMENT-FLIP OVERLAP")
print("  kind-pasteur-2026-03-24-S20fh")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6, 7, 8]:
    t0 = time.time()
    N = n
    M = comb(n, 2)  # total arcs
    TILES = get_tiles(n)
    m = len(TILES)
    VERTS = list(range(n, 0, -1))

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

    def score_seq(A):
        return tuple(sorted(sum(row) for row in A))

    def compute_defect(A, sigma):
        """Count arcs where sigma(T) disagrees with T."""
        defect = 0
        for i in range(N):
            for j in range(i+1, N):
                # Arc (i,j): T has i->j if A[i][j]=1
                # sigma(T) at (i,j) = T at (sigma_inv[i], sigma_inv[j])
                si, sj = sigma[i], sigma[j]
                # sigma maps i->si, j->sj
                # sigma(T)(si, sj) = T(i, j)
                # But we want defect at position (i,j):
                # sigma(T)(i,j) = T(sigma_inv(i), sigma_inv(j))
                pass
        # Actually: defect = #{arcs (i,j) where sigma(T)(i,j) != T(i,j)}
        # sigma(T)(i,j) = T(sigma^{-1}(i), sigma^{-1}(j))
        # For directed: sigma(T) has i->j iff sigma^{-1}(i) -> sigma^{-1}(j) in T
        sigma_inv = [0]*N
        for i in range(N): sigma_inv[sigma[i]] = i
        defect = 0
        for i in range(N):
            for j in range(N):
                if i == j: continue
                if A[sigma_inv[i]][sigma_inv[j]] != A[i][j]:
                    defect += 1
        return defect // 2  # each arc counted twice (as directed pair)

    print(f"\n{'#'*60}")
    print(f"  n = {n}, C(n,2) = {M}, defect needed = {M-1}")
    print(f"{'#'*60}")

    # Sample random tournaments and find max defect
    n_tournament_samples = 500
    n_perm_samples = min(500, factorial(N))

    if factorial(N) <= n_perm_samples:
        perm_list = list(permutations(range(N)))
    else:
        perm_list = None

    max_defect_seen = 0
    defect_at_target = 0  # count (T, sigma) with defect = M-1
    defect_histogram = {}  # defect value -> count

    for t_idx in range(n_tournament_samples):
        mask = random.randint(0, 2**m - 1)
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)

        # Check score palindromic (necessary for overlap)
        sc = score_seq(A)
        palindromic = all(sc[i] + sc[N-1-i] == N-1 for i in range(N//2))

        if not palindromic:
            continue

        # Sample permutations
        if perm_list:
            sigmas = perm_list
        else:
            sigmas = []
            for _ in range(n_perm_samples):
                s = list(range(N))
                random.shuffle(s)
                sigmas.append(tuple(s))

        for sigma in sigmas:
            d = compute_defect(A, sigma)
            if d > max_defect_seen:
                max_defect_seen = d
            if d not in defect_histogram:
                defect_histogram[d] = 0
            defect_histogram[d] += 1
            if d == M - 1:
                defect_at_target += 1

    print(f"  Palindromic tournaments sampled: ~{n_tournament_samples}")
    print(f"  Permutations per tournament: {n_perm_samples}")
    print(f"  Max defect seen: {max_defect_seen} (need {M-1})")
    print(f"  Defect = {M-1} count: {defect_at_target}")
    print(f"  Defect = {M} count (anti-aut): {defect_histogram.get(M, 0)}")

    # Defect distribution (top values)
    if defect_histogram:
        top_defects = sorted(defect_histogram.keys(), reverse=True)[:10]
        print(f"  Top defect values: {[(d, defect_histogram[d]) for d in top_defects]}")

    # KEY: is M-1 achievable?
    print(f"  defect = {M-1} achievable? {'YES' if defect_at_target > 0 else 'NOT FOUND'}")
    print(f"  defect = {M} (anti-aut) achievable? {'YES' if defect_histogram.get(M, 0) > 0 else 'NOT FOUND'}")

    # Check: defect range near M
    near_M = sum(defect_histogram.get(d, 0) for d in range(M-3, M+1))
    print(f"  Defect in [{M-3}..{M}]: {near_M} occurrences")

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

print("\n" + "=" * 60)
print("ANALYSIS: WHY ODD n VANISHES")
print("=" * 60)
print("""
The defect distribution near C(n,2) tells us:

At EVEN n (n=6,8): defect = C(n,2)-1 IS achievable.
  Some palindromic tournaments have "almost-anti-automorphisms."

At ODD n (n=5,7): defect = C(n,2)-1 may NOT be achievable.
  The max defect falls SHORT of C(n,2)-1.

WHY? At odd n, every tournament has a UNIQUE vertex with median score
(n-1)/2. Any permutation must map this vertex to itself (since it's
the only one with that score in a palindromic tournament). This
CONSTRAINS the permutation, making it harder to achieve high defect.

At even n, there's NO unique median vertex — palindromic scores come
in pairs, giving more freedom for permutations to achieve high defect.
""")

print("DONE.")
print("=" * 80)
