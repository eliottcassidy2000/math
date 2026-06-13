#!/usr/bin/env python3
"""
large_n_tests_s20fg.py — Extend key claims to n=8,9,10 via smart methods
kind-pasteur-2026-03-24-S20fg

CLAIM 1: Grid-symmetry pairing (PROVE theoretically — no computation)
CLAIM 2: Skip hierarchy (score-change rate as proxy for SL rate)
CLAIM 3: Complement-flip overlap sequence (sample-based search)

Score-change trick: flipping tile (x,y) changes the direction of arc (x,y).
This changes scores[x] and scores[y] by +1 and -1 (or vice versa).
The SORTED score sequence changes iff the resulting sorted scores differ.
Score-change is a NECESSARY condition for class-change (different scores
=> different class). So score-preservation rate is an UPPER BOUND on SL rate.
"""

import sys
import random
from math import factorial, comb
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EXTENDING CLAIMS TO n=8,9,10")
print("  kind-pasteur-2026-03-24-S20fg")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

# ================================================================
# CLAIM 1: Grid-symmetry pairing — THEORETICAL PROOF
# ================================================================
print("\n" + "=" * 60)
print("  CLAIM 1: GRID-SYMMETRY PAIRING (PROOF)")
print("=" * 60)
print("""
THEOREM: For all n >= 3, grid-symmetry paired wiggly classes
produce identical metagraph edge sets.

PROOF:
Let sigma_grid be the vertex permutation i -> n+1-i.
This maps tile (x,y) to tile (n+1-y, n+1-x).

For any tiling T (binary assignment to tiles):
  sigma_grid(T) is the tiling obtained by permuting tile values
  according to the grid symmetry map.

Key fact: sigma_grid is a vertex permutation, so it maps
isomorphic tournaments to isomorphic tournaments. I.e.,
[sigma_grid(T)] = [T] as iso classes.

Now: let X be a wiggly class (tile position), and X' = sigma_grid(X).
For any tiling T:
  - Flipping X in T gives T^X with class [T^X].
  - Flipping X' in sigma_grid(T) gives (sigma_grid(T))^{X'}.

Since sigma_grid maps tile X to X', we have:
  (sigma_grid(T))^{X'} = sigma_grid(T^X).

Therefore: [T^X] = [(sigma_grid(T))^{X'}] = [sigma_grid(T^X)].

The transition T -> T^X in class space is the SAME as
sigma_grid(T) -> (sigma_grid(T))^{X'} in class space.

Since sigma_grid is a bijection on tilings preserving iso classes,
the set of (source_class, target_class) pairs for wiggly class X
is IDENTICAL to that of wiggly class X'.

Therefore: edge_set(X) = edge_set(X'). QED.

This proof works for ALL n. No computation needed.
The grid-symmetry pairing is a THEOREM, not just a computation.
""")

# ================================================================
# CLAIM 2: Skip hierarchy via score-change rate
# ================================================================
print("=" * 60)
print("  CLAIM 2: SKIP HIERARCHY (score-change proxy)")
print("=" * 60)

for n in range(4, 11):
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))

    # For each tile: when flipped, which two vertices change score?
    # Tile (x,y): arc between vertex x and vertex y in VERTS labeling.
    # In adjacency: xi = VERTS.index(x), yi = VERTS.index(y).
    # Flipping changes: scores[xi] by +-1, scores[yi] by -+1.

    # Score-change rate: fraction of random tilings where flipping tile X
    # changes the SORTED score sequence.
    # This is an UPPER BOUND on the SL rate (score-preservation is necessary
    # but not sufficient for class-preservation).

    n_samples = min(100000, 2**m)  # sample size
    if 2**m <= n_samples:
        # Enumerate all
        masks = list(range(2**m))
        n_samples = len(masks)
    else:
        masks = [random.randint(0, 2**m - 1) for _ in range(n_samples)]

    # Precompute tile -> vertex indices
    tile_verts = []
    for i, (x, y) in enumerate(TILES):
        xi = VERTS.index(x)
        yi = VERTS.index(y)
        tile_verts.append((xi, yi))

    # For each tile: count score-preserving flips
    score_preserve = [0] * m

    for mask in masks:
        # Compute scores
        bits = [(mask >> k) & 1 for k in range(m)]

        # Build score array efficiently
        scores = [0] * N
        # Base path contributes: k -> k+1 for k=0..N-2
        for k in range(N-1):
            scores[k] += 1  # k beats k+1

        # Tiles contribute
        for i in range(m):
            xi, yi = tile_verts[i]
            if bits[i] == 0:
                scores[xi] += 1  # xi beats yi
            else:
                scores[yi] += 1  # yi beats xi

        sorted_scores = tuple(sorted(scores))

        # For each tile: flip and check score change
        for wi in range(m):
            xi, yi = tile_verts[wi]
            # Current: if bits[wi]==0, xi beats yi. Flip: yi beats xi.
            # Score change: scores[xi] -= 1, scores[yi] += 1 (if currently 0->1)
            # or: scores[xi] += 1, scores[yi] -= 1 (if currently 1->0)
            new_scores = list(scores)
            if bits[wi] == 0:
                new_scores[xi] -= 1
                new_scores[yi] += 1
            else:
                new_scores[xi] += 1
                new_scores[yi] -= 1

            new_sorted = tuple(sorted(new_scores))
            if new_sorted == sorted_scores:
                score_preserve[wi] += 1

    # Group by skip
    skip_groups = {}
    for i, (x, y) in enumerate(TILES):
        skip = x - y - 1
        if skip not in skip_groups:
            skip_groups[skip] = []
        skip_groups[skip].append(score_preserve[i] / n_samples * 100)

    print(f"\n  n = {n}, m = {m} tiles, samples = {n_samples}")
    print(f"  Score-preserve rate by skip (UPPER BOUND on SL rate):")
    print(f"  {'Skip':>4} {'Count':>5} {'Avg rate':>9} {'Range':>15}")
    for skip in sorted(skip_groups.keys()):
        rates = skip_groups[skip]
        print(f"  {skip:4d} {len(rates):5d} {sum(rates)/len(rates):8.2f}% {min(rates):6.2f}-{max(rates):6.2f}%")

    # Overall: skip=1 vs skip>=2
    skip1_rates = [r for s, rates in skip_groups.items() if s == 1 for r in rates]
    skip2plus_rates = [r for s, rates in skip_groups.items() if s >= 2 for r in rates]
    if skip1_rates and skip2plus_rates:
        avg1 = sum(skip1_rates) / len(skip1_rates)
        avg2 = sum(skip2plus_rates) / len(skip2plus_rates)
        gap = avg2 / avg1 if avg1 > 0 else float('inf')
        print(f"  Skip=1 avg: {avg1:.3f}%, Skip>=2 avg: {avg2:.3f}%, Gap: {gap:.1f}x")

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

# ================================================================
# CLAIM 3: Complement-flip overlap at n=8
# ================================================================
print("\n" + "=" * 60)
print("  CLAIM 3: COMPLEMENT-FLIP OVERLAP (sample search)")
print("=" * 60)

for n in [8]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
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

    def adj_complement(A):
        return [[1-A[i][j] if i!=j else 0 for j in range(N)] for i in range(N)]

    def score_seq(A):
        return tuple(sorted(sum(row) for row in A))

    # For complement-flip overlap: T and complement(T) must be isomorphic
    # after flipping one arc. This means: Hamming distance between T and
    # some relabeling of complement(T) is exactly 1.
    #
    # Necessary condition: T and complement(T) have the same score sequence
    # (since isomorphic tournaments have same scores).
    # complement(T) has scores (n-1-s_i) for each vertex i.
    # Sorted complement scores = (n-1-s_{max}, ..., n-1-s_{min}).
    # For sorted scores to match: s_i + s_{n-1-i} = n-1 for all i.
    # This is the PALINDROMIC score condition.
    #
    # So: only tournaments with palindromic scores CAN be in overlap pairs.
    # And even among those: the class must NOT be self-complementary,
    # AND some relabeling of the complement must be Hamming-1 away.

    print(f"\n  n = {n}, m = {m}")
    print(f"  Strategy: sample random palindromic-score tournaments,")
    print(f"  check if complement is Hamming-1 reachable (via defect counting)")

    n_palindromic = 0
    n_almost_anti_aut = 0
    n_samples = 50000

    for _ in range(n_samples):
        # Random tiling
        mask = random.randint(0, 2**m - 1)
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        sc = score_seq(A)

        # Check palindromic
        is_palindromic = all(sc[i] + sc[N-1-i] == N-1 for i in range(N//2))
        if not is_palindromic:
            continue
        n_palindromic += 1

        # Check if complement is single-flip reachable
        # Build complement
        Ac = adj_complement(A)

        # For each permutation sigma (EXPENSIVE — sample a subset)
        # Check: Hamming(A, sigma(Ac)) = ? Want it to be 2 (1 arc flip = 2 adj entries)
        # Sample random permutations
        import itertools
        perm_sample = random.sample(range(factorial(N)), min(200, factorial(N)))

        found = False
        for perm_idx in perm_sample:
            # Generate permutation from index (Lehmer code)
            remaining = list(range(N))
            sigma = []
            idx = perm_idx
            for k in range(N, 0, -1):
                q, idx = divmod(idx, factorial(k-1))
                if q >= len(remaining): q = len(remaining) - 1
                sigma.append(remaining[q])
                remaining.pop(q)

            # Compute Hamming distance
            hamming = 0
            for i in range(N):
                for j in range(N):
                    if i != j and A[i][j] != Ac[sigma[i]][sigma[j]]:
                        hamming += 1
            # Want hamming == 2 (one arc flip = 2 matrix entries)
            if hamming == 2:
                n_almost_anti_aut += 1
                found = True
                break

    print(f"  Samples: {n_samples}")
    print(f"  Palindromic score: {n_palindromic} ({n_palindromic/n_samples*100:.1f}%)")
    print(f"  Almost-anti-automorphisms found: {n_almost_anti_aut}")
    print(f"  (Each means: complement is single-flip reachable)")

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
