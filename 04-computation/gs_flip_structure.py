#!/usr/bin/env python3
"""
Grid-Symmetric (GS) Tiling Flip Structure.
Instance: opus-2026-03-06-S11

Key theorem (proved below): GS tilings are closed under flip because
flip (complementing all bits) preserves the pairwise equality condition
that defines grid symmetry.

This script investigates:
1. The GS-restricted flip pairing
2. How GS flip pairs distribute across isomorphism classes
3. The "GS skeleton" — the class-level structure within the GS subspace
4. Connection to self-converse (SC) tournaments
5. Cross-scale patterns in the GS structure
"""
from itertools import permutations
from collections import defaultdict, Counter
import sys

def build_gs_data(n):
    """Build GS-restricted tiling data."""
    verts = list(range(n, 0, -1))
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}
    flip_mask = (1 << m) - 1

    # Transpose map
    trans_map = []
    fixed_tiles = []  # tiles mapped to themselves
    paired_tiles = []  # pairs of tiles mapped to each other
    seen = set()
    for i, (x, y) in enumerate(tiles):
        nx, ny = n - y + 1, n - x + 1
        j = tile_idx.get((nx, ny), -1)
        trans_map.append(j)
        if j == i:
            fixed_tiles.append(i)
        elif j >= 0 and i not in seen:
            paired_tiles.append((i, j))
            seen.add(i)
            seen.add(j)

    gs_dim = len(fixed_tiles) + len(paired_tiles)
    print(f"  Tile structure: {m} total = {len(fixed_tiles)} fixed + {2*len(paired_tiles)} paired")
    print(f"  GS dimension: {gs_dim}, so #GS = 2^{gs_dim} = {2**gs_dim}")

    def is_grid_sym(bits):
        for i in range(m):
            j = trans_map[i]
            if j >= 0 and j != i and ((bits >> i) & 1) != ((bits >> j) & 1):
                return False
        return True

    def bits_to_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1):
            A[k][k+1] = 1
        for i, (xL, yL) in enumerate(tiles):
            xi = verts.index(xL)
            yi = verts.index(yL)
            if (bits >> i) & 1 == 0:
                A[xi][yi] = 1
            else:
                A[yi][xi] = 1
        return tuple(tuple(row) for row in A)

    def canonicalize(A):
        best = None
        for p in permutations(range(n)):
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def scores(A):
        return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))

    def hamiltonian_paths(A):
        dp = {}
        for v in range(n):
            dp[(1 << v, v)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or (mask, v) not in dp:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        key = (mask | (1 << u), u)
                        dp[key] = dp.get(key, 0) + dp[(mask, v)]
        full = (1 << n) - 1
        return sum(dp.get((full, v), 0) for v in range(n))

    # Build class data for ALL tilings
    classes = {}
    class_members = defaultdict(list)
    class_scores = {}
    class_H = {}
    tiling_class = {}
    cid_counter = 0

    for bits in range(1 << m):
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        if canon not in classes:
            classes[canon] = cid_counter
            class_scores[cid_counter] = scores(A)
            class_H[cid_counter] = hamiltonian_paths(A)
            cid_counter += 1
        cid = classes[canon]
        class_members[cid].append(bits)
        tiling_class[bits] = cid

    num_classes = cid_counter

    # Collect GS tilings
    gs_tilings = [b for b in range(1 << m) if is_grid_sym(b)]

    # GS flip pairs
    gs_flip_pairs = []
    gs_seen = set()
    for b in gs_tilings:
        if b in gs_seen:
            continue
        fb = b ^ flip_mask
        assert is_grid_sym(fb), f"Flip of GS tiling {b} is not GS!"
        if fb == b:
            # Self-flip GS tiling (impossible since m > 0)
            print(f"  WARNING: self-flip GS tiling found: {b}")
        else:
            gs_flip_pairs.append((b, fb))
            gs_seen.add(b)
            gs_seen.add(fb)

    print(f"  GS tilings: {len(gs_tilings)}, GS flip pairs: {len(gs_flip_pairs)}")

    # GS flip pairs across classes
    gs_cross_class = 0
    gs_same_class = 0
    gs_pair_class_map = []
    for (b1, b2) in gs_flip_pairs:
        c1 = tiling_class[b1]
        c2 = tiling_class[b2]
        gs_pair_class_map.append((c1, c2))
        if c1 == c2:
            gs_same_class += 1
        else:
            gs_cross_class += 1

    print(f"  GS pairs within same class: {gs_same_class}")
    print(f"  GS pairs crossing classes: {gs_cross_class}")

    # GS class structure
    gs_class_count = defaultdict(int)
    for b in gs_tilings:
        gs_class_count[tiling_class[b]] += 1

    gs_classes = sorted(gs_class_count.keys())
    print(f"  Classes containing GS tilings: {len(gs_classes)}")

    # The "GS skeleton": for each pair of GS-containing classes, how many GS pairs cross between them?
    gs_cross_matrix = defaultdict(int)
    gs_self_matrix = defaultdict(int)
    for (c1, c2) in gs_pair_class_map:
        if c1 == c2:
            gs_self_matrix[c1] += 1
        else:
            key = (min(c1, c2), max(c1, c2))
            gs_cross_matrix[key] += 1

    print(f"\n  --- GS-containing classes ---")
    print(f"  {'Class':>6} {'GS':>4} {'Size':>5} {'H':>4} {'SelfPairs':>10} {'Scores'}")
    for c in gs_classes:
        gc = gs_class_count[c]
        sz = len(class_members[c])
        H = class_H[c]
        sp = gs_self_matrix.get(c, 0)
        sc = class_scores[c]
        print(f"  {c:>6} {gc:>4} {sz:>5} {H:>4} {sp:>10} {sc}")

    if gs_cross_matrix:
        print(f"\n  --- GS cross-class flip pairs ---")
        for (c1, c2), count in sorted(gs_cross_matrix.items()):
            sc1 = class_scores[c1]
            sc2 = class_scores[c2]
            print(f"  Class {c1} <-> Class {c2}: {count} pair(s)")
            print(f"    {c1}: scores={sc1}, GS={gs_class_count[c1]}/{len(class_members[c1])}")
            print(f"    {c2}: scores={sc2}, GS={gs_class_count[c2]}/{len(class_members[c2])}")

    # Hamming weight distribution of GS tilings
    hw_dist = Counter()
    for b in gs_tilings:
        hw = bin(b).count('1')
        hw_dist[hw] += 1
    print(f"\n  Hamming weight distribution of GS tilings:")
    for hw in sorted(hw_dist):
        print(f"    hw={hw}: {hw_dist[hw]} tilings")

    # Flip distance from perpendicular midpoint
    # The "perpendicular midpoint" is at Hamming weight m/2
    # GS tilings at different distances from midpoint
    print(f"  (Midpoint at hw = {m/2})")

    return {
        'n': n, 'm': m, 'gs_dim': gs_dim,
        'gs_tilings': gs_tilings,
        'gs_flip_pairs': gs_flip_pairs,
        'gs_classes': gs_classes,
        'gs_class_count': dict(gs_class_count),
        'gs_same_class': gs_same_class,
        'gs_cross_class': gs_cross_class,
        'gs_cross_matrix': dict(gs_cross_matrix),
        'num_classes': num_classes,
        'class_H': class_H,
        'class_scores': class_scores,
        'class_members': dict(class_members),
    }


def main():
    results = []
    for n in range(3, 7):
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")
        r = build_gs_data(n)
        results.append(r)

    # Cross-scale summary
    print(f"\n{'='*60}")
    print("CROSS-SCALE GS SUMMARY")
    print(f"{'='*60}")
    print(f"{'n':>3} {'m':>3} {'GS-dim':>7} {'#GS':>5} {'#pairs':>7} {'same':>5} {'cross':>6} {'#GS-cls':>8}")
    for r in results:
        print(f"{r['n']:>3} {r['m']:>3} {r['gs_dim']:>7} {len(r['gs_tilings']):>5} "
              f"{len(r['gs_flip_pairs']):>7} {r['gs_same_class']:>5} {r['gs_cross_class']:>6} "
              f"{len(r['gs_classes']):>8}")

    # Pattern: GS-dim sequence
    print(f"\nGS dimension sequence: {[r['gs_dim'] for r in results]}")
    print(f"(Should be: floor(m/2) + fixed tiles)")

    # Key question: what fraction of GS pairs cross classes?
    print(f"\nFraction of GS pairs crossing classes:")
    for r in results:
        total = len(r['gs_flip_pairs'])
        if total > 0:
            frac = r['gs_cross_class'] / total
            print(f"  n={r['n']}: {r['gs_cross_class']}/{total} = {frac:.3f}")


if __name__ == '__main__':
    main()
