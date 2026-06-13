#!/usr/bin/env python3
"""
Tiling isomorphism analysis at n=7, optimized for speed.
Focus: SC+SF chain, H-maximizer symmetry, grid-symmetric distribution.

Uses nauty-like canonical form via sorted adjacency comparison.
Instance: opus-2026-03-06-S7
"""

from itertools import permutations
from collections import defaultdict, Counter
import time

def build_n7_data():
    n = 7
    verts = list(range(n, 0, -1))

    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)  # = 15
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}

    trans_map = []
    for (x, y) in tiles:
        nx, ny = n - y + 1, n - x + 1
        trans_map.append(tile_idx[(nx, ny)])

    def is_grid_sym(mask):
        bits = [(mask >> k) & 1 for k in range(m)]
        for i in range(m):
            if bits[i] != bits[trans_map[i]]:
                return False
        return True

    def mask_to_adj(mask):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1):
            A[k][k+1] = 1
        for i, (xL, yL) in enumerate(tiles):
            xi = verts.index(xL)
            yi = verts.index(yL)
            if (mask >> i) & 1 == 0:
                A[xi][yi] = 1
            else:
                A[yi][xi] = 1
        return A

    # Precompute vertex index lookups
    vert_to_idx = {v: i for i, v in enumerate(verts)}

    # Use score sequence as a fast pre-filter for canonicalization
    def scores(A):
        return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))

    def canon_key(A):
        """Canonical form: minimum adjacency string over all permutations."""
        best = None
        for p in permutations(range(n)):
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def count_ham_paths(A):
        count = 0
        for p in permutations(range(n)):
            valid = True
            for k in range(n-1):
                if A[p[k]][p[k+1]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
        return count

    print(f"n=7, m={m}, 2^m = {1<<m} tilings")
    print("Phase 1: Grouping tilings by score sequence...")
    t0 = time.time()

    # Group masks by score sequence (fast pre-filter)
    score_groups = defaultdict(list)
    gs_count = 0
    for mask in range(1 << m):
        A = mask_to_adj(mask)
        sc = scores(A)
        score_groups[sc].append(mask)
        if is_grid_sym(mask):
            gs_count += 1

    print(f"  {len(score_groups)} distinct score sequences, {gs_count} grid-symmetric tilings")
    print(f"  Time: {time.time()-t0:.1f}s")

    print("Phase 2: Canonical forms within each score group...")
    t1 = time.time()

    # Within each score group, canonicalize
    classes = defaultdict(list)  # canon_key -> list of masks
    processed = 0
    for sc, masks in sorted(score_groups.items(), key=lambda x: -len(x[1])):
        # For large groups, we need full canonicalization
        for mask in masks:
            A = mask_to_adj(mask)
            ck = canon_key(A)
            classes[ck].append(mask)
            processed += 1
            if processed % 5000 == 0:
                print(f"  Processed {processed}/{1<<m} ({100*processed/(1<<m):.1f}%)")

    print(f"  {len(classes)} isomorphism classes")
    print(f"  Time: {time.time()-t1:.1f}s")

    print("Phase 3: Computing class properties...")
    t2 = time.time()

    # Compute class-level data
    class_list = []
    for ci, (ck, masks) in enumerate(sorted(classes.items())):
        rep_A = mask_to_adj(masks[0])
        H = count_ham_paths(rep_A)
        sc = scores(rep_A)

        # Count automorphisms
        aut_count = sum(1 for p in permutations(range(n))
                        if all(rep_A[p[i]][p[j]] == rep_A[i][j]
                               for i in range(n) for j in range(n)))

        # Count anti-automorphisms
        anti_aut_count = sum(1 for p in permutations(range(n))
                             if all(rep_A[p[i]][p[j]] == rep_A[j][i]
                                    for i in range(n) for j in range(n)))

        # Grid-symmetric members
        gs_masks = [mask for mask in masks if is_grid_sym(mask)]

        # Self-converse check via canonical form of transpose
        A_op = [[rep_A[j][i] for j in range(n)] for i in range(n)]
        ck_op = canon_key(A_op)
        is_sc = ck_op == ck

        # Flip targets
        flip_masks_in_class = set()
        for mask in masks:
            flip_mask = mask ^ ((1 << m) - 1)
            flip_masks_in_class.add(flip_mask)
        # Check if any flip lands back in same class
        is_sf = bool(set(masks) & flip_masks_in_class)

        class_list.append({
            'ci': ci,
            'ck': ck,
            'masks': masks,
            'size': len(masks),
            'H': H,
            'scores': sc,
            'aut': aut_count,
            'anti_aut': anti_aut_count,
            'gs_count': len(gs_masks),
            'is_sc': is_sc,
            'is_sf': is_sf,
        })

        if ci % 50 == 0 and ci > 0:
            print(f"  Class {ci}...")

    print(f"  Time: {time.time()-t2:.1f}s")
    return class_list, m


def analyze_n7(class_list, m):
    n = 7
    print(f"\n{'='*70}")
    print(f"n=7 ANALYSIS: {len(class_list)} isomorphism classes")
    print(f"{'='*70}")

    # Summary counts
    sc_count = sum(1 for c in class_list if c['is_sc'])
    sf_count = sum(1 for c in class_list if c['is_sf'])
    both_count = sum(1 for c in class_list if c['is_sc'] and c['is_sf'])
    high_aut = sum(1 for c in class_list if c['aut'] > 1)
    has_anti = sum(1 for c in class_list if c['anti_aut'] > 0)

    print(f"\nSelf-converse: {sc_count}")
    print(f"Self-flip: {sf_count}")
    print(f"Both SC+SF: {both_count}")
    print(f"|Aut| > 1: {high_aut}")
    print(f"Has anti-automorphisms: {has_anti}")

    # SC+SF classes (the symmetry kernel)
    print(f"\n--- SC+SF classes (symmetry kernel) ---")
    for c in class_list:
        if c['is_sc'] and c['is_sf']:
            print(f"  #{c['ci']}: H={c['H']}, |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, "
                  f"scores={c['scores']}, size={c['size']}, GS={c['gs_count']}/{c['size']}")

    # H-maximizers
    max_h = max(c['H'] for c in class_list)
    print(f"\n--- H-maximizers (H={max_h}) ---")
    for c in class_list:
        if c['H'] == max_h:
            print(f"  #{c['ci']}: H={c['H']}, |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, "
                  f"SC={c['is_sc']}, SF={c['is_sf']}, scores={c['scores']}, "
                  f"size={c['size']}, GS={c['gs_count']}")

    # Top-10 by H
    print(f"\n--- Top 10 classes by H ---")
    top10 = sorted(class_list, key=lambda c: -c['H'])[:10]
    for c in top10:
        print(f"  #{c['ci']}: H={c['H']}, |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, "
              f"SC={c['is_sc']}, SF={c['is_sf']}, scores={c['scores']}, "
              f"size={c['size']}, GS={c['gs_count']}")

    # Grid-symmetric distribution
    total_gs = sum(c['gs_count'] for c in class_list)
    gs_classes = [c for c in class_list if c['gs_count'] > 0]
    print(f"\n--- Grid-symmetric distribution ---")
    print(f"Total GS tilings: {total_gs}")
    print(f"Classes with GS tilings: {len(gs_classes)}/{len(class_list)}")

    # Top GS classes
    print(f"\nTop 10 classes by GS count:")
    top_gs = sorted(gs_classes, key=lambda c: -c['gs_count'])[:10]
    for c in top_gs:
        print(f"  #{c['ci']}: GS={c['gs_count']}/{c['size']} ({c['gs_count']/c['size']:.3f}), "
              f"H={c['H']}, |Aut|={c['aut']}, SC={c['is_sc']}, SF={c['is_sf']}")

    # Classes with highest |Aut|
    print(f"\n--- Highest symmetry classes ---")
    top_aut = sorted(class_list, key=lambda c: -(c['aut'] + c['anti_aut']))[:10]
    for c in top_aut:
        print(f"  #{c['ci']}: |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, H={c['H']}, "
              f"SC={c['is_sc']}, SF={c['is_sf']}, scores={c['scores']}, size={c['size']}")

    # Score sequence distribution
    score_classes = defaultdict(list)
    for c in class_list:
        score_classes[c['scores']].append(c)

    print(f"\n--- Score sequence -> classes ---")
    for sc in sorted(score_classes.keys()):
        cls = score_classes[sc]
        h_vals = sorted(set(c['H'] for c in cls))
        sc_count = sum(1 for c in cls if c['is_sc'])
        print(f"  {sc}: {len(cls)} classes, H in {h_vals}, SC={sc_count}/{len(cls)}")


if __name__ == '__main__':
    class_list, m = build_n7_data()
    analyze_n7(class_list, m)
