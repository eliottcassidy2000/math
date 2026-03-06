#!/usr/bin/env python3
"""
Fast n=7 tiling isomorphism analysis using score-based hashing
and incremental canonicalization.

Key optimization: group by score sequence first, then use a hash
based on sorted neighborhood lists to partition further, then only
do full canon comparison within small groups.

Instance: opus-2026-03-06-S7
"""

from itertools import permutations
from collections import defaultdict
import time
import sys

def build_data_fast(n):
    verts = list(range(n, 0, -1))
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}

    trans_map = []
    for (x, y) in tiles:
        trans_map.append(tile_idx[(n - y + 1, n - x + 1)])

    def is_grid_sym(mask):
        for i in range(m):
            if ((mask >> i) & 1) != ((mask >> trans_map[i]) & 1):
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

    def adj_invariant(A):
        """Compute isomorphism-invariant hash of adjacency matrix.
        Uses sorted out-degree sequence, sorted in-degree sequence,
        and sorted (out_deg, in_deg, mutual_score) tuple sequence."""
        out_degs = tuple(sorted(sum(A[i]) for i in range(n)))
        # Neighborhood signature: for each vertex, sorted tuple of
        # (out-degree of successors, out-degree of predecessors)
        sigs = []
        for v in range(n):
            succ_degs = tuple(sorted(sum(A[u]) for u in range(n) if A[v][u]))
            pred_degs = tuple(sorted(sum(A[u]) for u in range(n) if A[u][v]))
            sigs.append((sum(A[v]), succ_degs, pred_degs))
        return (out_degs, tuple(sorted(sigs)))

    def are_isomorphic(A, B):
        """Check if two adjacency matrices are isomorphic."""
        # Quick checks
        if tuple(sorted(sum(A[i]) for i in range(n))) != tuple(sorted(sum(B[i]) for i in range(n))):
            return False
        for p in permutations(range(n)):
            if all(A[p[i]][p[j]] == B[i][j] for i in range(n) for j in range(n)):
                return True
        return False

    def count_ham_paths_dp(A):
        """Count Hamiltonian paths using DP (bitmask)."""
        # dp[mask][v] = number of paths visiting exactly the vertices in mask, ending at v
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if dp[mask][v] == 0:
                    continue
                if not (mask & (1 << v)):
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << n) - 1
        return sum(dp[full][v] for v in range(n))

    print(f"n={n}, m={m}, 2^m={1<<m}")
    t0 = time.time()

    # Phase 1: compute invariant hash for each mask
    print("Phase 1: Computing invariant hashes...")
    hash_groups = defaultdict(list)  # invariant -> list of (mask, A)
    gs_total = 0

    for mask in range(1 << m):
        A = mask_to_adj(mask)
        inv = adj_invariant(A)
        hash_groups[inv].append((mask, A))
        if is_grid_sym(mask):
            gs_total += 1

    print(f"  {len(hash_groups)} invariant groups, {gs_total} GS tilings")
    print(f"  Largest group: {max(len(v) for v in hash_groups.values())}")
    print(f"  Time: {time.time()-t0:.1f}s")

    # Phase 2: within each invariant group, find iso classes
    print("Phase 2: Isomorphism testing within groups...")
    t1 = time.time()

    classes = []  # list of (masks_list, rep_A)
    processed = 0
    total = 1 << m

    for inv, group in hash_groups.items():
        # Within this group, cluster by isomorphism
        local_classes = []  # list of (masks_list, rep_A)
        for mask, A in group:
            found = False
            for lc in local_classes:
                if are_isomorphic(A, lc[1]):
                    lc[0].append(mask)
                    found = True
                    break
            if not found:
                local_classes.append(([mask], A))
        classes.extend(local_classes)
        processed += len(group)
        if processed % 5000 < len(group):
            print(f"  {processed}/{total} ({100*processed/total:.0f}%)")

    print(f"  {len(classes)} isomorphism classes")
    print(f"  Time: {time.time()-t1:.1f}s")

    # Phase 3: compute class properties
    print("Phase 3: Computing class properties...")
    t2 = time.time()

    class_list = []
    mask_to_class = {}

    for ci, (masks, rep_A) in enumerate(sorted(classes, key=lambda x: x[1][0])):
        for mask in masks:
            mask_to_class[mask] = ci

        H = count_ham_paths_dp(rep_A)
        sc = tuple(sorted([sum(rep_A[i]) for i in range(n)], reverse=True))

        aut = sum(1 for p in permutations(range(n))
                  if all(rep_A[p[i]][p[j]] == rep_A[i][j] for i in range(n) for j in range(n)))

        anti_aut = sum(1 for p in permutations(range(n))
                       if all(rep_A[p[i]][p[j]] == rep_A[j][i] for i in range(n) for j in range(n)))

        gs_masks = [mask for mask in masks if is_grid_sym(mask)]

        A_op = [[rep_A[j][i] for j in range(n)] for i in range(n)]
        inv_op = adj_invariant(A_op)
        is_sc = are_isomorphic(rep_A, A_op)

        flip_set = set(mask ^ ((1 << m) - 1) for mask in masks)
        is_sf = bool(set(masks) & flip_set)

        class_list.append({
            'ci': ci, 'size': len(masks), 'H': H, 'scores': sc,
            'aut': aut, 'anti_aut': anti_aut,
            'gs_count': len(gs_masks), 'is_sc': is_sc, 'is_sf': is_sf,
            'masks': masks, 'rep_A': rep_A,
        })

    print(f"  Time: {time.time()-t2:.1f}s")
    print(f"  Total time: {time.time()-t0:.1f}s")
    return class_list, m, mask_to_class


def analyze(class_list, m):
    n_classes = len(class_list)
    max_h = max(c['H'] for c in class_list)

    sc_count = sum(1 for c in class_list if c['is_sc'])
    sf_count = sum(1 for c in class_list if c['is_sf'])
    both = sum(1 for c in class_list if c['is_sc'] and c['is_sf'])
    high_aut = sum(1 for c in class_list if c['aut'] > 1)
    has_anti = sum(1 for c in class_list if c['anti_aut'] > 0)

    print(f"\n{'='*70}")
    print(f"ANALYSIS: {n_classes} classes, max H={max_h}")
    print(f"{'='*70}")
    print(f"Self-converse: {sc_count}")
    print(f"Self-flip: {sf_count}")
    print(f"Both SC+SF: {both}")
    print(f"|Aut| > 1: {high_aut}")
    print(f"Has anti-aut: {has_anti}")

    # SC+SF kernel
    print(f"\n--- SC+SF KERNEL CLASSES ---")
    for c in class_list:
        if c['is_sc'] and c['is_sf']:
            is_max = c['H'] == max_h
            print(f"  #{c['ci']}: H={c['H']}{'*MAX*' if is_max else ''}, |Aut|={c['aut']}, "
                  f"|AntiAut|={c['anti_aut']}, GS={c['gs_count']}/{c['size']} "
                  f"({c['gs_count']/c['size']:.3f}), scores={c['scores']}")

    # SF-only
    print(f"\n--- SF-only classes ---")
    for c in class_list:
        if c['is_sf'] and not c['is_sc']:
            print(f"  #{c['ci']}: H={c['H']}, |Aut|={c['aut']}, GS={c['gs_count']}/{c['size']}, "
                  f"scores={c['scores']}")

    # H-maximizers
    print(f"\n--- H-maximizers (H={max_h}) ---")
    for c in class_list:
        if c['H'] == max_h:
            print(f"  #{c['ci']}: |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, SC={c['is_sc']}, "
                  f"SF={c['is_sf']}, scores={c['scores']}, size={c['size']}, GS={c['gs_count']}")

    # Top 15 by H
    print(f"\n--- Top 15 by H ---")
    for c in sorted(class_list, key=lambda c: -c['H'])[:15]:
        print(f"  #{c['ci']}: H={c['H']}, |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, "
              f"SC={c['is_sc']}, SF={c['is_sf']}, scores={c['scores']}, GS={c['gs_count']}/{c['size']}")

    # Highest symmetry
    print(f"\n--- Top 10 by symmetry (|Aut|+|AntiAut|) ---")
    for c in sorted(class_list, key=lambda c: -(c['aut']+c['anti_aut']))[:10]:
        print(f"  #{c['ci']}: |Aut|={c['aut']}, |AntiAut|={c['anti_aut']}, H={c['H']}, "
              f"SC={c['is_sc']}, SF={c['is_sf']}, scores={c['scores']}, GS={c['gs_count']}/{c['size']}")

    # GS distribution
    total_gs = sum(c['gs_count'] for c in class_list)
    gs_classes = [c for c in class_list if c['gs_count'] > 0]
    print(f"\n--- GS distribution ---")
    print(f"Total GS: {total_gs}")
    print(f"Classes with GS: {len(gs_classes)}/{n_classes}")
    for c in sorted(gs_classes, key=lambda c: -c['gs_count'])[:15]:
        print(f"  #{c['ci']}: GS={c['gs_count']}/{c['size']} ({c['gs_count']/c['size']:.3f}), "
              f"H={c['H']}, SC={c['is_sc']}, SF={c['is_sf']}, |Aut|={c['aut']}")

    # Score sequence analysis
    score_groups = defaultdict(list)
    for c in class_list:
        score_groups[c['scores']].append(c)
    print(f"\n--- Score sequences ---")
    for sc in sorted(score_groups.keys()):
        cls = score_groups[sc]
        h_vals = sorted(set(c['H'] for c in cls))
        sc_ct = sum(1 for c in cls if c['is_sc'])
        sf_ct = sum(1 for c in cls if c['is_sf'])
        print(f"  {sc}: {len(cls)} classes, H in {h_vals}, SC={sc_ct}, SF={sf_ct}")


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    class_list, m, mask_to_class = build_data_fast(n)
    analyze(class_list, m)
