#!/usr/bin/env python3
"""
Focused analysis of the "symmetry kernel" — the chain of SC+SF classes
across n values, and the GS tiling ratio within these classes.

Key questions:
1. Does the SC+SF chain always have exactly 2 classes at each n >= 5?
2. What is the GS/size ratio pattern?
3. Do the H values in the SC+SF chain relate to those of the H-maximizer?
4. What is the automorphism structure of these "kernel" classes?

Instance: opus-2026-03-06-S7
"""

from itertools import permutations
from collections import defaultdict
import sys

def build_data(n):
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

    def canon_key(A):
        best = None
        for p in permutations(range(n)):
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def count_ham_paths(A):
        count = 0
        for p in permutations(range(n)):
            valid = all(A[p[k]][p[k+1]] == 1 for k in range(n-1))
            if valid:
                count += 1
        return count

    classes = defaultdict(list)
    for mask in range(1 << m):
        A = mask_to_adj(mask)
        ck = canon_key(A)
        classes[ck].append(mask)

    class_list = []
    for ci, (ck, masks) in enumerate(sorted(classes.items())):
        rep_A = mask_to_adj(masks[0])
        H = count_ham_paths(rep_A)
        sc = tuple(sorted([sum(rep_A[i]) for i in range(n)], reverse=True))

        aut = sum(1 for p in permutations(range(n))
                  if all(rep_A[p[i]][p[j]] == rep_A[i][j] for i in range(n) for j in range(n)))

        anti_aut = sum(1 for p in permutations(range(n))
                       if all(rep_A[p[i]][p[j]] == rep_A[j][i] for i in range(n) for j in range(n)))

        gs_masks = [mask for mask in masks if is_grid_sym(mask)]

        A_op = [[rep_A[j][i] for j in range(n)] for i in range(n)]
        is_sc = canon_key(A_op) == ck

        flip_set = set(mask ^ ((1 << m) - 1) for mask in masks)
        is_sf = bool(set(masks) & flip_set)

        class_list.append({
            'ci': ci, 'size': len(masks), 'H': H, 'scores': sc,
            'aut': aut, 'anti_aut': anti_aut,
            'gs_count': len(gs_masks), 'is_sc': is_sc, 'is_sf': is_sf,
            'masks': masks, 'ck': ck,
        })

    return class_list, m, tiles, trans_map, verts, tile_idx


def find_parent(child_A, parent_classes, parent_n):
    """Find parent class at n-1 by deleting vertex 0 from child."""
    n = len(child_A)
    sub_A = [[child_A[i][j] for j in range(1, n)] for i in range(1, n)]

    best = None
    for p in permutations(range(parent_n)):
        s = tuple(sub_A[p[i]][p[j]] for i in range(parent_n) for j in range(parent_n))
        if best is None or s < best:
            best = s

    for c in parent_classes:
        if c['ck'] == best:
            return c['ci']
    return -1


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    all_data = {}
    for n in range(3, max_n + 1):
        print(f"Building n={n}...")
        cl, m, tiles, trans_map, verts, tile_idx = build_data(n)
        all_data[n] = cl

    print(f"\n{'='*70}")
    print(f"SYMMETRY KERNEL CHAIN ANALYSIS")
    print(f"{'='*70}")

    for n in sorted(all_data.keys()):
        cl = all_data[n]
        scsf = [c for c in cl if c['is_sc'] and c['is_sf']]
        sc_only = [c for c in cl if c['is_sc'] and not c['is_sf']]
        sf_only = [c for c in cl if c['is_sf'] and not c['is_sc']]
        max_h = max(c['H'] for c in cl)

        print(f"\nn={n}: {len(cl)} classes, max H={max_h}")
        print(f"  SC+SF: {len(scsf)}, SC only: {len(sc_only)}, SF only: {len(sf_only)}")

        for c in scsf:
            is_max = c['H'] == max_h
            print(f"  KERNEL #{c['ci']}: H={c['H']}{'*' if is_max else ''}, |Aut|={c['aut']}, "
                  f"|AntiAut|={c['anti_aut']}, GS={c['gs_count']}/{c['size']} "
                  f"({c['gs_count']/c['size']:.3f}), scores={c['scores']}")

        for c in sf_only:
            print(f"  SF-only #{c['ci']}: H={c['H']}, |Aut|={c['aut']}, "
                  f"GS={c['gs_count']}/{c['size']}, scores={c['scores']}")

    # Parent-child tracking of SC+SF
    print(f"\n{'='*70}")
    print(f"SC+SF PARENT-CHILD CHAIN")
    print(f"{'='*70}")

    for n in range(4, max_n + 1):
        if n not in all_data or n-1 not in all_data:
            continue

        cl_n = all_data[n]
        cl_prev = all_data[n-1]
        scsf_prev = {c['ci'] for c in cl_prev if c['is_sc'] and c['is_sf']}

        print(f"\nn={n-1} -> n={n}:")
        for c in cl_n:
            if c['is_sc'] or c['is_sf']:
                # Find parent
                from itertools import permutations as perms
                def mask_to_adj_n(mask, nn):
                    vv = list(range(nn, 0, -1))
                    tt = []
                    for y in range(1, nn-1):
                        for x in range(nn, y+1, -1):
                            tt.append((x, y))
                    A = [[0]*nn for _ in range(nn)]
                    for k in range(nn-1):
                        A[k][k+1] = 1
                    for i, (xL, yL) in enumerate(tt):
                        xi = vv.index(xL)
                        yi = vv.index(yL)
                        if (mask >> i) & 1 == 0:
                            A[xi][yi] = 1
                        else:
                            A[yi][xi] = 1
                    return A

                rep_A = mask_to_adj_n(c['masks'][0], n)
                parent_ci = find_parent(rep_A, cl_prev, n-1)
                parent_c = cl_prev[parent_ci] if parent_ci >= 0 else None

                p_sc = parent_c['is_sc'] if parent_c else '?'
                p_sf = parent_c['is_sf'] if parent_c else '?'
                from_kernel = parent_ci in scsf_prev

                tag = ""
                if c['is_sc'] and c['is_sf']:
                    tag = " <== KERNEL"
                elif c['is_sc']:
                    tag = " (SC only)"
                elif c['is_sf']:
                    tag = " (SF only)"

                print(f"  #{c['ci']}: H={c['H']}, SC={c['is_sc']}, SF={c['is_sf']}{tag}")
                print(f"    Parent #{parent_ci}: H={parent_c['H'] if parent_c else '?'}, "
                      f"SC={p_sc}, SF={p_sf}{' <== from kernel' if from_kernel else ''}")

    # GS ratio analysis
    print(f"\n{'='*70}")
    print(f"GS RATIO ANALYSIS")
    print(f"{'='*70}")

    for n in sorted(all_data.keys()):
        cl = all_data[n]
        m = n*(n-1)//2 - (n-1)  # = C(n-1,2)
        # fixed + orbit pairs
        # The transpose map has fixed_tiles tiles on the diagonal
        # For staircase delta_{n-2}, the "diagonal" (self-transpose tiles) count is floor((n-2)/2)

        scsf = [c for c in cl if c['is_sc'] and c['is_sf']]
        max_h = max(c['H'] for c in cl)
        max_classes = [c for c in cl if c['H'] == max_h]

        print(f"\nn={n}: m={m}")
        print(f"  H-max classes: {[c['ci'] for c in max_classes]}, H={max_h}")
        for c in max_classes:
            print(f"    #{c['ci']}: GS/size={c['gs_count']}/{c['size']}={c['gs_count']/c['size']:.4f}, "
                  f"|Aut|={c['aut']}, |AntiAut|={c['anti_aut']}")
        if scsf:
            print(f"  SC+SF kernel classes:")
            for c in scsf:
                print(f"    #{c['ci']}: GS/size={c['gs_count']}/{c['size']}={c['gs_count']/c['size']:.4f}, "
                      f"H={c['H']}")


if __name__ == '__main__':
    main()
