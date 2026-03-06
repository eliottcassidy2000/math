#!/usr/bin/env python3
"""
Skeletal Structure of Tiling Classes and Blue Line Pairs.
Instance: opus-2026-03-06-S11

Investigates:
1. The "skeleton" of unpaired isomorphism classes formed by blue/black flip pairs
2. Cross-scale self-similarity patterns between different n values
3. How the perpendicular structure (SC tournaments on the perpendicular hyperplane)
   constrains the class hierarchy
4. Connection between even cycle vanishing and the paired/unpaired structure
5. The Tribonacci constraint on T_full and its implications

Key concepts:
- BLUE LINE: flip pair where BOTH tilings are grid-symmetric (SC tournaments)
- BLACK LINE: flip pair where tilings are NOT grid-symmetric
- UNPAIRED CLASS: isomorphism class that is NOT self-flip (its flip goes to a different class)
- PAIRED CLASS: self-flip class (flip stays within the class)
- The "skeleton" = the structure of how classes are connected by flip lines
"""

from itertools import permutations
from collections import defaultdict, Counter
import sys

def build_tournament_data(n):
    """Build all tiling/tournament data for given n."""
    verts = list(range(n, 0, -1))  # [n, n-1, ..., 1]

    # Tiles: triangular grid positions
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))

    m = len(tiles)  # C(n-1, 2)

    # Tile index lookup
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}

    # Transpose map: (x,y) -> (n-y+1, n-x+1)
    trans_map = []
    for (x, y) in tiles:
        nx, ny = n - y + 1, n - x + 1
        if (nx, ny) in tile_idx:
            trans_map.append(tile_idx[(nx, ny)])
        else:
            trans_map.append(-1)

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

    # Build all tilings
    classes = {}  # canonical -> class_id
    class_members = defaultdict(list)
    class_canon = {}
    class_scores = {}
    class_H = {}
    class_gs_count = {}  # grid-symmetric members per class
    tiling_class = {}
    tiling_gs = {}

    class_id_counter = 0

    for bits in range(1 << m):
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        gs = is_grid_sym(bits)
        tiling_gs[bits] = gs

        if canon not in classes:
            classes[canon] = class_id_counter
            class_canon[class_id_counter] = canon
            class_scores[class_id_counter] = scores(A)
            class_H[class_id_counter] = hamiltonian_paths(A)
            class_gs_count[class_id_counter] = 0
            class_id_counter += 1

        cid = classes[canon]
        class_members[cid].append(bits)
        tiling_class[bits] = cid
        if gs:
            class_gs_count[cid] += 1

    num_classes = class_id_counter

    # Build flip map between classes
    flip_mask = (1 << m) - 1
    class_flip = {}  # class_id -> flip_class_id
    for cid in range(num_classes):
        rep_bits = class_members[cid][0]
        flip_bits = rep_bits ^ flip_mask
        flip_cid = tiling_class[flip_bits]
        class_flip[cid] = flip_cid

    # Categorize classes
    self_flip_classes = []
    flip_pairs = []
    seen = set()
    for cid in range(num_classes):
        if cid in seen:
            continue
        fcid = class_flip[cid]
        if fcid == cid:
            self_flip_classes.append(cid)
            seen.add(cid)
        else:
            if fcid not in seen:
                flip_pairs.append((cid, fcid))
                seen.add(cid)
                seen.add(fcid)

    # Determine blue vs black for self-flip classes
    blue_self = []
    black_self = []
    for cid in self_flip_classes:
        # Check if self-flip members are grid-symmetric
        sf_members = []
        for bits in class_members[cid]:
            fbits = bits ^ flip_mask
            if tiling_class[fbits] == cid:
                sf_members.append(bits)
        gs_in_sf = sum(1 for b in sf_members if tiling_gs[b])
        if gs_in_sf == len(sf_members) and gs_in_sf > 0:
            blue_self.append(cid)
        elif gs_in_sf == 0:
            black_self.append(cid)
        else:
            # Mixed - shouldn't happen by THM-022
            black_self.append(cid)

    # Check if class is self-converse
    def is_class_sc(cid):
        A = class_canon[cid]
        nn = n
        A_mat = [list(A[i*nn:(i+1)*nn]) for i in range(nn)]
        A_op = [[A_mat[j][i] for j in range(nn)] for i in range(nn)]
        A_op_tuple = tuple(tuple(row) for row in A_op)
        canon_op = canonicalize(A_op_tuple)
        op_cid = classes.get(canon_op, -1)
        return op_cid == cid

    class_sc = {cid: is_class_sc(cid) for cid in range(num_classes)}

    return {
        'n': n, 'm': m, 'num_classes': num_classes,
        'class_members': dict(class_members),
        'class_H': class_H,
        'class_scores': class_scores,
        'class_gs_count': class_gs_count,
        'class_flip': class_flip,
        'class_sc': class_sc,
        'self_flip_classes': self_flip_classes,
        'blue_self': blue_self,
        'black_self': black_self,
        'flip_pairs': flip_pairs,
    }


def analyze_skeleton(data):
    """Analyze the skeletal structure of unpaired classes."""
    n = data['n']
    num = data['num_classes']
    sf = data['self_flip_classes']
    fp = data['flip_pairs']
    blue = data['blue_self']
    black = data['black_self']

    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")
    print(f"Total tiling classes: {num}")
    print(f"Flip pairs (external lines): {len(fp)}")
    print(f"Self-flip classes: {len(sf)}")
    print(f"  Blue self: {len(blue)}")
    print(f"  Black self: {len(black)}")
    print(f"Unpaired classes (in external flip pairs): {2 * len(fp)}")

    # The "skeleton" = unpaired classes + self-flip classes
    print(f"\nSkeleton: {len(sf)} self-flip nodes + {len(fp)} pair-edges connecting {2*len(fp)} nodes")

    # H-values of self-flip vs non-self-flip classes
    sf_H = [data['class_H'][c] for c in sf]
    non_sf_H = []
    for (a, b) in fp:
        non_sf_H.extend([data['class_H'][a], data['class_H'][b]])

    if sf_H:
        print(f"\nSelf-flip H values: {sorted(sf_H)}")
    if non_sf_H:
        print(f"Non-self-flip H range: [{min(non_sf_H)}, {max(non_sf_H)}]")
        print(f"Non-self-flip H mean: {sum(non_sf_H)/len(non_sf_H):.1f}")

    # SC structure
    sc_count = sum(1 for c in range(num) if data['class_sc'][c])
    sf_sc = sum(1 for c in sf if data['class_sc'][c])
    print(f"\nSelf-converse classes: {sc_count}/{num}")
    print(f"Self-flip that are SC: {sf_sc}/{len(sf)}")

    # Grid-symmetric tiling distribution
    gs_total = sum(data['class_gs_count'][c] for c in range(num))
    gs_in_sf = sum(data['class_gs_count'][c] for c in sf)
    print(f"\nGrid-symmetric tilings: {gs_total} total")
    print(f"GS tilings in self-flip classes: {gs_in_sf}")

    # Full tiling class
    full_bits = (1 << data['m']) - 1
    # The "full tiling" is all bits set
    # The "transitive tiling" is all bits 0
    # They should be in a flip pair
    for c in range(num):
        members = data['class_members'][c]
        if 0 in members:
            trans_class = c
        if full_bits in members:
            full_class = c

    print(f"\nTransitive class: #{trans_class} (H={data['class_H'][trans_class]}, size={len(data['class_members'][trans_class])})")
    print(f"Full tiling class: #{full_class} (H={data['class_H'][full_class]}, size={len(data['class_members'][full_class])})")
    print(f"Transitive flip -> Full: {data['class_flip'][trans_class] == full_class}")

    # Detailed self-flip class info
    print(f"\n--- Self-flip classes detail ---")
    for c in sorted(sf, key=lambda x: data['class_H'][x]):
        H = data['class_H'][c]
        sc = data['class_sc'][c]
        gs = data['class_gs_count'][c]
        sz = len(data['class_members'][c])
        sc_str = "SC" if sc else "NSC"
        color = "BLUE" if c in blue else "BLACK"
        scores = data['class_scores'][c]
        print(f"  #{c:>3}: H={H:>4}, size={sz:>4}, GS={gs:>3}, {sc_str:>3}, {color:>5}, scores={scores}")

    # The perpendicular structure
    print(f"\n--- Perpendicular structure ---")
    # Classes on the perpendicular hyperplane = those with equal numbers of 0 and 1 bits
    mid_classes = []
    for c in range(num):
        bits_example = data['class_members'][c][0]
        hw = bin(bits_example).count('1')
        if hw == data['m'] // 2 or hw == (data['m'] + 1) // 2:
            mid_classes.append(c)

    print(f"Classes near Hamming midpoint: {len(mid_classes)}/{num}")
    sf_at_mid = [c for c in sf if c in mid_classes]
    print(f"Self-flip classes near midpoint: {len(sf_at_mid)}/{len(sf)}")

    return {
        'n': n, 'num_classes': num,
        'self_flip': sf, 'blue': blue, 'black': black,
        'flip_pairs': fp, 'trans_class': trans_class, 'full_class': full_class,
    }


def cross_scale_analysis(results):
    """Analyze patterns across different n values."""
    print(f"\n{'='*70}")
    print("CROSS-SCALE ANALYSIS")
    print(f"{'='*70}")

    print(f"\n{'n':>3} {'classes':>8} {'self-flip':>10} {'blue':>6} {'black':>7} {'pairs':>6} {'ratio_sf':>10}")
    for r in results:
        n = r['n']
        nc = r['num_classes']
        nsf = len(r['self_flip'])
        nb = len(r['blue'])
        nbl = len(r['black'])
        np = len(r['flip_pairs'])
        ratio = nsf / nc if nc > 0 else 0
        print(f"{n:>3} {nc:>8} {nsf:>10} {nb:>6} {nbl:>7} {np:>6} {ratio:>10.4f}")

    print(f"\nFull tiling H sequence (Tribonacci):")
    for r in results:
        n = r['n']
        # Get H of full class from the data
        # (we stored this in the analysis)
        print(f"  n={n}: full_class=#{r.get('full_class', '?')}")


def main():
    results = []
    for n in range(3, 8):
        print(f"\nBuilding data for n={n}...")
        if n <= 6:
            data = build_tournament_data(n)
            r = analyze_skeleton(data)
            results.append(r)
        else:
            print(f"  n={n} skipped (too slow for full canonicalization)")

    if results:
        cross_scale_analysis(results)

if __name__ == '__main__':
    main()
