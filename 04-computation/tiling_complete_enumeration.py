"""
tiling_complete_enumeration.py -- kind-pasteur-2026-03-14-S77
Complete enumeration of ALL tiling diagram aspects for n=3..7.

WHAT THE TOURNAMENT-TILING-EXPLORER SHOWS:
1. Triangular tile grid delta_{n-2} with m = C(n-1,2) tiles
2. Each tiling = binary word in {0,1}^m
3. Isomorphism classes: tilings grouped by tournament isomorphism
4. BLUE LINES: GS tiling whose flip lands in a different class
5. BLACK LINES: non-GS tiling whose flip lands in a different class
6. RED LINES: transpose-paired classes
7. Blueself: GS tiling whose flip lands in the SAME class
8. Blackself: non-GS tiling whose flip lands in the SAME class

GOAL: Fast formulas and complete census for each n.
For each isomorphism class, compute:
  - Size (# tilings in class)
  - H value
  - Score sequence
  - # GS tilings in class
  - # non-GS tilings in class
  - # blue lines FROM this class (GS flips to other classes)
  - # blue self-loops (blueself: GS flip stays in class)
  - # black lines FROM this class (non-GS flips to other classes)
  - # black self-loops (blackself: non-GS flip stays in class)
  - Transpose partner class
  - t3 (3-cycle count) for bipartition
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

def build_all(n):
    """Build complete tiling data for tournament size n."""
    # Vertices (1-indexed as in the paper)
    verts = list(range(n, 0, -1))  # [n, n-1, ..., 1]

    # Tiles: (a,b) with a >= b+2, a,b in {1,...,n}
    tiles = []
    for b in range(1, n-1):
        for a in range(b+2, n+1):
            tiles.append((a, b))
    m = len(tiles)

    # Tile index
    tile_idx = {t: i for i, t in enumerate(tiles)}

    # GS (transpose) map: arc (a,b) -> (n+1-b, n+1-a)
    trans_map = []
    for a, b in tiles:
        a2, b2 = n+1-b, n+1-a
        if a2 < b2:
            a2, b2 = b2, a2  # ensure a >= b+2
        trans_map.append(tile_idx.get((a2, b2), -1))

    def is_gs(bits):
        """Check if tiling bits are grid-symmetric."""
        for i in range(m):
            j = trans_map[i]
            if j != i and ((bits >> i) & 1) != ((bits >> j) & 1):
                return False
        return True

    def bits_to_adj(bits):
        """Convert tiling bits to adjacency matrix (0-indexed vertices)."""
        A = [[0]*n for _ in range(n)]
        # Backbone: vertex i beats vertex i-1 (0-indexed: position n-1-v in verts)
        # Base path: n -> n-1 -> ... -> 1 means A[vi][vj] = 1 when vi is before vj
        # In 0-indexed: vertex n-1 beats n-2, n-2 beats n-3, etc.
        for i in range(1, n):
            A[i][i-1] = 1  # i -> i-1

        for idx, (a, b) in enumerate(tiles):
            a0, b0 = a-1, b-1  # 0-indexed
            if (bits >> idx) & 1:
                A[b0][a0] = 1  # backward: b -> a
            else:
                A[a0][b0] = 1  # forward: a -> b
        return A

    def adj_to_str(A):
        return ''.join(str(A[i][j]) for i in range(n) for j in range(n))

    def canonicalize(A):
        """Canonical form under vertex permutation."""
        perms = list(permutations(range(n)))
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def scores(A):
        return tuple(sorted([sum(row) for row in A]))

    def count_c3(A):
        total = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        total += 1
        return total // 3

    # Enumerate all tilings
    all_tilings = []
    for mask in range(1 << m):
        A = bits_to_adj(mask)
        canon = canonicalize(A)
        gs = is_gs(mask)
        flip_mask = ((1 << m) - 1) ^ mask

        # Transpose mask
        trans_bits = 0
        for i in range(m):
            if (mask >> i) & 1:
                trans_bits |= (1 << trans_map[i])

        all_tilings.append({
            'mask': mask,
            'adj': A,
            'canon': canon,
            'gs': gs,
            'flip': flip_mask,
            'trans': trans_bits,
            'scores': scores(A),
            'c3': count_c3(A),
        })

    # Group by canonical form -> isomorphism classes
    groups = defaultdict(list)
    for t in all_tilings:
        groups[t['canon']].append(t)

    class_list = sorted(groups.keys())
    class_idx = {c: i for i, c in enumerate(class_list)}

    # Assign class index to each tiling
    for t in all_tilings:
        t['class'] = class_idx[t['canon']]

    # Build mask -> tiling lookup
    mask_to_tiling = {t['mask']: t for t in all_tilings}

    return {
        'n': n, 'm': m, 'tiles': tiles, 'trans_map': trans_map,
        'tilings': all_tilings, 'groups': groups, 'class_list': class_list,
        'class_idx': class_idx, 'mask_to_tiling': mask_to_tiling,
    }

def analyze_classes(data):
    """Compute all class-level statistics."""
    n = data['n']
    m = data['m']
    groups = data['groups']
    class_list = data['class_list']
    class_idx = data['class_idx']
    mask_to_tiling = data['mask_to_tiling']

    results = []
    for ci, canon in enumerate(class_list):
        members = groups[canon]
        rep = members[0]

        # Basic info
        size = len(members)
        H = sum(1 for perm in permutations(range(n))
                if all(rep['adj'][perm[i]][perm[i+1]] for i in range(n-1)))
        scores = rep['scores']
        c3 = rep['c3']

        # GS breakdown
        gs_members = [t for t in members if t['gs']]
        ngs_members = [t for t in members if not t['gs']]
        n_gs = len(gs_members)
        n_ngs = len(ngs_members)

        # Flip analysis for GS tilings (blue lines)
        blue_self = 0  # blueself: GS flip lands in SAME class
        blue_cross = Counter()  # blue_cross[cj] = # GS flips landing in class cj != ci
        for t in gs_members:
            flip_t = mask_to_tiling[t['flip']]
            if flip_t['class'] == ci:
                blue_self += 1
            else:
                blue_cross[flip_t['class']] += 1

        # Flip analysis for non-GS tilings (black lines)
        black_self = 0  # blackself
        black_cross = Counter()
        for t in ngs_members:
            flip_t = mask_to_tiling[t['flip']]
            if flip_t['class'] == ci:
                black_self += 1
            else:
                black_cross[flip_t['class']] += 1

        # Transpose partner
        trans_class = mask_to_tiling[rep['trans']]['class']

        results.append({
            'ci': ci,
            'size': size,
            'H': H,
            'scores': scores,
            'c3': c3,
            'n_gs': n_gs,
            'n_ngs': n_ngs,
            'blue_self': blue_self,
            'blue_cross': dict(blue_cross),
            'black_self': black_self,
            'black_cross': dict(black_cross),
            'trans_class': trans_class,
            'is_sc': (trans_class == ci),
        })

    return results

def main():
    print("=" * 70)
    print("COMPLETE TILING ENUMERATION — ALL DIAGRAM ASPECTS")
    print("kind-pasteur-2026-03-14-S77")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        data = build_all(n)
        m = data['m']
        N = 1 << m

        print(f"  m = {m} tiles, {N} tilings, {len(data['class_list'])} iso classes")

        # GS count
        n_gs = sum(1 for t in data['tilings'] if t['gs'])
        print(f"  GS tilings: {n_gs} ({100*n_gs/N:.1f}%)")

        classes = analyze_classes(data)

        # CLASS TABLE
        print(f"\n  {'CI':>3} {'Size':>5} {'H':>4} {'c3':>3} {'Scores':>20} {'GS':>3} {'nGS':>4} "
              f"{'BluS':>4} {'BluX':>4} {'BlkS':>4} {'BlkX':>4} {'SC':>3} {'Trans':>5}")
        print(f"  {'-'*85}")

        for c in classes:
            blue_cross_total = sum(c['blue_cross'].values())
            black_cross_total = sum(c['black_cross'].values())
            sc_str = 'SC' if c['is_sc'] else ''
            print(f"  {c['ci']:3d} {c['size']:5d} {c['H']:4d} {c['c3']:3d} {str(c['scores']):>20s} "
                  f"{c['n_gs']:3d} {c['n_ngs']:4d} {c['blue_self']:4d} {blue_cross_total:4d} "
                  f"{c['black_self']:4d} {black_cross_total:4d} {sc_str:>3s} {c['trans_class']:5d}")

        # SUMMARY STATISTICS
        total_blue_self = sum(c['blue_self'] for c in classes)
        total_blue_cross = sum(sum(c['blue_cross'].values()) for c in classes)
        total_black_self = sum(c['black_self'] for c in classes)
        total_black_cross = sum(sum(c['black_cross'].values()) for c in classes)
        total_sc = sum(1 for c in classes if c['is_sc'])
        total_nsc_pairs = sum(1 for c in classes if not c['is_sc']) // 2

        print(f"\n  SUMMARY:")
        print(f"    Total classes: {len(classes)}")
        print(f"    SC classes: {total_sc}")
        print(f"    NSC pairs: {total_nsc_pairs}")
        print(f"    GS tilings: {n_gs}")
        print(f"    Blueself (GS flip = same class): {total_blue_self}")
        print(f"    Blue cross (GS flip = diff class): {total_blue_cross}")
        print(f"    Blackself (non-GS flip = same class): {total_black_self}")
        print(f"    Black cross (non-GS flip = diff class): {total_black_cross}")
        print(f"    Check: blueself + blue_cross = 2 * GS_count? "
              f"{total_blue_self + total_blue_cross} = {2 * n_gs} = {total_blue_self + total_blue_cross == 2 * n_gs}")
        print(f"    Check: blackself + black_cross = 2 * (N - GS_count)? "
              f"{total_black_self + total_black_cross} = {2 * (N - n_gs)} = {total_black_self + total_black_cross == 2 * (N - n_gs)}")

        # BLUE LINE SKELETON
        blue_edges = set()
        for c in classes:
            for cj, count in c['blue_cross'].items():
                if cj != c['ci']:
                    edge = (min(c['ci'], cj), max(c['ci'], cj))
                    blue_edges.add(edge)

        print(f"\n  BLUE LINE SKELETON:")
        print(f"    Vertices (SC classes): {total_sc}")
        print(f"    Edges: {len(blue_edges)}")
        for e in sorted(blue_edges):
            c1, c2 = classes[e[0]], classes[e[1]]
            weight = c1['blue_cross'].get(e[1], 0) + c2['blue_cross'].get(e[0], 0)
            print(f"      {e[0]}(H={c1['H']},c3={c1['c3']}) --- {e[1]}(H={c2['H']},c3={c2['c3']}) [weight={weight}]")

        # Bipartite check (t3 parity)
        if n % 2 == 1:
            bipartite = all(
                classes[e[0]]['c3'] % 2 != classes[e[1]]['c3'] % 2
                for e in blue_edges
            )
            print(f"    Bipartite by t3 parity? {bipartite}")

        # BLACK LINE SKELETON
        black_edges = set()
        for c in classes:
            for cj, count in c['black_cross'].items():
                if cj != c['ci']:
                    edge = (min(c['ci'], cj), max(c['ci'], cj))
                    black_edges.add(edge)

        print(f"\n  BLACK LINE SKELETON:")
        print(f"    Edges: {len(black_edges)}")

        # FORMULAS
        print(f"\n  KEY FORMULAS:")
        print(f"    N = 2^m = 2^{m} = {N}")
        print(f"    #classes = {len(classes)}")
        print(f"    #GS = {n_gs}")

        # GS count formula: 2^(f + p) where f = fixed, p = paired
        gs_fixed = sum(1 for i in range(m) if data['trans_map'][i] == i)
        gs_pairs = (m - gs_fixed) // 2
        gs_dof = gs_fixed + gs_pairs
        print(f"    GS formula: 2^({gs_fixed} + {gs_pairs}) = 2^{gs_dof} = {2**gs_dof}")
        print(f"    Check: {n_gs} = {2**gs_dof}? {n_gs == 2**gs_dof}")

        # Total flip connections
        print(f"    Total flip connections: {N} (each tiling has exactly 1 flip)")
        print(f"    Of which GS: {n_gs} (landing in GS)")
        print(f"    Of which non-GS: {N - n_gs}")

        # Relationship: flip of GS tiling is also GS?
        gs_flip_gs = sum(1 for t in data['tilings'] if t['gs'] and data['mask_to_tiling'][t['flip']]['gs'])
        print(f"    GS flip -> GS: {gs_flip_gs}")
        print(f"    GS flip -> non-GS: {n_gs - gs_flip_gs}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
