#!/usr/bin/env python3
"""
Flip Graph Analysis: The true structure of tiling-level flips and classes.
Instance: opus-2026-03-06-S11

Key insight: Flip (bit complement) is NOT well-defined on isomorphism classes.
Different tilings in the same class can flip to different classes. This means
the "skeleton" must be analyzed at the tiling level, with class membership
as a coloring.

We build:
1. The flip graph on tilings (each tiling connected to its complement)
2. The "flip scatter matrix" F[i][j] = #{tilings in class i whose flip is in class j}
3. Grid-symmetric (GS) tiling analysis
4. Self-flip tilings within each class
"""
from itertools import permutations
from collections import defaultdict, Counter
import sys

def build_data(n):
    verts = list(range(n, 0, -1))
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}
    flip_mask = (1 << m) - 1

    # Transpose map for grid symmetry
    trans_map = []
    for (x, y) in tiles:
        nx, ny = n - y + 1, n - x + 1
        trans_map.append(tile_idx.get((nx, ny), -1))

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

    # Build classes
    classes = {}
    class_members = defaultdict(list)
    class_scores = {}
    class_H = {}
    tiling_class = {}
    tiling_gs = {}
    cid_counter = 0

    for bits in range(1 << m):
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        gs = is_grid_sym(bits)
        tiling_gs[bits] = gs

        if canon not in classes:
            classes[canon] = cid_counter
            class_scores[cid_counter] = scores(A)
            class_H[cid_counter] = hamiltonian_paths(A)
            cid_counter += 1

        cid = classes[canon]
        class_members[cid].append(bits)
        tiling_class[bits] = cid

    num_classes = cid_counter

    return {
        'n': n, 'm': m, 'num_classes': num_classes,
        'class_members': dict(class_members),
        'class_H': class_H, 'class_scores': class_scores,
        'tiling_class': tiling_class, 'tiling_gs': tiling_gs,
        'flip_mask': flip_mask,
    }


def analyze_flip_graph(data):
    n = data['n']
    m = data['m']
    nc = data['num_classes']
    fm = data['flip_mask']
    tc = data['tiling_class']
    gs = data['tiling_gs']
    cm = data['class_members']

    print(f"\n{'='*70}")
    print(f"n = {n}, m = {m}, classes = {nc}, tilings = {1 << m}")
    print(f"{'='*70}")

    # 1. Flip scatter matrix F[i][j]
    F = [[0]*nc for _ in range(nc)]
    for bits in range(1 << m):
        i = tc[bits]
        j = tc[bits ^ fm]
        F[i][j] += 1

    # Self-flip count per class: tilings in class i that flip to class i
    self_flip_count = [F[i][i] for i in range(nc)]

    # 2. GS tiling analysis
    gs_count = [0]*nc
    gs_flip_same_class = [0]*nc  # GS tilings whose flip is in same class
    gs_flip_gs = 0  # GS tilings whose flip is also GS
    for bits in range(1 << m):
        if gs[bits]:
            cid = tc[bits]
            gs_count[cid] += 1
            fbits = bits ^ fm
            if tc[fbits] == cid:
                gs_flip_same_class[cid] += 1
            if gs[fbits]:
                gs_flip_gs += 1

    total_gs = sum(gs_count)
    print(f"\nGrid-symmetric tilings: {total_gs}")
    print(f"GS tilings whose flip is also GS: {gs_flip_gs}")

    # 3. Class-level flip structure
    print(f"\n--- Class details ---")
    print(f"{'Cls':>4} {'Size':>5} {'H':>4} {'GS':>3} {'SelfF':>6} {'#Tgts':>6} {'Targets':>20} {'Scores'}")
    for cid in range(nc):
        sz = len(cm[cid])
        H = data['class_H'][cid]
        gs_c = gs_count[cid]
        sf = self_flip_count[cid]
        targets = {j: F[cid][j] for j in range(nc) if F[cid][j] > 0}
        n_targets = len(targets)
        scores = data['class_scores'][cid]
        tgt_str = str(dict(sorted(targets.items()))) if n_targets <= 5 else f"{n_targets} classes"
        print(f"  {cid:>3} {sz:>5} {H:>4} {gs_c:>3} {sf:>6} {n_targets:>6}  {tgt_str:<30} {scores}")

    # 4. "Pure flip" classes: classes where ALL members flip to the same target
    pure_flip = []
    for cid in range(nc):
        targets = [j for j in range(nc) if F[cid][j] > 0]
        if len(targets) == 1:
            pure_flip.append((cid, targets[0]))

    print(f"\n--- Pure flip classes (all members flip to single class) ---")
    if pure_flip:
        for (cid, tgt) in pure_flip:
            label = "SELF" if tgt == cid else f"-> {tgt}"
            print(f"  Class {cid} (size {len(cm[cid])}) {label}")
    else:
        print("  None!")

    # 5. Diagonal dominance: classes with highest self-flip ratio
    print(f"\n--- Self-flip ratio (fraction of members whose flip stays in same class) ---")
    ratios = []
    for cid in range(nc):
        sz = len(cm[cid])
        ratio = self_flip_count[cid] / sz
        ratios.append((ratio, cid))
    ratios.sort(reverse=True)
    for ratio, cid in ratios[:10]:
        sz = len(cm[cid])
        H = data['class_H'][cid]
        gs_c = gs_count[cid]
        sc = data['class_scores'][cid]
        print(f"  Class {cid}: {ratio:.3f} ({self_flip_count[cid]}/{sz}), H={H}, GS={gs_c}, scores={sc}")

    # 6. Transitive and Full classes
    trans_class = tc[0]
    full_class = tc[fm]
    print(f"\nTransitive class: {trans_class} (size={len(cm[trans_class])}, H={data['class_H'][trans_class]})")
    print(f"Full class: {full_class} (size={len(cm[full_class])}, H={data['class_H'][full_class]})")
    print(f"  Transitive tiling flips to class: {tc[fm]}")
    print(f"  Full tiling flips to class: {tc[0]}")

    # 7. Flip scatter matrix symmetry check
    # F[i][j] vs F[j][i]: are they equal?
    sym_violations = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if F[i][j] != F[j][i]:
                sym_violations += 1
    print(f"\nFlip matrix symmetry violations: {sym_violations}")

    # 8. GS tiling distribution
    print(f"\n--- GS tilings per class ---")
    gs_classes = [(cid, gs_count[cid]) for cid in range(nc) if gs_count[cid] > 0]
    for cid, gc in sorted(gs_classes, key=lambda x: -x[1]):
        sz = len(cm[cid])
        print(f"  Class {cid}: {gc}/{sz} GS, H={data['class_H'][cid]}, scores={data['class_scores'][cid]}")

    return {
        'n': n, 'F': F, 'self_flip_count': self_flip_count,
        'gs_count': gs_count, 'pure_flip': pure_flip,
        'trans_class': trans_class, 'full_class': full_class,
    }


def cross_scale(results):
    print(f"\n{'='*70}")
    print("CROSS-SCALE SUMMARY")
    print(f"{'='*70}")
    print(f"{'n':>3} {'cls':>5} {'pure':>5} {'totGS':>6} {'GSflipGS':>9} {'transH':>7} {'fullH':>6} {'fullSz':>7}")
    for r in results:
        n = r['n']
        data = r['data']
        nc = data['num_classes']
        pure = len(r['pure_flip'])
        total_gs = sum(r['gs_count'])
        # Recompute GS-flip-GS
        fm = data['flip_mask']
        gs = data['tiling_gs']
        tc = data['tiling_class']
        gfg = sum(1 for b in range(1 << data['m']) if gs[b] and gs[b ^ fm])
        tH = data['class_H'][r['trans_class']]
        fH = data['class_H'][r['full_class']]
        fSz = len(data['class_members'][r['full_class']])
        print(f"{n:>3} {nc:>5} {pure:>5} {total_gs:>6} {gfg:>9} {tH:>7} {fH:>6} {fSz:>7}")

    # Check the 1+2^(n-2) pattern for full class size
    print(f"\nFull class size vs 1+2^(n-2):")
    for r in results:
        n = r['n']
        fSz = len(r['data']['class_members'][r['full_class']])
        formula = 1 + 2**(n-2)
        print(f"  n={n}: actual={fSz}, 1+2^(n-2)={formula}, match={fSz == formula}")


def main():
    results = []
    for n in range(3, 7):  # n=7 too slow for full canonicalization
        print(f"\nBuilding n={n}...")
        data = build_data(n)
        r = analyze_flip_graph(data)
        r['data'] = data
        results.append(r)

    cross_scale(results)


if __name__ == '__main__':
    main()
