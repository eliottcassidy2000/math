#!/usr/bin/env python3
"""
Investigate the odd-n / even-n GS flip dichotomy.
Instance: opus-2026-03-06-S11

At odd n: ALL GS flip pairs cross classes (0% same-class)
At even n: some GS pairs stay within same class

Key question: WHY?

A GS tiling b and its flip b' = b XOR mask are BOTH grid-symmetric.
They stay in the same class iff T(b) is isomorphic to T(b').

T(b) is the tournament with path + non-path edges determined by b.
T(b') reverses ALL non-path edges.

If T(b) is self-converse (T isomorphic to T^op), and T(b') is also
self-converse, when can T(b) be isomorphic to T(b')?

T(b') is obtained from T(b) by reversing all non-path edges.
T(b)^op reverses ALL edges (including path edges).

So T(b') is NOT T(b)^op in general — T(b') reverses only non-path edges,
while T(b)^op reverses all edges.

The difference: T(b') has the SAME path edges as T(b) (both have i->i+1),
but reversed non-path edges. T(b)^op has REVERSED path edges (i+1->i) and
reversed non-path edges.

For T(b) ≅ T(b'), we need the tournament to be "self-flip" — invariant
under reversing only non-path edges.

Let's investigate when this happens and its connection to n parity.
"""
from itertools import permutations
from collections import defaultdict

def investigate_same_class_gs_flips(n):
    """Find GS tilings b where T(b) ≅ T(b XOR mask)."""
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

    # Build class data
    classes = {}
    tiling_class = {}
    class_scores = {}
    cid = 0
    for bits in range(1 << m):
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        if canon not in classes:
            classes[canon] = cid
            class_scores[cid] = scores(A)
            cid += 1
        tiling_class[bits] = classes[canon]

    # Find GS tilings
    gs_tilings = [b for b in range(1 << m) if is_grid_sym(b)]

    # Find same-class GS flip pairs
    same_class_pairs = []
    for b in gs_tilings:
        fb = b ^ flip_mask
        if tiling_class[b] == tiling_class[fb] and b < fb:
            same_class_pairs.append((b, fb))

    print(f"\nn={n}: {len(gs_tilings)} GS tilings, {len(same_class_pairs)} same-class flip pairs")

    if same_class_pairs:
        print(f"  Same-class GS flip pairs:")
        for (b, fb) in same_class_pairs:
            c = tiling_class[b]
            sc = class_scores[c]
            hw_b = bin(b).count('1')
            hw_fb = bin(fb).count('1')
            print(f"  b={b:0{m}b} (hw={hw_b}) <-> fb={fb:0{m}b} (hw={hw_fb})")
            print(f"    Class {c}, scores={sc}")
            # Show the tournaments
            Ab = bits_to_adj(b)
            Afb = bits_to_adj(fb)
            print(f"    T(b) edges: ", end="")
            for i in range(n):
                for j in range(n):
                    if Ab[i][j] and j != i+1:
                        print(f"{verts[i]}->{verts[j]} ", end="")
            print()
            print(f"    T(fb) edges: ", end="")
            for i in range(n):
                for j in range(n):
                    if Afb[i][j] and j != i+1:
                        print(f"{verts[i]}->{verts[j]} ", end="")
            print()
    else:
        print(f"  No same-class GS flip pairs (all cross classes)")

    # Hamming weight analysis of GS flip pairs
    hw_analysis = defaultdict(lambda: [0, 0])
    for b in gs_tilings:
        fb = b ^ flip_mask
        if b < fb:
            same = (tiling_class[b] == tiling_class[fb])
            hw = bin(b).count('1')
            hw_analysis[hw][0 if same else 1] += 1

    print(f"\n  Hamming weight of GS flip pairs (lower hw partner):")
    print(f"  {'hw':>4} {'same':>5} {'cross':>6}")
    for hw in sorted(hw_analysis):
        s, c = hw_analysis[hw]
        print(f"  {hw:>4} {s:>5} {c:>6}")

    return same_class_pairs

def main():
    for n in range(3, 7):
        investigate_same_class_gs_flips(n)

if __name__ == '__main__':
    main()
