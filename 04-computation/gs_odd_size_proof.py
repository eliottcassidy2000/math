#!/usr/bin/env python3
"""
Investigate WHY GS class sizes are always odd.
Instance: opus-2026-03-06-S11

Hypothesis: The perpendicular reflection sigma acts on the Hamiltonian path
embeddings of a self-converse tournament T. Since T = T^op, sigma maps each
Ham path to another Ham path. The number of GS tilings in a class equals
the number of Ham paths of T that are "symmetric" (invariant under sigma
composed with path reversal).

If sigma is an involution on Ham paths, then #fixed_points = #total mod 2.
Since #total = H(T), and if H(T) is always odd for SC tournaments...
but H(T) is NOT always odd (e.g., H=15 for the regular 5-tournament, 15 is odd,
but H=45 for a 6-tournament class, 45 is odd too).

Actually, H(T) might always be odd for SC tournaments! Let's check.
Then #GS = H(T) mod 2 (by involution argument) = odd.
And since #GS > 0 (T is SC, so at least one GS embedding), #GS is odd.

Wait, the GS count is the number of GS tilings, not all tilings.
|class| = H(T)/|Aut(T)| (approximately), and #GS = ?

Let me think differently. The GS tiling count for class C is the number of
Hamiltonian path orderings v_1,...,v_n of T such that the resulting tiling
is grid-symmetric. This means: the tiling is invariant under the perpendicular
reflection, which maps tile (x,y) to (n+1-y, n+1-x).

At the level of the tournament: the perpendicular reflection maps the
ordering v_1,...,v_n to some other ordering related to T^op.

Let's verify the odd-size conjecture and investigate the H-values of SC tournaments.
"""
from itertools import permutations
from collections import defaultdict

def check_H_parity(n):
    """Check if H(T) is always odd for SC tournaments at given n."""
    verts = list(range(n, 0, -1))
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}

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
    tiling_class = {}
    class_H = {}
    cid_counter = 0

    for bits in range(1 << m):
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        if canon not in classes:
            classes[canon] = cid_counter
            class_H[cid_counter] = hamiltonian_paths(A)
            cid_counter += 1
        cid = classes[canon]
        class_members[cid].append(bits)
        tiling_class[bits] = cid

    # Count GS tilings per class
    gs_count = defaultdict(int)
    for bits in range(1 << m):
        if is_grid_sym(bits):
            gs_count[tiling_class[bits]] += 1

    # SC classes = those with GS > 0
    sc_classes = [c for c in range(cid_counter) if gs_count[c] > 0]

    print(f"\nn={n}: {cid_counter} classes, {len(sc_classes)} SC")
    print(f"{'Class':>6} {'Size':>6} {'H':>6} {'GS':>4} {'H%2':>4} {'GS%2':>5} {'H/Size':>8}")

    all_H_odd = True
    all_GS_odd = True
    for c in sc_classes:
        sz = len(class_members[c])
        H = class_H[c]
        gc = gs_count[c]
        h_parity = H % 2
        g_parity = gc % 2
        aut = H // sz if sz > 0 else '?'
        if h_parity == 0:
            all_H_odd = False
        if g_parity == 0:
            all_GS_odd = False
        print(f"  {c:>4} {sz:>6} {H:>6} {gc:>4} {h_parity:>4} {g_parity:>5} {aut:>8}")

    print(f"  All H odd for SC? {all_H_odd}")
    print(f"  All GS odd? {all_GS_odd}")

    # Also check: is |Aut(T)| always odd for SC tournaments?
    print(f"  Aut orders: {[class_H[c] // len(class_members[c]) for c in sc_classes]}")

    return all_H_odd, all_GS_odd

def main():
    for n in range(3, 7):
        check_H_parity(n)

if __name__ == '__main__':
    main()
