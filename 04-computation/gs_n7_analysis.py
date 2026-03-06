#!/usr/bin/env python3
"""
GS structure at n=7: exploit the GS subspace to avoid full 2^15 enumeration.
Instance: opus-2026-03-06-S11

At n=7, m=C(6,2)=15 tiles, so 2^15=32768 total tilings.
GS dimension = #fixed + #pairs. For n=7:
  Tiles: y=1: (7,1)(6,1)(5,1)(4,1)(3,1)  [5]
         y=2: (7,2)(6,2)(5,2)(4,2)        [4]
         y=3: (7,3)(6,3)(5,3)             [3]
         y=4: (7,4)(6,4)                  [2]
         y=5: (7,5)                       [1]
  Trans: (x,y) -> (8-y, 8-x)
  (7,1)->(7,1) FIXED
  (6,1)->(7,2) pair
  (5,1)->(7,3) pair
  (4,1)->(7,4) pair
  (3,1)->(7,5) pair
  (7,2)->(6,1) (same pair)
  (6,2)->(6,2) FIXED
  (5,2)->(6,3) pair
  (4,2)->(6,4) pair
  (7,3)->(5,1) (same pair)
  (6,3)->(5,2) (same pair)
  (5,3)->(5,3) FIXED
  (7,4)->(4,1) (same pair)
  (6,4)->(4,2) (same pair)
  (7,5)->(3,1) (same pair)

  Fixed: (7,1), (6,2), (5,3) = 3 fixed tiles
  Pairs: (6,1)-(7,2), (5,1)-(7,3), (4,1)-(7,4), (3,1)-(7,5), (5,2)-(6,3), (4,2)-(6,4) = 6 pairs
  GS-dim = 3 + 6 = 9, so #GS = 2^9 = 512 (manageable!)
"""
from itertools import permutations
from collections import defaultdict, Counter
import sys
import time

def build_n7_gs():
    n = 7
    verts = list(range(n, 0, -1))  # [7,6,5,4,3,2,1]
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}
    flip_mask = (1 << m) - 1

    print(f"n={n}, m={m}, total tilings = {1<<m}")

    # Transpose map
    trans_map = []
    fixed_tiles = []
    paired_tiles = []
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
    print(f"Fixed tiles: {len(fixed_tiles)}, Paired: {len(paired_tiles)}, GS-dim: {gs_dim}")
    print(f"#GS tilings = 2^{gs_dim} = {2**gs_dim}")

    # Generate all GS tilings
    def generate_gs_tilings():
        """Generate all GS tilings by choosing free bits."""
        gs_list = []
        for code in range(1 << gs_dim):
            bits = 0
            # Set fixed tile bits
            for k, fix_idx in enumerate(fixed_tiles):
                if (code >> k) & 1:
                    bits |= (1 << fix_idx)
            # Set paired tile bits (both in pair get same value)
            for k, (i, j) in enumerate(paired_tiles):
                if (code >> (len(fixed_tiles) + k)) & 1:
                    bits |= (1 << i)
                    bits |= (1 << j)
            gs_list.append(bits)
        return gs_list

    gs_tilings = generate_gs_tilings()
    print(f"Generated {len(gs_tilings)} GS tilings")

    # Build adjacency and canonicalize for GS tilings
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

    # Canonicalize all GS tilings
    print("Canonicalizing GS tilings...")
    t0 = time.time()

    gs_classes = {}
    gs_class_members = defaultdict(list)
    gs_class_scores = {}
    gs_class_H = {}
    gs_tiling_class = {}
    cid = 0

    for idx, bits in enumerate(gs_tilings):
        if idx % 100 == 0 and idx > 0:
            elapsed = time.time() - t0
            print(f"  {idx}/{len(gs_tilings)} ({elapsed:.1f}s)")
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        if canon not in gs_classes:
            gs_classes[canon] = cid
            gs_class_scores[cid] = scores(A)
            gs_class_H[cid] = hamiltonian_paths(A)
            cid += 1
        c = gs_classes[canon]
        gs_class_members[c].append(bits)
        gs_tiling_class[bits] = c

    num_gs_classes = cid
    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s. Found {num_gs_classes} distinct GS classes (SC tournament classes)")

    # GS flip analysis
    gs_flip_pairs = []
    gs_seen = set()
    for bits in gs_tilings:
        if bits in gs_seen:
            continue
        fbits = bits ^ flip_mask
        # Verify fbits is also GS
        is_gs = True
        for i in range(m):
            j = trans_map[i]
            if j >= 0 and j != i and ((fbits >> i) & 1) != ((fbits >> j) & 1):
                is_gs = False
                break
        assert is_gs, "Flip of GS tiling is not GS!"
        gs_flip_pairs.append((bits, fbits))
        gs_seen.add(bits)
        gs_seen.add(fbits)

    gs_same_class = 0
    gs_cross_class = 0
    gs_cross_edges = defaultdict(int)  # (min_class, max_class) -> count
    for (b1, b2) in gs_flip_pairs:
        c1 = gs_tiling_class[b1]
        c2 = gs_tiling_class[b2]
        if c1 == c2:
            gs_same_class += 1
        else:
            gs_cross_class += 1
            key = (min(c1, c2), max(c1, c2))
            gs_cross_edges[key] += 1

    print(f"\nGS flip pairs: {len(gs_flip_pairs)}")
    print(f"  Same-class: {gs_same_class}")
    print(f"  Cross-class: {gs_cross_class}")
    print(f"  Distinct cross-class edges: {len(gs_cross_edges)}")

    # GS class details
    print(f"\n--- GS Classes at n=7 ---")
    print(f"{'Cls':>4} {'#GS':>5} {'H':>5} {'Scores'}")
    for c in range(num_gs_classes):
        gc = len(gs_class_members[c])
        H = gs_class_H[c]
        sc = gs_class_scores[c]
        print(f"  {c:>3} {gc:>5} {H:>5} {sc}")

    # GS cross-class graph structure
    print(f"\n--- GS cross-class flip graph ---")
    # Compute degree of each GS class in the cross-class graph
    degree = defaultdict(int)
    for (c1, c2), cnt in gs_cross_edges.items():
        degree[c1] += cnt
        degree[c2] += cnt

    print(f"{'Cls':>4} {'Degree':>7} {'Scores'}")
    for c in range(num_gs_classes):
        d = degree.get(c, 0)
        sc = gs_class_scores[c]
        print(f"  {c:>3} {d:>7} {sc}")

    # Check if the GS cross-class graph has any interesting topology
    # (trees? cycles? components?)
    from collections import deque
    adj = defaultdict(set)
    for (c1, c2) in gs_cross_edges:
        adj[c1].add(c2)
        adj[c2].add(c1)

    # BFS to find connected components
    visited = set()
    components = []
    for c in range(num_gs_classes):
        if c in visited:
            continue
        comp = []
        queue = deque([c])
        visited.add(c)
        while queue:
            v = queue.popleft()
            comp.append(v)
            for u in adj[v]:
                if u not in visited:
                    visited.add(u)
                    queue.append(u)
        components.append(comp)

    print(f"\nConnected components: {len(components)}")
    for i, comp in enumerate(components):
        E = sum(1 for (c1,c2) in gs_cross_edges if c1 in comp and c2 in comp)
        print(f"  Component {i}: {len(comp)} vertices, {E} edges, excess = {E - len(comp) + 1}")

    # Hamming weight distribution of GS tilings at n=7
    hw_dist = Counter()
    for bits in gs_tilings:
        hw = bin(bits).count('1')
        hw_dist[hw] += 1
    print(f"\nHamming weight distribution of GS tilings:")
    for hw in sorted(hw_dist):
        print(f"  hw={hw:>2}: {hw_dist[hw]:>3} tilings")


if __name__ == '__main__':
    build_n7_gs()
