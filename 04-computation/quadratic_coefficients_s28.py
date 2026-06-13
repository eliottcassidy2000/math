#!/usr/bin/env python3
"""
Order-2 coefficients of H as multilinear polynomial in tiling bits.
opus-2026-04-03-S28

The order-2 coefficient c_{ij} = H(flip i,j) - H(flip i) - H(flip j) + H(0)

where H(0) = 1 (transitive), H(flip i) = 1 + 2^(skip_i - 1).

So: c_{ij} = H(flip i,j) - 2^(skip_i - 1) - 2^(skip_j - 1) - 1

We need H(transitive + flip two arcs).

Tiles (x1,y1) and (x2,y2) can either:
  1. SHARE a vertex (adjacent in the conflict graph)
  2. Be DISJOINT (non-adjacent)

For shared-vertex pairs: the interaction should capture how the two
flipped arcs create/destroy Hamiltonian paths jointly.

For disjoint pairs: the new paths can use both flipped arcs independently,
and the interaction captures the "product" effect.

HYPOTHESIS: c_{ij} = -2^(skip_i + skip_j - 3) for shared-vertex pairs?
At n=4: all shared pairs had c=-2.
  AC: skip=(3,2), -2^(3+2-3) = -2^2 = -4? No, actual was -2.
  Hmm, that doesn't work.

At n=5: shared pairs all had c=-2.
  AB: skip=(4,3), actual=-2. But -2^(4+3-3) = -2^4 = -16. Way off.

So the coefficient for shared-vertex pairs might always be -2.
Let me check more carefully.
"""

def tournament_matrix(n, flips=None):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                T[i][j] = 1
    if flips:
        for (x, y) in flips:
            T[x][y] = 0
            T[y][x] = 1
    return T

def count_hp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def get_tiles(n):
    """1-indexed tiles (x,y) with x-y >= 2."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

for n in range(4, 9):
    tiles = get_tiles(n)
    m = len(tiles)
    print(f"\n{'='*70}")
    print(f"ORDER-2 COEFFICIENTS: n={n}, m={m}")
    print(f"{'='*70}")

    # Compute H(0) = 1
    H0 = 1

    # Compute H(flip i) for each tile
    H1 = {}
    for i in range(m):
        x, y = tiles[i]
        T = tournament_matrix(n, flips=[(x-1, y-1)])
        H1[i] = count_hp(T, n)

    # Compute H(flip i, flip j) for each pair
    results = []
    for i in range(m):
        for j in range(i+1, m):
            x1, y1 = tiles[i]
            x2, y2 = tiles[j]
            T = tournament_matrix(n, flips=[(x1-1, y1-1), (x2-1, y2-1)])
            H2 = count_hp(T, n)

            c2 = H2 - H1[i] - H1[j] + H0

            shared = set(tiles[i]) & set(tiles[j])
            skip_i = x1 - y1
            skip_j = x2 - y2

            results.append((i, j, c2, shared, skip_i, skip_j, tiles[i], tiles[j]))

    # Categorize
    shared_results = [(i,j,c,s,si,sj,ti,tj) for i,j,c,s,si,sj,ti,tj in results if s]
    disjoint_results = [(i,j,c,s,si,sj,ti,tj) for i,j,c,s,si,sj,ti,tj in results if not s]

    print(f"\nSHARED-VERTEX pairs ({len(shared_results)}):")
    # Check: are all c = -2?
    all_minus_2 = all(c == -2 for _,_,c,_,_,_,_,_ in shared_results)
    print(f"  ALL coefficients = -2? {all_minus_2}")
    if not all_minus_2:
        for i,j,c,s,si,sj,ti,tj in shared_results:
            if c != -2:
                print(f"  ({ti[0]},{ti[1]})×({tj[0]},{tj[1]}) shared={s} "
                      f"skip=({si},{sj}): c2={c}")

    # Show all shared-vertex results for small n
    if n <= 6:
        for i,j,c,s,si,sj,ti,tj in shared_results:
            # Which vertex is shared?
            shared_v = list(s)[0]
            # What is the relationship of the shared vertex to the arcs?
            # Tile (x1,y1): connects x1 and y1. Shared vertex v.
            # So v is either x1 or y1, and either x2 or y2.
            v = shared_v
            pos_i = 'x' if v == ti[0] else 'y'  # v is upper or lower end of tile i
            pos_j = 'x' if v == tj[0] else 'y'
            config = pos_i + pos_j  # e.g., "xy" means v is upper for i, lower for j

            print(f"  {chr(65+i)}{chr(65+j)} ({ti[0]},{ti[1]})×({tj[0]},{tj[1]}) "
                  f"shared v={v} config={config}: c2={c}")

    print(f"\nDISJOINT pairs ({len(disjoint_results)}):")
    if n <= 6:
        for i,j,c,s,si,sj,ti,tj in disjoint_results:
            print(f"  {chr(65+i)}{chr(65+j)} ({ti[0]},{ti[1]})×({tj[0]},{tj[1]}) "
                  f"skip=({si},{sj}): c2={c}")
    else:
        # Just show statistics
        c2_vals = [c for _,_,c,_,_,_,_,_ in disjoint_results]
        from collections import Counter
        print(f"  Value distribution: {dict(sorted(Counter(c2_vals).items()))}")

    # DISJOINT analysis: is c_{ij} related to the overlap of the "new path" sets?
    # When tiles are disjoint, flipping both creates paths that use BOTH new arcs.
    # The number of such "double-new" paths should be the interaction term.
    #
    # If tile i creates 2^(si-1) new paths and tile j creates 2^(sj-1) new paths,
    # and their between-vertex sets are disjoint (which happens iff the tiles are disjoint),
    # then the double-flip creates 2^(si-1)*2^(sj-1) = 2^(si+sj-2) "both-new" paths
    # PLUS the individual new paths.
    #
    # H(flip i,j) = 1 + 2^(si-1) + 2^(sj-1) + interaction
    # c2 = interaction

    # Let me check: for disjoint pairs, is c2 = 2^(si+sj-2)?
    if disjoint_results:
        print(f"\nDisjoint: is c2 = 2^(si+sj-2)?")
        for i,j,c,s,si,sj,ti,tj in disjoint_results[:20]:
            predicted = 2**(si+sj-2)
            label = f"{chr(65+i)}{chr(65+j)}" if i < 26 and j < 26 else f"({i},{j})"
            match = "✓" if c == predicted else "✗"
            print(f"  {label} ({ti[0]},{ti[1]})×({tj[0]},{tj[1]}) "
                  f"skip=({si},{sj}): c2={c}, 2^({si}+{sj}-2)={predicted} {match}")

    # SUMMARY for shared vertex: always c = -2?
    # This would mean: when two tiles share a vertex, flipping both
    # DESTROYS 2 paths compared to the sum of individual effects.
    # This is the "overcounting" from shared vertices.

    # Let me also check: for shared vertex, is the coefficient always
    # -2 regardless of n, skip values, or which vertex is shared?
    shared_all_c = [c for _,_,c,_,_,_,_,_ in shared_results]
    if shared_all_c:
        min_c = min(shared_all_c)
        max_c = max(shared_all_c)
        print(f"\nShared-vertex c2 range: [{min_c}, {max_c}]")
        if min_c == max_c == -2:
            print(f"  ALL SHARED-VERTEX COEFFICIENTS = -2")
        else:
            from collections import Counter
            print(f"  Distribution: {dict(sorted(Counter(shared_all_c).items()))}")
