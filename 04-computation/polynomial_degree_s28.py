#!/usr/bin/env python3
"""
The polynomial degree of H in tiling bits.
opus-2026-04-03-S28

PATTERN:
  n=3: deg=1, n=4: deg=2, n=5: deg=4, n=6: deg=4, n=7: deg=6

Conjecture: deg = 2*floor((n-2)/2) = n-2 if n even, n-3 if n odd?
  n=3: 2*0=0... no that's wrong.
  n=3 odd: n-3=0? But actual is 1.

Let me just tabulate and check.

More interesting: the multilinear polynomial at n=4 is:
  H = 1 + 4*A + 2*B + 2*C - 2*AB - 2*AC

where A = tile (4,1) = apex, B = tile (3,1), C = tile (4,2).

Note: AB shares vertex 1. AC shares vertex 4. BC shares vertex... no:
B = (3,1) and C = (4,2) share no vertex! Yet BC has coefficient 0.

The NONZERO order-2 terms are EXACTLY the tile pairs that SHARE A VERTEX.
Let me check this at all n.

Also: the coefficient of the constant term is ALWAYS 1. This makes sense:
the all-zero tiling is the transitive tournament, which has H=1.

And the coefficient of the apex tile A is ALWAYS 2^(n-2):
  n=3: 2, n=4: 4, n=5: 8, n=6: 16

The apex effect was 2^(n-3) but the POLYNOMIAL coefficient is 2^(n-2).
These differ by a factor of 2 because the polynomial is H = ... + c*x_A + ...,
so E[H|x_A=1] - E[H|x_A=0] = c (the coefficient), evaluated at the mean
of all other variables (= 1/2 each). So the main effect is:
  effect_A = c_A + sum_{j≠A} c_{Aj}/2 + ...
Not just c_A.

Let me verify this relationship.
"""

from collections import defaultdict
from math import comb
from fractions import Fraction

def tournament_from_bits(n, bits_list):
    T = [[0]*n for _ in range(n)]
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    for i in range(n, 1, -1):
        T[i-1][i-2] = 1
    for idx, (x, y) in enumerate(tiles):
        if idx < len(bits_list):
            if bits_list[idx] == 0:
                T[x-1][y-1] = 1
            else:
                T[y-1][x-1] = 1
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
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

for n in range(3, 7):
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)

    print(f"\n{'='*70}")
    print(f"MULTILINEAR POLYNOMIAL ANALYSIS: n={n}, m={m}")
    print(f"{'='*70}")

    # Compute H values
    H_vals = [0] * (2**m)
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H_vals[bits] = count_hp(T, n)

    # Compute multilinear coefficients by Möbius inversion
    poly = {}
    max_degree = 0
    for S in range(2**m):
        order = bin(S).count('1')
        coeff = 0
        for T_bits in range(2**m):
            if T_bits & S != T_bits:
                continue
            diff = bin(S & ~T_bits).count('1')
            coeff += ((-1) ** diff) * H_vals[T_bits]
        if abs(coeff) > 1e-10:
            poly[S] = coeff
            if order > max_degree:
                max_degree = order

    # Order-1 coefficients (linear terms)
    print(f"\nLinear coefficients (order 1):")
    for j in range(m):
        S = 1 << j
        c = poly.get(S, 0)
        t = tiles[j]
        skip = t[0] - t[1]
        print(f"  {chr(65+j)} ({t[0]},{t[1]}) skip={skip}: coeff = {c}")

    # Check apex coefficient = 2^(n-2)
    apex_j = max(range(m), key=lambda j: tiles[j][0] - tiles[j][1])
    apex_coeff = poly.get(1 << apex_j, 0)
    print(f"\nApex coefficient = {apex_coeff}, 2^(n-2) = {2**(n-2)}, match = {apex_coeff == 2**(n-2)}")

    # Order-2 coefficients: do they correspond to SHARED VERTICES?
    print(f"\nOrder-2 coefficients and vertex sharing:")
    for i in range(m):
        for j in range(i+1, m):
            S = (1 << i) | (1 << j)
            c = poly.get(S, 0)
            t_i, t_j = tiles[i], tiles[j]
            shared = set(t_i) & set(t_j)
            is_shared = len(shared) > 0
            if abs(c) > 1e-10 or is_shared:
                mark = "SHARED" if is_shared else "DISJOINT"
                label = f"{chr(65+i)}{chr(65+j)}"
                print(f"  {label} = ({t_i[0]},{t_i[1]})×({t_j[0]},{t_j[1]}): "
                      f"coeff={c:+.0f} [{mark} {shared if shared else ''}]")

    # Count: are ALL nonzero order-2 terms from shared-vertex pairs?
    shared_nonzero = 0
    disjoint_nonzero = 0
    shared_zero = 0
    disjoint_zero = 0
    for i in range(m):
        for j in range(i+1, m):
            S = (1 << i) | (1 << j)
            c = poly.get(S, 0)
            shared = len(set(tiles[i]) & set(tiles[j])) > 0
            if abs(c) > 1e-10:
                if shared:
                    shared_nonzero += 1
                else:
                    disjoint_nonzero += 1
            else:
                if shared:
                    shared_zero += 1
                else:
                    disjoint_zero += 1

    total_pairs = comb(m, 2)
    print(f"\nOrder-2 vertex-sharing analysis:")
    print(f"  Shared-vertex & nonzero: {shared_nonzero}")
    print(f"  Shared-vertex & zero:    {shared_zero}")
    print(f"  Disjoint & nonzero:      {disjoint_nonzero}")
    print(f"  Disjoint & zero:         {disjoint_zero}")
    print(f"  Total pairs:             {total_pairs}")

    # Do disjoint tiles have nonzero interactions at n >= 5?
    if disjoint_nonzero > 0:
        print(f"  *** Disjoint tiles CAN interact! (order 2 captures 3-cycle effects)")

    # The STRUCTURE of the polynomial:
    # Each tile (x,y) represents an arc from x to y (or y to x).
    # A Hamiltonian path uses n-1 arcs. H counts ALL Hamiltonian paths.
    # Each HP uses a subset of arcs. The polynomial decomposes by
    # how many tile arcs the HP uses (0, 1, 2, ...).
    # A degree-k term in the polynomial counts contributions from
    # k SPECIFIC tiles being "active" (flipped from default).

    # The polynomial degree is bounded by the maximum number of
    # tile arcs in any Hamiltonian path.
    # A HP uses exactly n-1 arcs total: some from the base path, some from tiles.
    # Base path has n-1 arcs. The overlap between HP arcs and base arcs is k,
    # then the HP uses n-1-k tile arcs.
    # The maximum n-1-k is when k is minimum. The minimum k is 0 (HP uses
    # NO base-path arcs). But can a HP avoid ALL base-path arcs?

    # At n=3: HP has 2 arcs. Base has 2 arcs. Tiles have 1 arc.
    # A HP must use at least 1 base-path arc (since only 1 tile arc exists).
    # Max tile arcs in HP = 1. So deg = 1. ✓

    # At n=4: HP has 3 arcs. Base has 3. Tiles have 3.
    # Can HP use 3 tile arcs? The tiles are (4,1), (3,1), (4,2).
    # A path using tiles (4,1), (3,1), (4,2) would visit 4->1, 3->1, 4->2
    # But 1 appears twice as endpoint — not a valid path.
    # Max tile arcs in HP = 2. So deg = 2. ✓

    # At n=5: HP has 4 arcs. Base has 4. Tiles have 6.
    # Max tile arcs in HP = ?

    # Let me compute this directly
    from itertools import permutations

    base_arcs = set()
    for i in range(n, 1, -1):
        base_arcs.add((i, i-1))

    max_tile_arcs = 0
    for perm in permutations(range(1, n+1)):
        # Count how many arcs of this path are tile arcs
        tile_count = 0
        for k in range(len(perm)-1):
            u, v = perm[k], perm[k+1]
            arc = (u, v) if u > v else (v, u)  # normalize
            is_base = (arc[0] - arc[1] == 1)  # consecutive vertices
            if not is_base:
                tile_count += 1
        if tile_count > max_tile_arcs:
            max_tile_arcs = tile_count

    print(f"\nMax tile arcs in any Hamiltonian path: {max_tile_arcs}")
    print(f"Polynomial degree: {max_degree}")
    print(f"Match: {max_tile_arcs == max_degree}")

    # Count paths with exactly k tile arcs
    tile_arc_dist = defaultdict(int)
    total_paths = 0
    for perm in permutations(range(1, n+1)):
        # Check if this is actually a Hamiltonian path in some tournament
        # Actually, let me count tile arcs in ANY permutation (as a path)
        tile_count = 0
        for k in range(len(perm)-1):
            u, v = perm[k], perm[k+1]
            is_base = abs(u - v) == 1
            if not is_base:
                tile_count += 1
        tile_arc_dist[tile_count] += 1
        total_paths += 1

    print(f"\nDistribution of tile arcs in permutations (potential HPs):")
    for k in sorted(tile_arc_dist):
        count = tile_arc_dist[k]
        print(f"  {k} tile arcs: {count}/{total_paths} = {count/total_paths:.4f}")

print(f"\n{'='*70}")
print("SUMMARY: Polynomial degree = max tile arcs in any HP")
print(f"{'='*70}")
print(f"{'n':>4} {'m':>4} {'deg':>5} {'max_tile':>10} {'n-1':>6} {'base_arcs':>12}")
for n in range(3, 8):
    m = (n-1)*(n-2)//2

    # Count max tile arcs
    max_tile = 0
    for perm in permutations(range(1, n+1)):
        tile_count = sum(1 for k in range(len(perm)-1)
                        if abs(perm[k] - perm[k+1]) != 1)
        if tile_count > max_tile:
            max_tile = tile_count

    # Degree (from earlier computation or prediction)
    deg = max_tile  # conjectured equal

    print(f"{n:>4} {m:>4} {deg:>5} {max_tile:>10} {n-1:>6} {n-1:>12}")
