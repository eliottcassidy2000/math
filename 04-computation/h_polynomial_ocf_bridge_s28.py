#!/usr/bin/env python3
"""
Bridging the multilinear polynomial H(x) with the OCF formula H = I(Ω,2).
opus-2026-04-03-S28

The multilinear polynomial: H(x) = Σ_S c_S · Π_{i∈S} x_i
The OCF formula: H(T) = I(Ω(T), 2) = Σ_k α_k(T) · 2^k

Both are representations of H. The multilinear polynomial is in terms of
TILING bits (which arcs to flip from transitive). The OCF is in terms of
the independence polynomial of the conflict graph.

QUESTION: How does the conflict graph Ω(T) change as we flip tile bits?

For the TRANSITIVE tournament (all bits = 0):
  - H = 1, so I(Ω, 2) = 1
  - This means α_k = 0 for all k ≥ 1
  - The transitive tournament has NO directed odd cycles!
  - So Ω is the empty graph.

When we flip tile (x,y) from 0 to 1:
  - We reverse the arc x→y to y→x
  - This can CREATE new directed cycles (that weren't there before)
  - H increases to 1 + 2^(skip-1)
  - So I(Ω_new, 2) = 1 + 2^(skip-1)
  - This means α_1(Ω_new) · 2 + ... = 2^(skip-1)
  - For skip=2: α_1 = 1 (one 3-cycle created)
  - For skip=3: α_1 · 2 = 4 → α_1 = 2 (two 3-cycles? or one 5-cycle giving α_1=0 + one independent set structure?)

Actually: I(Ω, 2) = 1 + α_1 · 2 + α_2 · 4 + ...
  1 + 2^(s-1) = 1 + 2^(s-1)

For skip=2: 2^1 = 2. So α_1·2 = 2, α_1 = 1. One odd cycle created.
For skip=3: 2^2 = 4. So α_1·2 + α_2·4 = 4.
  Options: α_1=2, α_2=0 (two vertex-disjoint cycles creating 4)
  Or α_1=0, α_2=1 (one pair of independent cycles creating 4)

We need α_1·2 + α_2·4 = 4. Since α's are non-negative integers:
  (α_1, α_2) ∈ {(2,0), (0,1)}. Both give 4.

Which is it? When we flip arc (x,y) with skip=3 from the transitive:
  Vertices x, x-1, y (3 vertices between base-path and flipped arc)
  This creates a 3-cycle: y → x → x-1 → ... → y?

Let me just compute the actual conflict graphs.
"""

from itertools import combinations

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

def find_directed_odd_cycles(T, n, max_length=None):
    """Find all directed odd cycles in tournament T."""
    if max_length is None:
        max_length = n
    cycles = []
    for length in range(3, max_length + 1, 2):  # odd lengths: 3, 5, 7, ...
        for combo in combinations(range(n), length):
            # Check all cyclic orderings
            from itertools import permutations
            for perm in permutations(combo):
                if all(T[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    # Normalize: smallest vertex first
                    min_idx = perm.index(min(perm))
                    normalized = perm[min_idx:] + perm[:min_idx]
                    if normalized not in cycles:
                        cycles.append(normalized)
    return cycles

def conflict_graph(cycles):
    """Build conflict graph: cycles adjacent iff they share a vertex."""
    n = len(cycles)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True
    return adj

def independence_polynomial(adj, n_nodes):
    """Compute independence polynomial coefficients."""
    from itertools import combinations
    alpha = [0] * (n_nodes + 1)
    alpha[0] = 1
    for k in range(1, n_nodes + 1):
        for combo in combinations(range(n_nodes), k):
            # Check if independent
            independent = True
            for i in range(len(combo)):
                for j in range(i+1, len(combo)):
                    if adj[combo[i]][combo[j]]:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                alpha[k] += 1
    return alpha

for n in range(3, 7):
    m = (n-1)*(n-2)//2
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    print(f"\n{'='*70}")
    print(f"MULTILINEAR ↔ OCF BRIDGE: n={n}")
    print(f"{'='*70}")

    # For each single-flip from transitive:
    for ti, (x, y) in enumerate(tiles):
        skip = x - y
        T = tournament_matrix(n, flips=[(x-1, y-1)])
        H = count_hp(T, n)

        cycles = find_directed_odd_cycles(T, n)
        adj = conflict_graph(cycles)
        alpha = independence_polynomial(adj, len(cycles))

        # Compute I(Ω, 2)
        I_val = sum(alpha[k] * (2**k) for k in range(len(alpha)))

        label = chr(65+ti) if ti < 26 else f'T{ti}'
        print(f"\n  Flip tile {label} ({x},{y}) skip={skip}:")
        print(f"    H = {H} (should be 1 + 2^{skip-1} = {1 + 2**(skip-1)})")
        print(f"    Odd cycles: {len(cycles)}")
        for i, c in enumerate(cycles):
            verts = [v+1 for v in c]
            length = len(c)
            print(f"      C{i}: {verts} (length {length})")
        print(f"    α = {alpha[:5]}")
        print(f"    I(Ω, 2) = {I_val}")
        assert H == I_val, f"OCF MISMATCH: H={H} ≠ I={I_val}"

    # Now for DOUBLE flips, check the OCF decomposition
    print(f"\n  --- DOUBLE FLIPS ---")
    for ti in range(min(m, 5)):
        for tj in range(ti+1, min(m, 5)):
            x1, y1 = tiles[ti]
            x2, y2 = tiles[tj]
            skip1, skip2 = x1-y1, x2-y2

            T = tournament_matrix(n, flips=[(x1-1, y1-1), (x2-1, y2-1)])
            H = count_hp(T, n)

            cycles = find_directed_odd_cycles(T, n)
            adj = conflict_graph(cycles)
            alpha = independence_polynomial(adj, len(cycles))

            I_val = sum(alpha[k] * (2**k) for k in range(len(alpha)))
            assert H == I_val

            shared = set(tiles[ti]) & set(tiles[tj])

            label = f"{chr(65+ti)}{chr(65+tj)}"
            print(f"\n  Flip {label} ({x1},{y1})×({x2},{y2}) skip=({skip1},{skip2}) "
                  f"{'shared='+str(shared) if shared else 'disjoint'}:")
            print(f"    H = {H}, cycles = {len(cycles)}, α = {alpha[:5]}")

            # The c₂ coefficient
            H_i = count_hp(tournament_matrix(n, [(x1-1, y1-1)]), n)
            H_j = count_hp(tournament_matrix(n, [(x2-1, y2-1)]), n)
            c2 = H - H_i - H_j + 1
            print(f"    c₂ = {c2}")
            print(f"    H_decomposition: 1 + {H_i-1} + {H_j-1} + {c2} = {H}")

    # For the COMPLEMENT of the transitive (all bits = 1):
    print(f"\n  --- ALL-ONE TILING (complement of transitive) ---")
    T = tournament_matrix(n, flips=[(x-1, y-1) for x, y in tiles])
    H = count_hp(T, n)
    cycles = find_directed_odd_cycles(T, n)
    adj = conflict_graph(cycles)
    alpha = independence_polynomial(adj, len(cycles))
    I_val = sum(alpha[k] * (2**k) for k in range(len(alpha)))
    print(f"  H = {H}, cycles = {len(cycles)}, α = {alpha}")
    print(f"  I(Ω, 2) = {I_val}")

    # This is the "anti-transitive" tournament: all non-base arcs reversed
    # Its score sequence is...
    scores = sorted([sum(T[i]) for i in range(n)])
    print(f"  Score sequence: {scores}")
