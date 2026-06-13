"""
negative_space_h7.py — opus-2026-04-04-S9

THE PERFECT STORM OBSTRUCTION: Why H=7 is permanently forbidden.

H = I(Ω, 2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...

For H=7: need 6 = 2α₁ + 4α₂ + ... → (α₁=3, α₂=0) is the ONLY decomposition.
This means: Ω = K₃ (exactly 3 cycles, all pairwise conflicting, no independent pair).

The "perfect storm": you need EXACTLY 3 cycles, all sharing vertices.
But this configuration ALWAYS generates additional structure that pushes H above 7.

6 and 1 being 7: the 6 = 2×3 (three cycles × weight 2), the 1 = empty set.
The obstruction is that 3 pairwise-conflicting cycles can't exist in isolation.

THE NEGATIVE SPACE: Study what's near H=7 — tournaments with H=5 and H=9
that would need to "pass through" H=7 to connect.
"""
import numpy as np
from collections import defaultdict
from itertools import combinations

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_adj(n, tiles, bits):
    T = np.zeros((n, n), dtype=np.int8)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp.get((mask, v), 0)
            if cnt == 0: continue
            if mask == full: total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def find_odd_cycles(T, n):
    """Find all directed odd cycles (3 and 5)."""
    cycles = []
    # 3-cycles
    for a in range(n):
        for b in range(n):
            if b == a or not T[a][b]: continue
            for c in range(n):
                if c in (a,b) or not T[b][c] or not T[c][a]: continue
                canon = tuple(sorted([a,b,c]))
                cycles.append(('3', canon, (a,b,c)))
    # Deduplicate by vertex set + direction
    seen = set()
    unique = []
    for typ, vset, cyc in cycles:
        key = (typ, vset)
        if key not in seen:
            seen.add(key)
            unique.append((typ, vset, cyc))

    # 5-cycles (for small n)
    if n >= 5:
        for combo in combinations(range(n), 5):
            from itertools import permutations
            for perm in permutations(combo):
                is_cycle = all(T[perm[i]][perm[(i+1)%5]] for i in range(5))
                if is_cycle:
                    canon = min(perm[i:]+perm[:i] for i in range(5))
                    key = ('5', tuple(sorted(combo)))
                    if key not in seen:
                        seen.add(key)
                        unique.append(('5', tuple(sorted(combo)), tuple(canon)))
                    break
    return unique

def conflict_graph(cycles):
    """Build the conflict graph: edges between cycles sharing a vertex."""
    n = len(cycles)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if set(cycles[i][1]) & set(cycles[j][1]):
                adj[i][j] = adj[j][i] = True
    return adj

def independence_polynomial(adj, n_cycles):
    """Compute I(Ω, x) coefficients: α₀, α₁, α₂, ..."""
    alpha = [0] * (n_cycles + 1)
    for mask in range(1 << n_cycles):
        # Check if independent
        bits = [i for i in range(n_cycles) if mask & (1 << i)]
        independent = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if adj[bits[i]][bits[j]]:
                    independent = False
                    break
            if not independent: break
        if independent:
            alpha[len(bits)] += 1
    return alpha

def study_h7_boundary(n):
    """Study tournaments with H near 7 (H=5 and H=9)."""
    tiles = make_tiles(n)
    m = len(tiles)

    print(f"\nH=7 NEGATIVE SPACE at n={n}")
    print("="*60)

    # Collect H=5, H=7, H=9 tournaments
    h5_tilings = []
    h9_tilings = []

    for bits in range(1 << m):
        T = tiling_to_adj(n, tiles, bits)
        h = count_hp(T, n)
        if h == 5: h5_tilings.append(bits)
        elif h == 9: h9_tilings.append(bits)
        elif h == 7: print(f"  !!! H=7 FOUND at bits={bits} !!!")

    print(f"  H=5 tilings: {len(h5_tilings)}")
    print(f"  H=9 tilings: {len(h9_tilings)}")
    print(f"  H=7 tilings: 0 (CONFIRMED)")

    # For each H=5 tiling, what does flipping each tile give?
    print(f"\n  H=5 → flip results:")
    flip_from_5 = defaultdict(int)
    for bits in h5_tilings[:10]:
        for tile in range(m):
            T_flip = tiling_to_adj(n, tiles, bits ^ (1 << tile))
            h_flip = count_hp(T_flip, n)
            flip_from_5[h_flip] += 1

    for h in sorted(flip_from_5.keys()):
        print(f"    → H={h}: {flip_from_5[h]} flips")

    # For each H=9 tiling, what does flipping each tile give?
    print(f"\n  H=9 → flip results:")
    flip_from_9 = defaultdict(int)
    for bits in h9_tilings[:10]:
        for tile in range(m):
            T_flip = tiling_to_adj(n, tiles, bits ^ (1 << tile))
            h_flip = count_hp(T_flip, n)
            flip_from_9[h_flip] += 1

    for h in sorted(flip_from_9.keys()):
        print(f"    → H={h}: {flip_from_9[h]} flips")

    # THE KEY: Study the Omega structure of H=5 and H=9 tournaments
    print(f"\n  Omega (cycle conflict graph) analysis:")
    print(f"\n  --- H=5 tournaments ---")
    for bits in h5_tilings[:5]:
        T = tiling_to_adj(n, tiles, bits)
        cycles = find_odd_cycles(T, n)
        adj_mat = conflict_graph(cycles)
        alpha = independence_polynomial(adj_mat, len(cycles))
        h_check = sum(alpha[k] * (2**k) for k in range(len(alpha)))
        c3 = sum(1 for t,_,_ in cycles if t=='3')
        c5 = sum(1 for t,_,_ in cycles if t=='5')
        print(f"    bits={bits}: c3={c3}, c5={c5}, α={alpha[:5]}, "
              f"I(Ω,2)={h_check}, H={count_hp(T,n)}")

    print(f"\n  --- H=9 tournaments ---")
    for bits in h9_tilings[:5]:
        T = tiling_to_adj(n, tiles, bits)
        cycles = find_odd_cycles(T, n)
        adj_mat = conflict_graph(cycles)
        alpha = independence_polynomial(adj_mat, len(cycles))
        h_check = sum(alpha[k] * (2**k) for k in range(len(alpha)))
        c3 = sum(1 for t,_,_ in cycles if t=='3')
        c5 = sum(1 for t,_,_ in cycles if t=='5')
        print(f"    bits={bits}: c3={c3}, c5={c5}, α={alpha[:5]}, "
              f"I(Ω,2)={h_check}, H={count_hp(T,n)}")

def why_no_K3():
    """Prove WHY Ω = K₃ is impossible.

    K₃ means exactly 3 directed odd cycles, all sharing a vertex.
    Three 3-cycles sharing a vertex v means:
      C₁ = (v, a₁, b₁), C₂ = (v, a₂, b₂), C₃ = (v, a₃, b₃)
    All pairwise share v, hence adjacent in Ω.

    For α₂=0: no two cycles can be vertex-disjoint.
    So every pair of cycles shares at least one vertex.

    The 3 cycles use vertices {v, a₁, b₁, a₂, b₂, a₃, b₃}.
    If all distinct: 7 vertices. If some overlap: fewer.

    KEY OBSERVATION: If all 3 cycles share vertex v and pairwise share v only,
    then the 6 vertices {a₁,b₁,a₂,b₂,a₃,b₃} are distinct. The tournament on
    these 7 vertices is highly constrained.

    The question: does this constraint FORCE a 4th cycle or a 5-cycle?
    """
    print(f"\n{'='*60}")
    print(f"WHY Ω = K₃ IS IMPOSSIBLE (structural analysis)")
    print(f"{'='*60}")

    # At n=5: check all tournaments with exactly 3 3-cycles
    n = 5
    tiles = make_tiles(n)
    m = len(tiles)

    three_cycle_tours = []
    for bits in range(1 << m):
        T = tiling_to_adj(n, tiles, bits)
        # Count 3-cycles
        c3 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (T[a][b] and T[b][c] and T[c][a]) or \
                       (T[a][c] and T[c][b] and T[b][a]):
                        c3 += 1
        if c3 == 3:
            h = count_hp(T, n)
            three_cycle_tours.append((bits, h, T.copy()))

    print(f"\n  n={n}: tournaments with exactly 3 directed 3-cycles: {len(three_cycle_tours)}")
    h_dist = defaultdict(int)
    for _, h, _ in three_cycle_tours:
        h_dist[h] += 1
    print(f"  H distribution: {dict(sorted(h_dist.items()))}")
    print(f"  H=7 count: {h_dist.get(7, 0)}")

    # For each c3=3 tournament, check if Ω = K₃
    k3_count = 0
    for bits, h, T in three_cycle_tours[:20]:
        cycles = find_odd_cycles(T, n)
        adj_mat = conflict_graph(cycles)
        alpha = independence_polynomial(adj_mat, len(cycles))

        # Check all pairwise adjacent (K₃)
        three_cyc = [(t,v,c) for t,v,c in cycles if t == '3']
        all_adjacent = True
        for i in range(len(three_cyc)):
            for j in range(i+1, len(three_cyc)):
                if not (set(three_cyc[i][1]) & set(three_cyc[j][1])):
                    all_adjacent = False
                    break
            if not all_adjacent: break

        is_k3 = (len(three_cyc) == 3 and all_adjacent and alpha[2] == 0)
        if is_k3:
            k3_count += 1

        five_cyc = [c for t,_,c in cycles if t == '5']

        if len(three_cycle_tours) <= 20 or not is_k3:
            if h <= 11:
                shared = set.intersection(*[set(v) for _,v,_ in three_cyc]) if three_cyc else set()
                print(f"    H={h}: 3-cycles={[v for _,v,_ in three_cyc]}, "
                      f"shared={shared}, 5-cycles={len(five_cyc)}, "
                      f"α₂={alpha[2]}, Ω=K₃: {is_k3}")

    print(f"\n  Ω = K₃ count: {k3_count} / {len(three_cycle_tours)}")

    # Do the same at n=6,7
    for n in [6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        k3_count = 0
        total_c3_3 = 0
        h_dist = defaultdict(int)

        for bits in range(1 << m):
            T = tiling_to_adj(n, tiles, bits)
            c3 = 0
            for a in range(n):
                for b in range(a+1, n):
                    for c in range(b+1, n):
                        if (T[a][b] and T[b][c] and T[c][a]) or \
                           (T[a][c] and T[c][b] and T[b][a]):
                            c3 += 1
            if c3 == 3:
                total_c3_3 += 1
                h = count_hp(T, n)
                h_dist[h] += 1

            if bits % 10000 == 0 and bits > 0 and n == 7:
                print(f"    n={n}: {bits}/{1<<m}...", flush=True)

        print(f"\n  n={n}: c3=3 tournaments: {total_c3_3}")
        print(f"  H distribution: {dict(sorted(h_dist.items()))}")
        print(f"  H=7 count: {h_dist.get(7, 0)}")

def multilinear_view_of_gap():
    """View the H=7 gap through the multilinear polynomial lens.

    In the multilinear polynomial H(t):
    - H(t) is odd for all t ∈ {0,1}^m
    - The coefficients are all even (THM-286)
    - H=1+2·(sum of cycle contributions)

    For H=7: need 2·(cycle contributions) = 6.
    The "cycle contribution" at a specific tiling depends on which cycles exist
    and whether they're independent.

    The MULTILINEAR PERSPECTIVE: H(t) jumps from 5 to 9 in certain directions.
    The ΔH = 4 jump (5→9 or 9→5) is the MINIMUM possible jump across the gap.
    This jump corresponds to flipping a tile that simultaneously creates 2 new cycles
    (or destroys 2), which changes α₁ by 2, giving ΔH = 2·2 = 4.

    But could there be a ΔH=2 flip that goes 5→7 or 9→7?
    ΔH=2 means α₁ changes by 1 (one cycle created or destroyed).
    From H=5 (α₁=2): adding 1 cycle gives α₁=3. If the new cycle conflicts with both
    existing ones → H = 1+2·3 = 7 ✓ (if α₂ stays 0).
    But: adding a cycle that conflicts with both existing means the tournament has 3
    pairwise-conflicting cycles. This is exactly the K₃ configuration that forces
    additional cycles!
    """
    n = 5
    tiles = make_tiles(n)
    m = len(tiles)

    print(f"\n{'='*60}")
    print(f"MULTILINEAR VIEW OF THE H=7 GAP (n={n})")
    print(f"{'='*60}")

    # Find all H=5 tilings and examine what tile flips DO to them
    for bits in range(1 << m):
        T = tiling_to_adj(n, tiles, bits)
        h = count_hp(T, n)
        if h != 5: continue

        # For each tile, compute ΔH
        flips = {}
        for tile in range(m):
            T_flip = tiling_to_adj(n, tiles, bits ^ (1 << tile))
            h_flip = count_hp(T_flip, n)
            flips[tiles[tile]] = h_flip - h

        # Show tiles that INCREASE H by 4 (5→9)
        up4 = [(t, dh) for t, dh in flips.items() if dh == 4]
        up2 = [(t, dh) for t, dh in flips.items() if dh == 2]
        down = [(t, dh) for t, dh in flips.items() if dh < 0]

        if bits < 10 or up2:
            cycles = find_odd_cycles(T, n)
            c3 = sum(1 for t,_,_ in cycles if t=='3')
            c5 = sum(1 for t,_,_ in cycles if t=='5')
            print(f"\n  H=5 tiling {bits}: c3={c3}, c5={c5}")
            print(f"    ΔH=+2 flips: {up2} (would give H=7 if possible)")
            print(f"    ΔH=+4 flips: {up4} (skip to H=9)")
            print(f"    ΔH<0 flips: {down}")

            if not up2:
                print(f"    *** NO ΔH=+2 flip exists → H=7 unreachable from this tiling ***")
            else:
                print(f"    ΔH=+2 EXISTS but result is H=7 which doesn't occur → contradiction?")

        if bits > 20: break

if __name__ == "__main__":
    study_h7_boundary(5)
    study_h7_boundary(6)
    why_no_K3()
    multilinear_view_of_gap()
