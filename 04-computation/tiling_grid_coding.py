"""
tiling_grid_coding.py -- kind-pasteur-2026-03-14-S76
Exploring the triangular pin grid, blueself/blackself, blue line skeleton,
and their connections to coding theory.

THE PIN GRID:
- Positions (r,c) with r >= 1, c >= 1, r+c <= n-1
- This is the staircase Young diagram delta_{n-2}
- m = C(n-1, 2) positions (tiles)
- Each tile gets a bit: 0=forward, 1=backward

THE SYMMETRIES:
- GS (Grid-Symmetric): 180° rotation of pin grid
  (r,c) -> (n-1-r-c, n-1-r-c)... actually (r,c) -> (n-1-c-r, n-1-c-r)?
  Let me compute it properly: transpose maps vertex v -> n-1-v,
  which maps arc (a,b) to (n-1-b, n-1-a). In pin grid coordinates:
  (r,c) where r=a-b-1, c=b, so a=r+c+1, b=c.
  Transpose: a'=n-1-c, b'=n-1-(r+c+1)=n-2-r-c, so
  r'=a'-b'-1 = (n-1-c)-(n-2-r-c)-1 = r, c'=b'=n-2-r-c.
  So GS: (r,c) -> (r, n-2-r-c). This is a REFLECTION along the anti-diagonal!

- Flip: reverse all non-backbone arcs. In tiling: bit -> 1-bit for every tile.

CODING THEORY CONNECTIONS:
1. The set of tilings = {0,1}^m = ALL binary words of length m
2. Each tournament has a unique tiling (via base path P_0)
3. H(T) is a function on {0,1}^m of degree d = 2*floor((n-1)/2) [Degree Drop]
4. So H is a codeword of Reed-Muller code RM(d, m)!
5. The GS tilings form a LINEAR subcode (intersection with symmetry constraint)
6. The blueself tilings = GS ∩ self-flip code

QUESTIONS:
1. What is the dimension of the GS subcode?
2. What are the weight enumerators of the H-fiber codes on the pin grid?
3. How does the triangular grid structure affect the code properties?
4. Is the blue line skeleton related to the Tanner graph of a code?
5. What error-correcting properties do tournament tilings have?
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def pin_grid_positions(n):
    """Return all pin grid positions (r, c) for delta_{n-2}."""
    positions = []
    for r in range(1, n-1):
        for c in range(1, n-1-r+1):
            positions.append((r, c))
    return positions

def position_to_arc(r, c):
    """Pin grid position (r,c) -> arc (a,b) where a = r+c+1, b = c."""
    a = r + c + 1
    b = c
    return (a, b)

def arc_to_position(a, b):
    """Arc (a,b) with a > b+1 -> pin grid position (r,c)."""
    r = a - b - 1
    c = b
    return (r, c)

def gs_map(r, c, n):
    """Grid-symmetric map: (r,c) -> (r, n-r-c).
    Derived from vertex transpose v -> n+1-v (1-indexed)."""
    return (r, n-r-c)

def bits_to_adj(bits, n):
    """Convert tiling bits to adjacency matrix.
    Base path P_0 (1-indexed): n -> n-1 -> ... -> 1.
    0-indexed: (n-1) -> (n-2) -> ... -> 0.
    Tiles use 1-indexed vertices: arc (a,b) with a >= b+2.
    Convert to 0-indexed: (a-1, b-1).
    """
    A = np.zeros((n, n), dtype=int)
    # Backbone: vertex i beats vertex i-1 (0-indexed)
    for i in range(1, n):
        A[i][i-1] = 1

    # Non-backbone arcs: tile (a,b) 1-indexed -> (a-1, b-1) 0-indexed
    positions = pin_grid_positions(n)
    for idx, (r, c) in enumerate(positions):
        a, b = position_to_arc(r, c)  # 1-indexed
        a0, b0 = a - 1, b - 1  # 0-indexed
        if bits & (1 << idx):
            A[b0][a0] = 1  # backward: b -> a
        else:
            A[a0][b0] = 1  # forward: a -> b
    return A

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("TRIANGULAR PIN GRID, TILING CODES, AND BLUE LINE SKELETON")
    print("kind-pasteur-2026-03-14-S76")
    print("=" * 70)

    # ========================================
    # PART 1: Pin Grid Structure
    # ========================================
    print(f"\n{'='*70}")
    print("PART 1: THE TRIANGULAR PIN GRID")
    print(f"{'='*70}")

    for n in [4, 5, 6, 7]:
        positions = pin_grid_positions(n)
        m = len(positions)
        print(f"\n  n={n}: m={m} tiles on delta_{{{n-2}}}")
        print(f"  Pin grid:")

        # Display the grid
        for r in range(1, n-1):
            row_pos = [(r, c) for c in range(1, n-1-r+1)]
            arcs = [position_to_arc(r, c) for r, c in row_pos]
            print(f"    r={r}: {['({},{})'.format(a,b) for a,b in arcs]}")

        # GS map: fixed points?
        gs_fixed = []
        gs_pairs = []
        seen = set()
        for r, c in positions:
            r2, c2 = gs_map(r, c, n)
            if (r, c) == (r2, c2):
                gs_fixed.append((r, c))
            elif (r2, c2) not in seen:
                gs_pairs.append(((r, c), (r2, c2)))
                seen.add((r, c))
                seen.add((r2, c2))

        print(f"\n  GS symmetry (r,c) -> (r, {n}-r-c):")
        print(f"    Fixed points: {gs_fixed}")
        print(f"    Paired positions: {len(gs_pairs)}")
        print(f"    GS degrees of freedom: {len(gs_fixed) + len(gs_pairs)}")
        print(f"    #GS tilings = 2^({len(gs_fixed) + len(gs_pairs)}) = {2**(len(gs_fixed) + len(gs_pairs))}")

    # ========================================
    # PART 2: GS Tilings as a Linear Code
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: GS TILINGS AS A LINEAR CODE")
    print("  GS constraint: tile(r,c) = tile(r, n-2-r-c)")
    print("  This is a linear constraint over F_2")
    print(f"{'='*70}")

    for n in [5, 7]:
        positions = pin_grid_positions(n)
        m = len(positions)
        pos_idx = {p: i for i, p in enumerate(positions)}

        # GS constraint matrix over F_2
        # For each paired position, tile(r,c) = tile(r', c') mod 2
        # This is the constraint tile(r,c) + tile(r',c') = 0 mod 2
        constraints = []
        for (p1, p2) in []:
            pass  # build below

        gs_fixed = []
        gs_pairs = []
        for r, c in positions:
            r2, c2 = gs_map(r, c, n)
            if (r, c) == (r2, c2):
                gs_fixed.append((r, c))
            elif pos_idx.get((r, c), -1) < pos_idx.get((r2, c2), -1):
                gs_pairs.append(((r, c), (r2, c2)))

        # GS code: dimension = #fixed + #pairs = m - #pairs (since each pair adds 1 constraint)
        gs_dim = len(gs_fixed) + len(gs_pairs)
        print(f"\n  n={n}: m={m}, GS dim = {gs_dim}, rate = {gs_dim/m:.4f}")
        print(f"    #GS codewords = {2**gs_dim}")

        # Enumerate GS tilings and compute H
        gs_H_values = []
        gs_tilings = []

        for gs_bits in range(2**gs_dim):
            # Expand gs_bits to full tiling
            full_bits = 0
            bit_idx = 0

            # Fixed positions: free bits
            for p in gs_fixed:
                if gs_bits & (1 << bit_idx):
                    full_bits |= (1 << pos_idx[p])
                bit_idx += 1

            # Paired positions: one free bit controls both
            for p1, p2 in gs_pairs:
                if gs_bits & (1 << bit_idx):
                    full_bits |= (1 << pos_idx[p1])
                    full_bits |= (1 << pos_idx[p2])
                bit_idx += 1

            A = bits_to_adj(full_bits, n)
            H = compute_H_dp(A, n)
            gs_H_values.append(H)
            gs_tilings.append(full_bits)

        H_dist = Counter(gs_H_values)
        print(f"\n  GS tiling H distribution:")
        for H in sorted(H_dist.keys()):
            print(f"    H={H}: {H_dist[H]} GS tilings")

        # Hamming weights of GS codewords
        weights = [bin(t).count('1') for t in gs_tilings]
        weight_dist = Counter(weights)
        print(f"\n  Weight distribution of GS tilings:")
        for w in sorted(weight_dist.keys()):
            print(f"    weight {w}: {weight_dist[w]}")

        # Minimum distance between different-H GS tilings
        H_groups = defaultdict(list)
        for i, t in enumerate(gs_tilings):
            H_groups[gs_H_values[i]].append(t)

        # Self-flip: complement of tiling
        flip_mask = (1 << m) - 1
        print(f"\n  Self-flip analysis:")
        blueself_count = 0
        for i, t in enumerate(gs_tilings):
            t_flip = t ^ flip_mask
            # Check if t_flip is also a GS tiling with same H
            if t_flip in set(gs_tilings):
                j = gs_tilings.index(t_flip)
                if gs_H_values[j] == gs_H_values[i]:
                    blueself_count += 1

        print(f"    Blueself (GS + self-flip same H): {blueself_count} tilings")
        print(f"    Expected: {0 if n%2==1 else 'some'} (THM-023: blueself requires even n)")

    # ========================================
    # PART 3: Reed-Muller Connection
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: REED-MULLER CODE CONNECTION")
    print("  H(T) has degree d = 2*floor((n-1)/2) [Degree Drop Theorem]")
    print("  So H is in RM(d, m) — the Reed-Muller code of order d on m variables")
    print("  RM(d,m) has dimension sum_{k=0}^{d} C(m,k)")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        m = n * (n-1) // 2  # This is C(n,2), the number of arc variables
        m_tiling = (n-1) * (n-2) // 2  # This is C(n-1,2), the number of TILING positions
        d = 2 * ((n-1) // 2)

        rm_dim = sum(math.comb(m, k) for k in range(d+1))
        rm_dim_tiling = sum(math.comb(m_tiling, k) for k in range(d+1))

        print(f"\n  n={n}: degree d={d}")
        print(f"    Arc variables: m={m}, RM({d},{m}) dim = {rm_dim}")
        print(f"    Tiling variables: m_tiling={m_tiling}, RM({d},{m_tiling}) dim = {rm_dim_tiling}")
        print(f"    Total tournaments: 2^{m} = {2**m}")
        print(f"    H values from {2**m} tournaments: at most {2**m} (in practice much fewer)")

        # H as a vector in F_2^{2^m}... actually H is integer-valued, not binary
        # But we can consider H mod 2, H mod 4, etc.
        # H mod 2 = 1 always (Redei). So H mod 2 is the constant function 1.
        # H mod 4 = 1 + 2*alpha_1 mod 4. Since alpha_1 can be even or odd, this is 1 or 3.
        # H mod 4 as a Boolean function has degree...

    # ========================================
    # PART 4: Tiling Weights and Strips
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: STRIP STRUCTURE AND TILING CODES")
    print("  Strip Str(k) = {(r,c) : r+c=k, r,c >= 1} has k-1 tiles")
    print("  Strips are the 'rows' of the triangular grid (anti-diagonals)")
    print(f"{'='*70}")

    for n in [5, 7]:
        positions = pin_grid_positions(n)
        m = len(positions)
        print(f"\n  n={n}:")

        # Organize by strip
        strips = defaultdict(list)
        for r, c in positions:
            strips[r+c].append((r, c))

        for k in sorted(strips.keys()):
            arcs = [position_to_arc(r, c) for r, c in strips[k]]
            print(f"    Strip {k} ({len(strips[k])} tiles): {arcs}")

        # For each tiling, compute the strip weight profile
        # (number of 1-bits in each strip)
        if n <= 5:
            strip_weight_profiles = Counter()
            for bits in range(2**m):
                profile = []
                for k in sorted(strips.keys()):
                    w = sum(1 for r, c in strips[k] if bits & (1 << positions.index((r, c))))
                    profile.append(w)
                strip_weight_profiles[tuple(profile)] += 1

                # Only count a few
                if bits > 100:
                    break

            # For GS tilings specifically
            print(f"\n  Strip weight profiles for GS tilings:")
            pos_idx = {p: i for i, p in enumerate(positions)}
            gs_fixed = []
            gs_pairs = []
            for r, c in positions:
                r2, c2 = gs_map(r, c, n)
                if (r, c) == (r2, c2):
                    gs_fixed.append((r, c))
                elif pos_idx[(r, c)] < pos_idx[(r2, c2)]:
                    gs_pairs.append(((r, c), (r2, c2)))

            gs_dim = len(gs_fixed) + len(gs_pairs)

            strip_H = defaultdict(list)
            for gs_bits in range(2**gs_dim):
                full_bits = 0
                bit_idx = 0
                for p in gs_fixed:
                    if gs_bits & (1 << bit_idx):
                        full_bits |= (1 << pos_idx[p])
                    bit_idx += 1
                for p1, p2 in gs_pairs:
                    if gs_bits & (1 << bit_idx):
                        full_bits |= (1 << pos_idx[p1])
                        full_bits |= (1 << pos_idx[p2])
                    bit_idx += 1

                A = bits_to_adj(full_bits, n)
                H = compute_H_dp(A, n)

                # Strip weights
                profile = []
                for k in sorted(strips.keys()):
                    w = sum(1 for r, c in strips[k] if full_bits & (1 << pos_idx[(r, c)]))
                    profile.append(w)

                strip_H[tuple(profile)].append(H)

            for profile in sorted(strip_H.keys()):
                H_vals = strip_H[profile]
                print(f"    strips={list(profile)}: H values = {sorted(set(H_vals))}")

    # ========================================
    # PART 5: Blue Line Skeleton as Tanner Graph
    # ========================================
    print(f"\n{'='*70}")
    print("PART 5: BLUE LINE SKELETON AND TANNER GRAPHS")
    print("  In coding theory, the Tanner graph connects variable nodes to")
    print("  check nodes. Does the skeleton have this structure?")
    print(f"{'='*70}")

    n = 5
    positions = pin_grid_positions(n)
    m = len(positions)
    pos_idx = {p: i for i, p in enumerate(positions)}

    # GS tilings
    gs_fixed = []
    gs_pairs = []
    for r, c in positions:
        r2, c2 = gs_map(r, c, n)
        if (r, c) == (r2, c2):
            gs_fixed.append((r, c))
        elif pos_idx[(r, c)] < pos_idx[(r2, c2)]:
            gs_pairs.append(((r, c), (r2, c2)))

    gs_dim = len(gs_fixed) + len(gs_pairs)

    # Build the skeleton: SC classes connected by GS flips
    # First compute all GS tilings with their H and isomorphism class
    gs_data = []
    for gs_bits in range(2**gs_dim):
        full_bits = 0
        bit_idx = 0
        for p in gs_fixed:
            if gs_bits & (1 << bit_idx):
                full_bits |= (1 << pos_idx[p])
            bit_idx += 1
        for p1, p2 in gs_pairs:
            if gs_bits & (1 << bit_idx):
                full_bits |= (1 << pos_idx[p1])
                full_bits |= (1 << pos_idx[p2])
            bit_idx += 1

        A = bits_to_adj(full_bits, n)
        H = compute_H_dp(A, n)

        # Simple isomorphism invariant: sorted score sequence + c3
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        c3 = int(np.trace(A @ A @ A)) // 3

        gs_data.append({
            'gs_bits': gs_bits,
            'full_bits': full_bits,
            'H': H,
            'scores': scores,
            'c3': c3,
        })

    # Group by (scores, c3, H) as proxy for isomorphism class
    class_map = defaultdict(list)
    for d in gs_data:
        key = (d['scores'], d['c3'], d['H'])
        class_map[key].append(d)

    print(f"\n  n={n}: {len(gs_data)} GS tilings -> {len(class_map)} classes")
    for key, members in sorted(class_map.items()):
        scores, c3, H = key
        print(f"    {scores}, c3={c3}, H={H}: {len(members)} GS tilings")

    # GS flip: bit flip corresponds to complementing
    # The skeleton connects class A to class B if some GS tiling in A
    # is the complement (flip all non-backbone arcs) of some GS tiling in B
    flip_mask = (1 << m) - 1
    skeleton_edges = set()
    for d in gs_data:
        t_flip = d['full_bits'] ^ flip_mask
        # Find which class t_flip belongs to
        A_flip = bits_to_adj(t_flip, n)
        H_flip = compute_H_dp(A_flip, n)
        scores_flip = tuple(sorted([sum(A_flip[i]) for i in range(n)]))
        c3_flip = int(np.trace(A_flip @ A_flip @ A_flip)) // 3
        key_orig = (d['scores'], d['c3'], d['H'])
        key_flip = (scores_flip, c3_flip, H_flip)

        if key_orig != key_flip:
            edge = tuple(sorted([key_orig, key_flip]))
            skeleton_edges.add(edge)

    print(f"\n  Skeleton edges (GS flip between different classes):")
    for e in sorted(skeleton_edges):
        print(f"    {e[0]} <-> {e[1]}")

    # Is it bipartite? Check t3 parity
    print(f"\n  t3 parity bipartition check:")
    for key in sorted(class_map.keys()):
        _, c3, H = key
        side = "EVEN" if c3 % 2 == 0 else "ODD"
        print(f"    {key}: t3={c3} -> side {side}")

    # Verify: all edges connect EVEN to ODD?
    bipartite = all(
        (e[0][1] % 2) != (e[1][1] % 2)
        for e in skeleton_edges
    )
    print(f"\n  Skeleton bipartite (t3 parity)? {bipartite}")

    # ========================================
    # PART 6: Hamming Weight on Pin Grid = "Backward Arc Count"
    # ========================================
    print(f"\n{'='*70}")
    print("PART 6: HAMMING WEIGHT = BACKWARD ARC COUNT")
    print("  The Hamming weight of a tiling = # backward (blue) arcs")
    print("  = # arcs pointing 'against' the base path direction")
    print("  This is related to the 'writhe' from knot theory (S69)")
    print(f"{'='*70}")

    for n in [5]:
        positions = pin_grid_positions(n)
        m = len(positions)

        # For each tiling, Hamming weight = # backward arcs
        hw_H = defaultdict(list)
        for bits in range(2**m):
            hw = bin(bits).count('1')
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            hw_H[hw].append(H)

        print(f"\n  n={n}: Hamming weight (backward arcs) vs H:")
        for hw in sorted(hw_H.keys()):
            vals = hw_H[hw]
            mean_H = np.mean(vals)
            H_set = sorted(set(vals))
            print(f"    hw={hw:2d}: mean H={mean_H:6.2f}, H range=[{min(vals)},{max(vals)}], "
                  f"count={len(vals)}")

        # The "perpendicular" structure: H vs Hamming weight
        # From earlier sessions: this is a symmetric bell curve at n=5
        print(f"\n  Mean H by Hamming weight (should be symmetric):")
        hw_means = {hw: np.mean(vals) for hw, vals in hw_H.items()}
        for hw in sorted(hw_means.keys()):
            bar = '#' * int(hw_means[hw])
            print(f"    hw={hw:2d}: mean H = {hw_means[hw]:6.2f} {bar}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
