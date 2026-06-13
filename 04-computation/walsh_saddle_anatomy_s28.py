#!/usr/bin/env python3
"""
Anatomy of the H-saddle through Walsh-Fourier lens.
opus-2026-04-03-S28

BREAKTHROUGH: The g(k,d) matrix shows that at d=m (complement):
  g(even_order, m) = 0   — EVEN Walsh orders contribute NOTHING
  g(odd_order, m) = 4    — ODD Walsh orders contribute MAXIMALLY

This means:
  E[delta²] at d=m = 4 * (order-1 energy + order-3 energy + order-5 energy + ...)

But we already know odd Walsh orders have ZERO energy for H!
(Because H is complement-invariant: H(T) = H(T^c))

Wait... that can't be right. H IS complement-invariant, so H(t) = H(comp(t)),
meaning delta = 0 ALWAYS at d=m... except that's NOT what we observed!

CORRECTION: H(T) = H(T^op) (opposite tournament, flip ALL arcs),
but the complement TILING (flip all TILE bits) does NOT give T^op!
The complement tiling preserves base-path arcs but flips all tile arcs.
T^op flips ALL arcs including base-path arcs.

So H(tiling) ≠ H(complement_tiling) in general!

The Walsh-Fourier analysis above used TILING bits, not tournament arcs.
H(T) = H(T^op) is a statement about ALL arcs.
The tiling complement only flips C(n-1,2) arcs but preserves n-1 base-path arcs.

THIS IS THE KEY: The tiling complement is NOT the tournament complement!
There are m = C(n-1,2) tile arcs and n-1 base arcs. Total = C(n,2).
Flipping all tiles = flipping m out of C(n,2) arcs.

Let's investigate this carefully.
"""

from collections import defaultdict
from fractions import Fraction
import sys

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

def tournament_opposite(T, n):
    """T^op: flip ALL arcs."""
    Top = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                Top[i][j] = T[j][i]
    return Top

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

for n in range(3, 8):
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)

    print(f"\n{'='*70}")
    print(f"COMPLEMENT ANATOMY: n={n}, m={m}")
    print(f"Total arcs = C({n},2) = {n*(n-1)//2}")
    print(f"Base-path arcs = {n-1}")
    print(f"Tile arcs = {m}")
    print(f"Fraction of arcs flipped by tiling complement = {m}/{n*(n-1)//2} = {m*100//(n*(n-1)//2)}%")
    print(f"{'='*70}")

    if 2**m > 40000:
        # Sample for n=7
        import random
        random.seed(42)
        sample = 5000
        same_count = 0
        comp_deltas = []
        op_deltas = []

        for _ in range(sample):
            bits = random.getrandbits(m) & ((1 << m) - 1)
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)

            # Tiling complement
            bc = [1 - x for x in b]
            Tc = tournament_from_bits(n, bc)
            Hc = count_hp(Tc, n)
            comp_deltas.append(Hc - H)

            # Tournament opposite
            Top = tournament_opposite(T, n)
            Hop = count_hp(Top, n)
            op_deltas.append(Hop - H)

            if H == Hc:
                same_count += 1

        print(f"\nH(tiling) vs H(complement_tiling): same {same_count}/{sample} = {same_count/sample:.3f}")
        print(f"H(T) vs H(T^op): always equal? ", end="")
        print("YES" if all(d == 0 for d in op_deltas) else "NO")
        print(f"Mean |comp delta| = {sum(abs(d) for d in comp_deltas)/len(comp_deltas):.2f}")
        print(f"Mean |op delta| = {sum(abs(d) for d in op_deltas)/len(op_deltas):.2f}")
    else:
        same_count = 0
        op_always_same = True
        comp_deltas = []

        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)

            # Tiling complement
            bc = [1 - x for x in b]
            Tc = tournament_from_bits(n, bc)
            Hc = count_hp(Tc, n)
            comp_deltas.append(Hc - H)

            if H == Hc:
                same_count += 1

            # Tournament opposite
            Top = tournament_opposite(T, n)
            Hop = count_hp(Top, n)
            if Hop != H:
                op_always_same = False

        print(f"\nH(tiling) = H(complement_tiling)? {same_count}/{2**m} = {same_count/2**m:.3f}")
        print(f"H(T) = H(T^op) always? {op_always_same}")
        if comp_deltas:
            mean_abs_comp = sum(abs(d) for d in comp_deltas) / len(comp_deltas)
            print(f"Mean |tiling complement delta| = {mean_abs_comp:.4f}")

    # THE KEY INSIGHT: What does tiling complement correspond to in tournament space?
    # Base path is n -> n-1 -> ... -> 1 (vertex i beats vertex j for i > j on base path)
    # Complement tiling: for each tile (x,y) with x-y >= 2, flip the arc direction
    # BUT keep base-path arcs (x, x-1) unchanged
    #
    # So the complement tiling is the tournament where:
    #   - Arc (i, i-1) for consecutive vertices: SAME as original
    #   - Arc (x, y) for non-consecutive: REVERSED
    #
    # This is like a "partial complement" — flipping only the non-consecutive arcs

    # Let's see what the GRID REFLECTION (THM-280) does
    # Grid reflection R maps tiling t to the tiling of the complement tournament T^op
    # But wait — THM-280 says grid reflection = tournament complement for ALL n
    # So grid reflection of tiling bits ≠ complement of tiling bits

    # Let's verify: grid-reflect each tiling and check if it gives T^op
    print(f"\n--- GRID REFLECTION vs TOURNAMENT COMPLEMENT ---")

    if 2**m <= 2000:
        # For each tiling, compute:
        # 1. The grid-reflected tiling
        # 2. The tournament of the grid-reflected tiling
        # 3. Check if it equals T^op

        tile_to_idx = {t: i for i, t in enumerate(tiles)}

        def grid_reflect_tiling(bits_list, n, tiles, tile_to_idx):
            """Grid-reflect a tiling: (x,y) -> (n+1-y, n+1-x), flip bit."""
            reflected = [0] * len(tiles)
            for i, (x, y) in enumerate(tiles):
                rx, ry = n + 1 - y, n + 1 - x
                if (rx, ry) in tile_to_idx:
                    j = tile_to_idx[(rx, ry)]
                    reflected[j] = 1 - bits_list[i]  # flip + remap
                else:
                    # This shouldn't happen if grid reflection preserves tiles
                    pass
            return reflected

        match_op = 0
        match_comp = 0
        match_neither = 0

        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)

            # Grid-reflected tiling
            br = grid_reflect_tiling(b, n, tiles, tile_to_idx)
            Tr = tournament_from_bits(n, br)
            Hr = count_hp(Tr, n)

            # Tournament opposite
            Top = tournament_opposite(T, n)
            Hop = count_hp(Top, n)

            # Tiling complement
            bc = [1 - x for x in b]
            Tc = tournament_from_bits(n, bc)
            Hc = count_hp(Tc, n)

            # Check if grid-reflected tournament = T^op
            is_op = all(Tr[i][j] == Top[i][j] for i in range(n) for j in range(n))
            # Check if grid-reflected tiling = complement tiling
            is_comp = (br == bc)

            if is_op:
                match_op += 1
            if is_comp:
                match_comp += 1
            if not is_op and not is_comp:
                match_neither += 1

        print(f"Grid-reflected tournament = T^op? {match_op}/{2**m} = {match_op/2**m:.3f}")
        print(f"Grid-reflected tiling = complement tiling? {match_comp}/{2**m} = {match_comp/2**m:.3f}")
        if match_neither > 0:
            print(f"Neither: {match_neither}/{2**m}")
    else:
        print(f"  (Skipping exhaustive for n={n})")

    # THE ACTUAL DUALITY: What IS the relationship between:
    # (a) complement_tiling(t) — flip all tile bits
    # (b) grid_reflect(t) — which gives T^op
    # (c) The composition complement_tiling ∘ grid_reflect

    # If grid_reflect gives T^op, and T^op has H(T^op) = H(T),
    # then grid_reflect preserves H.
    # But complement_tiling does NOT generally preserve H.
    # So: complement_tiling = grid_reflect ∘ (some other operation)
    # That "other operation" is the one that changes H.

    if 2**m <= 2000:
        print(f"\n--- THE THREE OPERATIONS ---")
        # For a sample tiling, trace all three operations
        import random
        random.seed(n)

        for trial in range(min(5, 2**m)):
            bits = random.getrandbits(m) & ((1 << m) - 1) if 2**m > 5 else trial
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)

            bc = [1 - x for x in b]
            Hc = count_hp(tournament_from_bits(n, bc), n)

            br = grid_reflect_tiling(b, n, tiles, tile_to_idx)
            Hr = count_hp(tournament_from_bits(n, br), n)

            # Composition: complement then reflect
            bcr = grid_reflect_tiling(bc, n, tiles, tile_to_idx)
            Hcr = count_hp(tournament_from_bits(n, bcr), n)

            # Composition: reflect then complement
            brc = [1 - x for x in br]
            Hrc = count_hp(tournament_from_bits(n, brc), n)

            hw = sum(b)
            print(f"  t={bits:0{m}b} (w={hw}): H={H}  comp={Hc}  refl={Hr}  comp∘refl={Hcr}  refl∘comp={Hrc}")

        # Check if reflect∘complement = complement∘reflect (do they commute?)
        commute_count = 0
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            bc = [1 - x for x in b]
            br = grid_reflect_tiling(b, n, tiles, tile_to_idx)
            bcr = grid_reflect_tiling(bc, n, tiles, tile_to_idx)
            brc = [1 - x for x in br]
            if bcr == brc:
                commute_count += 1

        print(f"\nComplement and grid-reflect COMMUTE? {commute_count}/{2**m} = {commute_count/2**m:.3f}")

        # The composition complement∘reflect: what IS it?
        # It flips all tiles AND remaps coordinates.
        # In tournament terms: it's "reflect the staircase then invert all tiles"
        # = "rotate the staircase 180° (with bit-flip from reflection) then invert all bits again"
        # = "rotate 180° with NO bit flip" = pure geometric rotation!

        # Check: is complement∘reflect an involution?
        invol_count = 0
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            bc = [1 - x for x in b]
            bcr = grid_reflect_tiling(bc, n, tiles, tile_to_idx)
            # Apply again
            bcrc = [1 - x for x in bcr]
            bcrcr = grid_reflect_tiling(bcrc, n, tiles, tile_to_idx)
            if bcrcr == b:
                invol_count += 1

        print(f"Complement∘reflect is involution? {invol_count}/{2**m} = {invol_count/2**m:.3f}")

        # Check: what is complement∘reflect in terms of base tournament operations?
        # Is it vertex relabeling? Score reversal?
        print(f"\nComplement∘reflect: H preservation?")
        cr_preserves = 0
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)
            bc = [1 - x for x in b]
            bcr = grid_reflect_tiling(bc, n, tiles, tile_to_idx)
            Tcr = tournament_from_bits(n, bcr)
            Hcr = count_hp(Tcr, n)
            if H == Hcr:
                cr_preserves += 1

        print(f"  H(t) = H(comp∘refl(t))? {cr_preserves}/{2**m} = {cr_preserves/2**m:.3f}")

for n in range(3, 8):
    m = (n-1)*(n-2)//2
    print(f"\nn={n}: m={m}, base-path arcs={n-1}, total arcs={n*(n-1)//2}")
    print(f"  Tiles flipped by complement: {m}/{n*(n-1)//2} = {m*2/(n*(n-1)):.4f}")
    print(f"  This is {m}/{n*(n-1)//2} = (n-1)(n-2)/2 / n(n-1)/2 = (n-2)/n")
    print(f"  (n-2)/n = {(n-2)/n:.4f}")
