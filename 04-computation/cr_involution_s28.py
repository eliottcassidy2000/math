#!/usr/bin/env python3
"""
The CR involution: complement∘reflect as a tournament operation.
opus-2026-04-03-S28

DISCOVERY: complement∘reflect is an H-preserving involution on tilings.
What tournament operation does it correspond to?

Hypotheses:
1. It might be vertex relabeling (some permutation sigma of [n])
2. It might be "base-path reversal" (reverse the base path, keep other arcs)
3. It might be something entirely new

The tiling complement flips all tile arcs but keeps base-path arcs.
The grid reflection remaps tile (x,y) -> (n+1-y, n+1-x) with bit flip.

Composition: for each tile position, the CR involution:
  - Takes tile (x,y) with bit b
  - Maps it to tile (n+1-y, n+1-x) with bit (1 - (1 - b)) = b
  Wait — reflect flips bits AND remaps. Complement flips bits.
  So complement∘reflect: remap coordinates (x,y) -> (n+1-y, n+1-x) and keep bit.

  That's PURE COORDINATE REMAP without bit change!

  In tournament terms: the arc between vertices x and y (x-y >= 2) gets
  mapped to the arc between vertices n+1-y and n+1-x.
  The direction is PRESERVED (same bit value).

  But what about base-path arcs? The base path is unchanged in the tiling model
  (both complement and reflect preserve it... or do they?).

  Let me verify this algebraically and check what sigma it corresponds to.
"""

from collections import defaultdict
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

def apply_cr(b, n, tiles, tile_to_idx):
    """Apply complement∘reflect to a tiling."""
    m = len(tiles)
    # Complement: flip all bits
    bc = [1 - x for x in b]
    # Then reflect: (x,y) -> (n+1-y, n+1-x) with bit flip
    result = [0] * m
    for i, (x, y) in enumerate(tiles):
        rx, ry = n + 1 - y, n + 1 - x
        j = tile_to_idx[(rx, ry)]
        result[j] = 1 - bc[i]  # reflect flips the bit
    return result

def apply_pure_remap(b, n, tiles, tile_to_idx):
    """Apply pure coordinate remap (x,y) -> (n+1-y, n+1-x) WITHOUT bit flip."""
    m = len(tiles)
    result = [0] * m
    for i, (x, y) in enumerate(tiles):
        rx, ry = n + 1 - y, n + 1 - x
        j = tile_to_idx[(rx, ry)]
        result[j] = b[i]  # NO bit flip
    return result

for n in range(3, 8):
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    print(f"\n{'='*70}")
    print(f"CR INVOLUTION ANATOMY: n={n}, m={m}")
    print(f"{'='*70}")

    # First verify that CR = pure remap
    if 2**m <= 40000:
        cr_eq_remap = 0
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            cr = apply_cr(b, n, tiles, tile_to_idx)
            remap = apply_pure_remap(b, n, tiles, tile_to_idx)
            if cr == remap:
                cr_eq_remap += 1

        print(f"CR = pure remap? {cr_eq_remap}/{2**m} = {cr_eq_remap/2**m:.3f}")

    # What vertex permutation does the coordinate remap correspond to?
    # Tile (x,y) maps to tile (n+1-y, n+1-x).
    # In the tiling model, tile (x,y) represents the arc between vertices x and y.
    # After remap: it represents the arc between vertices n+1-y and n+1-x.
    #
    # So vertex v maps to vertex n+1-v. That's the "reversal" permutation!
    # sigma(v) = n+1-v
    #
    # Let's verify: does the CR tournament equal the tournament obtained by
    # applying the permutation sigma(v) = n+1-v to all vertices?

    print(f"\nVertex permutation sigma(v) = {n+1}-v:")
    for v in range(1, n+1):
        print(f"  {v} -> {n+1-v}")

    # Apply sigma to a tournament: T_sigma(sigma(u), sigma(v)) = T(u, v)
    # In matrix terms: T_sigma[n-v][n-u] = T[u-1][v-1] (0-indexed)
    # Or more carefully: T_sigma[sigma(u)-1][sigma(v)-1] = T[u-1][v-1]
    # where sigma(u) = n+1-u, so sigma(u)-1 = n-u

    if 2**m <= 40000:
        sigma_match = 0
        for bits in range(min(2**m, 10000)):
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)

            # CR tournament
            cr = apply_cr(b, n, tiles, tile_to_idx)
            Tcr = tournament_from_bits(n, cr)

            # Sigma-permuted tournament
            Ts = [[0]*n for _ in range(n)]
            for u in range(n):
                for v in range(n):
                    if u != v:
                        su = n - 1 - u  # sigma(u+1)-1 = n-(u+1) = n-1-u
                        sv = n - 1 - v
                        Ts[su][sv] = T[u][v]

            if all(Tcr[i][j] == Ts[i][j] for i in range(n) for j in range(n)):
                sigma_match += 1

        total_checked = min(2**m, 10000)
        print(f"\nCR tournament = sigma-permuted tournament? {sigma_match}/{total_checked} = {sigma_match/total_checked:.3f}")

    # WAIT: the remap only touches TILE arcs. Base-path arcs are a separate thing.
    # Under the vertex permutation sigma(v) = n+1-v:
    #   Base path n -> n-1 -> ... -> 1 maps to 1 -> 2 -> ... -> n
    #   Which is the REVERSE base path!
    # But in the tiling model, the base path is FIXED.
    # So sigma relabeling can't be equivalent to pure tile remap,
    # because it also changes the base path direction.

    # Let me check directly what base-path arcs look like in CR vs sigma
    if 2**m <= 40000:
        # Take a specific tiling and compare arc by arc
        bits = 0  # All-zero tiling = transitive tournament
        b = [(bits >> i) & 1 for i in range(m)]
        T_orig = tournament_from_bits(n, b)

        cr_bits = apply_cr(b, n, tiles, tile_to_idx)
        T_cr = tournament_from_bits(n, cr_bits)

        # sigma-permuted
        Ts = [[0]*n for _ in range(n)]
        for u in range(n):
            for v in range(n):
                if u != v:
                    su = n - 1 - u
                    sv = n - 1 - v
                    Ts[su][sv] = T_orig[u][v]

        print(f"\nTransitive tournament (bits=0):")
        print(f"Base-path arcs in original:  ", end="")
        for i in range(n-1):
            print(f"  T[{i}][{i+1}]={'>' if T_orig[i][i+1] else '<'}{'>' if T_orig[i+1][i] else '<'}", end="")
        print()

        print(f"Base-path arcs in CR:        ", end="")
        for i in range(n-1):
            print(f"  T[{i}][{i+1}]={'>' if T_cr[i][i+1] else '<'}{'>' if T_cr[i+1][i] else '<'}", end="")
        print()

        print(f"Base-path arcs in sigma:     ", end="")
        for i in range(n-1):
            print(f"  T[{i}][{i+1}]={'>' if Ts[i][i+1] else '<'}{'>' if Ts[i+1][i] else '<'}", end="")
        print()

        print(f"\nDifference CR vs sigma (all arcs):")
        diff_count = 0
        for u in range(n):
            for v in range(u+1, n):
                if T_cr[u][v] != Ts[u][v]:
                    diff_count += 1
                    print(f"  Arc ({u+1},{v+1}): CR={'→' if T_cr[u][v] else '←'} vs sigma={'→' if Ts[u][v] else '←'}")
        if diff_count == 0:
            print(f"  NONE — they are identical!")
        else:
            print(f"  Total differences: {diff_count}")

    # Now let's understand what CR does to the TOURNAMENT
    # without reference to vertex relabeling.
    # CR preserves H. What other invariants does it preserve or flip?

    if 2**m <= 2000:
        print(f"\n--- CR INVARIANTS ---")

        # For each tiling, compute: scores, cycle count, etc
        from itertools import combinations

        def scores(T, n):
            return tuple(sorted([sum(T[i]) for i in range(n)]))

        def count_3cycles(T, n):
            count = 0
            for i, j, k in combinations(range(n), 3):
                if T[i][j] and T[j][k] and T[k][i]:
                    count += 1
                if T[i][k] and T[k][j] and T[j][i]:
                    count += 1
            return count

        score_preserved = 0
        c3_preserved = 0

        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            cr = apply_cr(b, n, tiles, tile_to_idx)
            Tcr = tournament_from_bits(n, cr)

            s1 = scores(T, n)
            s2 = scores(Tcr, n)
            if s1 == s2:
                score_preserved += 1

            c1 = count_3cycles(T, n)
            c2 = count_3cycles(Tcr, n)
            if c1 == c2:
                c3_preserved += 1

        print(f"Score sequence preserved? {score_preserved}/{2**m} = {score_preserved/2**m:.3f}")
        print(f"3-cycle count preserved? {c3_preserved}/{2**m} = {c3_preserved/2**m:.3f}")

        # Check: is CR the same as the ISOMORPHISM that identifies
        # T with its complement T^c (which is the grid-reflected tournament)?
        # THM-280 says grid reflection gives the complement tournament T^c (= T^op for tournaments).
        # So: reflect(t) -> T^op and complement(t) -> flip tile arcs only.
        # CR(t) = complement(reflect(t)) = complement of T^op-tiling.
        #
        # Actually let me think about this differently.
        # reflect(t) has H = H(T^op) = H(T).
        # complement(reflect(t)) has some H. We showed this H equals H(T).
        # So CR preserves H, but what IS the resulting tournament?

        # Let me just print a few examples
        print(f"\nExamples of T and CR(T):")
        import random
        random.seed(n)
        for trial in range(min(5, 2**m)):
            bits = random.getrandbits(m) & ((1 << m) - 1) if 2**m > 5 else trial
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H1 = count_hp(T, n)
            s1 = scores(T, n)

            cr = apply_cr(b, n, tiles, tile_to_idx)
            Tcr = tournament_from_bits(n, cr)
            H2 = count_hp(Tcr, n)
            s2 = scores(Tcr, n)

            # Check if T and CR(T) are isomorphic
            # (Same score and same H is necessary but not sufficient)
            same_scores = (s1 == s2)
            same_H = (H1 == H2)

            # Count fixed points
            is_fixed = (b == cr)

            print(f"  bits={bits:0{m}b}: scores={s1} H={H1} | CR scores={s2} H={H2} "
                  f"{'SAME' if same_scores else 'DIFF'} {'FIXED' if is_fixed else ''}")

        # Count fixed points of CR
        fixed_count = 0
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            cr = apply_cr(b, n, tiles, tile_to_idx)
            if b == cr:
                fixed_count += 1

        print(f"\nFixed points of CR: {fixed_count}/{2**m}")
        print(f"Fixed = grid-symmetric tilings: 2^(fixed tiles in reflection)")

        # Count self-reflected tiles (fixed points of the coordinate remap)
        self_reflected = sum(1 for i, (x, y) in enumerate(tiles)
                          if (n + 1 - y, n + 1 - x) == (x, y))
        print(f"Self-reflected tiles: {self_reflected}")
        print(f"Predicted fixed points: 2^{self_reflected} = {2**self_reflected}")
