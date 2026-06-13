#!/usr/bin/env python3
"""
Deep dive into the H-saddle structure and wiggly/complement duality.
opus-2026-04-03-S28

KEY OBSERVATIONS from first exploration:
1. The skip-max tile (A at n=5, A at n=6) has the HIGHEST |mean delta| —
   it controls H most strongly. This is the APEX tile.
2. Distance layers show PALINDROMIC R² values at n=6:
   dist 4,5,6,7,8,9,10 -> R² 0.038, 0.021, 0.013, 0.068, 0.013, 0.021, 0.038
   The middle layer (dist=7, which contains the apex tile) explains the most!
3. Skip layers: skip-max (apex) always explains the most H-variance.
4. Complement deltas are SYMMETRIC around 0 (as expected from H(T)=H(T^c)).
5. Wiggly-complement correlation is NEAR ZERO — these are truly orthogonal operations!

NEW QUESTIONS:
- Is the palindromic R² a general pattern?
- What is the JOINT distribution of (hyp_bits, apex_bit, H)?
- How does the saddle relate to the Walsh-Fourier spectrum?
- Is there a "saddle coordinate" that explains more variance than any single layer?
"""

from itertools import product
from collections import defaultdict
import sys

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    for i in range(n, 1, -1):
        T[i-1][i-2] = 1

    for idx, (x, y) in enumerate(tiles):
        if idx < len(bits):
            if bits[idx] == 0:
                T[x-1][y-1] = 1
            else:
                T[y-1][x-1] = 1

    return T, tiles

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

def analyze_saddle_structure(n):
    """Deep analysis of the saddle structure."""
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)
    assert len(tiles) == m

    print(f"\n{'='*70}")
    print(f"DEEP SADDLE ANALYSIS: n={n}, m={m}")
    print(f"{'='*70}")

    # Compute tile properties
    tile_skip = [x - y for x, y in tiles]
    tile_dist = [x + y for x, y in tiles]
    max_skip = max(tile_skip)
    max_dist = max(tile_dist)
    min_dist = min(tile_dist)

    apex_idx = [i for i in range(m) if tile_skip[i] == max_skip]
    hyp_idx = [i for i in range(m) if tile_dist[i] == max_dist]
    base_idx = [i for i in range(m) if tile_dist[i] == min_dist]

    print(f"Apex tile(s) (max skip={max_skip}): {[(tiles[i], chr(65+i)) for i in apex_idx]}")
    print(f"Hypotenuse tile(s) (max dist={max_dist}): {[(tiles[i], chr(65+i)) for i in hyp_idx]}")
    print(f"Base tile(s) (min dist={min_dist}): {[(tiles[i], chr(65+i)) for i in base_idx]}")

    # Enumerate
    data = []
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T, _ = tournament_from_bits(n, b)
        H = count_hp(T, n)
        data.append((H, b))

    # 1. APEX vs HYPOTENUSE as saddle coordinates
    print(f"\n--- APEX BIT vs HYPOTENUSE WEIGHT as saddle coordinates ---")

    # For single apex tile (when unique)
    if len(apex_idx) == 1:
        ai = apex_idx[0]
        apex_vs_rest = defaultdict(list)
        for H, b in data:
            apex_bit = b[ai]
            rest_weight = sum(b[i] for i in range(m) if i != ai)
            apex_vs_rest[(apex_bit, rest_weight)].append(H)

        print(f"\nMean H by (apex_bit, rest_weight):")
        max_rest = m - 1
        print(f"{'':>12}", end='')
        for rw in range(max_rest + 1):
            print(f"  rw={rw:>2}", end='')
        print()
        for ab in [0, 1]:
            print(f"apex={ab:>2}:", end='')
            for rw in range(max_rest + 1):
                vals = apex_vs_rest.get((ab, rw), [])
                if vals:
                    print(f"  {sum(vals)/len(vals):>6.1f}", end='')
                else:
                    print(f"  {'---':>6}", end='')
            print()

    # 2. PRODUCT COORDINATE: apex_bit * (m - total_weight)
    # Idea: the saddle should be captured by something like apex*(m-w) - (1-apex)*w
    print(f"\n--- SADDLE COORDINATE SEARCH ---")

    # Try various linear combinations of tile bits as H predictors
    overall_H = [H for H, b in data]
    mean_H = sum(overall_H) / len(overall_H)
    var_H = sum((h - mean_H)**2 for h in overall_H) / len(overall_H)

    # For each pair of tiles, compute the R² of the interaction term b_i * b_j
    if m <= 12:
        print(f"\nPairwise interaction R² (top 10):")
        pair_r2 = []
        for i in range(m):
            for j in range(i+1, m):
                # Compute mean H for each (b_i, b_j) combination
                groups = defaultdict(list)
                for H, b in data:
                    groups[(b[i], b[j])].append(H)
                between = sum(len(g) * (sum(g)/len(g) - mean_H)**2 for g in groups.values()) / len(data)
                r2 = between / var_H if var_H > 0 else 0
                pair_r2.append((r2, i, j))

        pair_r2.sort(reverse=True)
        for r2, i, j in pair_r2[:10]:
            ti, tj = tiles[i], tiles[j]
            si, sj = tile_skip[i], tile_skip[j]
            print(f"  {chr(65+i)}{chr(65+j)} = ({ti[0]},{ti[1]})×({tj[0]},{tj[1]}) "
                  f"skip=({si},{sj}) dist=({tile_dist[i]},{tile_dist[j]}): R²={r2:.4f}")

    # 3. Walsh-like interaction: sum of pairwise products
    print(f"\n--- PARITY AND INTERACTION EFFECTS ---")

    # Total parity = sum(b) mod 2
    even_H = [H for H, b in data if sum(b) % 2 == 0]
    odd_H = [H for H, b in data if sum(b) % 2 == 1]
    print(f"Mean H (even parity): {sum(even_H)/len(even_H):.4f} ({len(even_H)} tilings)")
    print(f"Mean H (odd parity):  {sum(odd_H)/len(odd_H):.4f} ({len(odd_H)} tilings)")

    # Apex parity
    if len(apex_idx) == 1:
        ai = apex_idx[0]
        a0_H = [H for H, b in data if b[ai] == 0]
        a1_H = [H for H, b in data if b[ai] == 1]
        print(f"Mean H (apex=0): {sum(a0_H)/len(a0_H):.4f}")
        print(f"Mean H (apex=1): {sum(a1_H)/len(a1_H):.4f}")
        print(f"Apex effect: {sum(a1_H)/len(a1_H) - sum(a0_H)/len(a0_H):.4f}")

    # 4. The COMPLEMENT as a REFLECTION
    print(f"\n--- COMPLEMENT REFLECTION ANALYSIS ---")

    # For each tiling, pair it with its complement
    # H(t) vs H(complement(t))
    complement_pairs = []
    for bits in range(2**(m-1)):  # Only need half (avoid double-counting)
        b = [(bits >> i) & 1 for i in range(m)]
        bc = [1 - x for x in b]
        bits_c = sum(bc[i] << i for i in range(m))
        if bits < bits_c:  # Avoid self-complement and double-counting
            T1, _ = tournament_from_bits(n, b)
            H1 = count_hp(T1, n)
            T2, _ = tournament_from_bits(n, bc)
            H2 = count_hp(T2, n)
            complement_pairs.append((H1, H2, sum(b)))

    print(f"Number of complement pairs: {len(complement_pairs)}")
    if complement_pairs:
        diffs = [abs(h1 - h2) for h1, h2, _ in complement_pairs]
        print(f"Mean |H(t) - H(comp(t))| = {sum(diffs)/len(diffs):.4f}")
        print(f"Max |H(t) - H(comp(t))| = {max(diffs)}")

        # Does the weight of the tiling predict the complement gap?
        print(f"\nComplement gap by Hamming weight:")
        gap_by_weight = defaultdict(list)
        for h1, h2, w in complement_pairs:
            gap_by_weight[w].append(abs(h1 - h2))
        for w in sorted(gap_by_weight):
            vals = gap_by_weight[w]
            print(f"  w={w}: mean |gap| = {sum(vals)/len(vals):.2f}, max = {max(vals)}, n = {len(vals)}")

    # 5. WIGGLY vs COMPLEMENT: detailed comparison at tile level
    print(f"\n--- TILE-LEVEL WIGGLY vs COMPLEMENT EFFECTS ---")

    # For each tile j, compute:
    #   - The wiggly delta when that tile alone is flipped (d=1)
    #   - The complement delta "contribution" — how much does tile j contribute
    #     to the total complement delta?
    # Decompose complement delta = sum of marginal contributions?

    # Actually, let's look at it differently:
    # complement_delta = H(complement) - H(original)
    # wiggly_delta_j = H(flip j) - H(original)
    # Is complement_delta ≈ sum of wiggly_delta_j? (No, because interactions!)
    # The INTERACTION = complement_delta - sum(wiggly_delta_j) measures higher-order effects

    print(f"\nDecomposing complement delta into wiggly sum + interaction:")
    interactions = []
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T0, _ = tournament_from_bits(n, b)
        H0 = count_hp(T0, n)

        # Complement
        bc = [1 - x for x in b]
        Tc, _ = tournament_from_bits(n, bc)
        Hc = count_hp(Tc, n)
        comp_delta = Hc - H0

        # Sum of individual wiggly deltas
        wig_sum = 0
        for j in range(m):
            bj = b.copy()
            bj[j] = 1 - bj[j]
            Tj, _ = tournament_from_bits(n, bj)
            Hj = count_hp(Tj, n)
            wig_sum += (Hj - H0)

        interaction = comp_delta - wig_sum
        interactions.append(interaction)

    mean_int = sum(interactions) / len(interactions)
    mean_abs_int = sum(abs(x) for x in interactions) / len(interactions)
    var_int = sum((x - mean_int)**2 for x in interactions) / len(interactions)
    print(f"  Mean interaction = {mean_int:.4f}")
    print(f"  Mean |interaction| = {mean_abs_int:.4f}")
    print(f"  Var(interaction) = {var_int:.4f}")
    print(f"  Var(comp_delta) = {sum((Hc - H0)**2 for H0, Hc in [(count_hp(tournament_from_bits(n, [(bits>>i)&1 for i in range(m)])[0], n), count_hp(tournament_from_bits(n, [1-((bits>>i)&1) for i in range(m)])[0], n)) for bits in range(min(2**m, 100))]) / min(2**m, 100):.4f} (sampled)")

    # Distribution of interactions
    int_hist = defaultdict(int)
    for x in interactions:
        int_hist[x] += 1
    print(f"  Interaction distribution: {dict(sorted(int_hist.items()))}")

    # 6. The GRADIENT of H — which direction increases H fastest?
    print(f"\n--- H GRADIENT ANALYSIS ---")

    # For each tiling, find the tile flip that increases H the most
    best_flip_hist = defaultdict(int)
    worst_flip_hist = defaultdict(int)
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T0, _ = tournament_from_bits(n, b)
        H0 = count_hp(T0, n)

        deltas = []
        for j in range(m):
            bj = b.copy()
            bj[j] = 1 - bj[j]
            Tj, _ = tournament_from_bits(n, bj)
            Hj = count_hp(Tj, n)
            deltas.append(Hj - H0)

        best_j = max(range(m), key=lambda j: deltas[j])
        worst_j = min(range(m), key=lambda j: deltas[j])
        best_flip_hist[best_j] += 1
        worst_flip_hist[worst_j] += 1

    print(f"Which tile is MOST OFTEN the best flip (steepest ascent)?")
    for j in sorted(best_flip_hist, key=best_flip_hist.get, reverse=True):
        t = tiles[j]
        pct = best_flip_hist[j] / 2**m * 100
        print(f"  {chr(65+j)} ({t[0]},{t[1]}) skip={tile_skip[j]}: {best_flip_hist[j]} times ({pct:.1f}%)")

    print(f"Which tile is MOST OFTEN the worst flip (steepest descent)?")
    for j in sorted(worst_flip_hist, key=worst_flip_hist.get, reverse=True):
        t = tiles[j]
        pct = worst_flip_hist[j] / 2**m * 100
        print(f"  {chr(65+j)} ({t[0]},{t[1]}) skip={tile_skip[j]}: {worst_flip_hist[j]} times ({pct:.1f}%)")

    return data

# Run for n=3,4,5,6
for n in [3, 4, 5]:
    analyze_saddle_structure(n)

# n=6 is larger, run it too
analyze_saddle_structure(6)
