#!/usr/bin/env python3
"""
Effect symmetry under grid reflection and the n=7 breaking.
opus-2026-04-03-S28

KEY DISCOVERY: At n=6, the tile effects show a beautiful reflection symmetry:
  D(3,1) = J(6,4) = +4.0   (base corners, reflected pair)
  G(4,2) = I(5,3) = 0.0    (center tiles, reflected pair)
  C(4,1) = H(6,3) = +3.0   (reflected pair)
  B(5,1) = E(6,2) = +3.25  (reflected pair)
  F(5,2) = -0.5             (self-reflected? or has a partner?)

This is EXACTLY the grid reflection R: (x,y) -> (n+1-y, n+1-x)!
Under grid reflection, tile (x,y) maps to tile (n+1-y, n+1-x).

THEOREM to verify: effect(tile (x,y)) = effect(tile (n+1-y, n+1-x))

Also: the apex effect = 2^(n-3) BREAKS at n=7. Why?
At n=7: effect = 159/8 = 19.875, not 16.
BUT 159/8 = 19.875 is close to 20. And 20 = C(6,3). Hmm.
Actually 159/8. Let's see what fraction 159/8 represents more cleanly.
159 = 3 * 53. Not obviously nice.

Let's also check: is the effect a RATIONAL number with denominator 2^k?
At n=7: 159/8 = 159/2^3. Denominator = 2^3 = 8.
At n=6: 8 = 8/1 = 2^3/2^0.
At n=5: 4 = 2^2.
Pattern: denominator grows?

The key question: is there a CORRECTION TERM that explains the departure?
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

def grid_reflect(x, y, n):
    """Grid reflection R: (x,y) -> (n+1-y, n+1-x)."""
    return (n + 1 - y, n + 1 - x)

for n in range(3, 8):
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)

    print(f"\n{'='*70}")
    print(f"n={n}, m={m}: EFFECT SYMMETRY UNDER GRID REFLECTION")
    print(f"{'='*70}")

    # Verify grid reflection maps tiles to tiles
    tile_set = set(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    print(f"\nGrid reflection R: (x,y) -> (n+1-y, n+1-x) = ({n+1}-y, {n+1}-x)")
    print(f"{'Tile':>8} {'R(tile)':>10} {'Self-paired?':>14}")
    for i, (x, y) in enumerate(tiles):
        rx, ry = grid_reflect(x, y, n)
        is_self = (rx, ry) == (x, y)
        if (rx, ry) in tile_set:
            j = tile_to_idx[(rx, ry)]
            print(f"  {chr(65+i)}({x},{y}) -> {chr(65+j)}({rx},{ry})  {'SELF' if is_self else ''}")
        else:
            print(f"  {chr(65+i)}({x},{y}) -> ({rx},{ry}) NOT A TILE!")

    # Compute exact effects (exhaustive for small n)
    if 2**m > 100000:
        print(f"  (Skipping exhaustive — too large)")
        continue

    total_H = {j: {0: 0, 1: 0} for j in range(m)}
    count_per = {j: {0: 0, 1: 0} for j in range(m)}

    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H = count_hp(T, n)
        for j in range(m):
            total_H[j][b[j]] += H
            count_per[j][b[j]] += 1

    effects = {}
    for j in range(m):
        e = Fraction(total_H[j][1], count_per[j][1]) - Fraction(total_H[j][0], count_per[j][0])
        effects[j] = e

    # Check reflection symmetry
    print(f"\nEffect values and reflection symmetry:")
    for i, (x, y) in enumerate(tiles):
        rx, ry = grid_reflect(x, y, n)
        if (rx, ry) in tile_to_idx:
            j = tile_to_idx[(rx, ry)]
            match = "=" if effects[i] == effects[j] else "≠"
            print(f"  {chr(65+i)}({x},{y}) effect={float(effects[i]):+.4f} "
                  f" {match}  {chr(65+j)}({rx},{ry}) effect={float(effects[j]):+.4f}")
        else:
            print(f"  {chr(65+i)}({x},{y}) effect={float(effects[i]):+.4f}  (NO PAIR)")

    # Group by reflection orbits
    seen = set()
    print(f"\nReflection orbit effects:")
    for i, (x, y) in enumerate(tiles):
        if i in seen:
            continue
        rx, ry = grid_reflect(x, y, n)
        if (rx, ry) in tile_to_idx:
            j = tile_to_idx[(rx, ry)]
            if i == j:
                print(f"  Fixed point: {chr(65+i)}({x},{y}) effect={float(effects[i]):+.4f}")
                seen.add(i)
            else:
                avg = (effects[i] + effects[j]) / 2
                diff = effects[i] - effects[j]
                print(f"  Pair: {chr(65+i)}({x},{y}) & {chr(65+j)}({rx},{ry}) "
                      f" effects={float(effects[i]):+.4f},{float(effects[j]):+.4f} "
                      f" avg={float(avg):+.4f} diff={float(diff):+.4f}")
                seen.add(i)
                seen.add(j)

    # How much of total H-variance is explained by each tile's effect?
    all_H = []
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H = count_hp(T, n)
        all_H.append(H)

    mean_H = Fraction(sum(all_H), len(all_H))
    var_H = sum((h - mean_H)**2 for h in all_H) / len(all_H)
    print(f"\nVar(H) = {float(var_H):.4f}")

    # The linear model: H_predicted = mean + sum(effect_j * (b_j - 1/2))
    # R² of the linear model
    ss_pred = sum(effects[j]**2 for j in range(m)) / 4  # Each b_j has var 1/4
    r2_linear = float(ss_pred / var_H) if var_H > 0 else 0
    print(f"R² of linear (main effects only) model = {r2_linear:.6f}")
    print(f"=> Pairwise+ interactions explain {1-r2_linear:.6f} of variance")

    # Does the APEX ALONE explain more than any other single tile?
    skips = [x - y for x, y in tiles]
    for j in range(m):
        r2_j = float(effects[j]**2 / (4 * var_H)) if var_H > 0 else 0
        print(f"  {chr(65+j)}({tiles[j][0]},{tiles[j][1]}) skip={skips[j]}: "
              f"effect={float(effects[j]):+.4f} R²={r2_j:.6f}")

print(f"\n{'='*70}")
print("KEY PATTERN: Effects by (skip, position-on-boundary)")
print(f"{'='*70}")

# At each n, classify tiles as:
#   - BOUNDARY: y=1 (bottom row) or x=n (right column)
#   - INTERIOR: all others
# The boundary tiles connect to the base path endpoints

for n in range(3, 8):
    if (n-1)*(n-2)//2 > 17:
        continue
    tiles = get_tiles(n)
    m = len(tiles)
    print(f"\nn={n}:")
    for i, (x, y) in enumerate(tiles):
        is_bottom = (y == 1)
        is_right = (x == n)
        boundary = "BOTTOM" if is_bottom else ("RIGHT" if is_right else "INTERIOR")
        print(f"  {chr(65+i)}({x},{y}) skip={x-y} {boundary}")
