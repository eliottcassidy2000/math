#!/usr/bin/env python3
"""
Tensor product structure of the H-polynomial.
opus-2026-04-03-S28

INSIGHT: The multilinear polynomial H(x₁,...,xₘ) has a natural
tensor factorization:
  - The m tile bits live on {0,1}^m
  - But the EFFECTIVE dimension of H is much lower
  - H decomposes along the "skip layers" of the staircase

The skip layers partition tiles by x-y value:
  skip=2: tiles (3,1), (4,2), (5,3), (6,4), ...
  skip=3: tiles (4,1), (5,2), (6,3), ...
  skip=4: tiles (5,1), (6,2), ...
  ...

QUESTION: Does H factor as a PRODUCT of functions on these layers?
  H = f_2(skip-2 bits) × f_3(skip-3 bits) × ... × f_{n-1}(apex bit)

If so, the multilinear polynomial would have EXPONENTIALLY fewer
effective terms than the full product.

Let me check: does Var(H | skip-2 bits, skip-3 bits, ...) decompose
into a sum of per-layer variances (independence) or not?

Also: the SCORE is a linear function of the tile bits.
  score(v) = (incoming arcs to v) = ... involves specific tiles.
  How does the score relate to the skip decomposition?
"""

from collections import defaultdict
from math import comb

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

for n in [5, 6]:
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)

    print(f"\n{'='*70}")
    print(f"TENSOR STRUCTURE ANALYSIS: n={n}, m={m}")
    print(f"{'='*70}")

    # Skip layers
    skip_layers = defaultdict(list)
    for i, (x, y) in enumerate(tiles):
        skip_layers[x-y].append(i)

    print(f"\nSkip layers:")
    for s in sorted(skip_layers):
        idxs = skip_layers[s]
        tile_labels = [chr(65+i) for i in idxs]
        print(f"  skip={s}: {tile_labels} ({len(idxs)} tiles)")

    # Compute all H values
    H_vals = {}
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H_vals[bits] = count_hp(T, n)

    mean_H = sum(H_vals.values()) / len(H_vals)
    var_H = sum((h - mean_H)**2 for h in H_vals.values()) / len(H_vals)

    # Test: does H factor across skip layers?
    # If H = f(skip-2 bits) * g(other bits), then
    # H / E[H | skip-2 bits] should be independent of skip-2 bits.

    # More precisely, compute the ANOVA decomposition across skip layers
    # R² for each layer individually
    print(f"\nIndividual layer R² (variance explained):")
    total_r2 = 0
    for s in sorted(skip_layers):
        idxs = skip_layers[s]
        layer_size = len(idxs)

        # Group by this layer's bit pattern
        group_H = defaultdict(list)
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            layer_pattern = tuple(b[i] for i in idxs)
            group_H[layer_pattern].append(H_vals[bits])

        between_var = sum(len(g) * (sum(g)/len(g) - mean_H)**2
                        for g in group_H.values()) / len(H_vals)
        r2 = between_var / var_H if var_H > 0 else 0
        total_r2 += r2
        print(f"  skip={s} ({layer_size} tiles): R²={r2:.6f}")

    print(f"  Sum of individual R²: {total_r2:.6f}")
    print(f"  (If layers are independent, this would equal the full model R²)")

    # Now compute R² for the FULL model (all tiles jointly)
    # This is just 1 (since the bits determine H exactly)
    print(f"  Full model R² = 1.000000 (bits determine H)")
    print(f"  Interaction contribution = {1 - total_r2:.6f}")

    # TEST MULTIPLICATIVITY: Is H ≈ product of per-layer functions?
    # Compute per-layer conditional means
    print(f"\n--- MULTIPLICATIVITY TEST ---")

    layer_means = {}
    for s in sorted(skip_layers):
        idxs = skip_layers[s]
        group_H = defaultdict(list)
        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            layer_pattern = tuple(b[i] for i in idxs)
            group_H[layer_pattern].append(H_vals[bits])

        layer_means[s] = {k: sum(v)/len(v) for k, v in group_H.items()}

    # For multiplicativity: H ≈ mean_H * prod_s (E[H|layer_s] / mean_H)
    # Check: does this product predict H well?
    product_predicted = {}
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        pred = mean_H
        for s in sorted(skip_layers):
            idxs = skip_layers[s]
            layer_pattern = tuple(b[i] for i in idxs)
            pred *= layer_means[s][layer_pattern] / mean_H
        product_predicted[bits] = pred

    # R² of the product model
    ss_res = sum((H_vals[bits] - product_predicted[bits])**2 for bits in range(2**m))
    ss_tot = sum((H_vals[bits] - mean_H)**2 for bits in range(2**m))
    product_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    print(f"Product model R² = {product_r2:.6f}")
    print(f"(1.0 = perfect multiplicative factorization)")

    # Show some predictions vs actual
    import random
    random.seed(n)
    print(f"\nSample predictions:")
    for _ in range(10):
        bits = random.getrandbits(m) & ((1 << m) - 1)
        actual = H_vals[bits]
        pred = product_predicted[bits]
        print(f"  bits={bits:0{m}b}: actual={actual}, predicted={pred:.2f}, error={actual-pred:.2f}")

    # ADDITIVE test: H ≈ mean_H + sum_s (E[H|layer_s] - mean_H)
    additive_predicted = {}
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        pred = mean_H
        for s in sorted(skip_layers):
            idxs = skip_layers[s]
            layer_pattern = tuple(b[i] for i in idxs)
            pred += layer_means[s][layer_pattern] - mean_H
        additive_predicted[bits] = pred

    ss_res_add = sum((H_vals[bits] - additive_predicted[bits])**2 for bits in range(2**m))
    additive_r2 = 1 - ss_res_add / ss_tot if ss_tot > 0 else 0
    print(f"\nAdditive model R² = {additive_r2:.6f}")
    print(f"(This equals sum of individual R²: {total_r2:.6f})")

    # So the question is: which is better, additive or multiplicative?
    print(f"\nModel comparison:")
    print(f"  Additive:       R² = {additive_r2:.6f}")
    print(f"  Multiplicative: R² = {product_r2:.6f}")
    print(f"  Winner: {'MULTIPLICATIVE' if product_r2 > additive_r2 else 'ADDITIVE'}")
