"""
fast_n8_recursion.py — opus-2026-04-04-S5

Speed up n=8 computation using the n→n+2 recursion (THM-291).

IDEA: H_8(inner, boundary) where inner = tiles of 6-tournament (10 tiles),
boundary = 2·6-1 = 11 tiles. Total = 21 tiles = C(7,2).

For EACH of the 2^10 = 1024 inner tilings, we know the 6-tournament's HP
count from precomputation. Then for each of the 2^11 = 2048 boundary
configurations, we compute H_8 using the DP.

OPTIMIZATION 1: Cache the DP state at the "inner only" level.
For the 6 inner vertices {2,...,7} (0-indexed: {1,...,6}), compute the DP
up to the mask including all inner vertices. Then extend to include
vertices 0 and 7 (the boundary endpoints) using boundary tile settings.

OPTIMIZATION 2: Exploit the recursion structure.
H_8(inner, 0) = H_6(inner) — skip the DP for boundary=0 entirely!

GOAL: Compute:
  1. E[H_8] = mean H over all 2^21 tilings
  2. P(8) = number of nonzero multilinear coefficients
  3. The H-spectrum at n=8 (set of achievable H values)
  4. The multilinear polynomial summary statistics
"""
import numpy as np
import time
from collections import defaultdict

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, bits):
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

def count_hp_fast(T, n):
    """Optimized HP count using int arrays."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            if mask == full:
                total += cnt
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def compute_n8_fast():
    """Compute H for all 2^21 tilings at n=8."""
    n = 8
    tiles = make_tiles(n)
    m = len(tiles)  # 21

    print(f"n={n}, m={m} tiles, 2^m = {1 << m} tilings")

    # Compute H for ALL tilings using optimized DP
    # 2^21 = 2,097,152 tilings × ~n^2 2^n DP each = significant
    # At n=8: 2^8 * 8^2 = 16384 ops per tiling × 2M tilings ≈ 34 billion ops
    # With Python overhead: ~hours. Let's sample first.

    # Strategy: compute in batches, track statistics
    h_spectrum = defaultdict(int)
    h_sum = 0
    h_count = 0
    h_min = float('inf')
    h_max = 0

    batch_size = 100000
    total = 1 << m

    t0 = time.time()

    # PHASE 1: Verify recursion at n=8 (small sample)
    print("\nPhase 1: Verify H_8(inner, 0) = H_6(inner) (first 100 inner tilings)...")
    tiles_6 = make_tiles(6)
    m_6 = len(tiles_6)
    tile_idx_6 = {t: i for i, t in enumerate(tiles_6)}
    tile_idx_8 = {t: i for i, t in enumerate(tiles)}

    # Map inner tiles: 6-tournament tile (x,y) → 8-tournament tile (x+1, y+1)
    inner_map = {}
    for i6, (x, y) in enumerate(tiles_6):
        t8 = (x+1, y+1)
        if t8 in tile_idx_8:
            inner_map[i6] = tile_idx_8[t8]

    recursion_ok = 0
    for inner_bits_6 in range(min(100, 1 << m_6)):
        # H_6
        T6 = tiling_to_tournament(6, tiles_6, inner_bits_6)
        h6 = count_hp_fast(T6, 6)

        # H_8 with boundary=0
        outer_bits = 0
        for i6, i8 in inner_map.items():
            if inner_bits_6 & (1 << i6):
                outer_bits |= (1 << i8)
        T8 = tiling_to_tournament(8, tiles, outer_bits)
        h8 = count_hp_fast(T8, 8)

        if h6 == h8:
            recursion_ok += 1

    print(f"  Recursion check: {recursion_ok}/100 match")

    # PHASE 2: Compute summary statistics over a SAMPLE
    print(f"\nPhase 2: Sample {min(200000, total)} random tilings for statistics...")
    import random
    random.seed(42)
    sample_size = min(200000, total)

    for trial in range(sample_size):
        if trial % 50000 == 0 and trial > 0:
            elapsed = time.time() - t0
            rate = trial / elapsed
            print(f"  {trial}/{sample_size} ({rate:.0f}/s, "
                  f"elapsed {elapsed:.0f}s)...", flush=True)

        bits = random.randint(0, total - 1)
        T = tiling_to_tournament(n, tiles, bits)
        h = count_hp_fast(T, n)

        h_spectrum[h] += 1
        h_sum += h
        h_count += 1
        h_min = min(h_min, h)
        h_max = max(h_max, h)

    elapsed = time.time() - t0
    mean_h = h_sum / h_count

    print(f"\n{'='*60}")
    print(f"n=8 STATISTICS ({sample_size} samples)")
    print(f"{'='*60}")
    print(f"H range: [{h_min}, {h_max}]")
    print(f"Mean H: {mean_h:.2f}")
    print(f"Distinct H values: {len(h_spectrum)}")
    print(f"Time: {elapsed:.1f}s")

    # E[H] from linear approximation
    mean_linear = 1 + sum((n - s) * (2 ** (s - 2)) for s in range(2, n))
    print(f"\nLinear approx E_lin[H] = {mean_linear}")
    print(f"Correction: {mean_h - mean_linear:.2f}")

    # H-spectrum density
    all_h_odd = sorted(h_spectrum.keys())
    print(f"\nH values (first 20): {all_h_odd[:20]}")
    print(f"H values (last 10): {all_h_odd[-10:]}")

    # Check: is H always odd? (Redei)
    all_odd = all(h % 2 == 1 for h in h_spectrum.keys())
    print(f"All H odd (Redei): {all_odd}")

    # PHASE 3: Use recursion to compute E[H_8] EXACTLY
    print(f"\nPhase 3: Exact E[H_8] via recursion...")
    # E[H_8] = (1/2^21) Σ_{all tilings} H_8(t)
    # = (1/2^11) Σ_{boundary b} (1/2^10) Σ_{inner i} H_8(i, b)
    # For b=0: (1/2^10) Σ_i H_8(i,0) = (1/2^10) Σ_i H_6(i) = E[H_6]

    # So E[H_8] = (1/2^11) [E[H_6]·2^10 + Σ_{b≠0} Σ_i H_8(i,b)/2^10]
    # = E[H_6]/2 + (1/2^11) Σ_{b≠0} E_i[H_8(·,b)]

    # Compute E[H_6] exactly
    print("  Computing E[H_6] exactly...")
    H6_all = np.zeros(1 << m_6, dtype=np.int64)
    for bits6 in range(1 << m_6):
        T6 = tiling_to_tournament(6, tiles_6, bits6)
        H6_all[bits6] = count_hp_fast(T6, 6)

    E_H6 = float(np.mean(H6_all))
    sum_H6 = int(np.sum(H6_all))
    print(f"  E[H_6] = {E_H6:.4f}")
    print(f"  Σ H_6 = {sum_H6}")

    # The exact E[H_8] = E[H_6]/2 + boundary_correction/2
    # The boundary correction = (1/2^11) Σ_{b=1}^{2^11-1} E_inner[H_8(·,b)]
    # Each E_inner requires summing over 1024 inner configs — too slow for all 2048 boundary configs.

    # Instead, use RANDOM SAMPLING for boundary:
    print("  Sampling 1000 boundary configs...")
    boundary_idx = sorted(i for i in range(m) if i not in set(inner_map.values()))
    n_bnd = len(boundary_idx)
    print(f"  Boundary tiles: {n_bnd}")

    boundary_means = []
    for trial in range(1000):
        bnd_bits_local = random.randint(0, (1 << n_bnd) - 1)
        # Convert to full boundary bits
        bnd_bits_full = 0
        for j, b_idx in enumerate(boundary_idx):
            if bnd_bits_local & (1 << j):
                bnd_bits_full |= (1 << b_idx)

        # For this boundary, compute E_inner[H_8]
        # Sample 100 inner configs
        inner_sum = 0
        for _ in range(100):
            inner_bits_6 = random.randint(0, (1 << m_6) - 1)
            outer_bits = bnd_bits_full
            for i6, i8 in inner_map.items():
                if inner_bits_6 & (1 << i6):
                    outer_bits |= (1 << i8)
            T8 = tiling_to_tournament(n, tiles, outer_bits)
            inner_sum += count_hp_fast(T8, n)
        boundary_means.append(inner_sum / 100)

    overall_mean = np.mean(boundary_means)
    print(f"  E[H_8] ≈ {overall_mean:.2f} (from boundary sampling)")
    print(f"  E[H_8] from direct sampling: {mean_h:.2f}")

    return h_spectrum

if __name__ == "__main__":
    compute_n8_fast()
