#!/usr/bin/env python3
"""
The Waggly Spectrum: H-effects as a function of Hamming distance d.
opus-2026-04-03-S28

The waggly layers d=1,...,m connect tilings at different Hamming distances.
d=1 is wiggly (flip 1 tile), d=m is complement (flip all).

QUESTION: How does the mean |delta H| vary with d?
Is it monotone? U-shaped? What's the duality?

Also: Is there a FOURIER DUALITY between d=1 and d=m effects?
In the Walsh-Fourier picture:
  - Order-k Walsh coefficients contribute to d-flip effects via inclusion-exclusion
  - d=1 flips see ALL Walsh orders (each contributes proportionally)
  - d=m flips see only EVEN orders (by complement invariance)
  - So d=1 and d=m are seeing DIFFERENT subsets of the spectrum!

This is the waggly version of the Uncertainty Principle:
  d=1 is "localized" in tile space (flip one tile) but "spread" in Walsh space
  d=m is "spread" in tile space (flip all tiles) but "localized" in Walsh space
"""

from collections import defaultdict
from math import comb
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

def hamming_weight(x, m):
    return bin(x & ((1 << m) - 1)).count('1')

def analyze_waggly_spectrum(n):
    m = (n-1)*(n-2)//2
    print(f"\n{'='*70}")
    print(f"WAGGLY SPECTRUM ANALYSIS: n={n}, m={m}")
    print(f"{'='*70}")

    # Precompute all H values
    H_vals = {}
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H_vals[bits] = count_hp(T, n)

    mean_H = sum(H_vals.values()) / len(H_vals)
    var_H = sum((h - mean_H)**2 for h in H_vals.values()) / len(H_vals)
    print(f"Mean H = {mean_H:.4f}, Var H = {var_H:.4f}")

    # For each Hamming distance d, compute statistics of delta-H
    print(f"\n--- MEAN |delta H| BY HAMMING DISTANCE ---")
    print(f"{'d':>3} {'C(m,d)':>8} {'mean_delta':>12} {'mean_|delta|':>14} {'var_delta':>12} {'neutral_%':>10}")

    delta_by_d = defaultdict(list)

    # For small m, do exhaustive; for large m, sample
    if 2**m <= 2048:
        for bits_a in range(2**m):
            H_a = H_vals[bits_a]
            # Instead of all pairs, sample neighbors at each distance
            for bits_b in range(bits_a + 1, 2**m):
                d = hamming_weight(bits_a ^ bits_b, m)
                H_b = H_vals[bits_b]
                delta = H_b - H_a
                delta_by_d[d].append(delta)
    else:
        # Sample: for each tiling, sample random neighbors at each distance
        import random
        random.seed(42)
        sample_per_d = min(5000, 2**m)

        for _ in range(sample_per_d):
            bits_a = random.getrandbits(m) & ((1 << m) - 1)
            H_a = H_vals[bits_a]

            for d in range(1, m+1):
                # Generate a random neighbor at distance d
                # Pick d random bit positions to flip
                positions = random.sample(range(m), d)
                bits_b = bits_a
                for p in positions:
                    bits_b ^= (1 << p)
                H_b = H_vals[bits_b]
                delta = H_b - H_a
                delta_by_d[d].append(delta)

    for d in sorted(delta_by_d):
        deltas = delta_by_d[d]
        n_pairs = len(deltas)
        mean_d = sum(deltas) / n_pairs
        mean_abs = sum(abs(x) for x in deltas) / n_pairs
        var_d = sum((x - mean_d)**2 for x in deltas) / n_pairs
        neutral = sum(1 for x in deltas if x == 0) / n_pairs * 100
        print(f"{d:>3} {comb(m,d):>8} {mean_d:>12.4f} {mean_abs:>14.4f} {var_d:>12.4f} {neutral:>9.2f}%")

    # KEY: Is mean |delta| a U-shaped function of d?
    # Or monotone increasing? Or what?
    print(f"\n--- SHAPE OF THE WAGGLY SPECTRUM ---")
    if delta_by_d:
        d_values = sorted(delta_by_d.keys())
        abs_means = [sum(abs(x) for x in delta_by_d[d]) / len(delta_by_d[d]) for d in d_values]

        # Normalize by 2*var_H for comparison
        if var_H > 0:
            norm_means = [x / (2*var_H)**0.5 for x in abs_means]
        else:
            norm_means = abs_means

        print(f"d vs normalized mean |delta H|:")
        for d, nm in zip(d_values, norm_means):
            bar = '#' * int(nm * 40)
            print(f"  d={d:>2}: {nm:.4f} {bar}")

    # WALSH-FOURIER DECOMPOSITION
    print(f"\n--- WALSH-FOURIER SPECTRUM OF H ---")
    # H_hat(S) = (1/2^m) * sum_t (-1)^{<S,t>} * H(t)
    # where <S,t> = sum of S_j * t_j mod 2

    walsh = {}
    for S in range(2**m):
        total = 0
        for bits in range(2**m):
            # (-1)^{popcount(S & bits)}
            sign = (-1) ** bin(S & bits).count('1')
            total += sign * H_vals[bits]
        walsh[S] = total / (2**m)

    # Group by Walsh order (= popcount of S)
    walsh_by_order = defaultdict(float)
    walsh_energy_by_order = defaultdict(float)
    for S, coeff in walsh.items():
        order = bin(S).count('1')
        walsh_energy_by_order[order] += coeff**2

    total_energy = sum(walsh_energy_by_order.values())
    print(f"Walsh spectrum (energy by order):")
    for order in sorted(walsh_energy_by_order):
        energy = walsh_energy_by_order[order]
        pct = energy / total_energy * 100 if total_energy > 0 else 0
        print(f"  Order {order}: energy = {energy:.4f} ({pct:.2f}%)")

    # KEY QUESTION: How does each Walsh order contribute to delta at distance d?
    # Theory: for a d-flip (flipping d specific tiles),
    #   delta_H = 2 * sum_{S: S intersects flipped set non-trivially} walsh(S) * ...
    # More precisely, the expected |delta|^2 at distance d involves
    #   sum_{S} |walsh(S)|^2 * (something depending on |S| and d)

    # The "something" is: E[(-1)^{popcount(S∩F)} - 1]^2 where F is the set of d flipped tiles
    # This equals 4 * sin^2(pi * overlap / 2) averaged over all C(m,d) choices of F
    # For order-k coefficient with d flips:
    #   E[((-1)^{popcount(S∩F)} - 1)^2] = 2 * (1 - prod over the m-d non-flipped positions)
    # Actually, let me compute this directly.

    print(f"\n--- WALSH ORDER CONTRIBUTION TO DISTANCE-d VARIANCE ---")
    # Var(delta_d) = sum_S |walsh_hat(S)|^2 * g(|S|, d)
    # where g(k, d) = E_{|F|=d} [((-1)^{|S∩F|} - 1)^2]
    # = 2(1 - E[(-1)^{|S∩F|}])
    # E[(-1)^{|S∩F|}] for |S|=k, |F|=d:
    #   = sum_{j=0}^{min(k,d)} C(k,j)*C(m-k,d-j)/C(m,d) * (-1)^j
    #   = Krawtchouk polynomial K_d(k; m) / C(m,d)?
    # Actually this is the Krawtchouk evaluation:
    #   E[(-1)^{|S∩F|}] = K_k(d; m) / C(m,k) ... no.
    # Let me just compute it directly.

    def expected_minus_one_overlap(k, d, m):
        """E[(-1)^{|S∩F|}] where |S|=k, |F|=d, universe size m."""
        total = 0
        for j in range(min(k, d) + 1):
            c = comb(k, j) * comb(m - k, d - j)
            total += c * ((-1) ** j)
        return total / comb(m, d)

    print(f"Expected (-1)^overlap for Walsh order k, flip distance d:")
    hdr = "k\\d"
    print(f"{hdr:>4}", end='')
    for d in range(1, m+1):
        print(f"  d={d:>2}", end='')
    print()
    for k in range(m+1):
        if walsh_energy_by_order.get(k, 0) < 1e-10:
            continue
        print(f"k={k:>2}:", end='')
        for d in range(1, m+1):
            val = expected_minus_one_overlap(k, d, m)
            print(f"  {val:>5.3f}", end='')
        print()

    # So g(k,d) = 2*(1 - E[(-1)^overlap])
    print(f"\ng(k,d) = 2*(1 - E[(-1)^overlap]) — sensitivity of order-k to distance-d flips:")
    print(f"{hdr:>4}", end='')
    for d in range(1, m+1):
        print(f"  d={d:>2}", end='')
    print()
    for k in range(m+1):
        if walsh_energy_by_order.get(k, 0) < 1e-10:
            continue
        print(f"k={k:>2}:", end='')
        for d in range(1, m+1):
            val = expected_minus_one_overlap(k, d, m)
            g = 2 * (1 - val)
            print(f"  {g:>5.3f}", end='')
        print()

    # The predicted variance at distance d:
    print(f"\nPredicted vs actual Var(delta) at each distance d:")
    for d in sorted(delta_by_d):
        predicted = 0
        for k in range(m+1):
            energy = walsh_energy_by_order.get(k, 0)
            if energy < 1e-10:
                continue
            val = expected_minus_one_overlap(k, d, m)
            g = 2 * (1 - val)
            predicted += energy * g

        actual_deltas = delta_by_d[d]
        actual_mean = sum(actual_deltas) / len(actual_deltas)
        actual_var = sum((x - actual_mean)**2 for x in actual_deltas) / len(actual_deltas)
        # Note: var(delta) = E[delta^2] since E[delta]=0
        actual_e_sq = sum(x**2 for x in actual_deltas) / len(actual_deltas)

        print(f"  d={d}: predicted E[delta²] = {predicted:.4f}, actual E[delta²] = {actual_e_sq:.4f}, "
              f"ratio = {actual_e_sq/predicted:.4f}" if predicted > 0 else f"  d={d}: predicted = 0")

    return walsh_energy_by_order, delta_by_d

# Run
for n in [4, 5, 6]:
    analyze_waggly_spectrum(n)
