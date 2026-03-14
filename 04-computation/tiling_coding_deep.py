"""
tiling_coding_deep.py -- kind-pasteur-2026-03-14-S76
Deep coding theory analysis of the tiling model.

KEY DISCOVERY: GS weight distribution at n=7 is [1,3,9,19,33,51,65,75,75,65,51,33,19,9,3,1]
This is SYMMETRIC and unimodal. Is it a known sequence?

QUESTIONS:
1. Is the GS weight distribution a BINOMIAL convolution?
2. Does the strip structure create a PRODUCT CODE?
3. What is the dual code of the GS code?
4. Is there a connection to Hamming codes or BCH codes?
5. How does the weight distribution relate to H-distribution?

ALSO: The tiling model is a CODE on the triangular lattice.
The "codewords" = tournaments with specific H values.
The "parity checks" = the OCF constraint H = I(Omega, 2).
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import sys, math
from functools import reduce

sys.stdout.reconfigure(encoding='utf-8')

def pin_grid_positions(n):
    positions = []
    for r in range(1, n-1):
        for c in range(1, n-r):
            positions.append((r, c))
    return positions

def position_to_arc(r, c):
    return (r + c + 1, c)  # 1-indexed

def gs_map(r, c, n):
    return (r, n - r - c)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    for i in range(1, n):
        A[i][i-1] = 1
    positions = pin_grid_positions(n)
    for idx, (r, c) in enumerate(positions):
        a, b = position_to_arc(r, c)
        a0, b0 = a - 1, b - 1
        if bits & (1 << idx):
            A[b0][a0] = 1
        else:
            A[a0][b0] = 1
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

def get_gs_tilings(n):
    """Generate all GS tilings for tournament size n."""
    positions = pin_grid_positions(n)
    m = len(positions)
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
    tilings = []

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
        tilings.append(full_bits)

    return tilings, m, gs_dim, gs_fixed, gs_pairs

def main():
    print("=" * 70)
    print("DEEP CODING THEORY OF TOURNAMENT TILINGS")
    print("kind-pasteur-2026-03-14-S76")
    print("=" * 70)

    # ========================================
    # PART 1: Weight enumerator analysis
    # ========================================
    print(f"\n{'='*70}")
    print("PART 1: WEIGHT ENUMERATOR OF GS CODE")
    print(f"{'='*70}")

    for n in [5, 7, 9]:
        tilings, m, gs_dim, gs_fixed, gs_pairs = get_gs_tilings(n)
        N = len(tilings)

        weights = [bin(t).count('1') for t in tilings]
        W = Counter(weights)

        print(f"\n  n={n}: m={m}, gs_dim={gs_dim}, N={N}")
        print(f"  Weight enumerator W(z) = sum A_w z^w:")
        W_list = [W.get(w, 0) for w in range(m+1)]
        print(f"    {W_list}")
        print(f"    Sum = {sum(W_list)}")
        print(f"    Symmetric? {W_list == W_list[::-1]}")

        # Check: is this C(gs_dim, k) * something?
        # Or: product of (1+z) factors?
        # GS code = {x in F_2^m : x_i = x_{gs(i)} for all paired i}
        # Each fixed position contributes factor (1+z) to weight enumerator
        # Each paired position contributes factor (1+z^2) to weight enumerator
        # W(z) = (1+z)^{#fixed} * (1+z^2)^{#pairs}

        predicted = [0] * (m+1)
        # Expand (1+z)^f * (1+z^2)^p where f = #fixed, p = #pairs
        f = len(gs_fixed)
        p = len(gs_pairs)
        print(f"\n  Predicted: (1+z)^{f} * (1+z^2)^{p}")

        # Compute this product
        poly = [1]
        for _ in range(f):
            # multiply by (1+z)
            new_poly = [0] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new_poly[i] += c
                new_poly[i+1] += c
            poly = new_poly

        for _ in range(p):
            # multiply by (1+z^2)
            new_poly = [0] * (len(poly) + 2)
            for i, c in enumerate(poly):
                new_poly[i] += c
                new_poly[i+2] += c
            poly = new_poly

        # Pad to length m+1
        while len(poly) < m+1:
            poly.append(0)
        poly = poly[:m+1]

        print(f"    Predicted: {poly}")
        print(f"    Actual:    {W_list}")
        print(f"    Match? {poly == W_list}")

    # ========================================
    # PART 2: GS + H joint structure
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: GS TILINGS — JOINT (WEIGHT, H) DISTRIBUTION")
    print(f"{'='*70}")

    for n in [5, 7]:
        tilings, m, gs_dim, _, _ = get_gs_tilings(n)

        wh_dist = defaultdict(int)
        H_at_weight = defaultdict(list)

        for t in tilings:
            w = bin(t).count('1')
            A = bits_to_adj(t, n)
            H = compute_H_dp(A, n)
            wh_dist[(w, H)] += 1
            H_at_weight[w].append(H)

        print(f"\n  n={n}:")
        for w in sorted(H_at_weight.keys()):
            vals = H_at_weight[w]
            H_set = sorted(set(vals))
            mean_H = np.mean(vals)
            print(f"    weight {w:2d}: mean H={mean_H:7.2f}, "
                  f"H in {H_set[:8]}{'...' if len(H_set)>8 else ''}")

        # Is H determined by weight for GS tilings?
        det = all(len(set(v)) == 1 for v in H_at_weight.values())
        print(f"  Weight determines H? {det}")

        # Correlation
        ws = []
        hs = []
        for t in tilings:
            ws.append(bin(t).count('1'))
            A = bits_to_adj(t, n)
            hs.append(compute_H_dp(A, n))
        if len(set(ws)) > 1:
            print(f"  Corr(weight, H) = {np.corrcoef(ws, hs)[0,1]:.6f}")

    # ========================================
    # PART 3: The FLIP as a code operation
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: THE FLIP AS COMPLEMENT OPERATION")
    print("  Flip = bitwise complement of tiling")
    print("  In coding theory: the complement code")
    print(f"{'='*70}")

    for n in [5, 7]:
        tilings, m, gs_dim, _, _ = get_gs_tilings(n)
        flip_mask = (1 << m) - 1

        # For each GS tiling, what is H(flip(T))?
        flip_pairs = []
        for t in tilings:
            A = bits_to_adj(t, n)
            H_orig = compute_H_dp(A, n)

            t_flip = t ^ flip_mask
            A_flip = bits_to_adj(t_flip, n)
            H_flip = compute_H_dp(A_flip, n)

            flip_pairs.append((H_orig, H_flip))

        print(f"\n  n={n}: (H_original, H_flipped) distribution:")
        pair_dist = Counter(flip_pairs)
        for (h1, h2), cnt in sorted(pair_dist.items()):
            print(f"    ({h1:3d}, {h2:3d}): {cnt}")

        # Key: H + H_flip is always what?
        sums = [h1 + h2 for h1, h2 in flip_pairs]
        print(f"\n  H + H_flip values: {sorted(set(sums))}")
        print(f"  H + H_flip constant? {len(set(sums)) == 1}")

    # ========================================
    # PART 4: Blueself as self-dual code
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: BLUESELF = SELF-COMPLEMENTARY GS TILING")
    print("  A tiling is blueself iff T = flip(T) as tournaments (up to iso)")
    print("  This is like a SELF-DUAL codeword")
    print(f"{'='*70}")

    for n in [4, 5, 6, 7, 8]:
        if n > 7:
            print(f"\n  n={n}: skipping (too large)")
            continue

        tilings, m, gs_dim, _, _ = get_gs_tilings(n)
        flip_mask = (1 << m) - 1

        blueself_count = 0
        blueself_tilings = []

        for t in tilings:
            A_orig = bits_to_adj(t, n)
            H_orig = compute_H_dp(A_orig, n)
            scores_orig = tuple(sorted([sum(A_orig[i]) for i in range(n)]))

            t_flip = t ^ flip_mask
            A_flip = bits_to_adj(t_flip, n)
            H_flip = compute_H_dp(A_flip, n)
            scores_flip = tuple(sorted([sum(A_flip[i]) for i in range(n)]))

            # Blueself: same H AND same score sequence (necessary for isomorphism)
            if H_orig == H_flip and scores_orig == scores_flip:
                blueself_count += 1
                blueself_tilings.append((t, H_orig, scores_orig, bin(t).count('1')))

        print(f"\n  n={n}: {blueself_count} blueself GS tilings out of {len(tilings)}")
        for t, H, scores, w in blueself_tilings[:10]:
            print(f"    bits={t:0{m}b}, weight={w}, H={H}, scores={scores}")

        # Blueself weight = m/2? (complement of itself has same weight iff weight = m/2)
        if blueself_tilings:
            bs_weights = [w for _, _, _, w in blueself_tilings]
            print(f"  Blueself weights: {sorted(set(bs_weights))}")
            print(f"  Expected weight = m/2 = {m/2}: {all(w == m/2 for w in bs_weights)}")

    # ========================================
    # PART 5: H distribution from tiling perspective
    # ========================================
    print(f"\n{'='*70}")
    print("PART 5: H DISTRIBUTION — TILING vs ARC PERSPECTIVE")
    print("  Tiling model: only C(n-1,2) = m_tiling variables (not all C(n,2) arcs)")
    print("  Because the backbone arcs are fixed!")
    print("  Q: Does the tiling perspective give different Fourier structure?")
    print(f"{'='*70}")

    for n in [5]:
        positions = pin_grid_positions(n)
        m = len(positions)
        N = 2**m

        print(f"\n  n={n}: m_tiling={m} positions")

        # Compute H for all tilings
        H_values = np.zeros(N)
        for bits in range(N):
            A = bits_to_adj(bits, n)
            H_values[bits] = compute_H_dp(A, n)

        # Fourier transform on tiling space (m variables, not C(n,2))
        H_hat = H_values.copy()
        for i in range(m):
            step = 1 << (i + 1)
            half = 1 << i
            for j in range(0, N, step):
                for k in range(half):
                    u, v = H_hat[j+k], H_hat[j+k+half]
                    H_hat[j+k], H_hat[j+k+half] = u+v, u-v
        H_hat /= N

        # Energy by level
        level_energy = defaultdict(float)
        level_count = defaultdict(int)
        for S in range(N):
            level = bin(S).count('1')
            level_energy[level] += H_hat[S]**2
            if abs(H_hat[S]) > 1e-10:
                level_count[level] += 1

        total_E = sum(level_energy.values())
        print(f"\n  TILING-SPACE FOURIER SPECTRUM (m={m} variables):")
        for level in sorted(level_energy.keys()):
            pct = 100 * level_energy[level] / total_E
            nz = level_count[level]
            print(f"    level {level}: energy={pct:6.2f}%, {nz} nonzero")

        # Compare with arc-space Fourier (m_arc = C(n,2) = 10 variables)
        # From S73: levels 0 (75.95%), 2 (22.78%), 4 (1.27%)
        print(f"\n  Compare arc-space (m=10): levels 0 (75.95%), 2 (22.78%), 4 (1.27%)")
        print(f"  Tiling-space has DIFFERENT Fourier structure because:")
        print(f"  - Fewer variables (6 vs 10)")
        print(f"  - Different basis (tiling positions vs arc indicators)")
        print(f"  - Backbone arcs are NOT free variables in tiling space")

        # KEY: odd levels in tiling space?
        odd_energy = sum(level_energy[k] for k in level_energy if k % 2 == 1)
        print(f"\n  Odd-level energy in tiling space: {100*odd_energy/total_E:.4f}%")
        print(f"  (In arc space: exactly 0% because H(T)=H(T^op))")
        print(f"  In tiling space: the T^op symmetry maps to DIFFERENT tiling,")
        print(f"  so the odd-level vanishing does NOT hold in tiling coordinates!")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
