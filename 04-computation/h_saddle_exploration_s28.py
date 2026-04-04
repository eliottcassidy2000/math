#!/usr/bin/env python3
"""
H-Saddle Exploration in (hypotenuse, perpendicular) coordinates.
opus-2026-04-03-S28

The staircase triangle has two natural perpendicular coordinate families:
  - x+y = c  (parallel to hypotenuse, "distance from origin corner")
  - x-y = c  (perpendicular to hypotenuse, "skip value")

where tile (x,y) has x=upper vertex, y=lower vertex, skip = x-y >= 2.

HYPOTENUSE tiles: those on the anti-diagonal (highest x+y values)
PERIPHERY tiles: those near the base path (lowest x+y values, skip=1 edges)

The H-saddle conjecture: H is maximized at OPPOSITE corners of the
(hyp_weight, per_weight) plane, creating a saddle structure.

This script maps H exhaustively across these coordinates for n=3..7.
"""

from itertools import product, permutations
from collections import defaultdict
import sys

def tournament_from_bits(n, bits):
    """Convert m-bit tiling to tournament adjacency matrix."""
    T = [[0]*n for _ in range(n)]
    # Base path: n-1 -> n-2 -> ... -> 0  (i.e., i beats j if i > j for base path)
    # In our convention: vertex labels 1..n, base path n->n-1->...->1
    # Tile (x,y) with x-y >= 2 means arc between vertices x and y
    # bit=0: x->y (forward/default), bit=1: y->x (backward/flipped)

    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    # Set base path arcs (not tiles, always present)
    for i in range(n, 1, -1):
        T[i-1][i-2] = 1  # i -> i-1 (0-indexed: T[i-1][i-2] = 1)

    # Set tile arcs
    for idx, (x, y) in enumerate(tiles):
        if idx < len(bits):
            if bits[idx] == 0:
                T[x-1][y-1] = 1  # x -> y
            else:
                T[y-1][x-1] = 1  # y -> x

    return T, tiles

def count_hamiltonian_paths(T, n):
    """Count Hamiltonian paths using DP on bitmask."""
    # dp[mask][v] = number of Hamiltonian paths ending at vertex v
    # using exactly the vertices in mask
    dp = [[0]*n for _ in range(1 << n)]

    for v in range(n):
        dp[1 << v][v] = 1

    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:  # v -> u
                    dp[mask | (1 << u)][u] += dp[mask][v]

    return sum(dp[full])

def get_tile_coordinates(n):
    """Get tiles with their (x,y) coordinates and derived quantities."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    tile_info = []
    for x, y in tiles:
        skip = x - y          # perpendicular coordinate (x - y)
        dist = x + y          # parallel to hypotenuse (x + y)
        tile_info.append({
            'x': x, 'y': y,
            'skip': skip,
            'dist': dist,
            'is_hypotenuse': (x + y == n),   # on the main anti-diagonal
            'is_periphery': (x - y == 2),      # skip-2 = minimum non-base
        })
    return tiles, tile_info

def classify_tiles(tile_info, n):
    """Classify tiles into hypotenuse vs periphery vs interior."""
    max_dist = max(t['dist'] for t in tile_info)
    min_skip = min(t['skip'] for t in tile_info)

    hyp_indices = [i for i, t in enumerate(tile_info) if t['dist'] == max_dist]
    per_indices = [i for i, t in enumerate(tile_info) if t['skip'] == min_skip]
    int_indices = [i for i, t in enumerate(tile_info)
                   if t['dist'] != max_dist and t['skip'] != min_skip]

    return hyp_indices, per_indices, int_indices

def analyze_n(n, verbose=True):
    """Full H-saddle analysis for tournament size n."""
    m = (n-1)*(n-2)//2
    tiles, tile_info = get_tile_coordinates(n)
    assert len(tiles) == m

    if verbose:
        print(f"\n{'='*70}")
        print(f"n = {n}, m = {m} tiles")
        print(f"{'='*70}")
        print(f"\nTile layout:")
        for i, t in enumerate(tile_info):
            label = chr(65 + i) if i < 26 else f"T{i}"
            flags = []
            if t['is_hypotenuse']: flags.append("HYP")
            if t['is_periphery']: flags.append("PER")
            if not flags: flags.append("INT")
            print(f"  {label}: ({t['x']},{t['y']}) skip={t['skip']} dist={t['dist']} [{','.join(flags)}]")

    # Distance layers (x+y = c)
    dist_layers = defaultdict(list)
    for i, t in enumerate(tile_info):
        dist_layers[t['dist']].append(i)

    # Skip layers (x-y = s)
    skip_layers = defaultdict(list)
    for i, t in enumerate(tile_info):
        skip_layers[t['skip']].append(i)

    if verbose:
        print(f"\nDistance layers (x+y=c, parallel to hypotenuse):")
        for c in sorted(dist_layers):
            indices = dist_layers[c]
            print(f"  x+y={c}: {len(indices)} tiles — {[chr(65+i) for i in indices]}")

        print(f"\nSkip layers (x-y=s, perpendicular to hypotenuse):")
        for s in sorted(skip_layers):
            indices = skip_layers[s]
            print(f"  skip={s}: {len(indices)} tiles — {[chr(65+i) for i in indices]}")

    # Enumerate all 2^m tilings and compute H
    print(f"\nEnumerating {2**m} tilings...")

    # For each tiling, compute:
    #   - H value
    #   - hyp_weight = sum of bits on hypotenuse tiles
    #   - per_weight = sum of bits on periphery (skip-2) tiles
    #   - total_weight = Hamming weight
    #   - dist_profile = tuple of sums per distance layer
    #   - skip_profile = tuple of sums per skip layer

    max_dist = max(dist_layers.keys())
    min_skip = min(skip_layers.keys())
    hyp_idx = [i for i, t in enumerate(tile_info) if t['dist'] == max_dist]
    per_idx = [i for i, t in enumerate(tile_info) if t['skip'] == min_skip]

    # Data collection
    H_by_hyp_per = defaultdict(list)  # (hyp_w, per_w) -> [H values]
    H_by_dist_profile = defaultdict(list)
    H_by_skip_profile = defaultdict(list)
    H_all = []

    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T, _ = tournament_from_bits(n, b)
        H = count_hamiltonian_paths(T, n)

        hyp_w = sum(b[i] for i in hyp_idx)
        per_w = sum(b[i] for i in per_idx)
        total_w = sum(b)

        dist_sorted = sorted(dist_layers.keys())
        dist_prof = tuple(sum(b[i] for i in dist_layers[c]) for c in dist_sorted)

        skip_sorted = sorted(skip_layers.keys())
        skip_prof = tuple(sum(b[i] for i in skip_layers[s]) for s in skip_sorted)

        H_by_hyp_per[(hyp_w, per_w)].append(H)
        H_by_dist_profile[dist_prof].append(H)
        H_by_skip_profile[skip_prof].append(H)
        H_all.append((H, hyp_w, per_w, total_w, dist_prof, skip_prof, bits))

    # Saddle analysis
    print(f"\n--- H-SADDLE ANALYSIS (hypotenuse weight, periphery weight) ---")
    print(f"Hypotenuse tiles ({len(hyp_idx)}): {[chr(65+i) for i in hyp_idx]}")
    print(f"Periphery tiles ({len(per_idx)}): {[chr(65+i) for i in per_idx]}")

    max_hyp = len(hyp_idx)
    max_per = len(per_idx)

    print(f"\nMean H in (hyp_weight, per_weight) grid:")
    print(f"{'':>8}", end='')
    for pw in range(max_per + 1):
        print(f"  per={pw:>2}", end='')
    print()

    saddle_grid = {}
    for hw in range(max_hyp + 1):
        print(f"hyp={hw:>2}:", end='')
        for pw in range(max_per + 1):
            vals = H_by_hyp_per.get((hw, pw), [])
            if vals:
                mean_h = sum(vals) / len(vals)
                saddle_grid[(hw, pw)] = mean_h
                print(f"  {mean_h:>6.1f}", end='')
            else:
                print(f"  {'---':>6}", end='')
        print()

    print(f"\nMax H in (hyp_weight, per_weight) grid:")
    print(f"{'':>8}", end='')
    for pw in range(max_per + 1):
        print(f"  per={pw:>2}", end='')
    print()

    for hw in range(max_hyp + 1):
        print(f"hyp={hw:>2}:", end='')
        for pw in range(max_per + 1):
            vals = H_by_hyp_per.get((hw, pw), [])
            if vals:
                print(f"  {max(vals):>6}", end='')
            else:
                print(f"  {'---':>6}", end='')
        print()

    print(f"\nMin H in (hyp_weight, per_weight) grid:")
    print(f"{'':>8}", end='')
    for pw in range(max_per + 1):
        print(f"  per={pw:>2}", end='')
    print()

    for hw in range(max_hyp + 1):
        print(f"hyp={hw:>2}:", end='')
        for pw in range(max_per + 1):
            vals = H_by_hyp_per.get((hw, pw), [])
            if vals:
                print(f"  {min(vals):>6}", end='')
            else:
                print(f"  {'---':>6}", end='')
        print()

    # The KEY question: is there a saddle?
    # A saddle means: H increases along one diagonal, decreases along the other
    # Check: moving from (0,0) -> (max,max) vs (0,max) -> (max,0)
    corners = {}
    for hw, pw in [(0,0), (0,max_per), (max_hyp,0), (max_hyp,max_per)]:
        vals = H_by_hyp_per.get((hw, pw), [])
        if vals:
            corners[(hw,pw)] = (min(vals), sum(vals)/len(vals), max(vals), len(vals))

    print(f"\nCorner analysis:")
    for (hw, pw), (mn, avg, mx, cnt) in sorted(corners.items()):
        print(f"  ({hw},{pw}): min={mn}, mean={avg:.1f}, max={mx}, count={cnt}")

    # Diagonal comparison
    if (0,0) in corners and (max_hyp, max_per) in corners:
        diag1_start = corners[(0,0)][1]
        diag1_end = corners[(max_hyp, max_per)][1]
        print(f"\nMain diagonal (0,0)->(max,max): {diag1_start:.1f} -> {diag1_end:.1f}")
    if (0, max_per) in corners and (max_hyp, 0) in corners:
        diag2_start = corners[(0, max_per)][1]
        diag2_end = corners[(max_hyp, 0)][1]
        print(f"Anti diagonal (0,max)->(max,0): {diag2_start:.1f} -> {diag2_end:.1f}")

    # Wiggly/complement duality analysis
    print(f"\n--- WIGGLY / COMPLEMENT DUALITY ---")
    # For each tiling, compute delta-H for each wiggly flip (d=1) and the complement flip (d=m)

    wiggly_deltas_by_tile = defaultdict(list)
    complement_delta = []

    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T, _ = tournament_from_bits(n, b)
        H = count_hamiltonian_paths(T, n)

        # Wiggly: flip each tile individually
        for j in range(m):
            b2 = b.copy()
            b2[j] = 1 - b2[j]
            T2, _ = tournament_from_bits(n, b2)
            H2 = count_hamiltonian_paths(T2, n)
            delta = H2 - H
            wiggly_deltas_by_tile[j].append(delta)

        # Complement: flip ALL tiles
        b_comp = [1 - x for x in b]
        T_comp, _ = tournament_from_bits(n, b_comp)
        H_comp = count_hamiltonian_paths(T_comp, n)
        complement_delta.append(H_comp - H)

    print(f"\nWiggly (d=1) delta-H statistics per tile:")
    for j in range(m):
        deltas = wiggly_deltas_by_tile[j]
        t = tile_info[j]
        label = chr(65+j)
        mean_d = sum(deltas) / len(deltas)
        mean_abs = sum(abs(d) for d in deltas) / len(deltas)
        # How often is delta = 0 (neutral flip)?
        neutral_rate = sum(1 for d in deltas if d == 0) / len(deltas)
        print(f"  {label} ({t['x']},{t['y']}) skip={t['skip']} dist={t['dist']}: "
              f"mean={mean_d:+.2f} |mean|={mean_abs:.2f} neutral={neutral_rate:.3f}")

    print(f"\nComplement (d=m) delta-H statistics:")
    mean_comp = sum(complement_delta) / len(complement_delta)
    mean_abs_comp = sum(abs(d) for d in complement_delta) / len(complement_delta)
    neutral_comp = sum(1 for d in complement_delta if d == 0) / len(complement_delta)
    print(f"  mean={mean_comp:+.2f} |mean|={mean_abs_comp:.2f} neutral={neutral_comp:.3f}")

    # Distribution of complement deltas
    comp_hist = defaultdict(int)
    for d in complement_delta:
        comp_hist[d] += 1
    print(f"  Distribution: {dict(sorted(comp_hist.items()))}")

    # KEY: correlation between wiggly and complement effects
    # For each tiling, the average wiggly |delta| vs the complement |delta|
    print(f"\n--- WIGGLY-COMPLEMENT CORRELATION ---")
    wig_means = []
    comp_vals = []
    for idx in range(2**m):
        avg_wig = sum(abs(wiggly_deltas_by_tile[j][idx]) for j in range(m)) / m
        wig_means.append(avg_wig)
        comp_vals.append(abs(complement_delta[idx]))

    # Pearson correlation
    n_pts = len(wig_means)
    mean_w = sum(wig_means) / n_pts
    mean_c = sum(comp_vals) / n_pts
    cov = sum((wig_means[i] - mean_w) * (comp_vals[i] - mean_c) for i in range(n_pts)) / n_pts
    var_w = sum((wig_means[i] - mean_w)**2 for i in range(n_pts)) / n_pts
    var_c = sum((comp_vals[i] - mean_c)**2 for i in range(n_pts)) / n_pts
    corr = cov / (var_w * var_c)**0.5 if var_w > 0 and var_c > 0 else 0
    print(f"  Pearson corr(avg |wiggly delta|, |complement delta|) = {corr:.4f}")

    # Distance-layer decomposition of H
    print(f"\n--- DISTANCE-LAYER DECOMPOSITION ---")
    print(f"Distance layers: {dict((c, len(dist_layers[c])) for c in sorted(dist_layers))}")

    # For each distance layer, what fraction of H-variance does it explain?
    overall_mean_H = sum(h for h,_,_,_,_,_,_ in H_all) / len(H_all)
    overall_var_H = sum((h - overall_mean_H)**2 for h,_,_,_,_,_,_ in H_all) / len(H_all)
    print(f"Overall: mean(H) = {overall_mean_H:.2f}, var(H) = {overall_var_H:.2f}")

    for layer_type, layers in [("dist (x+y)", dist_layers), ("skip (x-y)", skip_layers)]:
        sorted_keys = sorted(layers.keys())
        # Conditional mean of H given this layer's weight
        for c in sorted_keys:
            layer_idx = layers[c]
            layer_size = len(layer_idx)
            weight_to_H = defaultdict(list)
            for h, hw, pw, tw, dp, sp, bits in H_all:
                b = [(bits >> i) & 1 for i in range(m)]
                w = sum(b[i] for i in layer_idx)
                weight_to_H[w].append(h)

            # Explained variance by this layer
            between_var = sum(len(weight_to_H[w]) * (sum(weight_to_H[w])/len(weight_to_H[w]) - overall_mean_H)**2
                             for w in weight_to_H) / len(H_all)
            r2 = between_var / overall_var_H if overall_var_H > 0 else 0

            mean_per_w = {w: sum(weight_to_H[w])/len(weight_to_H[w]) for w in sorted(weight_to_H)}
            print(f"  {layer_type}={c} ({layer_size} tiles): R²={r2:.4f}  means={dict(sorted(mean_per_w.items()))}")

    return H_all, saddle_grid, tile_info

if __name__ == '__main__':
    for n in [3, 4, 5, 6]:
        analyze_n(n)
        print()
