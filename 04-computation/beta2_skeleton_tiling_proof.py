#!/usr/bin/env python3
"""
beta2_skeleton_tiling_proof.py — Connect skeleton structure to β₂=0 proof

KEY IDEA: The β₂=0 theorem is equivalent to surplus := dim(Ω₃) - dim(Z₂) ≥ 0.
The arc-flip identity (THM-100) gives:
  δ|A₃| = (n-3) · δ|A₂| under any single arc flip.

The tiling model provides a BACKBONE PATH p₀→p₁→...→p_{n-1}.
Each "tile" (non-backbone arc) can be independently flipped.
The GS (grid-symmetric) tilings form a subcube closed under flip.

STRATEGY:
1. Track how surplus changes along paths in the tiling cube
2. Start from transitive tournament (surplus = C(n-1,4) >> 0)
3. Show surplus stays ≥ 0 via the arc-flip formulas

In the tiling model, a single tile flip = reversing one non-backbone arc.
Kind-pasteur proved:
  δ|A₂| = 2(d_u - d_v - 1)  where u→v is the flipped arc
  δ|A₃| = (n-3) · δ|A₂|

QUESTION: Does δ(surplus) have a similar formula? Can we show
surplus(T') ≥ 0 whenever surplus(T) ≥ 0?

surplus = dim(Ω₃) - dim(Z₂) = dim(Ω₃) - [dim(Ω₂) - rk(∂₂|Ω₂)]
        = dim(Ω₃) - dim(Ω₂) + rk(∂₂|Ω₂)

This is NOT just a function of |A₂|,|A₃|. The dimensions of Ω₂,Ω₃
involve the STRUCTURE of which paths are TT/DT, not just the counts.

However, the arc-flip formulas for |A₂|,|A₃| hint at deeper structure.
Let's track dim(Ω₂), dim(Ω₃), dim(Z₂), surplus along tile flips.

Author: opus-2026-03-08-S44
"""
import sys, time, os, random
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def tournament_data(A, n, max_p=4):
    """Compute chain complex dimensions and surplus."""
    allowed = {}
    omega_dim = {}
    ranks = {}

    for p in range(max_p + 1):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega_dim[p] = n
        elif allowed[p]:
            om = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
            omega_dim[p] = om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0
        else:
            omega_dim[p] = 0

    # Compute ranks of boundary maps
    for p in range(1, max_p + 1):
        if omega_dim[p] == 0 or omega_dim[p-1] == 0:
            ranks[p] = 0
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1])
        om_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        om_prev = compute_omega_basis(A, n, p-1, allowed[p-1],
                                       allowed[p-2] if p >= 2 else [])
        bd_om = bd @ om_p
        # Project into Ω_{p-1} coordinates
        coords, _, _, _ = np.linalg.lstsq(om_prev, bd_om, rcond=None)
        S = np.linalg.svd(coords, compute_uv=False)
        ranks[p] = int(sum(s > 1e-8 for s in S))

    # Z₂ = ker(∂₂|Ω₂)
    z2 = omega_dim[2] - ranks[2]
    surplus = omega_dim[3] - z2 if omega_dim.get(3, 0) > 0 else -z2

    # Count TT triples and DT 4-paths
    tt_count = sum(1 for p in allowed[2] if A[p[0]][p[2]] == 1)
    nt_count = len(allowed[2]) - tt_count
    dt_count = sum(1 for p in allowed[3] if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)

    scores = sorted([sum(A[i]) for i in range(n)])
    t3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
             if i < j < k and A[i][j] and A[j][k] and A[k][i])
    t3 += sum(1 for i in range(n) for j in range(n) for k in range(n)
              if i < j < k and A[j][i] and A[i][k] and A[k][j])

    return {
        'A2': len(allowed[2]),
        'A3': len(allowed[3]),
        'Om2': omega_dim[2],
        'Om3': omega_dim[3],
        'Z2': z2,
        'surplus': surplus,
        'rk2': ranks[2],
        'rk3': ranks[3],
        'TT': tt_count,
        'NT': nt_count,
        'DT': dt_count,
        'scores': tuple(scores),
        't3': t3,
        'beta1': omega_dim[1] - ranks[1] - ranks[2],
        'beta2': z2 - ranks[3],
    }

def tiling_to_tournament(n, tiling_bits):
    """Convert a tiling (backbone path 0→1→...→n-1 plus tile bits) to adjacency matrix."""
    A = [[0]*n for _ in range(n)]
    # Backbone edges
    for i in range(n-1):
        A[i][i+1] = 1
    # Non-backbone edges (tiles)
    tile_idx = 0
    for i in range(n):
        for j in range(i+2, n):  # skip adjacent (backbone) edges
            if tiling_bits & (1 << tile_idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            tile_idx += 1
    return A

def count_tiles(n):
    """Number of non-backbone arc pairs."""
    return n*(n-1)//2 - (n-1)  # C(n,2) - (n-1)

def gs_free_bits(n):
    """Get the indices of GS-free tile bits (one per transpose-pair + fixed)."""
    tile_idx = 0
    tiles = []
    for i in range(n):
        for j in range(i+2, n):
            tiles.append((i, j))
            tile_idx += 1

    # Transpose map: (i,j) → (n-1-j, n-1-i)
    free_indices = []
    seen = set()
    for idx, (i, j) in enumerate(tiles):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if tj - ti < 2:  # backbone edge, skip
            continue
        partner_idx = tiles.index((ti, tj)) if (ti, tj) in tiles else None
        if partner_idx == idx:
            # Fixed tile (self-transpose)
            free_indices.append(idx)
            seen.add(idx)
        elif partner_idx is not None:
            free_indices.append(idx)
            seen.add(idx)
            seen.add(partner_idx)
    return free_indices

print("=" * 70)
print("TILING MODEL: SURPLUS TRACKING ALONG TILE FLIPS")
print("=" * 70)

for n in [5, 6]:
    m = count_tiles(n)
    total = 1 << m
    print(f"\n{'='*70}")
    print(f"n = {n}, m = {m} tiles, 2^m = {total} tilings")
    print(f"{'='*70}")

    # Compute data for ALL tilings
    all_data = {}
    t0 = time.time()
    for bits in range(total):
        if bits % 10000 == 0 and bits > 0:
            print(f"  ... {bits}/{total} ({time.time()-t0:.0f}s)")
        A = tiling_to_tournament(n, bits)
        all_data[bits] = tournament_data(A, n)

    print(f"  Computed all {total} tilings in {time.time()-t0:.1f}s")

    # === SURPLUS DISTRIBUTION ===
    surplus_dist = Counter()
    for bits, d in all_data.items():
        surplus_dist[d['surplus']] += 1
    print(f"\n  Surplus distribution:")
    for s in sorted(surplus_dist):
        print(f"    surplus={s}: {surplus_dist[s]} ({100*surplus_dist[s]/total:.1f}%)")
    print(f"  Min surplus: {min(surplus_dist)}")

    # === SINGLE TILE FLIP: δ(surplus) ===
    print(f"\n  Single tile flip: δ(surplus) distribution:")
    delta_surplus_dist = Counter()
    delta_surplus_by_dir = defaultdict(list)

    for bits in range(total):
        d = all_data[bits]
        for tile in range(m):
            bits2 = bits ^ (1 << tile)
            d2 = all_data[bits2]
            ds = d2['surplus'] - d['surplus']
            delta_surplus_dist[ds] += 1
            delta_surplus_by_dir[(d2['A2'] - d['A2'], d2['A3'] - d['A3'])].append(ds)

    for ds in sorted(delta_surplus_dist):
        print(f"    δ(surplus)={ds}: {delta_surplus_dist[ds]}")

    # === CORRELATION: δ(surplus) vs δ|A₂| ===
    print(f"\n  δ(surplus) vs δ|A₂| statistics:")
    delta_a2_vs_surplus = defaultdict(list)
    for bits in range(total):
        d = all_data[bits]
        for tile in range(m):
            bits2 = bits ^ (1 << tile)
            d2 = all_data[bits2]
            da2 = d2['A2'] - d['A2']
            ds = d2['surplus'] - d['surplus']
            delta_a2_vs_surplus[da2].append(ds)

    for da2 in sorted(delta_a2_vs_surplus):
        vals = delta_a2_vs_surplus[da2]
        print(f"    δ|A₂|={da2}: δ(surplus) range [{min(vals)}, {max(vals)}], "
              f"mean={np.mean(vals):.2f}, count={len(vals)}")

    # === GS TILING SUBNET ===
    gs_bits = gs_free_bits(n)
    n_gs = len(gs_bits)
    print(f"\n  GS tilings: {1 << n_gs} tilings, {n_gs} free bits")

    gs_surplus_dist = Counter()
    gs_tilings = []
    for gs_mask in range(1 << n_gs):
        # Build full tiling from GS free bits
        bits = 0
        for idx, fb in enumerate(gs_bits):
            if gs_mask & (1 << idx):
                bits |= (1 << fb)
        # Fill paired bits by symmetry
        tile_idx = 0
        tiles = []
        for i in range(n):
            for j in range(i+2, n):
                tiles.append((i, j))
                tile_idx += 1

        full_bits = 0
        for idx, (i, j) in enumerate(tiles):
            if bits & (1 << idx):
                full_bits |= (1 << idx)
            else:
                # Check if transpose partner is set
                ti, tj = n-1-j, n-1-i
                if ti > tj:
                    ti, tj = tj, ti
                if tj - ti >= 2:
                    partner_idx = tiles.index((ti, tj))
                    if bits & (1 << partner_idx):
                        full_bits |= (1 << idx)

        gs_tilings.append(full_bits)
        if full_bits in all_data:
            gs_surplus_dist[all_data[full_bits]['surplus']] += 1

    print(f"  GS surplus distribution:")
    for s in sorted(gs_surplus_dist):
        print(f"    surplus={s}: {gs_surplus_dist[s]}")

    # === TRANSITIVE TILING ===
    # Transitive tournament: i→j for all i<j. In tiling model, all tiles = 1
    trans_bits = (1 << m) - 1
    d_trans = all_data[trans_bits]
    print(f"\n  Transitive tournament (all tiles=1):")
    print(f"    surplus={d_trans['surplus']}, Ω₂={d_trans['Om2']}, Ω₃={d_trans['Om3']}, "
          f"Z₂={d_trans['Z2']}, |TT|={d_trans['TT']}, |DT|={d_trans['DT']}")

    # === PATH FROM TRANSITIVE TO WORST SURPLUS ===
    min_surplus_bits = min(all_data, key=lambda b: all_data[b]['surplus'])
    d_min = all_data[min_surplus_bits]
    print(f"\n  Minimum surplus tournament:")
    print(f"    bits={bin(min_surplus_bits)}, surplus={d_min['surplus']}, "
          f"scores={d_min['scores']}, t3={d_min['t3']}")

    # Find Hamming path from transitive to min-surplus
    current = trans_bits
    path = [current]
    target = min_surplus_bits
    for step in range(m):
        diff = current ^ target
        if diff == 0:
            break
        # Flip the first differing bit
        bit = diff & (-diff)  # lowest set bit
        tile_pos = bit.bit_length() - 1
        current ^= bit
        path.append(current)

    print(f"\n  Path from transitive to min-surplus ({len(path)-1} steps):")
    print(f"  {'step':>4}  {'surplus':>7}  {'δsurp':>6}  {'Ω₂':>3}  {'Ω₃':>3}  {'Z₂':>3}  {'δA₂':>4}  {'δA₃':>4}  {'TT':>3}  {'DT':>3}")
    prev_d = all_data[path[0]]
    for step, bits in enumerate(path):
        d = all_data[bits]
        ds = d['surplus'] - prev_d['surplus'] if step > 0 else 0
        da2 = d['A2'] - prev_d['A2'] if step > 0 else 0
        da3 = d['A3'] - prev_d['A3'] if step > 0 else 0
        print(f"  {step:>4}  {d['surplus']:>7}  {ds:>+6}  {d['Om2']:>3}  {d['Om3']:>3}  "
              f"{d['Z2']:>3}  {da2:>+4}  {da3:>+4}  {d['TT']:>3}  {d['DT']:>3}")
        prev_d = d

# === N=5 DEEP ANALYSIS: TILE FLIP EFFECT ON Ω₂, Ω₃ ===
print(f"\n{'='*70}")
print("N=5 DEEP: HOW TILE FLIPS CHANGE THE CHAIN COMPLEX")
print("="*70)

n = 5
m = count_tiles(n)
total = 1 << m

# For each tile flip, categorize the change in (δΩ₂, δΩ₃, δZ₂, δsurplus)
flip_effects = Counter()
flip_by_tile = defaultdict(lambda: Counter())

for bits in range(total):
    d = all_data[bits]
    for tile in range(m):
        bits2 = bits ^ (1 << tile)
        d2 = all_data[bits2]
        effect = (d2['Om2']-d['Om2'], d2['Om3']-d['Om3'],
                  d2['Z2']-d['Z2'], d2['surplus']-d['surplus'])
        flip_effects[effect] += 1
        flip_by_tile[tile][effect] += 1

print(f"\n  Flip effect (δΩ₂, δΩ₃, δZ₂, δsurplus) distribution:")
for effect in sorted(flip_effects, key=lambda e: (-flip_effects[e])):
    count = flip_effects[effect]
    pct = 100 * count / (total * m)
    print(f"    {effect}: {count} ({pct:.1f}%)")

# === KEY IDENTITY CHECK: δ|Ω₂| formula ===
print(f"\n  Checking if δ|Ω₂| depends only on δ|A₂| and δ(|TT|)...")
delta_om2_by_delta_tt_a2 = defaultdict(set)
for bits in range(total):
    d = all_data[bits]
    for tile in range(m):
        bits2 = bits ^ (1 << tile)
        d2 = all_data[bits2]
        dtt = d2['TT'] - d['TT']
        da2 = d2['A2'] - d['A2']
        dom2 = d2['Om2'] - d['Om2']
        delta_om2_by_delta_tt_a2[(dtt, da2)].add(dom2)

print(f"  (δTT, δA₂) → set of δΩ₂ values:")
for key in sorted(delta_om2_by_delta_tt_a2):
    vals = delta_om2_by_delta_tt_a2[key]
    print(f"    {key}: {sorted(vals)}")

# === CHECK: δΩ₂ = δ(TT) + δ(NT_cancel) ===
# Since Ω₂ = TT ⊕ NT_cancel, we need to understand how NT_cancel changes
print(f"\n  Checking δΩ₂ = δ(TT) + δ(NT_cancel_dim)...")
tt_nt_cancel_check = 0
for bits in range(total):
    d = all_data[bits]
    nt_cancel = d['Om2'] - d['TT']
    if nt_cancel < 0:
        tt_nt_cancel_check += 1

print(f"  Violations of Ω₂ ≥ TT: {tt_nt_cancel_check}")

# === ARC-FLIP vs TILE-FLIP RELATIONSHIP ===
print(f"\n{'='*70}")
print("TILE FLIP = ARC FLIP: RELATING δ|A₂|, δ|A₃| TO δΩ₂, δΩ₃")
print("="*70)

# For each flip, compute (u,v) = the arc being flipped
# Then check kind-pasteur's formulas
n = 5
tiles = []
for i in range(n):
    for j in range(i+2, n):
        tiles.append((i, j))

print(f"\n  Tiles (non-backbone arcs): {tiles}")

formula_check = 0
formula_total = 0
for bits in range(total):
    A = tiling_to_tournament(n, bits)
    d = all_data[bits]

    for tile_idx, (i, j) in enumerate(tiles):
        bits2 = bits ^ (1 << tile_idx)
        d2 = all_data[bits2]

        # Which direction is the arc?
        if A[i][j] == 1:
            u, v = i, j  # flip u→v to v→u
        else:
            u, v = j, i  # flip u→v to v→u

        du = sum(A[u])  # out-degree of u
        dv = sum(A[v])  # out-degree of v

        # Kind-pasteur formula: δ|A₂| = 2(d_u - d_v - 1)
        expected_da2 = 2 * (du - dv - 1)
        actual_da2 = d2['A2'] - d['A2']

        # THM-100: δ|A₃| = (n-3) * δ|A₂|
        expected_da3 = (n - 3) * actual_da2
        actual_da3 = d2['A3'] - d['A3']

        formula_total += 1
        if expected_da2 == actual_da2 and expected_da3 == actual_da3:
            formula_check += 1

print(f"\n  THM-100 verification: {formula_check}/{formula_total} correct")

# === CRITICAL: δΩ₂ FORMULA ===
print(f"\n  Seeking δΩ₂ formula in terms of (d_u, d_v, n)...")
delta_om2_by_degrees = defaultdict(list)
delta_om3_by_degrees = defaultdict(list)
delta_z2_by_degrees = defaultdict(list)
delta_surplus_by_degrees = defaultdict(list)

for bits in range(total):
    A = tiling_to_tournament(n, bits)
    d = all_data[bits]

    for tile_idx, (i, j) in enumerate(tiles):
        bits2 = bits ^ (1 << tile_idx)
        d2 = all_data[bits2]

        if A[i][j] == 1:
            u, v = i, j
        else:
            u, v = j, i

        du = sum(A[u])
        dv = sum(A[v])

        key = (du, dv)
        delta_om2_by_degrees[key].append(d2['Om2'] - d['Om2'])
        delta_om3_by_degrees[key].append(d2['Om3'] - d['Om3'])
        delta_z2_by_degrees[key].append(d2['Z2'] - d['Z2'])
        delta_surplus_by_degrees[key].append(d2['surplus'] - d['surplus'])

print(f"\n  δΩ₂ by (d_u, d_v) [flipping arc u→v]:")
print(f"  {'(du,dv)':>8}  {'mean':>6}  {'min':>4}  {'max':>4}  {'#unique':>7}  {'values':}")
for key in sorted(delta_om2_by_degrees):
    vals = delta_om2_by_degrees[key]
    uniqs = sorted(set(vals))
    print(f"  {str(key):>8}  {np.mean(vals):>6.2f}  {min(vals):>4}  {max(vals):>4}  {len(uniqs):>7}  {uniqs}")

print(f"\n  δΩ₃ by (d_u, d_v):")
print(f"  {'(du,dv)':>8}  {'mean':>6}  {'min':>4}  {'max':>4}  {'#unique':>7}  {'values':}")
for key in sorted(delta_om3_by_degrees):
    vals = delta_om3_by_degrees[key]
    uniqs = sorted(set(vals))
    print(f"  {str(key):>8}  {np.mean(vals):>6.2f}  {min(vals):>4}  {max(vals):>4}  {len(uniqs):>7}  {uniqs}")

print(f"\n  δZ₂ by (d_u, d_v):")
print(f"  {'(du,dv)':>8}  {'mean':>6}  {'min':>4}  {'max':>4}  {'#unique':>7}")
for key in sorted(delta_z2_by_degrees):
    vals = delta_z2_by_degrees[key]
    uniqs = sorted(set(vals))
    print(f"  {str(key):>8}  {np.mean(vals):>6.2f}  {min(vals):>4}  {max(vals):>4}  {len(uniqs):>7}")

print(f"\n  δ(surplus) by (d_u, d_v):")
print(f"  {'(du,dv)':>8}  {'mean':>6}  {'min':>4}  {'max':>4}  {'#unique':>7}")
for key in sorted(delta_surplus_by_degrees):
    vals = delta_surplus_by_degrees[key]
    uniqs = sorted(set(vals))
    print(f"  {str(key):>8}  {np.mean(vals):>6.2f}  {min(vals):>4}  {max(vals):>4}  {len(uniqs):>7}")

# === CHECK: IS SURPLUS MONOTONE FROM TRANSITIVE? ===
print(f"\n{'='*70}")
print("MONOTONICITY FROM TRANSITIVE TOURNAMENT")
print("="*70)
print("Can we reach every tournament from transitive by flips that never go negative?")

n = 5
m = count_tiles(n)
trans_bits = (1 << m) - 1

# BFS from transitive, only following non-negative surplus transitions
from collections import deque
reachable = {trans_bits}
queue = deque([trans_bits])
while queue:
    bits = queue.popleft()
    d = all_data[bits]
    for tile in range(m):
        bits2 = bits ^ (1 << tile)
        if bits2 in reachable:
            continue
        d2 = all_data[bits2]
        if d2['surplus'] >= 0:
            reachable.add(bits2)
            queue.append(bits2)

print(f"  n={n}: Reachable from transitive (surplus≥0 at each step): "
      f"{len(reachable)}/{total} ({100*len(reachable)/total:.1f}%)")
if len(reachable) == total:
    print(f"  *** ALL tilings reachable! surplus≥0 is CONNECTED. ***")

# === BLUESELF/BLACKSELF OVERLAY ===
print(f"\n{'='*70}")
print("BLUESELF/BLACKSELF STRUCTURE AND CHAIN COMPLEX")
print("="*70)

n = 5
m = count_tiles(n)
total = 1 << m

# GS flip = complement ALL tile bits
complement_bits = (1 << m) - 1

# Classify tilings
blueself_count = 0
blackself_count = 0
gs_count = 0

for bits in range(total):
    gs_bits_val = bits ^ complement_bits  # complement
    # Check if this is a GS tiling (grid-symmetric)
    A = tiling_to_tournament(n, bits)

    # Grid symmetry check: T ≅ T^{op} via reverse
    A_rev = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A_rev[n-1-i][n-1-j] = A[j][i]  # transpose + reverse labeling

    # Check if A_rev = A (tiling level, not isomorphism)
    is_gs = all(A[i][j] == A_rev[i][j] for i in range(n) for j in range(n))
    if is_gs:
        gs_count += 1

    # Blueself: complement = same tiling
    if gs_bits_val == bits:
        blueself_count += 1

print(f"  n={n}: GS tilings = {gs_count}, Blueself = {blueself_count}")

# Now: for GS tilings, compare their chain complex with their GS-complement
print(f"\n  GS flip pairs: chain complex comparison")
seen_gs_pairs = set()
for bits in range(total):
    A = tiling_to_tournament(n, bits)
    A_rev = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A_rev[n-1-i][n-1-j] = A[j][i]

    is_gs = all(A[i][j] == A_rev[i][j] for i in range(n) for j in range(n))
    if not is_gs:
        continue

    comp = bits ^ complement_bits
    pair = (min(bits, comp), max(bits, comp))
    if pair in seen_gs_pairs:
        continue
    seen_gs_pairs.add(pair)

    d1 = all_data[bits]
    d2 = all_data[comp]

    # How many tiles differ?
    diff = bin(bits ^ comp).count('1')

    if diff > 0:
        print(f"    {bin(bits):>10} ↔ {bin(comp):>10}: "
              f"Δtiles={diff}, "
              f"surplus {d1['surplus']}→{d2['surplus']}, "
              f"β₁ {d1['beta1']}→{d2['beta1']}, "
              f"Ω₂ {d1['Om2']}→{d2['Om2']}, "
              f"Ω₃ {d1['Om3']}→{d2['Om3']}")

# === KEY ANALYSIS: SURPLUS LOWER BOUND ===
print(f"\n{'='*70}")
print("SURPLUS LOWER BOUND ANALYSIS")
print("="*70)

for n_test in [5, 6]:
    m_test = count_tiles(n_test)
    total_test = 1 << m_test

    if n_test > 5:
        # n=6: already computed above
        pass

    # For each tournament, compute surplus and check:
    # surplus = Ω₃ - Z₂ = Ω₃ - (Ω₂ - rk(∂₂))
    # Known: Ω₂ = TT + NT_cancel
    # For transitive: Ω₂ = C(n,3) (all triples are TT), Ω₃ = C(n,4) (all DT)
    # surplus_trans = C(n,4) - C(n-1,3) = C(n-1,4) = ...

    # Find the tightest bound: surplus vs (n, |TT|, |DT|, |NT|)
    print(f"\n  n={n_test}: surplus as function of TT/NT/DT counts")

    surp_by_tt_dt = defaultdict(list)
    for bits, d in all_data.items():
        if not hasattr(d, '__getitem__'):
            continue
        surp_by_tt_dt[(d['TT'], d['DT'])].append(d['surplus'])

    # Find (TT,DT) pairs where surplus is minimized
    print(f"  (TT, DT) pairs with smallest min surplus:")
    min_pairs = []
    for key, surps in surp_by_tt_dt.items():
        min_pairs.append((min(surps), key, len(surps)))
    min_pairs.sort()
    for s, key, cnt in min_pairs[:10]:
        print(f"    TT={key[0]:>3}, DT={key[1]:>3}: min_surplus={s}, count={cnt}")

print("\nDone.")
