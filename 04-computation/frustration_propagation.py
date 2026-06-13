"""
frustration_propagation.py — opus-2026-04-04-S19

THE FRUSTRATION PROPAGATION FORMULA.

When we go n → n+1:
- OLD tiles are preserved (THM-299): their couplings unchanged
- NEW tiles (n-1 tiles involving vertex n+1) create:
  (a) NEW-NEW frustration: same upper vertex n+1 → all pairs are same-end → c = -2
  (b) NEW-OLD frustration: new tile shares lower vertex with old tiles → same-end → negative
  (c) NEW-OLD synergy: new tile is cross-end or disjoint with old tiles → positive

The TOTAL FRUSTRATION = sum of |negative c_{ij}| = old frustration + new frustration.

RECURSIVE FORMULA:
  F(n+1) = F(n) + F_new-new(n+1) + F_new-old(n+1)

where:
  F_new-new = sum of |c_{ij}| over pairs of new tiles (all same-end)
  F_new-old = sum of |c_{ij}| over (new, old) pairs that are same-end
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from math import comb
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def tile_relation(t1, t2):
    x1, y1 = t1; x2, y2 = t2
    if x1 == x2 or y1 == y2: return "same-end"
    if x1 == y2 or y1 == x2: return "cross-end"
    return "disjoint"

def frustration_decomposition(n):
    """Decompose the total quadratic frustration by generation (old vs new tiles)."""
    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # New tiles = those with x = n (upper vertex = n)
    new_indices = [i for i, (x,y) in enumerate(tiles) if x == n]
    old_indices = [i for i, (x,y) in enumerate(tiles) if x < n]

    # Decompose quadratic coefficients
    F_old_old = 0   # frustration between old tiles (inherited from n-1)
    F_new_new = 0   # frustration between new tiles
    F_new_old_neg = 0  # negative (frustration) between new and old
    F_new_old_pos = 0  # positive (synergy) between new and old

    P_old_old = 0   # positive coupling between old tiles
    P_new_new = 0
    P_new_old = 0

    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        i, j = bits

        is_i_new = i in new_indices
        is_j_new = j in new_indices

        if not is_i_new and not is_j_new:
            # old-old
            if c < 0:
                F_old_old += abs(c)
            else:
                P_old_old += c
        elif is_i_new and is_j_new:
            # new-new
            if c < 0:
                F_new_new += abs(c)
            else:
                P_new_new += c
        else:
            # new-old
            if c < 0:
                F_new_old_neg += abs(c)
            else:
                F_new_old_pos += c

    total_frustration = F_old_old + F_new_new + F_new_old_neg
    total_synergy = P_old_old + P_new_new + F_new_old_pos

    return {
        'n': n, 'm': m,
        'n_old': len(old_indices), 'n_new': len(new_indices),
        'F_old_old': F_old_old, 'F_new_new': F_new_new,
        'F_new_old': F_new_old_neg, 'P_new_old': F_new_old_pos,
        'P_old_old': P_old_old,
        'total_F': total_frustration, 'total_P': total_synergy,
    }

print("FRUSTRATION PROPAGATION: RECURSIVE DECOMPOSITION")
print("=" * 70)

prev_F = 0
for n in range(4, 8):
    d = frustration_decomposition(n)

    print(f"\n  n={d['n']}: m={d['m']} tiles ({d['n_old']} old + {d['n_new']} new)")
    print(f"    FRUSTRATION (negative couplings):")
    print(f"      Old-Old:  {d['F_old_old']:6d}  (inherited from n-1)")
    print(f"      New-New:  {d['F_new_new']:6d}  (among new tiles)")
    print(f"      New-Old:  {d['F_new_old']:6d}  (new frustrating old)")
    print(f"      TOTAL:    {d['total_F']:6d}")
    print(f"    SYNERGY (positive couplings):")
    print(f"      Old-Old:  {d['P_old_old']:6d}")
    print(f"      New-Old:  {d['P_new_old']:6d}")
    print(f"    RECURSIVE CHECK:")
    print(f"      F(n-1) inherited = {prev_F}")
    print(f"      F(n) - F(n-1) = {d['total_F'] - prev_F} "
          f"(= new-new + new-old = {d['F_new_new'] + d['F_new_old']})")
    match = (d['total_F'] - prev_F) == (d['F_new_new'] + d['F_new_old'])
    print(f"      Recursion exact? {match}")

    prev_F = d['total_F']

# THE FORMULA for new-new frustration
print(f"\n{'='*70}")
print(f"NEW-NEW FRUSTRATION FORMULA")
print(f"{'='*70}")

print(f"""
  New tiles at n+1: (n+1, 1), (n+1, 2), ..., (n+1, n-1).
  All share upper vertex n+1 → ALL pairs are same-end.
  Number of pairs: C(n-1, 2).

  For same-end pairs, c_{{ij}} = -2^{{max(1, |s_i - s_j| - 1)}} (THM-302).
  Tiles (n+1, y) have skip = n+1-y.

  F_new-new = sum over pairs (i,j) of 2^{{max(1, |s_i-s_j|-1)}}.
""")

for n in range(4, 9):
    # Compute F_new-new from formula
    new_tiles = [(n, y) for y in range(1, n-1)]
    total = 0
    for i in range(len(new_tiles)):
        for j in range(i+1, len(new_tiles)):
            s1 = new_tiles[i][0] - new_tiles[i][1]
            s2 = new_tiles[j][0] - new_tiles[j][1]
            c = 2**max(1, abs(s1-s2)-1)
            total += c
    print(f"  n={n}: {len(new_tiles)} new tiles, "
          f"C({len(new_tiles)},2)={comb(len(new_tiles),2)} pairs, "
          f"F_new-new = {total}")

# THE FORMULA for new-old frustration
print(f"\n{'='*70}")
print(f"NEW-OLD FRUSTRATION FORMULA")
print(f"{'='*70}")

print(f"""
  New tile (n+1, y) interacts with old tile (x, z):
  - Same lower vertex: y = z → same-end → c < 0
  - Same upper vertex: n+1 = x → impossible (x ≤ n < n+1)
  - Cross-end: n+1 = z or y = x → cross-end → c > 0
  - Disjoint: otherwise → c > 0

  So new-old FRUSTRATION comes ONLY from shared LOWER vertex.
  For new tile (n+1, y): old tiles sharing lower vertex y are
  (x, y) for x = y+2, ..., n. That's n-y-1 old tiles.

  For each pair (new, old) with same lower y:
  c = -2^{{max(1, |s_new - s_old| - 1)}} where s_new = n+1-y, s_old = x-y.
""")

for n in range(4, 8):
    total_new_old = 0
    for y in range(1, n):  # new tile (n+1, y), but in our notation n goes to n
        # Actually, at step n, new tiles have x = n
        s_new = n - y  # skip of new tile (n, y)
        # Old tiles sharing lower vertex y: (x, y) for x = y+2, ..., n-1
        for x in range(y+2, n):
            s_old = x - y
            c = 2**max(1, abs(s_new - s_old) - 1)
            total_new_old += c
    print(f"  n={n}: F_new-old (formula) = {total_new_old}")

# RECURSION TABLE
print(f"\n{'='*70}")
print(f"COMPLETE FRUSTRATION RECURSION TABLE")
print(f"{'='*70}")
print(f"  {'n':>3s} {'F(n)':>8s} {'F_nn':>8s} {'F_no':>8s} {'ΔF':>8s} {'F/m':>8s}")

prev_F = 0
for n in range(4, 8):
    d = frustration_decomposition(n)
    delta_F = d['total_F'] - prev_F
    f_per_tile = d['total_F'] / d['m'] if d['m'] > 0 else 0
    print(f"  {n:3d} {d['total_F']:8d} {d['F_new_new']:8d} {d['F_new_old']:8d} "
          f"{delta_F:8d} {f_per_tile:8.2f}")
    prev_F = d['total_F']

# THE FRUSTRATION THRESHOLD revisited
print(f"\n{'='*70}")
print(f"FRUSTRATION THRESHOLD vs LAYER")
print(f"{'='*70}")

for n in [5, 6, 7]:
    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # For each tile k, compute its "frustration load" = sum of |c_{kj}| for same-end j
    print(f"\n  n={n}:")
    for k in range(m):
        tile = tiles[k]
        skip = tile[0] - tile[1]
        frust_load = 0
        synergy_load = 0
        for S, c in coeffs.items():
            if bin(S).count('1') != 2: continue
            bits = [i for i in range(m) if S & (1 << i)]
            if k not in bits: continue
            if c < 0:
                frust_load += abs(c)
            else:
                synergy_load += c

        c_k = 2**(skip-1)
        # Threshold: the gradient flips when frust_load > c_k (roughly)
        threshold_ratio = frust_load / c_k if c_k > 0 else 0
        print(f"    tile {tile} skip={skip}: c_k={c_k:4d}, "
              f"frust_load={frust_load:5d}, synergy={synergy_load:5d}, "
              f"frust/c_k={threshold_ratio:.2f}")
