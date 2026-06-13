"""
sum_H_formula.py — opus-2026-04-04-S5

DERIVE: Σ_{t} H(t) = Σ_{valid σ} 2^{m - tile_arcs(σ)}

where "valid σ" means: no ascending-by-1 step in the HP ordering.
(Because v→(v+1) is a reverse base-path arc, which never exists.)

A valid HP σ = (v₁,...,v_n) has:
  - bp(σ) = base-path arcs = descents by 1 (v_i = v_{i+1}+1)
  - tile_arcs(σ) = n-1-bp(σ) = non-consecutive arcs
  - 0 reverse base-path arcs = no ascending by 1

The count of valid HPs with bp base-path arcs = N(n, bp)
And Σ H = Σ_{bp=0}^{n-1} N(n,bp) · 2^{m-n+1+bp}

N(n, bp) counts permutations of {0,...,n-1} with exactly bp descents-by-1
and ZERO ascents-by-1. These are "succession-free with bp selected descents."
"""
from itertools import permutations
from math import comb
import numpy as np

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def count_valid_hp_by_bp(n):
    """Count valid HPs (no ascending-by-1) by base-path arc count."""
    bp_count = {}
    for perm in permutations(range(n)):
        bp = 0
        valid = True
        for i in range(n-1):
            u, v = perm[i], perm[i+1]
            if v == u + 1:  # ascending by 1 = reverse base path
                valid = False
                break
            if u == v + 1:  # descending by 1 = base path arc
                bp += 1
        if valid:
            bp_count[bp] = bp_count.get(bp, 0) + 1
    return bp_count

def verify_sum_H_formula():
    print("Σ H FORMULA VERIFICATION")
    print("=" * 60)

    for n in range(3, 8):
        tiles = make_tiles(n)
        m = len(tiles)

        # Count valid HPs by bp
        bp_dist = count_valid_hp_by_bp(n)
        total_valid = sum(bp_dist.values())

        # Compute Σ H from formula
        formula_sum = 0
        for bp, count in sorted(bp_dist.items()):
            contrib = count * (2 ** (m - n + 1 + bp))
            formula_sum += contrib

        # Compute actual Σ H
        actual_sum = 0
        for bits in range(1 << m):
            T = np.zeros((n, n), dtype=np.int8)
            for k in range(n-1):
                T[k+1][k] = 1
            for i, (x, y) in enumerate(tiles):
                a, b = x-1, y-1
                if bits & (1 << i):
                    T[b][a] = 1
                else:
                    T[a][b] = 1
            # Count HP
            dp = {}
            for v in range(n):
                dp[(1 << v, v)] = 1
            full = (1 << n) - 1
            h = 0
            for mask in range(1, 1 << n):
                for v in range(n):
                    cnt = dp.get((mask, v), 0)
                    if cnt == 0: continue
                    if mask == full: h += cnt; continue
                    for u in range(n):
                        if mask & (1 << u): continue
                        if T[v][u]:
                            key = (mask | (1 << u), u)
                            dp[key] = dp.get(key, 0) + cnt
            actual_sum += h

        match = formula_sum == actual_sum

        print(f"\nn={n}: m={m}")
        print(f"  Valid HPs: {total_valid} (out of {n}! = {int(np.prod(range(1,n+1)))})")
        print(f"  BP distribution: {dict(sorted(bp_dist.items()))}")
        print(f"  Formula: Σ N(bp)·2^(m-n+1+bp) = {formula_sum}")
        print(f"  Actual Σ H = {actual_sum}")
        print(f"  {'✓ MATCH' if match else '�� MISMATCH'}")

        # Mean H
        mean_h = actual_sum / (1 << m)
        print(f"  E[H] = {mean_h:.6f}")

    # Now extend the formula to compute Σ H at larger n (without enumerating tilings!)
    print(f"\n{'='*60}")
    print(f"EXTENDED Σ H via BP FORMULA")
    print(f"{'='*60}")

    for n in range(3, 12):
        if n <= 10:  # Permutation enumeration feasible up to ~10
            bp_dist = count_valid_hp_by_bp(n)
        else:
            bp_dist = None

        m = comb(n-1, 2)

        if bp_dist:
            total = sum(count * (2 ** (m - n + 1 + bp))
                       for bp, count in bp_dist.items())
            mean_h = total / (2 ** m)
            total_valid = sum(bp_dist.values())
            print(f"  n={n}: Σ H = {total}, E[H] = {mean_h:.4f}, "
                  f"valid HPs = {total_valid}")
        else:
            print(f"  n={n}: (too large for permutation enumeration)")

    # The VALID HP COUNTS sequence: total succession-free permutations
    print(f"\nSuccession-free permutations (= valid HP count):")
    for n in range(3, 12):
        if n <= 10:
            bp_dist = count_valid_hp_by_bp(n)
            total = sum(bp_dist.values())
            print(f"  n={n}: {total}")

    # The Σ H sequence
    print(f"\nΣ H sequence (total HP count over all tilings):")
    for n in range(3, 11):
        bp_dist = count_valid_hp_by_bp(n)
        m = comb(n-1, 2)
        total = sum(count * (2 ** (m - n + 1 + bp))
                   for bp, count in bp_dist.items())
        print(f"  n={n}: {total}")

    # E[H] sequence
    print(f"\nE[H] = mean HP count with fixed base path:")
    for n in range(3, 11):
        bp_dist = count_valid_hp_by_bp(n)
        m = comb(n-1, 2)
        total = sum(count * (2 ** (m - n + 1 + bp))
                   for bp, count in bp_dist.items())
        mean = total / (2 ** m)
        print(f"  n={n}: E[H] = {mean}")

if __name__ == "__main__":
    verify_sum_H_formula()
