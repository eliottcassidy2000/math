"""
h_multilinear_analysis.py — opus-2026-04-04-S1

Deep analysis of the multilinear coefficients:
1. Apex effect decomposition
2. Disjoint pair formula search
3. Order-k coefficient structure
"""
import numpy as np
from collections import defaultdict
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def apex_effect_decomposition(n_values=[5, 6, 7]):
    """Decompose the apex effect E[H(apex=1) - H(apex=0)] by order."""
    print("=" * 60)
    print("APEX EFFECT DECOMPOSITION")
    print("=" * 60)

    for n in n_values:
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        # Find apex tile index
        apex_idx = None
        for i, (x,y) in enumerate(tiles):
            if x == n and y == 1:
                apex_idx = i
                break

        if apex_idx is None:
            print(f"n={n}: no apex tile")
            continue

        # Decompose apex effect by order
        total = 0
        print(f"\nn={n}: apex tile ({n},1), skip={n-1}, linear coeff = 2^{n-2} = {2**(n-2)}")
        for order in range(m+1):
            order_sum = sum(c for S, c in coeffs.items()
                          if (S & (1 << apex_idx)) and bin(S).count('1') == order)
            contribution = order_sum / (2 ** (order - 1)) if order > 0 else 0
            total += contribution
            if order_sum != 0:
                print(f"  Order {order}: sum(c_S) = {order_sum}, "
                      f"contribution = {order_sum}/{2**(order-1)} = {contribution}")

        expected = 2**(n-3)
        print(f"  TOTAL apex effect = {total}")
        print(f"  Expected 2^(n-3) = {expected}")
        print(f"  Correction = {total - expected} = {total - expected} ")
        if total == expected:
            print(f"  EXACT MATCH!")
        else:
            # Express as fraction
            from fractions import Fraction
            f = Fraction(total).limit_denominator(1000)
            print(f"  As fraction: {f}")

def disjoint_pair_analysis(n=7):
    """Analyze the formula for disjoint tile pair coefficients."""
    print("\n" + "=" * 60)
    print(f"DISJOINT PAIR COEFFICIENT ANALYSIS (n={n})")
    print("=" * 60)

    H, tiles, m = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # Extract order-2 disjoint pairs
    pairs = []
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = []
        for i in range(m):
            if S & (1 << i):
                bits.append(i)
        i1, i2 = bits
        t1, t2 = tiles[i1], tiles[i2]
        x1, y1 = t1
        x2, y2 = t2

        # Check if vertex-disjoint
        verts1 = {x1, y1}
        verts2 = {x2, y2}
        if verts1 & verts2:
            continue  # shared vertex, skip

        s1, s2 = x1-y1, x2-y2
        # Geometric properties
        # Position on staircase: (x-y-1, y) for tile (x,y) → row=skip-1, col=y
        r1, c1_pos = s1-1, y1
        r2, c2_pos = s2-1, y2

        # Manhattan distance between tiles on staircase
        manhattan = abs(r1-r2) + abs(c1_pos-c2_pos)

        # "Sum of skips" and "product"
        sum_s = s1 + s2
        prod_s = s1 * s2
        expected_prod = 2**(s1-1) * 2**(s2-1)  # product of linear coefficients

        pairs.append({
            'tiles': (t1, t2),
            'skips': (s1, s2),
            'c2': c,
            'sum_s': sum_s,
            'prod_s': prod_s,
            'expected_prod': expected_prod,
            'manhattan': manhattan,
            'staircase_pos': ((r1,c1_pos), (r2,c2_pos)),
            'ratio': c / expected_prod if expected_prod else 0,
        })

    # Sort by coefficient value
    pairs.sort(key=lambda p: p['c2'])

    print(f"\nTotal disjoint pairs: {len(pairs)}")
    print(f"\n{'Tiles':<25} {'Skips':<10} {'c₂':>6} {'2^(s1-1)·2^(s2-1)':>18} {'Ratio':>8} {'Manhattan':>9} {'Staircase pos'}")
    print("-" * 110)
    for p in pairs:
        t1, t2 = p['tiles']
        print(f"({t1[0]},{t1[1]})-({t2[0]},{t2[1]})  "
              f"{str(p['skips']):<10} {p['c2']:>6} {p['expected_prod']:>18} "
              f"{p['ratio']:>8.4f} {p['manhattan']:>9} {p['staircase_pos']}")

    # Try to find the formula
    print("\n\n=== FORMULA SEARCH ===")

    # Hypothesis 1: c₂ = f(s1, s2, some geometric property)
    # Group by (s1, s2) sorted
    by_skips = defaultdict(list)
    for p in pairs:
        key = tuple(sorted(p['skips']))
        by_skips[key].append(p)

    print("\nGrouped by skip pair:")
    for key in sorted(by_skips.keys()):
        entries = by_skips[key]
        values = [e['c2'] for e in entries]
        print(f"  skips {key}: values {sorted(set(values))}, count {len(entries)}")
        if len(set(values)) > 1:
            # Multiple values — what distinguishes them?
            for e in sorted(entries, key=lambda x: x['c2']):
                t1, t2 = e['tiles']
                # Check if tiles are "adjacent" on staircase
                (r1,c1), (r2,c2) = e['staircase_pos']
                adj = abs(r1-r2) <= 1 and abs(c1-c2) <= 1
                # Check various metrics
                overlap = min(t1[0], t2[0]) - max(t1[1], t2[1])
                print(f"    ({t1[0]},{t1[1]})-({t2[0]},{t2[1]}): c₂={e['c2']}, "
                      f"overlap={overlap}, manhattan={e['manhattan']}")

    # Hypothesis 2: c₂ depends on the "overlap" = min(x1,x2) - max(y1,y2)
    # (how many vertices are "between" the two arcs)
    print("\n\nGrouped by overlap (min(x1,x2) - max(y1,y2)):")
    by_overlap = defaultdict(list)
    for p in pairs:
        (x1,y1), (x2,y2) = p['tiles']
        overlap = min(x1,x2) - max(y1,y2)
        by_overlap[overlap].append(p)

    for ov in sorted(by_overlap.keys()):
        entries = by_overlap[ov]
        values = [e['c2'] for e in entries]
        skips = [(e['skips']) for e in entries]
        print(f"  overlap={ov}: values {sorted(set(values))}, count={len(entries)}")

    # Hypothesis 3: c₂ = 2^(s1+s2-2) * f(overlap)
    # where overlap = min(x1,x2) - max(y1,y2)
    print("\n\nRatio c₂ / 2^(s1+s2-2):")
    for p in sorted(pairs, key=lambda p: (p['skips'], p['c2'])):
        (x1,y1), (x2,y2) = p['tiles']
        overlap = min(x1,x2) - max(y1,y2)
        base = 2**(p['skips'][0] + p['skips'][1] - 2)
        ratio = p['c2'] / base
        from fractions import Fraction
        frac = Fraction(p['c2'], base)
        print(f"  ({x1},{y1})-({x2},{y2}) s={p['skips']} ov={overlap}: "
              f"c₂={p['c2']}, 2^(s1+s2-2)={base}, ratio={frac}")

def order_k_analysis(n=7):
    """Analyze the structure of order-k coefficients."""
    print("\n" + "=" * 60)
    print(f"ORDER-k COEFFICIENT STRUCTURE (n={n})")
    print("=" * 60)

    H, tiles, m = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    max_order = max(bin(S).count('1') for S in coeffs.keys())

    for order in range(max_order+1):
        entries = [(S, c) for S, c in coeffs.items() if bin(S).count('1') == order]
        if not entries:
            continue

        values = [c for _, c in entries]
        print(f"\nOrder {order}: {len(entries)} nonzero, values in [{min(values)}, {max(values)}]")

        # Check divisibility by 2
        min_v2 = min(v2(abs(c)) for _, c in entries)
        print(f"  Min 2-adic valuation: {min_v2}")
        print(f"  All divisible by 2^{min_v2} = {2**min_v2}")

        # Distribution of values
        val_counts = defaultdict(int)
        for _, c in entries:
            val_counts[c] += 1
        print(f"  Value distribution: {dict(sorted(val_counts.items()))}")

def v2(n):
    """2-adic valuation."""
    if n == 0: return float('inf')
    n = abs(n)
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v

if __name__ == "__main__":
    apex_effect_decomposition()
    disjoint_pair_analysis(n=7)
    order_k_analysis(n=7)
