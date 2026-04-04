"""
disjoint_and_cubic_deep.py — opus-2026-04-04-S3

Deep investigation of:
1. The disjoint pair formula (ratio c_{ij} / 2^{s1+s2-2})
2. Cubic coefficient structure via OCF
3. The "interaction range" — how far apart can interacting tiles be?
4. Degree-sum polynomial H(p) = sum_{d} S_d p^d — roots and factorization
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

# =================================================================
# PART 1: Disjoint pair formula — identify the key parameter
# =================================================================
def disjoint_formula_crack(n):
    """
    For disjoint pairs, we know ratio = c / 2^{s1+s2-2} is NOT always 1/4.
    The anomalous values: 1/2, 5/16, 3/8, 17/64, 13/64.

    Try: ratio depends on the "interaction strength" = number of odd cycles
    passing through both tiles with POSITIVE chi contribution.

    From THM-287: c_{ij} = 2*A(i,j) + 4*B(i,j)
    A = single-cycle contribution, B = disjoint-pair contribution.
    For disjoint tiles: A(i,j) is the number of single cycles through both tiles
    (with appropriate sign), and B(i,j) counts independent pairs.

    Let's compute A and B separately and see what determines the ratio.
    """
    print(f"\n{'='*60}")
    print(f"DISJOINT PAIR FORMULA ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # For each disjoint pair, compute the number of directed odd cycles
    # that pass through BOTH tiles' arcs

    # First, enumerate all directed odd cycles (3-cycles and 5-cycles)
    def enumerate_cycles_through_pair(tile1, tile2, n):
        """Count directed odd cycles using BOTH tile arcs (forward direction)."""
        x1, y1 = tile1
        x2, y2 = tile2
        # Convert to 0-indexed
        a1, b1 = x1-1, y1-1  # tile 1 endpoints
        a2, b2 = x2-1, y2-1  # tile 2 endpoints

        cycles_3 = 0
        cycles_5 = 0

        # 3-cycles through both tiles
        # A 3-cycle on vertices {u,v,w} uses 3 arcs.
        # For it to use both tiles, both (a1,b1) and (a2,b2) must be among the arcs.
        # 3-cycle on vertices of both tiles: need 3 vertices from {a1,b1,a2,b2}.
        # If tiles are disjoint (4 distinct vertices), can't have a 3-cycle through both.
        # UNLESS the 3-cycle shares some extra vertex... but tiles are disjoint,
        # so {a1,b1} and {a2,b2} are disjoint (4 distinct vertices).
        # A 3-cycle uses 3 vertices, so can't use all 4. Can use at most 3.
        # To use both tile arcs, need a1->b1 and a2->b2 (or reversed) in the cycle.
        # This requires a1,b1,a2,b2 to have at most 3 distinct values...
        # but they're 4 distinct values (disjoint tiles). So NO 3-cycle through both.

        # Actually wait: tiles are "disjoint" in the sense of no shared vertex label.
        # But a tile (x,y) covers the ARC from x to y, not all intermediate vertices.
        # Disjoint tiles: {x1,y1} ∩ {x2,y2} = ∅.
        # A 3-cycle uses 3 arcs on 3 vertices. To use arcs (a1,b1) and (a2,b2),
        # the 3-cycle needs a1,b1,a2,b2 to be at most 3 distinct vertices.
        # Since tiles are disjoint (4 distinct), this is impossible.
        # So: ALL c_{ij} for disjoint pairs come from 5-cycles or higher!

        # 5-cycles through both tiles
        # A 5-cycle uses 5 arcs on 5 vertices.
        # It must include arcs involving (a1,b1) and (a2,b2).
        # The cycle visits a1, b1, a2, b2 plus one more vertex.
        for v in range(n):
            if v in {a1, b1, a2, b2}:
                continue
            # Try all 5-cycles on {a1, b1, a2, b2, v} that use both tile arcs
            verts = [a1, b1, a2, b2, v]
            for perm in permutations(verts):
                # Check if this is a valid directed cycle using both tile arcs
                # The cycle is perm[0]→perm[1]→...→perm[4]→perm[0]
                arcs = [(perm[i], perm[(i+1)%5]) for i in range(5)]
                # Check both tile arcs are present
                has_t1 = (a1, b1) in arcs or (b1, a1) in arcs
                has_t2 = (a2, b2) in arcs or (b2, a2) in arcs
                if has_t1 and has_t2:
                    # Check canonical (lexicographically smallest rotation)
                    rotations = [perm[i:]+perm[:i] for i in range(5)]
                    if perm == min(rotations):
                        cycles_5 += 1

        return cycles_3, cycles_5

    # This is slow for n=7, so only do it for n≤6
    if n <= 6:
        print(f"  Computing cycle anatomy for disjoint pairs...")
        for S, c in sorted(coeffs.items(), key=lambda x: abs(x[1])):
            if bin(S).count('1') != 2: continue
            bits = [i for i in range(m) if S & (1 << i)]
            t1, t2 = tiles[bits[0]], tiles[bits[1]]
            if {t1[0],t1[1]} & {t2[0],t2[1]}: continue  # skip non-disjoint

            s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
            c3, c5 = enumerate_cycles_through_pair(t1, t2, n)

            # The "range overlap" of the tiles
            range1 = set(range(t1[1], t1[0]+1))
            range2 = set(range(t2[1], t2[0]+1))
            range_ov = len(range1 & range2)

            print(f"    {t1}-{t2} s=({s1},{s2}): c={c}, "
                  f"3-cyc={c3}, 5-cyc={c5}, "
                  f"range_ov={range_ov}, ratio={Fraction(c, 2**(s1+s2-2))}")

# =================================================================
# PART 2: The interaction range
# =================================================================
def interaction_range(n):
    """What is the maximum 'distance' between interacting tiles?"""
    print(f"\n{'='*60}")
    print(f"INTERACTION RANGE n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    max_deg = max(bin(S).count('1') for S in coeffs.keys())

    for deg in range(2, max_deg + 1):
        deg_coeffs = {S: c for S, c in coeffs.items() if bin(S).count('1') == deg}
        if not deg_coeffs:
            continue

        # For each degree-d coefficient, compute the "span" of the tile set
        # Span = max(x) - min(y) over tiles in the set
        spans = []
        for S in deg_coeffs:
            bits = [i for i in range(m) if S & (1 << i)]
            tile_set = [tiles[i] for i in bits]
            max_x = max(t[0] for t in tile_set)
            min_y = min(t[1] for t in tile_set)
            span = max_x - min_y
            spans.append(span)

        print(f"  Degree {deg}: {len(deg_coeffs)} coefficients, "
              f"span range [{min(spans)}, {max(spans)}]")

# =================================================================
# PART 3: Degree-sum polynomial deep analysis
# =================================================================
def degree_sum_polynomial():
    """
    The univariate H(p) = sum_d S_d p^d where S_d = sum_{|S|=d} c_S.

    Known: S_0 = 1, S_1 = sum 2^{skip-1}.

    Can we find formulas for S_2, S_3, ...?
    """
    print(f"\n{'='*60}")
    print(f"DEGREE-SUM POLYNOMIAL ANALYSIS")
    print(f"{'='*60}")

    degree_sums = {}

    for n in [3, 4, 5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        if m == 0:
            continue
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        by_deg = defaultdict(int)
        for S, c in coeffs.items():
            by_deg[bin(S).count('1')] += c

        max_d = max(by_deg.keys()) if by_deg else 0
        sums = [by_deg.get(d, 0) for d in range(max_d + 1)]
        degree_sums[n] = sums

    print(f"  Degree sums S_d by n:")
    for n, sums in sorted(degree_sums.items()):
        print(f"    n={n}: {sums}")

    # Analyze each degree separately
    print(f"\n  S_0 = 1 always (= H(transitive))")

    S1 = {n: sums[1] if len(sums) > 1 else 0 for n, sums in degree_sums.items()}
    print(f"\n  S_1 values: {S1}")
    # S_1 = sum of 2^{skip-1} over all tiles
    # = sum_{y=1}^{n-2} sum_{x=y+2}^{n} 2^{x-y-1}
    # = sum_{s=2}^{n-1} (n-s) * 2^{s-1}   where s = skip = x-y
    for n in range(3, 10):
        s1_formula = sum((n-s) * 2**(s-1) for s in range(2, n))
        print(f"    n={n}: S_1 = {s1_formula} = {sum((n-s)*2**(s-1) for s in range(2,n))}")

    # S_1 = sum_{s=2}^{n-1} (n-s)*2^{s-1}
    # = n * sum 2^{s-1} - sum s*2^{s-1}  for s=2..n-1
    # = n * (2^{n-1} - 2) - [sum s*2^{s-1} for s=2..n-1]
    # sum s*2^{s-1} for s=1..N = (N-1)*2^N + 1
    # So for s=2..n-1: [(n-2)*2^{n-1}+1] - [1*2^0+1] = (n-2)*2^{n-1} + 1 - 2 = (n-2)*2^{n-1} - 1
    # S_1 = n*(2^{n-1}-2) - (n-2)*2^{n-1} + 1 = n*2^{n-1} - 2n - (n-2)*2^{n-1} + 1
    # = 2*2^{n-1} - 2n + 1 = 2^n - 2n + 1
    print(f"\n  Formula: S_1 = 2^n - 2n + 1")
    for n in range(3, 10):
        print(f"    n={n}: 2^n-2n+1 = {2**n - 2*n + 1}, actual S_1 = {sum((n-s)*2**(s-1) for s in range(2,n))}")

    # S_2 values
    S2 = {n: sums[2] if len(sums) > 2 else 0 for n, sums in degree_sums.items()}
    print(f"\n  S_2 values: {S2}")
    # n=4: -4, n=5: 0, n=6: 92, n=7: 616
    # Differences: 4, 92, 524
    # Second diff: 88, 432
    # Not obviously polynomial. Let me compute for n=8.

    # Check: S_2 = sum of all quadratic coefficients
    # = sum of same-end + sum of cross-end + sum of disjoint
    for n in [5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        same_end_sum = 0
        cross_end_sum = 0
        disjoint_sum = 0
        for S, c in coeffs.items():
            if bin(S).count('1') != 2: continue
            bits = [i for i in range(m) if S & (1 << i)]
            t1, t2 = tiles[bits[0]], tiles[bits[1]]
            x1,y1 = t1; x2,y2 = t2
            if x1==x2 or y1==y2:
                same_end_sum += c
            elif x1==y2 or y1==x2:
                cross_end_sum += c
            else:
                disjoint_sum += c

        print(f"    n={n}: S_2={same_end_sum+cross_end_sum+disjoint_sum}, "
              f"same-end={same_end_sum}, cross-end={cross_end_sum}, disjoint={disjoint_sum}")

    # The same-end sum should be computable from the formula c = -2^{max(1,|s1-s2|-1)}
    # For each pair of tiles sharing a vertex:
    for n in [5, 6, 7, 8]:
        tiles_n = make_tiles(n)
        se_sum = 0
        for i in range(len(tiles_n)):
            for j in range(i+1, len(tiles_n)):
                t1, t2 = tiles_n[i], tiles_n[j]
                x1,y1 = t1; x2,y2 = t2
                if x1==x2 or y1==y2:
                    s1, s2 = x1-y1, x2-y2
                    c = -(2**max(1, abs(s1-s2)-1))
                    se_sum += c
        print(f"    n={n}: same-end sum from formula = {se_sum}")

# =================================================================
# PART 4: The H polynomial as generating function for something?
# =================================================================
def polynomial_as_gf():
    """
    Can we interpret H(t) as a generating function?

    H(t) evaluated at t_k = x_k gives a multivariate generating function.

    Key observations:
    - H(0,...,0) = 1 (transitive tournament, 1 HP)
    - H(1,...,1) = H(all-backward tiling)
    - Coefficients are even (except constant)
    - All max-degree coefficients are ±2

    What if we write P(t) = (H(t) - 1) / 2? Then P has integer coefficients.
    P is the "excess HP count" (over the transitive minimum) divided by 2.
    """
    print(f"\n{'='*60}")
    print(f"REDUCED POLYNOMIAL P(t) = (H(t)-1)/2")
    print(f"{'='*60}")

    for n in [5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        P_coeffs = {}
        for S, c in coeffs.items():
            if S == 0: continue
            P_coeffs[S] = c // 2  # All c_S are even for |S|>=1

        by_deg = defaultdict(list)
        for S, c in P_coeffs.items():
            by_deg[bin(S).count('1')].append(c)

        print(f"\n  n={n}: P(t) = (H(t)-1)/2")
        for d in sorted(by_deg.keys()):
            vals = by_deg[d]
            print(f"    d={d}: {len(vals)} terms, sum={sum(vals)}, "
                  f"values={sorted(set(vals))[:10]}")

        # At max degree: P coefficients are all ±1
        max_d = max(by_deg.keys())
        max_vals = set(by_deg[max_d])
        print(f"    Max degree {max_d}: P coefficients = {max_vals}")
        if max_vals == {1, -1} or max_vals == {1} or max_vals == {-1}:
            print(f"    *** P has ±1 coefficients at max degree! ***")

if __name__ == "__main__":
    for n in [5, 6]:
        disjoint_formula_crack(n)

    for n in [5, 6, 7]:
        interaction_range(n)

    degree_sum_polynomial()
    polynomial_as_gf()
