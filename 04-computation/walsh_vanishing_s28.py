#!/usr/bin/env python3
"""
Why do high Walsh orders vanish for H on tiling space?
opus-2026-04-03-S28

PUZZLE: The Walsh-Fourier spectrum of H(tiling) shows:
  n=4: orders 0,1,2 only (order 3 = 0)
  n=5: orders 0,1,2,3,4 (orders 5,6 = 0)
  n=6: orders 0,1,2,3,4 (orders 5-10 = 0)

WAIT: At n=5, orders 1 and 3 are nonzero! These are ODD orders.
But the earlier info-geometry reflection said odd orders are zero...

Let me recheck. The info-geometry was about the FULL arc space {0,1}^{C(n,2)}.
On the full arc space, T -> T^op flips ALL bits, so H (being complement-invariant)
must have zero odd-order Walsh coefficients.

But on TILING space {0,1}^m, the tiling complement (flip all tile bits) is
NOT tournament complement. So H on tiling space CAN have odd Walsh orders.

The question: why does H max out at Walsh order 4 (on tiling space)?

Hypothesis: H is a polynomial of degree <= 4 in the tiling bits.

From the OCF: H(T) = I(Omega(T), 2) = sum of 2^k * alpha_k(T)
where alpha_k = number of independent sets of size k in the conflict graph.
The independent sets of odd cycles are determined by the arc directions.

Each arc is a LINEAR function of the tiling bit (or the base-path arc, which is constant).
A 3-cycle involves 3 arcs, each of which is degree 1 in some tile bit.
The indicator of a 3-cycle being directed is a PRODUCT of 3 arc indicators,
which is degree 3 in the tiling bits.
Two vertex-disjoint 3-cycles: product of two degree-3 indicators = degree 6.

But H involves the independence polynomial, which includes pairs of
vertex-disjoint cycles. At n=5 you can have two disjoint 3-cycles?
No — 5 vertices, two disjoint 3-cycles need 6 vertices. So at n=5,
alpha_2 can only come from pairs involving 5-cycles, which are rare.

Actually, a single 3-cycle: its indicator function in terms of arcs is
  C_directed = T(a,b) * T(b,c) * T(c,a) = product of 3 arc bits.
But each arc is either a tile bit or a base-path constant.
If all 3 arcs are tiles, the indicator is degree 3 in tile bits.
If some arcs are base-path, the degree is reduced.

Let me compute the degree of H as a polynomial in the tile bits.
"""

from collections import defaultdict
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

for n in range(3, 8):
    m = (n-1)*(n-2)//2

    print(f"\n{'='*70}")
    print(f"WALSH SPECTRUM DEEP ANALYSIS: n={n}, m={m}")
    print(f"{'='*70}")

    if 2**m > 40000:
        print(f"Too large for exhaustive Walsh computation")
        continue

    # Compute all H values
    H_vals = [0] * (2**m)
    for bits in range(2**m):
        b = [(bits >> i) & 1 for i in range(m)]
        T = tournament_from_bits(n, b)
        H_vals[bits] = count_hp(T, n)

    # Full Walsh-Fourier transform
    walsh = {}
    for S in range(2**m):
        total = 0
        for bits in range(2**m):
            sign = (-1) ** bin(S & bits).count('1')
            total += sign * H_vals[bits]
        walsh[S] = total / (2**m)

    # Group by order
    walsh_by_order = defaultdict(list)
    for S, coeff in walsh.items():
        order = bin(S).count('1')
        if abs(coeff) > 1e-10:
            walsh_by_order[order].append((S, coeff))

    print(f"\nWalsh coefficients by order:")
    for order in sorted(walsh_by_order):
        coeffs = walsh_by_order[order]
        print(f"\n  Order {order}: {len(coeffs)} nonzero coefficients")
        for S, coeff in sorted(coeffs, key=lambda x: -abs(x[1]))[:10]:
            bits_str = f"{S:0{m}b}"
            support = [i for i in range(m) if (S >> i) & 1]
            print(f"    S={bits_str} support={support} coeff={coeff:.6f}")
        if len(coeffs) > 10:
            print(f"    ... and {len(coeffs) - 10} more")

    # Maximum nonzero order
    max_order = max(walsh_by_order.keys()) if walsh_by_order else 0
    print(f"\n  Maximum nonzero Walsh order: {max_order}")
    print(f"  Number of tile bits: {m}")
    print(f"  Polynomial degree of H in tile bits: {max_order}")

    # Check: is there a pattern to WHICH order-k subsets have nonzero coefficients?
    # At n=5, m=6: which specific pairs/triples/quads have nonzero Walsh?
    if order_2_coeffs := walsh_by_order.get(2, []):
        from itertools import combinations
        tiles = []
        for y in range(1, n-1):
            for x in range(n, y+1, -1):
                if x - y >= 2:
                    tiles.append((x, y))

        print(f"\n  Order-2 analysis (pairwise interactions):")
        for S, coeff in sorted(order_2_coeffs, key=lambda x: -abs(x[1])):
            support = [i for i in range(m) if (S >> i) & 1]
            t1, t2 = tiles[support[0]], tiles[support[1]]
            skip1, skip2 = t1[0]-t1[1], t2[0]-t2[1]
            # Do these tiles share a vertex?
            shared = set(t1) & set(t2)
            share_str = f"share {shared}" if shared else "disjoint"
            print(f"    {chr(65+support[0])}{chr(65+support[1])} = "
                  f"({t1[0]},{t1[1]})×({t2[0]},{t2[1]}) "
                  f"skip=({skip1},{skip2}) {share_str}: coeff={coeff:+.6f}")

    # KEY INSIGHT: Can we express H as an EXPLICIT multilinear polynomial?
    # H = sum_{S} walsh(S) * prod_{i in S} (2*x_i - 1)
    # = sum_{S} walsh(S) * (-1)^{|S|} * prod_{i in S} (1 - 2*x_i)
    # In terms of x_i (not chi_i = 2x_i - 1):
    # Using inclusion-exclusion, this gives the unique multilinear polynomial.

    # Let me compute the multilinear polynomial coefficients directly
    # f(x_1, ..., x_m) = sum_{S} c_S * prod_{i in S} x_i
    # where c_S = sum_{T subset S} (-1)^{|S|-|T|} * H(T)  (Möbius inversion)

    # Actually, let me use the standard method:
    # The multilinear polynomial is uniquely determined by its values at {0,1}^m.
    # The coefficient of x_S = prod_{i in S} x_i is:
    # c_S = sum_{T subset of S} (-1)^{|S\T|} * f(e_T)
    # where e_T is the indicator of T.

    # For small m, just do it by Möbius:
    if m <= 10:
        print(f"\n  Multilinear polynomial H(x_1,...,x_{m}) coefficients:")
        print(f"  (Only showing nonzero, up to 20)")

        nonzero_polys = []
        for S in range(2**m):
            order = bin(S).count('1')
            if order > max_order + 1:
                continue
            # c_S = sum_{T subset S} (-1)^{|S\T|} H(e_T)
            coeff = 0
            for T_bits in range(2**m):
                if T_bits & S != T_bits:  # T not subset of S
                    continue
                diff = bin(S & ~T_bits).count('1')
                coeff += ((-1) ** diff) * H_vals[T_bits]
            if abs(coeff) > 1e-10:
                nonzero_polys.append((S, coeff, order))

        nonzero_polys.sort(key=lambda x: (x[2], -abs(x[1])))
        for S, coeff, order in nonzero_polys[:30]:
            support = [i for i in range(m) if (S >> i) & 1]
            tile_labels = [chr(65+i) for i in support]
            print(f"    {''.join(tile_labels) if tile_labels else '1'} (order {order}): {coeff:+.0f}")

        print(f"\n  Total nonzero terms: {len(nonzero_polys)}")
        by_order_count = defaultdict(int)
        for _, _, order in nonzero_polys:
            by_order_count[order] += 1
        for order in sorted(by_order_count):
            count = by_order_count[order]
            total_possible = 1
            from math import comb
            total_possible = comb(m, order)
            print(f"    Order {order}: {count}/{total_possible} terms nonzero")
