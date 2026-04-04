"""
transfer_matrix_explore.py — opus-2026-04-04-S3

The RECURSIVE PRESERVATION theorem says all multilinear coefficients involving
only "old" tiles (those not involving the highest vertex) are exactly preserved
when n increases. This suggests a LAYERED structure.

IDEA 1: H_n(t) = product over "layers"?
  Each layer adds one new vertex and its associated tiles.
  If H factors as a product of layer contributions, we'd have
  H_n(t) = prod_{v=3}^{n} L_v(t_{new_v}, t_{old_v})

IDEA 2: H_n(t) as a determinant?
  det(I + M(t)) where M is a matrix built from tile variables.
  This would explain the multiplicative structure and the OCF (independence polynomial
  of the conflict graph).

IDEA 3: Connection to the OCF.
  H(T) = I(Omega(T), 2) = sum_{I independent} 2^|I|.
  This IS the partition function of the hard-core lattice gas at fugacity 2.
  Partition functions often admit transfer matrix representations.

IDEA 4: H(t) as graph polynomial on the conflict graph.
  The multilinear structure of H(t) comes from products of chi_C(t).
  Each chi_C is a product of tile factors. The OCF decomposes H into
  cycle contributions. Can we organize this by the "layer" (vertex order)
  in which each cycle lives?

Let's explore computationally.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def compute_layer_contribution(n):
    """
    Compute the "layer contribution" at vertex n.

    If H_n(t) = H_{n-1}(t_old) + Delta_n(t), what is Delta_n?
    Delta_n = sum of all coefficients c_S where S involves at least one new tile.
    New tiles at n: (n, y) for y = 1, ..., n-2.
    """
    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    new_tile_indices = {i for i, (x, y) in enumerate(tiles) if x == n}

    # Separate old and new coefficients
    old_coeffs = {}
    new_coeffs = {}
    for S, c in coeffs.items():
        bits = {i for i in range(m) if S & (1 << i)}
        if bits & new_tile_indices:
            new_coeffs[S] = c
        else:
            old_coeffs[S] = c

    return old_coeffs, new_coeffs, tiles, m, coeffs

def test_multiplicative_structure():
    """
    Test if H_n = H_{n-1} * (1 + something involving new tiles).

    If H_n(t) = H_{n-1}(t_old) * F_n(t_new, t_old), then the coefficients
    of the new terms would be products of old coefficients and F coefficients.
    """
    print(f"{'='*70}")
    print(f"MULTIPLICATIVE STRUCTURE TEST")
    print(f"{'='*70}")

    for n in [5, 6, 7]:
        old_c, new_c, tiles, m, full_c = compute_layer_contribution(n)
        print(f"\n  n={n}:")
        print(f"    Old coefficients (from n-1): {len(old_c)}")
        print(f"    New coefficients (involve vertex n): {len(new_c)}")

        # Group new coefficients by how many new tiles they involve
        by_new_count = defaultdict(int)
        for S, c in new_c.items():
            new_tile_indices = {i for i, (x, y) in enumerate(tiles) if x == n}
            bits = {i for i in range(m) if S & (1 << i)}
            new_count = len(bits & new_tile_indices)
            by_new_count[new_count] += 1

        print(f"    By number of new tiles involved: {dict(sorted(by_new_count.items()))}")

        # If multiplicative: H_n(t) / H_{n-1}(t_old) = F_n(t)
        # At the special point t_old = 0 (all old tiles forward):
        #   H_n(0_old, t_new) = H_{n-1}(0) * F_n(t_new, 0)
        #   H_n(0_old, t_new) = 1 * F_n(t_new, 0) = F_n(t_new, 0)
        # So F at old=0 is just H with old tiles set to 0.

        # Compute H_n with all old tiles at 0, varying new tiles
        H, _, _ = compute_all_H(n)
        new_indices = sorted(i for i, (x, y) in enumerate(tiles) if x == n)
        num_new = len(new_indices)

        print(f"    New tiles: {[tiles[i] for i in new_indices]}")
        print(f"    F_n(t_new, t_old=0) = H_n(t_old=0, t_new):")

        for bits in range(1 << num_new):
            tiling = 0
            for k, idx in enumerate(new_indices):
                if bits & (1 << k):
                    tiling |= (1 << idx)
            h = H[tiling]
            new_labels = []
            for k in range(num_new):
                if bits & (1 << k):
                    new_labels.append(f"1")
                else:
                    new_labels.append(f"0")
            if num_new <= 5:
                print(f"      t_new=({','.join(new_labels)}): H = {h}")

        # Now test multiplicative hypothesis more directly
        # If H_n = H_{n-1} * F_n, then for each (t_old, t_new):
        #   H_n(t_old, t_new) = H_{n-1}(t_old) * F_n(t_new, t_old)
        # This means F_n = H_n / H_{n-1} at each point.
        # But H_{n-1}(t_old) can be zero... well, H is always odd, so H_{n-1} is always >= 1.
        # Actually H >= 1 always (Redei).

        # Compute H_{n-1}(t_old) for each old tiling
        tiles_prev = make_tiles(n-1)
        m_prev = len(tiles_prev)
        H_prev, _, _ = compute_all_H(n-1)

        # Map old tiles at n to tiles at n-1
        old_indices = sorted(i for i, (x, y) in enumerate(tiles) if x != n)
        old_to_prev = {}
        for j, idx in enumerate(old_indices):
            t = tiles[idx]
            for k, tp in enumerate(tiles_prev):
                if tp == t:
                    old_to_prev[j] = k
                    break

        # For each tiling at n, decompose into old part and new part
        # Check if H_n / H_{n-1}(old) depends only on new tiles
        print(f"\n    Testing multiplicative factorization:")
        if num_new <= 4 and m <= 15:
            # For each new tiling pattern, check if H_n/H_{n-1}(old) is constant
            ratios_by_new = defaultdict(set)
            for t in range(1 << m):
                # Extract old part
                old_bits = 0
                for j, idx in enumerate(old_indices):
                    if t & (1 << idx):
                        old_bits |= (1 << old_to_prev[j])

                # Extract new part
                new_bits = 0
                for k, idx in enumerate(new_indices):
                    if t & (1 << idx):
                        new_bits |= (1 << k)

                h_n = int(H[t])
                h_prev = int(H_prev[old_bits])

                ratio = Fraction(h_n, h_prev)
                ratios_by_new[new_bits].add(ratio)

            # If multiplicative, each new_bits pattern should give a single ratio
            is_multiplicative = all(len(v) == 1 for v in ratios_by_new.values())
            print(f"      Multiplicative? {is_multiplicative}")
            if not is_multiplicative:
                for nb, ratios in sorted(ratios_by_new.items()):
                    if len(ratios) > 1:
                        new_labels = ''.join(str((nb >> k) & 1) for k in range(num_new))
                        print(f"        new={new_labels}: {len(ratios)} distinct ratios "
                              f"(range [{min(ratios)}, {max(ratios)}])")
            else:
                print(f"      Layer factors:")
                for nb in range(1 << num_new):
                    ratio = list(ratios_by_new[nb])[0]
                    new_labels = ''.join(str((nb >> k) & 1) for k in range(num_new))
                    print(f"        F_n({new_labels}) = {ratio}")

def test_additive_decomposition():
    """
    Test if H_n = H_{n-1} + G_n where G_n is a "layer polynomial" that
    depends on new and old tiles but is additive (not multiplicative).
    """
    print(f"\n{'='*70}")
    print(f"ADDITIVE LAYER DECOMPOSITION")
    print(f"{'='*70}")

    for n in [5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)

        tiles_prev = make_tiles(n-1)
        m_prev = len(tiles_prev)
        H_prev, _, _ = compute_all_H(n-1)

        old_indices = sorted(i for i, (x, y) in enumerate(tiles) if x != n)
        new_indices = sorted(i for i, (x, y) in enumerate(tiles) if x == n)
        num_new = len(new_indices)

        old_to_prev = {}
        for j, idx in enumerate(old_indices):
            t = tiles[idx]
            for k, tp in enumerate(tiles_prev):
                if tp == t:
                    old_to_prev[j] = k
                    break

        # Compute G_n(t) = H_n(t) - H_{n-1}(t_old)
        G = np.zeros(1 << m, dtype=np.int64)
        for t in range(1 << m):
            old_bits = 0
            for j, idx in enumerate(old_indices):
                if t & (1 << idx):
                    old_bits |= (1 << old_to_prev[j])
            G[t] = H[t] - H_prev[old_bits]

        # Extract multilinear coefficients of G
        coeffs_G = mobius_inversion(G, m)

        print(f"\n  n={n}: G_n = H_n - H_{{n-1}} (layer increment)")
        print(f"    G range: [{G.min()}, {G.max()}]")
        print(f"    G(0) = {G[0]} (should be 0, since H_n(0) = 1 = H_{{n-1}}(0))")
        print(f"    Nonzero G coefficients: {len(coeffs_G)}")

        # Does G only involve new tiles? (It shouldn't — new tiles interact with old.)
        only_new = sum(1 for S in coeffs_G if all(i in set(new_indices) for i in range(m) if S & (1<<i)))
        mixed = sum(1 for S in coeffs_G if any(i in set(new_indices) for i in range(m) if S & (1<<i))
                    and not all(i in set(new_indices) for i in range(m) if S & (1<<i)))
        print(f"    Pure-new coefficients: {only_new}")
        print(f"    Mixed (new+old) coefficients: {mixed}")
        print(f"    Pure-old coefficients: {len(coeffs_G) - only_new - mixed}")

        # Analyze G's structure by degree
        by_deg = defaultdict(list)
        for S, c in coeffs_G.items():
            by_deg[bin(S).count('1')].append(c)

        for d in sorted(by_deg.keys()):
            vals = by_deg[d]
            print(f"    d={d}: {len(vals)} terms, sum={sum(vals)}, range=[{min(vals)},{max(vals)}]")

def test_determinant_representation():
    """
    Test if H(t) can be expressed as det(I + M(t)) for some matrix M.

    If H = det(I + M), then log(H) = tr(log(I+M)) = tr(M - M²/2 + M³/3 - ...)
    This would connect the multilinear polynomial to matrix structure.

    Key insight: H = I(Omega, 2) = sum_{I indep} 2^|I|.
    This is the independence polynomial at x=2.
    For line graphs: I(L(G), x) relates to matching polynomials.
    For general graphs: no simple determinant formula exists.

    BUT: what if we look at H(t) - 1 and check if it equals det(I + M(t)) - 1?
    """
    print(f"\n{'='*70}")
    print(f"DETERMINANT REPRESENTATION TEST")
    print(f"{'='*70}")

    for n in [4, 5]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        print(f"\n  n={n}: Trying H(t) = det(I_k + M(t)) for various k")

        # For k=1: det(1 + m(t)) = 1 + m(t). This gives linear polynomial only. Skip.

        # For k=2: det([[1+a, b],[c, 1+d]]) = (1+a)(1+d) - bc
        # = 1 + a + d + ad - bc
        # Has degree ≤ 2. Need to check if n=4 (max degree 2) fits.

        if n == 4:
            # n=4, m=3, max degree 2. Coefficients:
            # d=0: 1
            # d=1: 8, -4 → tiles (4,1)=4, (3,1)=2, (4,2)=2 but 8=4+2+2=8, hmm
            # Actually at n=4: 3 tiles, sums: 1, 8, -4
            # Can 1 + 8p - 4p² = det(I_2 + M(p)) for some 2×2 M(p)?
            # det = 1 + tr(M)·p + (tr(M)² - tr(M²))/2 · p² [for M = p·A]
            # Need tr(A)=8, (tr(A)²-tr(A²))/2 = -4
            # (64-tr(A²))/2 = -4, tr(A²)=72.
            # Eigenvalues a,b of A: a+b=8, a²+b²=72. (a+b)²=64, 2ab=-8, ab=-4.
            # a,b = roots of x²-8x-4=0 = (8±sqrt(80))/2 = 4±2√5.
            # So A has eigenvalues 4+2√5, 4-2√5. Real! det(I+A·p) = (1+(4+2√5)p)(1+(4-2√5)p)
            # = 1 + 8p + (16-20)p² = 1 + 8p - 4p². YES!

            print(f"    n=4: univariate H(p) = 1 + 8p - 4p²")
            print(f"    = det(I_2 + p·A) where A has eigenvalues 4+2sqrt(5), 4-2sqrt(5)")
            print(f"    = (1 + (4+2√5)p)(1 + (4-2√5)p)")
            print(f"    This works for the UNIVARIATE version H(p,...,p)!")

        if n == 5:
            # n=5: univariate H(p) = 1 + 22p + 0p² - 20p³ + 6p⁴
            # = det(I_4 + p·A) for some 4×4 A?
            # Coefficients of det(I + pA) = 1 + e₁(A)p + e₂(A)p² + e₃(A)p³ + e₄(A)p⁴
            # where e_k are elementary symmetric polynomials of eigenvalues.
            # e₁ = 22, e₂ = 0, e₃ = -20, e₄ = 6.
            # e₂=0 means the sum of products of pairs of eigenvalues = 0.
            # The characteristic polynomial of A: x⁴ - 22x³ + 0x² + 20x - 6 = 0.

            from numpy.polynomial import polynomial as P
            # Roots of x⁴ - 22x³ + 20x - 6
            roots = np.roots([1, -22, 0, 20, -6])
            print(f"\n    n=5: univariate H(p) = 1 + 22p + 0p² - 20p³ + 6p⁴")
            print(f"    If H(p) = det(I₄ + pA), eigenvalues of A: {roots}")
            print(f"    e₁={sum(roots):.4f} (expect 22)")
            print(f"    e₂={sum(roots[i]*roots[j] for i in range(4) for j in range(i+1,4)):.4f} (expect 0)")
            print(f"    e₃={sum(roots[i]*roots[j]*roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4)):.4f} (expect -20)")
            print(f"    e₄={np.prod(roots):.4f} (expect 6)")
            print(f"    All real? {all(abs(r.imag) < 1e-10 for r in roots)}")

        # n=6: H(p) = 1 + 52p + 92p² - 208p³ + 80p⁴ (degree 4)
        # det(I₄ + pA): e₁=52, e₂=92, e₃=-208, e₄=80
        # Char poly: x⁴ - 52x³ + 92x² + 208x - 80

        if n == 5:
            print(f"\n    n=6: H(p) = 1 + 52p + 92p² - 208p³ + 80p⁴")
            roots6 = np.roots([1, -52, 92, 208, -80])
            print(f"    Eigenvalues: {roots6}")
            print(f"    All real? {all(abs(r.imag) < 1e-10 for r in roots6)}")

            print(f"\n    n=7: H(p) = 1 + 114p + 616p² - 948p³ - 138p⁴ + 532p⁵ - 146p⁶")
            roots7 = np.roots([1, -114, 616, 948, -138, -532, -146])
            # Wait, sign convention: det(I + pA) = sum_k e_k(A) p^k
            # where e_k are elementary symmetric polynomials of eigenvalues.
            # For degree 6: need 6×6 matrix.
            # Char poly of A: prod(x - lambda_i) = x^6 - e1*x^5 + e2*x^4 - e3*x^3 + ...
            # = x^6 - 114x^5 + 616x^4 + 948x^3 - 138x^2 - 532x + 146
            # Wait no: e3 = -(-948) = 948? Let me be more careful.
            # det(I + pA) = sum_{k=0}^{6} e_k p^k where:
            # e_0 = 1
            # e_1 = 114 (trace)
            # e_2 = 616
            # e_3 = -948  (note: negative!)
            # e_4 = -138
            # e_5 = 532
            # e_6 = -146 (determinant)
            # Char poly: det(xI - A) = x^6 - e1*x^5 + e2*x^4 - e3*x^3 + e4*x^2 - e5*x + e6
            # = x^6 - 114x^5 + 616x^4 + 948x^3 - 138x^2 - 532x - 146
            roots7 = np.roots([1, -114, 616, 948, -138, -532, -146])
            print(f"    Eigenvalues: {roots7}")
            print(f"    All real? {all(abs(r.imag) < 1e-10 for r in roots7)}")

            # Check: are the eigenvalues nice?
            for r in sorted(roots7, key=lambda x: -x.real):
                if abs(r.imag) < 1e-10:
                    print(f"      lambda = {r.real:.6f}, log2 = {np.log2(abs(r.real)):.4f}")

def test_layer_factorization_multivariate():
    """
    Test if the multilinear H(t) factors as a product over layers.

    H_n(t) = prod_{v=3}^{n} L_v(t_{tiles of v}, t_{tiles below v})

    This is a STRONGER claim than multiplicative — each layer only depends
    on its own tiles and the tiles below.
    """
    print(f"\n{'='*70}")
    print(f"MULTIVARIATE LAYER FACTORIZATION")
    print(f"{'='*70}")

    n = 5
    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)

    # Layer 3: tiles involving vertex 3: (3,1) [index 2 at n=5]
    # Layer 4: tiles involving vertex 4 but not higher: (4,1), (4,2) [indices 1, 4]
    # Layer 5: tiles involving vertex 5: (5,1), (5,2), (5,3) [indices 0, 3, 5]

    # The "layer decomposition" would be:
    # H_5 = L_3(t_{(3,1)}) * L_4(t_{(4,1)}, t_{(4,2)}, t_{(3,1)}) * L_5(...)
    # L_3 is the contribution from just having vertex 3 (with base path 5→4→3→2→1)

    # Actually, the recursive structure says:
    # H_3 = polynomial in t_{(3,1)}
    # H_4 = H_3(old) + new terms involving (4,1), (4,2)
    # H_5 = H_4(old) + new terms involving (5,1), (5,2), (5,3)

    # Can we factor H_5 = H_3 * F_4 * F_5?
    # This requires H_n / H_{n-1} to be independent of how we evaluate old tiles
    # We already tested this above (it's NOT multiplicative)

    # Instead, test a WEAKER claim: does H factor as H_3 * (something)?
    H3, _, _ = compute_all_H(3)
    # H_3 has 1 tile: (3,1), m3=1. H_3(0) = 1, H_3(1) = 3.
    # So H_3(t) = 1 + 2*t_{(3,1)}

    # Check: does H_5 = (1 + 2*t_{(3,1)}) * Q(t)?
    # At t_{(3,1)} = 0: H_5 = Q(t). At t_{(3,1)} = 1: H_5 = 3*Q(t).
    # So Q(t) = H_5(t | t_{(3,1)}=0), and 3*Q(t) should equal H_5(t | t_{(3,1)}=1).

    idx_31 = next(i for i, t in enumerate(tiles) if t == (3,1))

    print(f"\n  n=5: Testing H_5 = (1 + 2*t_{{(3,1)}}) * Q(t_{{'other'}})")
    matches = 0
    total = 0
    for t in range(1 << m):
        if t & (1 << idx_31):
            continue  # Only check with t_{(3,1)} = 0

        h0 = int(H[t])  # H with (3,1)=0
        h1 = int(H[t | (1 << idx_31)])  # H with (3,1)=1

        total += 1
        if h1 == 3 * h0:
            matches += 1
        elif total <= 5:
            print(f"    t={t:06b}: H(0)={h0}, H(1)={h1}, ratio={Fraction(h1,h0)}")

    print(f"  H_5(t|(3,1)=1) = 3 * H_5(t|(3,1)=0)? {matches}/{total} match")

    # More generally: for the skip-2 tile (k+1, k-1), which is the "minimal" tile,
    # does flipping it multiply H by a constant?
    # Skip-2 tiles: (3,1), (4,2), (5,3)
    for skip2_tile in [(3,1), (4,2), (5,3)]:
        idx = next(i for i, t in enumerate(tiles) if t == skip2_tile)
        ratios = set()
        for t in range(1 << m):
            if t & (1 << idx):
                continue
            h0 = int(H[t])
            h1 = int(H[t | (1 << idx)])
            ratios.add(Fraction(h1, h0))

        if len(ratios) == 1:
            print(f"  Tile {skip2_tile}: H(1)/H(0) = {list(ratios)[0]} (constant!)")
        else:
            print(f"  Tile {skip2_tile}: H(1)/H(0) takes {len(ratios)} values: "
                  f"{sorted(ratios)[:5]}...")

if __name__ == "__main__":
    test_multiplicative_structure()
    test_additive_decomposition()
    test_determinant_representation()
    test_layer_factorization_multivariate()
