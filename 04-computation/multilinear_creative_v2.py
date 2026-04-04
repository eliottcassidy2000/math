"""
multilinear_creative_v2.py — opus-2026-04-04-S3

Creative exploration round 2:
1. Factorization attempts for H(t)
2. Recursive n -> n+1 structure
3. Same-end formula proof: c_{ij} = -2^{max(1,|s1-s2|-1)} for same-end
4. Negative eigenspace geometry — connection to staircase
5. Polynomial as a function on the lattice path cube
6. Multilinear Hessian at special points (transitive, regular, Paley)
7. Boolean derivative / influence of each tile
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def v2(n):
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0: n //= 2; v += 1
    return v

def tile_relation(t1, t2):
    x1, y1 = t1; x2, y2 = t2
    if x1 == x2 or y1 == y2: return "same-end"
    if x1 == y2 or y1 == x2: return "cross-end"
    return "disjoint"

# =================================================================
# INVESTIGATION 1: Same-end formula proof
# =================================================================
def investigate_same_end_formula(n):
    """Test c_{ij} = -2^{max(1, |s_i - s_j| - 1)} for same-end pairs."""
    print(f"\n{'='*60}")
    print(f"SAME-END FORMULA TEST n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    matches = 0
    tests = 0
    for S, c in coeffs.items():
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        t1, t2 = tiles[bits[0]], tiles[bits[1]]
        if tile_relation(t1, t2) != "same-end": continue

        s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
        expected = -(2 ** max(1, abs(s1-s2)-1))
        tests += 1
        if c == expected:
            matches += 1
        else:
            print(f"  MISMATCH: {t1}-{t2} skips=({s1},{s2}): c={c}, expected={expected}")

    print(f"  {matches}/{tests} match formula c = -2^max(1,|s1-s2|-1)")

# =================================================================
# INVESTIGATION 2: Disjoint pair formula search
# =================================================================
def investigate_disjoint_formula(n):
    """Search for formula for disjoint pairs."""
    print(f"\n{'='*60}")
    print(f"DISJOINT PAIR FORMULA SEARCH n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    for S, c in sorted(coeffs.items(), key=lambda x: abs(x[1])):
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        t1, t2 = tiles[bits[0]], tiles[bits[1]]
        if tile_relation(t1, t2) != "disjoint": continue

        s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
        x1, y1 = t1; x2, y2 = t2

        # Staircase coordinates: row r = x-y-1, col c = y-1
        r1, c1 = x1-y1-1, y1-1
        r2, c2 = x2-y2-1, y2-1

        # Overlap = min(x1,x2) - max(y1,y2) (how many internal vertices shared)
        overlap = min(x1,x2) - max(y1,y2)

        # Manhattan distance in staircase
        manhattan = abs(r1-r2) + abs(c1-c2)

        # Various potential formulas
        formula_A = 2**(s1+s2-2)  # product of linear coefficients / 2
        formula_B = 2**(s1-1) * 2**(s2-1)  # same as A
        formula_C = 2**(min(s1,s2)-1) * 2**(max(s1,s2)-1)  # same as A

        ratio = Fraction(c, formula_A) if formula_A != 0 else None

        print(f"  {t1}-{t2} s=({s1},{s2}) ov={overlap} man={manhattan}: "
              f"c={c}, 2^(s1+s2-2)={formula_A}, ratio={ratio}")

# =================================================================
# INVESTIGATION 3: Recursive structure n -> n+1
# =================================================================
def investigate_recursive(n_small, n_big):
    """How do multilinear coefficients at n relate to n+1?"""
    print(f"\n{'='*60}")
    print(f"RECURSIVE STRUCTURE n={n_small} -> n={n_big}")
    print(f"{'='*60}")

    tiles_s = make_tiles(n_small)
    tiles_b = make_tiles(n_big)
    m_s = len(tiles_s)
    m_b = len(tiles_b)

    H_s, _, _ = compute_all_H(n_small)
    H_b, _, _ = compute_all_H(n_big)
    coeffs_s = mobius_inversion(H_s, m_s)
    coeffs_b = mobius_inversion(H_b, m_b)

    # Tiles at n+1 = tiles at n UNION new tiles in the bottom row
    # New tiles at n+1: those involving vertex n (0-indexed: n_big-1)
    # In 1-indexed: tile (x,y) with x = n_big or y = n_big... wait.
    # Actually: new tiles are those with x = n_big (the new top vertex in 1-indexed).
    # Wait, base path is n -> n-1 -> ... -> 1.
    # At n_big, new vertex is n_big (top of base path).
    # New tiles: (n_big, y) for y = 1, ..., n_big-2. That's n_big-2 = n_small-1 new tiles.

    new_tiles = [(n_big, y) for y in range(1, n_big-1)]
    old_tiles = [t for t in tiles_b if t not in new_tiles]

    print(f"  Tiles at n={n_small}: {m_s}")
    print(f"  Tiles at n={n_big}: {m_b}")
    print(f"  New tiles at n={n_big}: {new_tiles}")
    print(f"  Old tiles (should match n={n_small}): {old_tiles}")

    # Do old tiles at n+1 match tiles at n?
    # The tile (x,y) at n_big corresponds to (x,y) at n_small... but only if x <= n_small.
    # Actually tiles at n_small have x <= n_small, y >= 1, x-y >= 2.
    # Old tiles at n_big: tiles NOT involving the new vertex n_big in upper position.
    # But there's also a new base-path arc n_big -> n_big-1. This doesn't affect tiles.

    # Map old tiles at n_big to tiles at n_small
    tile_map = {}  # index in big -> index in small
    for j, t in enumerate(tiles_b):
        if t in new_tiles:
            continue
        for i, t_s in enumerate(tiles_s):
            if t_s == t:
                tile_map[j] = i
                break

    print(f"  Tile mapping (big->small): {len(tile_map)} mapped")

    # Check: are the OLD linear coefficients at n_big the same as at n_small?
    print(f"\n  Linear coefficient comparison:")
    for j in range(m_b):
        t = tiles_b[j]
        c_big = coeffs_b.get(1 << j, 0)
        if j in tile_map:
            i = tile_map[j]
            c_small = coeffs_s.get(1 << i, 0)
            match = "=" if c_big == c_small else f"!= (small={c_small})"
            print(f"    {t}: big={c_big}, small {match}")
        else:
            print(f"    {t}: NEW, c={c_big}")

    # Check: are OLD quadratic coefficients preserved?
    print(f"\n  Quadratic coefficient comparison (old-old pairs):")
    preserved = 0
    changed = 0
    for j1 in range(m_b):
        for j2 in range(j1+1, m_b):
            if j1 not in tile_map or j2 not in tile_map:
                continue
            S_big = (1 << j1) | (1 << j2)
            c_big = coeffs_b.get(S_big, 0)

            i1, i2 = tile_map[j1], tile_map[j2]
            S_small = (1 << i1) | (1 << i2)
            c_small = coeffs_s.get(S_small, 0)

            if c_big == c_small:
                preserved += 1
            else:
                changed += 1
                if changed <= 10:
                    print(f"    {tiles_b[j1]}-{tiles_b[j2]}: big={c_big}, small={c_small}, diff={c_big-c_small}")

    print(f"  Preserved: {preserved}, Changed: {changed}")

# =================================================================
# INVESTIGATION 4: Negative eigenspace and staircase geometry
# =================================================================
def investigate_eigenspace_geometry(n):
    """Study the negative eigenspace of Q and its relation to the staircase."""
    print(f"\n{'='*60}")
    print(f"EIGENSPACE GEOMETRY n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        Q[bits[0]][bits[1]] = c
        Q[bits[1]][bits[0]] = c

    eigvals, eigvecs = np.linalg.eigh(Q)

    # Staircase coordinates for each tile
    staircase = []
    for t in tiles:
        x, y = t
        r = x - y - 1  # row (0-indexed from bottom)
        c = y - 1       # column (0-indexed from left)
        staircase.append((r, c))

    print(f"  Staircase coords: {list(zip(tiles, staircase))}")

    # Negative eigenvectors: decompose in staircase coordinates
    print(f"\n  NEGATIVE EIGENVECTORS in staircase coords:")
    for idx in range(m):
        e = eigvals[idx]
        if e > -1e-10:
            continue
        v = eigvecs[:, idx]

        print(f"\n  lambda = {e:.4f}:")
        # Print as staircase grid
        max_r = max(r for r, c in staircase)
        max_c = max(c for r, c in staircase)
        grid = np.zeros((max_r+1, max_c+1))
        for i, (r, c) in enumerate(staircase):
            grid[r][c] = v[i]

        for r in range(max_r, -1, -1):
            row_str = f"    r={r}: "
            for c in range(max_c+1):
                if r + c < n - 2:  # valid staircase position
                    row_str += f"{grid[r][c]:7.3f}"
                else:
                    row_str += "       "
            print(row_str)

    # Null vector (if exists)
    zero_idx = [i for i in range(m) if abs(eigvals[i]) < 1e-10]
    if zero_idx:
        print(f"\n  NULL VECTOR in staircase coords:")
        for idx in zero_idx:
            v = eigvecs[:, idx]
            max_r = max(r for r, c in staircase)
            max_c = max(c for r, c in staircase)
            grid = np.zeros((max_r+1, max_c+1))
            for i, (r, c) in enumerate(staircase):
                grid[r][c] = v[i]
            for r in range(max_r, -1, -1):
                row_str = f"    r={r}: "
                for c in range(max_c+1):
                    if r + c < n - 2:
                        row_str += f"{grid[r][c]:7.3f}"
                    else:
                        row_str += "       "
                print(row_str)

# =================================================================
# INVESTIGATION 5: Boolean influence of each tile
# =================================================================
def investigate_influence(n):
    """Compute Boolean influence of each tile on H(t)."""
    print(f"\n{'='*60}")
    print(f"BOOLEAN INFLUENCE ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)

    # Influence of tile k = Pr_t[H(t) != H(t XOR e_k)] where t is uniform on {0,1}^m
    # But H takes integer values, so always != unless flipping k doesn't change H.
    # Better: compute E[|H(t) - H(t XOR e_k)|] and E[(H(t) - H(t XOR e_k))^2]
    # These measure the "sensitivity" of H to each tile.

    N = 1 << m
    for k in range(m):
        diff_sum = 0
        diff_sq_sum = 0
        unchanged_count = 0
        for t in range(N):
            t_flipped = t ^ (1 << k)
            diff = int(H[t_flipped]) - int(H[t])
            diff_sum += diff
            diff_sq_sum += diff * diff
            if diff == 0:
                unchanged_count += 1

        mean_diff = diff_sum / N
        mean_sq = diff_sq_sum / N
        variance = mean_sq - mean_diff**2
        unchanged_frac = unchanged_count / N

        tile = tiles[k]
        skip = tile[0] - tile[1]
        linear_coeff = 2**(skip-1)

        print(f"  tile {tile} skip={skip}: "
              f"E[dH]={mean_diff:.3f} (=c_k*E[deriv]?), "
              f"Var(dH)={variance:.3f}, "
              f"E[dH^2]={mean_sq:.3f}, "
              f"unchanged={unchanged_frac:.4f}")

# =================================================================
# INVESTIGATION 6: Factorization of H(t) - 1
# =================================================================
def investigate_factorization(n):
    """Test if H(t) - 1 or H(t) has useful factorizations."""
    print(f"\n{'='*60}")
    print(f"FACTORIZATION INVESTIGATION n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # H(t) - 1: remove constant term
    coeffs_centered = {S: c for S, c in coeffs.items() if S != 0}

    # Check: can H(t) - 1 be written as a product of two multilinear polynomials?
    # If H - 1 = P * Q with deg(P) >= 1 and deg(Q) >= 1, then
    # H(t) - 1 = 0 when P(t) = 0 or Q(t) = 0.
    # H(t) = 1 iff tournament is transitive. How many tilings give H=1?
    transitive_count = sum(1 for h in H if h == 1)
    print(f"  H(t) = 1 (transitive) for {transitive_count}/{1<<m} tilings")

    # GCD of all coefficients (excluding constant)
    from math import gcd
    vals = [abs(c) for c in coeffs_centered.values() if c != 0]
    if vals:
        g = vals[0]
        for v in vals[1:]:
            g = gcd(g, v)
        print(f"  GCD of all nonzero non-constant coefficients: {g}")
        # If g > 1, then H(t) - 1 = g * P(t) for some integer polynomial P
        if g > 1:
            print(f"  H(t) = 1 + {g} * P(t) where P(t) is a multilinear polynomial")

    # H(t) - 1 is always even (since H is odd, H-1 is even)
    # Dividing by 2: (H(t)-1)/2 is a multilinear polynomial
    half_coeffs = {S: c//2 for S, c in coeffs_centered.items() if c != 0}
    print(f"  (H(t)-1)/2 coefficients: {len(half_coeffs)} nonzero")
    half_vals = [c for c in half_coeffs.values()]
    if half_vals:
        half_g = half_vals[0]
        for v in half_vals[1:]:
            half_g = gcd(abs(half_g), abs(v))
        print(f"  GCD of (H-1)/2 coefficients: {half_g}")

    # Check if H(t) - 1 is divisible by any single-variable factor (1 + a*t_k)
    # H(t) - 1 = (1 + a*t_k) * Q(t) would mean H(t_k=0) - 1 | H(t_k=1) - 1 always
    print(f"\n  Linear factor check:")
    for k in range(m):
        # H(t) restricted to t_k=0: sum of coefficients for S not containing k
        h_at_0 = sum(c for S, c in coeffs.items() if not (S & (1 << k)))
        h_at_1 = sum(c for S, c in coeffs.items())  # This isn't right...

        # Actually, H at t_k=0 means the polynomial with t_k set to 0
        # = sum_{S: k not in S} c_S * prod_{i in S} t_i
        # This is still a polynomial in other variables!
        # To check if (1-t_k) | (H-1), we need H(t_k=1) = 1 for all other t.
        # That means flipping tile k to backward always gives H=1. Obviously not true.

        # Instead, check a different kind of factor.
        # Can we write H(t) - 1 = t_k * P(t) + (1-t_k) * Q(t)?
        # That's trivially true (P = terms with k, Q = terms without k).
        # More interesting: is (H(t)-1) / sum_{S ni k} c_S * prod t_i a polynomial?

        # Just check: does sum of coefficients involving tile k equal 0?
        # That would mean the net effect of tile k averaged over all other tiles is 0.
        s_k = sum(c for S, c in coeffs_centered.items() if S & (1 << k))
        if k < 6:
            print(f"    tile {tiles[k]}: sum of coefficients involving k = {s_k}")

    # Check H(t) mod 4
    print(f"\n  H(t) mod 4 distribution:")
    mod4_dist = defaultdict(int)
    for h in H:
        mod4_dist[int(h) % 4] += 1
    print(f"    {dict(sorted(mod4_dist.items()))}")

    # Check if there's a SUBSET of tiles that factors
    # Meaning: H(t) = P(t_{S1}) * Q(t_{S2}) for some partition S1 U S2 of tiles
    # This would mean the tile flips in S1 are independent of those in S2
    # Test: does H(t) = H(t_row1) * H(t_row2) for rows of the staircase?
    print(f"\n  Row independence check:")
    # Group tiles by row in staircase (row = y-1)
    rows = defaultdict(list)
    for i, (x, y) in enumerate(tiles):
        rows[y-1].append(i)

    # For each pair of rows, check if cross-row coefficients are all zero
    row_keys = sorted(rows.keys())
    for r1_idx, r1 in enumerate(row_keys):
        for r2 in row_keys[r1_idx+1:]:
            cross_coeffs = 0
            for i1 in rows[r1]:
                for i2 in rows[r2]:
                    S = (1 << i1) | (1 << i2)
                    if S in coeffs and coeffs[S] != 0:
                        cross_coeffs += 1
            if cross_coeffs == 0:
                print(f"    Rows {r1},{r2}: INDEPENDENT (no cross-row quadratic coefficients)")
            else:
                print(f"    Rows {r1},{r2}: {cross_coeffs} cross-row interactions")

# =================================================================
# INVESTIGATION 7: H evaluated at special points
# =================================================================
def investigate_special_points(n):
    """Evaluate H at special points of the hypercube and its interior."""
    print(f"\n{'='*60}")
    print(f"SPECIAL POINT EVALUATIONS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # H at (1/2, 1/2, ..., 1/2) = average over all tilings
    avg = sum(int(H[t]) for t in range(1 << m)) / (1 << m)
    print(f"  H(1/2,...,1/2) = average = {avg}")
    print(f"  = {Fraction(sum(int(H[t]) for t in range(1 << m)), 1 << m)}")

    # H at (p, p, ..., p) for various p
    print(f"\n  H(p, ..., p) as function of p:")
    for p_num in range(0, 11):
        p = Fraction(p_num, 10)
        val = Fraction(0)
        for S, c in coeffs.items():
            deg = bin(S).count('1')
            val += c * p**deg
        print(f"    p={float(p):.1f}: H = {float(val):.4f}")

    # The function H(p,...,p) is a UNIVARIATE polynomial in p!
    # Its coefficients are the sums of multilinear coefficients by degree.
    print(f"\n  Univariate H(p) = sum_d (sum_{{|S|=d}} c_S) * p^d:")
    by_deg = defaultdict(int)
    for S, c in coeffs.items():
        by_deg[bin(S).count('1')] += c

    max_d = max(by_deg.keys())
    poly_str = " + ".join(f"{by_deg[d]}*p^{d}" for d in range(max_d+1) if by_deg[d] != 0)
    print(f"    H(p) = {poly_str}")

    # Roots of H(p) - 1 = 0 (beyond p=0)
    # Since H(0) = 1, p=0 is a root of H(p)-1.
    poly_coeffs = [by_deg.get(d, 0) for d in range(max_d+1)]
    poly_coeffs[0] -= 1  # H(p) - 1
    # Factor out p: (H(p)-1)/p = sum_d c_d * p^(d-1) for d >= 1
    reduced = poly_coeffs[1:]
    roots = np.roots(reduced[::-1])  # numpy wants highest degree first
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-10]
    print(f"    Roots of (H(p)-1)/p = 0: {sorted(real_roots)}")

    # The critical points of H(p,...,p)
    # dH/dp = sum_d d * (sum_{{|S|=d}} c_S) * p^(d-1)
    deriv_coeffs = [d * by_deg.get(d, 0) for d in range(1, max_d+1)]
    if deriv_coeffs:
        crit_roots = np.roots(deriv_coeffs[::-1])
        real_crits = sorted([r.real for r in crit_roots if abs(r.imag) < 1e-10])
        print(f"    Critical points of H(p): {real_crits}")

    # H at the "all-equal" point t_k = 1 for long-range, t_k = 0 for short-range
    # (threshold at various skip values)
    print(f"\n  H with skip-threshold tilings:")
    for threshold in range(2, n):
        tiling = 0
        for i, (x, y) in enumerate(tiles):
            if x - y >= threshold:
                tiling |= (1 << i)
        print(f"    skip >= {threshold}: H = {H[tiling]}")

# =================================================================
# INVESTIGATION 8: Coefficient sign pattern at max degree
# =================================================================
def investigate_max_degree_signs(n):
    """Analyze the sign pattern of max-degree coefficients."""
    print(f"\n{'='*60}")
    print(f"MAX-DEGREE SIGN PATTERN n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    max_deg = max(bin(S).count('1') for S in coeffs.keys())
    max_coeffs = {S: c for S, c in coeffs.items() if bin(S).count('1') == max_deg}

    print(f"  Max degree: {max_deg}")
    print(f"  Number of max-degree coefficients: {len(max_coeffs)}")

    pos_count = sum(1 for c in max_coeffs.values() if c > 0)
    neg_count = sum(1 for c in max_coeffs.values() if c < 0)
    print(f"  Positive: {pos_count}, Negative: {neg_count}")

    # For each max-degree coefficient, which tiles are involved?
    # Look at the staircase structure
    for S, c in sorted(max_coeffs.items(), key=lambda x: x[0]):
        tile_set = [tiles[i] for i in range(m) if S & (1 << i)]
        # Staircase coordinates
        stair = [(x-y-1, y-1) for x, y in tile_set]
        # Which rows and columns are covered?
        rows = sorted(set(r for r, _ in stair))
        cols = sorted(set(c for _, c in stair))
        sign = "+" if c > 0 else "-"
        if len(max_coeffs) <= 50:
            print(f"  {sign} tiles={tile_set} stair={stair} rows={rows} cols={cols}")

    # Check: is the sign determined by the parity of some function of the tile set?
    print(f"\n  SIGN PARITY ANALYSIS:")
    # Hypothesis: sign = (-1)^{sum of row indices} or similar
    for func_name, func in [
        ("sum(rows)", lambda S: sum(tiles[i][0]-tiles[i][1]-1 for i in range(m) if S & (1<<i))),
        ("sum(cols)", lambda S: sum(tiles[i][1]-1 for i in range(m) if S & (1<<i))),
        ("sum(rows+cols)", lambda S: sum(tiles[i][0]-2 for i in range(m) if S & (1<<i))),
        ("num_tiles_with_even_skip", lambda S: sum(1 for i in range(m) if S & (1<<i) and (tiles[i][0]-tiles[i][1]) % 2 == 0)),
    ]:
        correct = 0
        total = len(max_coeffs)
        for S, c in max_coeffs.items():
            val = func(S)
            predicted_sign = (-1)**val
            actual_sign = 1 if c > 0 else -1
            if predicted_sign == actual_sign:
                correct += 1
        print(f"    sign = (-1)^{func_name}: {correct}/{total} correct")

# =================================================================
# MAIN
# =================================================================
print("MULTILINEAR CREATIVE EXPLORATION v2")
print("=" * 70)

for n in [5, 6, 7]:
    investigate_same_end_formula(n)
    investigate_disjoint_formula(n)

for n in [5, 6]:
    investigate_recursive(n, n+1)

for n in [5, 6, 7]:
    investigate_eigenspace_geometry(n)

for n in [5, 6]:
    investigate_influence(n)

for n in [5, 6, 7]:
    investigate_factorization(n)

for n in [5, 6, 7]:
    investigate_special_points(n)

for n in [5, 6, 7]:
    investigate_max_degree_signs(n)
