"""
max_degree_sign_and_betti.py — opus-2026-04-04-S3

1. Deep analysis of the sign pattern at max degree
2. Connection between multilinear polynomial and path homology Betti numbers
3. The H(t) landscape: local maxima, saddle points, gradient flow
4. Generating function approach: H_n(p) as function of n
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles, tiling_to_tournament, count_hamiltonian_paths

# =================================================================
# PART 1: Max-degree sign pattern — what determines +1 vs -1 in P?
# =================================================================
def max_degree_sign_deep(n):
    """
    At odd n, P(t) = (H(t)-1)/2 has ±1 at max degree d = 2*floor((n-1)/2).
    323 such terms at n=7: 125 positive, 198 negative.

    What determines the sign? Test many hypotheses.
    """
    print(f"\n{'='*60}")
    print(f"MAX-DEGREE SIGN ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    max_deg = max(bin(S).count('1') for S in coeffs.keys())
    max_coeffs = [(S, c//2) for S, c in coeffs.items() if bin(S).count('1') == max_deg]

    print(f"  Max degree: {max_deg}, {len(max_coeffs)} terms")
    print(f"  P values: {sorted(set(c for _,c in max_coeffs))}")

    pos = [(S,c) for S,c in max_coeffs if c > 0]
    neg = [(S,c) for S,c in max_coeffs if c < 0]
    print(f"  Positive: {len(pos)}, Negative: {len(neg)}")

    # For each max-degree subset S, compute various features
    features = []
    for S, c in max_coeffs:
        bits = [i for i in range(m) if S & (1 << i)]
        tile_set = [tiles[i] for i in bits]

        # Staircase coordinates: r = x-y-1 (row), col = y-1
        stair = [(t[0]-t[1]-1, t[1]-1) for t in tile_set]
        rows = [r for r,c in stair]
        cols = [c for r,c in stair]

        # Feature 1: sum of rows
        sum_r = sum(rows)
        # Feature 2: sum of columns
        sum_c = sum(cols)
        # Feature 3: sum of skip values
        sum_skip = sum(t[0]-t[1] for t in tile_set)
        # Feature 4: number of tiles with odd skip
        odd_skip = sum(1 for t in tile_set if (t[0]-t[1]) % 2 == 1)
        # Feature 5: number of tiles in left column (y=1)
        left_col = sum(1 for t in tile_set if t[1] == 1)
        # Feature 6: number of tiles in top row (x=n)
        top_row = sum(1 for t in tile_set if t[0] == n)
        # Feature 7: row multiset
        row_multiset = tuple(sorted(rows))
        # Feature 8: col multiset
        col_multiset = tuple(sorted(cols))
        # Feature 9: parity of product of skips
        prod_skip_parity = sum(t[0]-t[1] for t in tile_set) % 2
        # Feature 10: number of tiles below anti-diagonal
        below_diag = sum(1 for r,c in stair if r < c)
        above_diag = sum(1 for r,c in stair if r > c)
        on_diag = sum(1 for r,c in stair if r == c)

        features.append({
            'S': S, 'c': c, 'tiles': tile_set, 'stair': stair,
            'sum_r': sum_r, 'sum_c': sum_c, 'sum_skip': sum_skip,
            'odd_skip': odd_skip, 'left_col': left_col, 'top_row': top_row,
            'row_ms': row_multiset, 'col_ms': col_multiset,
            'below_diag': below_diag, 'above_diag': above_diag,
            'on_diag': on_diag
        })

    # Test each hypothesis
    for name, key in [
        ("sum_r parity", lambda f: f['sum_r'] % 2),
        ("sum_c parity", lambda f: f['sum_c'] % 2),
        ("sum_skip parity", lambda f: f['sum_skip'] % 2),
        ("odd_skip parity", lambda f: f['odd_skip'] % 2),
        ("left_col parity", lambda f: f['left_col'] % 2),
        ("top_row parity", lambda f: f['top_row'] % 2),
        ("(sum_r+sum_c) parity", lambda f: (f['sum_r']+f['sum_c']) % 2),
        ("(sum_r-sum_c) parity", lambda f: (f['sum_r']-f['sum_c']) % 2),
        ("below_diag > above_diag", lambda f: 1 if f['below_diag'] > f['above_diag'] else 0),
        ("(above - below) parity", lambda f: (f['above_diag'] - f['below_diag']) % 2),
        ("sum_r mod 3", lambda f: f['sum_r'] % 3),
        ("sum(r*c) parity", lambda f: sum(r*c for r,c in f['stair']) % 2),
        ("sum(r+c) parity", lambda f: sum(r+c for r,c in f['stair']) % 2),
        ("num tiles with r > c", lambda f: sum(1 for r,c in f['stair'] if r > c)),
    ]:
        correct = 0
        for f in features:
            predicted = (-1)**key(f)
            actual = 1 if f['c'] > 0 else -1
            if predicted == actual:
                correct += 1
        rate = correct / len(features)
        if rate > 0.8 or rate < 0.2:
            print(f"  *** {name}: {correct}/{len(features)} = {rate:.3f}")
        elif rate > 0.6 or rate < 0.4:
            print(f"    {name}: {correct}/{len(features)} = {rate:.3f}")

    # Try a different approach: row/column COVERAGE
    # At max degree, which rows and columns are covered?
    print(f"\n  Row/column coverage analysis:")
    row_coverage = defaultdict(int)
    col_coverage = defaultdict(int)
    for f in features:
        covered_rows = frozenset(r for r,c in f['stair'])
        covered_cols = frozenset(c for r,c in f['stair'])
        row_coverage[(covered_rows, f['c'] > 0)] += 1
        col_coverage[(covered_cols, f['c'] > 0)] += 1

    # Check: does the sign depend on which rows/columns are covered?
    # Group by row set
    by_rows = defaultdict(lambda: [0, 0])
    for f in features:
        key = frozenset(r for r,c in f['stair'])
        if f['c'] > 0:
            by_rows[key][0] += 1
        else:
            by_rows[key][1] += 1

    mixed_row_sets = sum(1 for k,v in by_rows.items() if v[0] > 0 and v[1] > 0)
    pure_row_sets = sum(1 for k,v in by_rows.items() if v[0] == 0 or v[1] == 0)
    print(f"  Row-set sign determination: {pure_row_sets} pure, {mixed_row_sets} mixed")

    by_cols = defaultdict(lambda: [0, 0])
    for f in features:
        key = frozenset(c for r,c in f['stair'])
        if f['c'] > 0:
            by_cols[key][0] += 1
        else:
            by_cols[key][1] += 1

    mixed_col_sets = sum(1 for k,v in by_cols.items() if v[0] > 0 and v[1] > 0)
    pure_col_sets = sum(1 for k,v in by_cols.items() if v[0] == 0 or v[1] == 0)
    print(f"  Col-set sign determination: {pure_col_sets} pure, {mixed_col_sets} mixed")

# =================================================================
# PART 2: Connection to Betti numbers
# =================================================================
def multilinear_vs_betti(n):
    """
    For each tiling, compute H(t) and related quantities.
    Check if the multilinear coefficients encode Betti numbers in any way.
    """
    print(f"\n{'='*60}")
    print(f"MULTILINEAR vs BETTI NUMBERS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # H mod 4 encodes beta_1 structure (from OCF: H mod 4 = 1 + 2*alpha_1 mod 4)
    # Can we extract beta_1 or beta_3 from the multilinear structure?

    # Recall: H(T) is odd always (Redei). H mod 4 ∈ {1, 3}.
    # H ≡ 1 (mod 4) iff alpha_1(T) is even (= even number of 3-cycles through any vertex)
    # The multilinear polynomial gives H(t) = 1 + 2P(t) where P is integer-valued.
    # H mod 4 = (1 + 2P) mod 4 = {1 if P even, 3 if P odd}

    # So P(t) mod 2 determines H(t) mod 4.
    # P(t) = sum_{S≠∅} (c_S/2) * prod t_i
    # P(t) mod 2 is a polynomial over GF(2)!

    # Compute P(t) mod 2 as a GF(2) multilinear polynomial
    P_mod2 = {}
    for S, c in coeffs.items():
        if S == 0: continue
        p_coeff = (c // 2) % 2
        if p_coeff != 0:
            P_mod2[S] = 1

    print(f"  P(t) mod 2: {len(P_mod2)} nonzero terms out of {2**m - 1}")

    # P mod 2 by degree
    by_deg = defaultdict(int)
    for S in P_mod2:
        by_deg[bin(S).count('1')] += 1
    print(f"  By degree: {dict(sorted(by_deg.items()))}")

    # This GF(2) polynomial determines which tilings give H ≡ 1 vs H ≡ 3 mod 4
    h_mod4 = defaultdict(int)
    for t in range(1 << m):
        h_mod4[int(H[t]) % 4] += 1
    print(f"  H mod 4 distribution: {dict(sorted(h_mod4.items()))}")

    # How many tilings have P even (H ≡ 1 mod 4)?
    p_even = h_mod4.get(1, 0)
    p_odd = h_mod4.get(3, 0)
    print(f"  P even (H≡1): {p_even}, P odd (H≡3): {p_odd}")

    # Check: is P mod 2 related to beta_1?
    # beta_1 = 0 or 1. beta_1 = 1 iff the tournament has a certain path homology structure.
    # Can we determine beta_1 from H mod 4? Not directly, but the correlation might be strong.

# =================================================================
# PART 3: The n-generating function
# =================================================================
def n_generating_function():
    """
    Study H_n(p) = sum_d S_d(n) p^d as a function of BOTH n and p.

    Known:
    - S_0(n) = 1 always
    - S_1(n) = 2^n - 2n (proved)
    - S_2(n) = ?

    Can we find formulas for S_d(n)?
    """
    print(f"\n{'='*60}")
    print(f"n-GENERATING FUNCTION FOR DEGREE SUMS")
    print(f"{'='*60}")

    degree_data = {}  # n -> list of S_d

    for n in [3, 4, 5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        if m == 0: continue
        H, _, _ = compute_all_H(n)
        coeffs = mobius_inversion(H, m)
        by_deg = defaultdict(int)
        for S, c in coeffs.items():
            by_deg[bin(S).count('1')] += c
        max_d = max(by_deg.keys())
        degree_data[n] = [by_deg.get(d, 0) for d in range(max_d + 1)]

    # Print table
    max_d = max(len(v) for v in degree_data.values())
    print(f"\n  {'n':>3s}", end="")
    for d in range(max_d):
        print(f"  {'S_'+str(d):>10s}", end="")
    print(f"  {'H(1)':>10s}")

    for n in sorted(degree_data.keys()):
        sums = degree_data[n]
        print(f"  {n:3d}", end="")
        for d in range(max_d):
            if d < len(sums):
                print(f"  {sums[d]:10d}", end="")
            else:
                print(f"  {'':>10s}", end="")
        print(f"  {sum(sums):10d}")

    # S_1(n) = 2^n - 2n. Check.
    print(f"\n  S_1(n) = 2^n - 2n:")
    for n in sorted(degree_data.keys()):
        actual = degree_data[n][1] if len(degree_data[n]) > 1 else 0
        formula = 2**n - 2*n
        print(f"    n={n}: formula={formula}, actual={actual}, match={formula==actual}")

    # S_2(n): values are 0, -4, 0, 92, 616
    # wait: n=3 has only degree 1, so S_2 = 0
    # n=4: S_2 = -4
    # n=5: S_2 = 0
    # n=6: S_2 = 92
    # n=7: S_2 = 616
    print(f"\n  S_2(n) values:")
    for n in sorted(degree_data.keys()):
        s2 = degree_data[n][2] if len(degree_data[n]) > 2 else 0
        print(f"    n={n}: S_2 = {s2}")

    # Try: S_2(n) as polynomial in n or 2^n
    s2_vals = [(n, degree_data[n][2] if len(degree_data[n]) > 2 else 0)
               for n in sorted(degree_data.keys())]

    # Check if S_2(n) = a*4^n + b*2^n + c*n^2 + d*n + e
    # We have 5 data points. Solve linear system.
    from numpy.linalg import lstsq
    ns = np.array([n for n, _ in s2_vals], dtype=float)
    s2s = np.array([s for _, s in s2_vals], dtype=float)

    # Try: S_2(n) = a*4^n + b*2^n + c*n^2 + d*n + e
    A = np.column_stack([4**ns, 2**ns, ns**2, ns, np.ones_like(ns)])
    x, res, _, _ = lstsq(A, s2s, rcond=None)
    print(f"\n  Fit S_2(n) = a*4^n + b*2^n + c*n^2 + d*n + e:")
    print(f"    a={x[0]:.6f}, b={x[1]:.6f}, c={x[2]:.6f}, d={x[3]:.6f}, e={x[4]:.6f}")
    for n, s2 in s2_vals:
        pred = x[0]*4**n + x[1]*2**n + x[2]*n**2 + x[3]*n + x[4]
        print(f"    n={n}: actual={s2}, predicted={pred:.1f}")

    # Try simpler: S_2(n) = a*4^n + b*2^n + c*n + d
    A2 = np.column_stack([4**ns, 2**ns, ns, np.ones_like(ns)])
    x2, _, _, _ = lstsq(A2, s2s, rcond=None)
    print(f"\n  Fit S_2(n) = a*4^n + b*2^n + c*n + d:")
    print(f"    a={x2[0]:.6f}, b={x2[1]:.6f}, c={x2[2]:.6f}, d={x2[3]:.6f}")
    for n, s2 in s2_vals:
        pred = x2[0]*4**n + x2[1]*2**n + x2[2]*n + x2[3]
        print(f"    n={n}: actual={s2}, predicted={pred:.1f}")

    # H(1) = sum of all degree sums = H(all backward)
    print(f"\n  H(1,...,1) sequence:")
    for n in sorted(degree_data.keys()):
        h1 = sum(degree_data[n])
        print(f"    n={n}: H(1) = {h1}")

    # H(1) = 3, 5, 9, 17, 31 for n=3..7
    # Check: 2H(1) - 1 = 5, 9, 17, 33, 61
    # Differences: 4, 8, 16, 28
    # These are NOT powers of 2 (28 = 32-4)

# =================================================================
# PART 4: Multilinear polynomial on random tilings — landscape
# =================================================================
def landscape_analysis(n):
    """Analyze the landscape of H(t) on the hypercube."""
    print(f"\n{'='*60}")
    print(f"H(t) LANDSCAPE ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)

    N = 1 << m

    # Find local maxima and minima on the hypercube
    # A tiling t is a local maximum iff H(t) >= H(t XOR e_k) for all k
    local_max = []
    local_min = []
    for t in range(N):
        h = int(H[t])
        is_max = True
        is_min = True
        for k in range(m):
            h_flip = int(H[t ^ (1 << k)])
            if h_flip > h:
                is_max = False
            if h_flip < h:
                is_min = False
        if is_max:
            local_max.append((t, h))
        if is_min:
            local_min.append((t, h))

    print(f"  Local maxima: {len(local_max)} (H values: {sorted(set(h for _,h in local_max), reverse=True)[:10]})")
    print(f"  Local minima: {len(local_min)} (H values: {sorted(set(h for _,h in local_min))[:10]})")
    print(f"  Global max: H = {max(H)}")
    print(f"  Global min: H = {min(H)}")

    # What fraction of tilings are local maxima?
    print(f"  Fraction local max: {len(local_max)}/{N} = {len(local_max)/N:.4f}")
    print(f"  Fraction local min: {len(local_min)}/{N} = {len(local_min)/N:.4f}")

    # The unique minimum is always the transitive (H=1)?
    min_tilings = [t for t, h in local_min if h == min(H)]
    print(f"  Tilings achieving global min (H={min(H)}): {len(min_tilings)}")

    # Number of Hamming distance-1 "uphill" neighbors
    uphill_counts = defaultdict(int)
    for t in range(N):
        h = int(H[t])
        up = sum(1 for k in range(m) if int(H[t ^ (1 << k)]) > h)
        uphill_counts[up] += 1

    print(f"\n  Uphill neighbor distribution:")
    for up in sorted(uphill_counts.keys()):
        print(f"    {up} uphill: {uphill_counts[up]} tilings ({uphill_counts[up]/N*100:.1f}%)")

if __name__ == "__main__":
    for n in [5, 7]:
        max_degree_sign_deep(n)

    for n in [5, 6]:
        multilinear_vs_betti(n)

    n_generating_function()

    for n in [5, 6, 7]:
        landscape_analysis(n)
