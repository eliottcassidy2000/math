#!/usr/bin/env python3
"""
thm477_adversarial_recheck_0611.py — adversarial subagent re-check of THM-477.

FRESH code, deliberately independent of rm_digit_ladder_macmini_0611s1.py:
  - tile enumeration: lexicographic x-major (x=3..n, y=1..x-2), NOT explorer order
  - kernel/image/dual computed by explicit GF(2) Gaussian elimination on the
    (1+sigma) matrix, NOT by orbit counting
  - isGridSym transcribed from tournament-tiling-explorer.html (bit EQUALITY at
    sigma-paired positions, fixed tiles unconstrained, NO negation) and compared
    against ker(1+sigma) membership by brute force over all 2^m tilings
  - affine variant (negation on pairs) tested to rule out "affine code" reading
  - sigma well-definedness WITHOUT sorting; proof that the sort in
    sigma_perm() of rm_digit_ladder is dead code
  - leg-swap / apex / interior-shift checks for the Mode-B recursion
"""
import sys
from math import comb

sys.stdout.reconfigure(line_buffering=True)

def tiles_lex(n):
    """My OWN enumeration: x-major lexicographic."""
    return [(x, y) for x in range(3, n + 1) for y in range(1, x - 1)]

def rref_gf2(rows, m):
    """rows: list of int bitmasks (rows of a matrix with m columns).
    Returns (rank, list of pivot rows in RREF, pivot column list)."""
    rows = [r for r in rows]
    pivots = []
    rr = []
    for col in range(m):
        piv = None
        for i, r in enumerate(rows):
            if (r >> col) & 1:
                piv = i
                break
        if piv is None:
            continue
        p = rows.pop(piv)
        for i in range(len(rows)):
            if (rows[i] >> col) & 1:
                rows[i] ^= p
        for i in range(len(rr)):
            if (rr[i] >> col) & 1:
                rr[i] ^= p
        rr.append(p)
        pivots.append(col)
    return len(rr), rr, pivots

def kernel_basis_gf2(rows, m):
    """Kernel of the linear map x -> (row_i . x)_i. rows = matrix rows as masks."""
    rank, rr, pivots = rref_gf2(rows, m)
    pivset = set(pivots)
    free = [c for c in range(m) if c not in pivset]
    basis = []
    for fcol in free:
        v = 1 << fcol
        # back-substitute: each pivot row determines pivot bit
        for prow, pcol in zip(rr, pivots):
            # value of row.v restricted to free cols already set
            if bin(prow & v).count('1') & 1:
                v ^= 1 << pcol
        basis.append(v)
    return basis

def span_membership_closure(basis):
    """Return set of all vectors in span(basis). Only for small dims."""
    S = {0}
    for b in basis:
        S |= {s ^ b for s in S}
    return S

print('=' * 78)
print('ADVERSARIAL RE-CHECK OF THM-477 (fresh code, independent tile order)')
print('=' * 78)

canon_exponents = {3: 0, 4: -1, 5: -2, 6: -4, 7: -6, 8: -9, 9: -12}
b_dim = {}
all_ok = True

for n in range(3, 10):
    TILES = tiles_lex(n)
    m = len(TILES)
    assert m == comb(n - 1, 2), 'tile count'
    idx = {t: i for i, t in enumerate(TILES)}

    # --- (1) sigma well-defined WITHOUT any sorting -----------------------
    sigma = []
    for (x, y) in TILES:
        sx, sy = n - y + 1, n - x + 1
        # claim: image is ALWAYS a legal tile in (big,small) orientation
        assert sx - sy == x - y, 'sigma preserves skip'
        assert sx - sy >= 2, f'sort would trigger at {(x,y)} n={n}!'
        assert 1 <= sy and sx <= n, 'image leaves staircase'
        assert (sx, sy) in idx, 'image not a tile'
        sigma.append(idx[(sx, sy)])
    assert all(sigma[sigma[i]] == i for i in range(m)), 'not an involution'

    # --- (2) fixed tiles = anti-diagonal, count floor((n-1)/2) ------------
    fixed = [i for i in range(m) if sigma[i] == i]
    f = len(fixed)
    assert f == (n - 1) // 2, f'fixed count {f}'
    assert all(TILES[i][0] + TILES[i][1] == n + 1 for i in fixed)
    antidiag = [i for i in range(m) if TILES[i][0] + TILES[i][1] == n + 1]
    assert set(antidiag) == set(fixed), 'fixed set != anti-diagonal'

    # --- (3) full linear algebra on A = 1 + sigma over GF(2) --------------
    # rows of A: row i = e_i + e_{sigma(i)}
    A_rows = [(1 << i) ^ (1 << sigma[i]) for i in range(m)]
    rankA, _, _ = rref_gf2(A_rows, m)
    kerB = kernel_basis_gf2(A_rows, m)
    dim_ker = len(kerB)
    # image basis: columns of A = same as rows (A symmetric since sigma^T=sigma)
    # but compute image independently as span of {x + sigma(x)} over basis e_i
    im_gens = [(1 << i) ^ (1 << sigma[i]) for i in range(m)]
    rank_im, im_rr, _ = rref_gf2(im_gens, m)
    dim_im = rank_im
    assert dim_ker == (m + f) // 2, f'dim ker {dim_ker} != (m+f)/2'
    assert dim_im == (m - f) // 2, f'dim im {dim_im} != (m-f)/2'
    assert dim_ker - dim_im == f, 'defect != f'
    assert rankA == dim_im and rankA + dim_ker == m, 'rank-nullity'

    # im(1+sigma) subset ker(1+sigma): apply A to every image basis vector
    def apply_A(v):
        out = 0
        for i in range(m):
            if (v >> i) & 1:
                out ^= A_rows[i]
        return out
    assert all(apply_A(g) == 0 for g in im_rr), 'im not inside ker'

    # ker^perp == im : compute ker^perp as kernel of matrix whose rows = kerB
    kperp = kernel_basis_gf2(kerB, m)
    rank_kp, _, _ = rref_gf2(kperp, m)
    # equality of spans: same dim + each kperp basis vector reduces to 0 against im rows
    r_join, _, _ = rref_gf2(kperp + im_rr, m)
    assert rank_kp == dim_im and r_join == dim_im, 'ker^perp != im over F_2'

    # all-ones vector: in ker, not in im (f>=1)
    ones = (1 << m) - 1
    assert apply_A(ones) == 0, '1 not in ker'
    r_with_ones, _, _ = rref_gf2(im_rr + [ones], m)
    if f >= 1:
        assert r_with_ones == dim_im + 1, '1 unexpectedly in im'

    # --- (4) leg swap + apex + interior shift (Mode B) ---------------------
    if n >= 5:
        # bottom row (x,1), 3<=x<=n-1 maps to right column (n, n-x+1)
        for x in range(3, n):
            assert sigma[idx[(x, 1)]] == idx[(n, n - x + 1)], 'leg swap fails'
        assert sigma[idx[(n, 1)]] == idx[(n, 1)], 'apex not fixed'
        # interior tiles (2<=y, x<=n-1) shift to sigma_{n-2}
        SUB = tiles_lex(n - 2)
        sidx = {t: i for i, t in enumerate(SUB)}
        for (x, y) in TILES:
            if y >= 2 and x <= n - 1:
                sx, sy = TILES[sigma[idx[(x, y)]]]
                assert sx <= n - 1 and sy >= 2, 'interior not sigma-closed'
                # shifted: (x-1,y-1) under sigma_{n-2} -> ((n-2)-(y-1)+1, (n-2)-(x-1)+1)
                ex, ey = (n - 2) - (y - 1) + 1, (n - 2) - (x - 1) + 1
                assert (sx - 1, sy - 1) == (ex, ey), 'interior shift mismatch'
    b_dim[n] = dim_ker
    rec = ''
    if n - 2 in b_dim:
        ok = b_dim[n] == b_dim[n - 2] + (n - 2)
        rec = f'  b({n})=b({n-2})+{n-2}: {ok}'
        assert ok, 'recursion fails'

    # --- (5) explorer isGridSym vs ker membership, brute force ------------
    exp = canon_exponents[n]
    assert (f - m) // 2 == exp and (f - m) % 2 == 0, f'exponent {(f-m)/2} != canon {exp}'
    bf = ''
    if n <= 8:
        ker_set = span_membership_closure(kerB) if dim_ker <= 20 else None
        n_blue = 0
        n_affine = 0
        mismatch = 0
        pairs = [(i, sigma[i]) for i in range(m) if i < sigma[i]]
        for bits in range(1 << m):
            # explorer isGridSym: equality at sigma pairs, no negation
            gs = all(((bits >> i) & 1) == ((bits >> j) & 1) for i, j in pairs)
            # affine variant: negation on pairs
            af = all(((bits >> i) & 1) == 1 - ((bits >> j) & 1) for i, j in pairs)
            if gs:
                n_blue += 1
            if af:
                n_affine += 1
            if gs != (bits in ker_set):
                mismatch += 1
        assert mismatch == 0, f'isGridSym != ker membership ({mismatch} mismatches)'
        assert n_blue == 1 << dim_ker, 'blue count != 2^dim ker'
        assert n_blue * 2 ** (-exp) == 1 << m, 'blue fraction != canon exponent'
        # affine variant gives the same COUNT but a DIFFERENT set (coset) when f<m
        affine_is_same_set = (n_affine == n_blue and all(
            all(((b >> i) & 1) == 1 - ((b >> j) & 1) for i, j in pairs) == (b in ker_set)
            for b in range(1 << m))) if pairs else True
        bf = (f'  blue={n_blue}/{1 << m}=2^{exp}  isGridSym==ker: EXACT  '
              f'affine-set==linear-set: {affine_is_same_set}')
    print(f'n={n}: m={m:2d} f={f} dimB={dim_ker:2d} dimB_perp={dim_im:2d} '
          f'defect={dim_ker - dim_im}{rec}{bf}')

print()
print('ALL ASSERTIONS PASSED — THM-477 claims (1)-(4) verified with independent code.')
print('Note: the coordinate sort in rm_digit_ladder sigma_perm() is DEAD CODE')
print('(sx - sy = x - y >= 2 always; asserted above for every tile, n=3..9).')
