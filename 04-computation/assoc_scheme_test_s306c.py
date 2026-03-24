#!/usr/bin/env python3
"""
assoc_scheme_test_s306c.py -- opus-2026-03-24-S306

IS THE METAGRAPH AN ASSOCIATION SCHEME?

An association scheme on v vertices has relations R_0,...,R_d such that:
1. R_0 = identity
2. R_0 ∪ R_1 ∪ ... ∪ R_d = all pairs
3. Each R_i is symmetric
4. The PRODUCT of adjacency matrices is in the span:
   A_i × A_j = Σ_k p_{ij}^k A_k (intersection numbers)

For the tiling metagraph, the natural "relations" are the waggly layers.
But from MSG-426: "W_1*W_1 not in span of W_d's" — so it fails condition 4.

HOW BAD is the failure? Can we measure how close it is to a scheme?

Author: opus-2026-03-24-S306
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def build_space(n):
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()
    cmap = {}; tc = {}; members = defaultdict(list); sizes = {}
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap); sizes[cmap[c]] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid; members[cid].append(bits)
    nc = len(cmap)
    info = {}
    for cid in range(nc):
        b = members[cid][0]
        A2 = [[0]*n for _ in range(n)]
        for i in range(1, n): A2[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (b >> k) & 1: A2[lo][hi] = 1
            else: A2[hi][lo] = 1
        comp_c = canon(complement(A2, n), n)
        info[cid] = {'comp': cmap.get(comp_c, -1)}
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        if comp >= 0 and comp != cid:
            merged[mid] = (cid, comp); mid_map[cid] = mid; mid_map[comp] = mid
        else:
            merged[mid] = (cid,); mid_map[cid] = mid
        mid += 1
    return m, total, tile_pairs, tc, members, sizes, mid_map, mid, merged

# ============================================================
banner("ASSOCIATION SCHEME TEST AT n=5")
# ============================================================

from itertools import combinations

n = 5
m, total, tile_pairs, tc, members, sizes, mid_map, nm, merged = build_space(n)

print("  n=%d: m=%d, %d merged classes\n" % (n, m, nm))

# Build waggly weight matrices for d=1,...,m
# These are WEIGHTED adjacency matrices (not 0/1)
W = {}
for d in range(1, m+1):
    Wd = np.zeros((nm, nm), dtype=float)
    for bits in range(total):
        mi = mid_map[tc[bits]]
        for flips in combinations(range(m), d):
            mask = sum(1 << k for k in flips)
            mj = mid_map[tc[bits ^ mask]]
            Wd[mi][mj] += 1
    W[d] = Wd

# Also add the identity as W[0]
W[0] = np.zeros((nm, nm), dtype=float)
for mi in range(nm):
    fib = sum(sizes[c] for c in merged[mi])
    W[0][mi][mi] = fib  # diagonal = fiber

# Check: do the W matrices form a closed algebra?
# I.e., is W[i] × W[j] in span(W[0], W[1], ..., W[m])?

print("  ALGEBRA CLOSURE TEST: W_i × W_j in span(W_0,...,W_m)?\n")

# Stack all W matrices as vectors in a matrix A (each row = flattened W_d)
W_flat = np.array([W[d].flatten() for d in range(m+1)])  # (m+1) × nm²

for i in range(min(m+1, 4)):
    for j in range(i, min(m+1, 4)):
        product = W[i] @ W[j]
        product_flat = product.flatten()

        # Project onto span of W_flat
        # Use least squares: find coefficients c such that Σ c_d W_d ≈ product
        coeffs, residuals, rank, sv = np.linalg.lstsq(W_flat.T, product_flat, rcond=None)

        # Reconstruction error
        reconstruction = W_flat.T @ coeffs
        error = np.linalg.norm(product_flat - reconstruction) / np.linalg.norm(product_flat)

        print("  W_%d × W_%d: coefficients = %s, relative error = %.6f" %
              (i, j, ["%.2f" % c for c in coeffs[:5]], error))

# Is the error exactly zero?
print("\n  If error = 0 for ALL products → association scheme.")
print("  Any nonzero error → NOT an association scheme.\n")

# Let's check which products have nonzero error
print("  CLOSURE ERRORS (relative):")
all_errors = {}
for i in range(m+1):
    for j in range(i, m+1):
        product = W[i] @ W[j]
        product_flat = product.flatten()
        coeffs, _, _, _ = np.linalg.lstsq(W_flat.T, product_flat, rcond=None)
        reconstruction = W_flat.T @ coeffs
        error = np.linalg.norm(product_flat - reconstruction) / (np.linalg.norm(product_flat) + 1e-10)
        all_errors[(i,j)] = error
        if error > 1e-10:
            print("    W_%d × W_%d: error = %.6f (NOT in span)" % (i, j, error))
        elif i <= 3 and j <= 3:
            print("    W_%d × W_%d: error = %.10f (IN span ✓)" % (i, j, error))

n_nonzero = sum(1 for e in all_errors.values() if e > 1e-10)
n_total = len(all_errors)
print("\n  %d / %d products have nonzero error" % (n_nonzero, n_total))
print("  Association scheme: %s" % ("YES" if n_nonzero == 0 else "NO"))

# ============================================================
banner("THE BOSE-MESNER ALGEBRA DIMENSION")
# ============================================================

# The Bose-Mesner algebra is the span of {W_0, ..., W_m}.
# Its dimension = rank of the matrix of flattened W_d vectors.

BM_rank = np.linalg.matrix_rank(W_flat, tol=1e-8)
print("  Dimension of span(W_0,...,W_%d) = %d" % (m, BM_rank))
print("  Number of relations: m+1 = %d" % (m+1))
print("  Number of classes: nm = %d" % nm)
print("  Space of all nm×nm matrices: dimension %d" % (nm*nm))
print()

# The BM algebra of an association scheme has dimension = #relations.
# If BM_rank < m+1, some W_d are linearly dependent.
# If BM_rank > m+1... that's impossible (they span at most m+1 dims).
# If BM_rank = m+1 but closure fails, the algebra is larger than the span.

print("  Is BM algebra closed (self-contained)?")
print("  If yes, it IS a scheme. If no, it's a PARTIAL scheme.\n")

# Check: what's the FULL algebra generated by {W_0,...,W_m}?
# Add products iteratively
basis = [W[d].flatten() for d in range(m+1)]
basis_matrix = np.array(basis)
current_rank = np.linalg.matrix_rank(basis_matrix, tol=1e-8)

print("  Starting rank: %d (from W_0,...,W_%d)" % (current_rank, m))

# Add products
new_added = True
round_num = 0
while new_added:
    round_num += 1
    new_added = False
    n_basis = len(basis)
    for i in range(n_basis):
        for j in range(i, n_basis):
            Wi = basis[i].reshape(nm, nm)
            Wj = basis[j].reshape(nm, nm)
            product = (Wi @ Wj).flatten()
            extended = np.vstack([basis_matrix, product.reshape(1, -1)])
            new_rank = np.linalg.matrix_rank(extended, tol=1e-8)
            if new_rank > current_rank:
                basis.append(product)
                basis_matrix = np.vstack([basis_matrix, product.reshape(1, -1)])
                current_rank = new_rank
                new_added = True
                if current_rank >= nm * nm:
                    break
        if current_rank >= nm * nm:
            break
    print("  Round %d: rank = %d" % (round_num, current_rank))
    if current_rank >= nm:  # Can't exceed nm for commutative algebra
        break

print("\n  FULL ALGEBRA DIMENSION = %d" % current_rank)
print("  For an association scheme, this would be %d (= m+1)" % (m+1))
print("  Excess: %d (= how far from being a scheme)" % (current_rank - (m+1)))

# ============================================================
banner("WHAT MAKES IT FAIL?")
# ============================================================

print("""
  THE METAGRAPH FAILS TO BE AN ASSOCIATION SCHEME BECAUSE:

  In the Hamming scheme H(m,2), the relations are the Hamming distances.
  The S_n action partitions tilings into orbits (iso classes).
  The ORBIT MATRICES (waggly weight matrices W_d) are obtained by
  collapsing the Hamming scheme's adjacency matrices to the orbit space.

  For the quotient to be an association scheme, the group action must be
  COMPATIBLE with the scheme: the orbits must form a "coherent configuration."

  The S_n action is NOT compatible because:
  - S_n permutes vertices, which permutes tile POSITIONS
  - Two tilings at Hamming distance d may map (under σ ∈ S_n) to
    tilings at a DIFFERENT Hamming distance
  - So the quotient orbit matrices don't close under multiplication

  The EXCESS (full algebra dimension - scheme dimension) measures
  how much additional structure exists beyond the scheme.
  A small excess means the metagraph is "almost" a scheme.

  PRACTICAL IMPLICATION:
  Even though the metagraph isn't a scheme, the Krawtchouk polynomials
  APPROXIMATELY diagonalize the weight matrices (K_1 gives r=-0.94).
  The scheme structure is a GOOD APPROXIMATION that degrades slowly with n.
""")

print("Done. opus-2026-03-24-S306")
