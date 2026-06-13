#!/usr/bin/env python3
"""
trace_F_analysis_s215.py -- opus-2026-03-23-S215

Analyze Tr(F) = sum_C F[C][C] = total self-reversible labeled pairs.
Compare with SL to understand the non-identity Burnside corrections.

SL = (1/n!) * [Tr(F) + non-identity corrections]

At n≤4: corrections = 0 (Tr(F)/n! = SL exactly)
At n≥5: corrections > 0

Also: compute Tr(F) at n=6 where we have the full F matrix.

Author: opus-2026-03-23-S215
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

def canonical(A, n):
    if n <= 5:
        best = None
        for perm in permutations(range(n)):
            form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or form < best: best = form
        return best
    else:
        scores = [sum(A[i]) for i in range(n)]
        score_groups = defaultdict(list)
        for v in range(n): score_groups[scores[v]].append(v)
        groups = [score_groups[s] for s in sorted(score_groups.keys())]
        best = None
        def gen_perms(grps):
            if not grps:
                yield []
                return
            for p in permutations(grps[0]):
                for rest in gen_perms(grps[1:]):
                    yield list(p) + rest
        for perm in gen_perms(groups):
            form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or form < best: best = form
        return best

banner("Tr(F) AND SL AT n=3..6")

SL_known = {3:2, 4:6, 5:16, 6:58}

for n in [3, 4, 5, 6]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    cmap = {}; sizes = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap:
            cmap[c] = len(cmap); sizes[len(cmap)-1] = 0
        sizes[cmap[c]] += 1
        tc[bits] = cmap[c]
    nc = len(cmap)

    F_diag = np.zeros(nc, dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            bits2 = bits ^ (1<<k)
            cj = tc[bits2]
            if ci == cj:
                F_diag[ci] += 1

    trF = int(np.sum(F_diag))
    sl = SL_known[n]
    nfact = factorial(n)

    correction = sl * nfact - trF

    # SL_trivial = sum of F[C][C]/size(C) for |Aut|=1 classes
    sl_trivial = 0
    sl_nontrivial = 0
    for cid in range(nc):
        aut_size = nfact // sizes[cid]
        self_arcs = F_diag[cid] // sizes[cid]  # arcs per T targeting same class
        if aut_size == 1:
            sl_trivial += self_arcs  # each arc = one orbit
        # For nontrivial: orbits = fewer than arcs

    print("  n=%d:" % n)
    print("    nc=%d, Tr(F)=%d, n!=%d" % (nc, trF, nfact))
    print("    Tr(F)/n! = %d/%d = %.4f" % (trF, nfact, trF/nfact))
    print("    SL = %d" % sl)
    print("    Non-id correction = SL*n! - Tr(F) = %d*%d - %d = %d" %
          (sl, nfact, trF, correction))
    print("    SL - Tr(F)/n! = %.4f" % (sl - trF/nfact))

    # Decompose: which classes have high F[C][C]?
    top_self = sorted(range(nc), key=lambda c: F_diag[c]//sizes[c], reverse=True)[:5]
    print("    Top self-loop classes:")
    for cid in top_self:
        arcs_self = F_diag[cid] // sizes[cid]
        aut = nfact // sizes[cid]
        if arcs_self > 0:
            print("      class %d: arcs_to_self=%d, |Aut|=%d, size=%d" %
                  (cid, arcs_self, aut, sizes[cid]))
    print()

# ============================================================
banner("SELF-REVERSIBLE ARC FRACTION")
# ============================================================

# What fraction of arcs in a random tournament are self-reversible?
# Fix(id) / (2^m * m) = Tr(F) / (2^m * m)

print("  Self-reversible fraction = Tr(F) / (2^m * m):\n")
TrF_known = {}
for n in [3, 4, 5, 6]:
    m = comb(n, 2)
    total_pairs = (2**m) * m
    # Recompute Tr(F)
    pairs_n = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs_n):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]

    trF = 0
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            if tc[bits ^ (1<<k)] == ci:
                trF += 1

    TrF_known[n] = trF
    frac = trF / total_pairs
    print("  n=%d: Tr(F)=%d, total=(2^%d)*%d=%d, fraction=%.6f" %
          (n, trF, m, m, total_pairs, frac))

# ============================================================
banner("Tr(F) SEQUENCE")
# ============================================================

print("  Tr(F) values: %s" % [TrF_known[n] for n in [3,4,5,6]])

# Tr(F) = total self-reversible (T, e) pairs (labeled)
# = sum_C F[C][C] = sum_C size(C) * arcs_to_self_per_T(C)

# As Burnside: Tr(F) = n! * SL_approx, where SL_approx ignores
# non-identity corrections and orbit grouping for nontrivial Aut.

# The KEY identity: SL = Tr(F)/n! + (corrections from non-id and Aut)
# At n=3,4: corrections = 0.
# At n=5,6: corrections = 16 - 1760/120 = 1.33, 58 - 37920/720 = 5.33.

# These corrections * n! = 160, 3840.

# Can we understand the corrections?
for n in [3, 4, 5, 6]:
    nfact = factorial(n)
    corr = SL_known[n] * nfact - TrF_known[n]
    print("  n=%d: correction = %d (= SL*n! - Tr(F))" % (n, corr))
    if corr > 0:
        # This comes from permutations sigma of ORDER 3+ with ≥2 fixed points
        # that fix tournaments and arcs where the arc is self-reversible.
        # At n=5: only cycle types (1,1,3) can contribute (2 fixed pts, 3-cycle).
        # There are C(5,2)*2 = 20 such permutations.
        # Correction = 160, so each contributes 160/20 = 8 on average.
        pass

print()

# ============================================================
banner("THE NATURE OF THE CORRECTION")
# ============================================================

print("""
  The Burnside correction comes from non-identity sigma with ≥2 fixed points
  AND odd cycle lengths (cycle types containing only odd cycles and ≥2 fixed points).

  Even cycles CANNOT contribute because sigma-swapping 2 vertices means
  the edge between them has no consistent orientation.

  So the only contributing cycle types are those with:
  - All odd cycles (including 1-cycles for fixed points)
  - At least 2 fixed points (1-cycles)

  At n=3: only type (1,1,1) = identity. No corrections.
  At n=4: type (1,1,1,1) = identity. Type (1,1,2) has a 2-cycle, excluded.
          Type (1,3) has 1 fixed point, not enough. No corrections.
  At n=5: types with all odd cycles AND ≥2 fixed pts:
          (1,1,1,1,1) = identity (already counted)
          (1,1,3) = 2 fixed + 3-cycle. CONTRIBUTES.
          Correction = 160.
  At n=6: (1,1,1,1,1,1) = identity
          (1,1,1,3) = 3 fixed + 3-cycle
          (1,1,1,1,... but all cycle lengths odd)
          (1,1,1,3): yes, CONTRIBUTES.
          (3,3): 0 fixed, excluded.
          (1,5): 1 fixed, excluded.
          Correction = 3840.

  So the correction at n comes from cycle types (1^k, c_1, c_2, ...)
  where k ≥ 2 and all c_i are ODD and > 1.

  At n=5: only (1,1,3). 20 permutations. Average contribution: 160/20 = 8.
  At n=6: only (1,1,1,3). 40 permutations. Average: 3840/40 = 96.

  CONJECTURE: The correction is computable from the number of
  "sigma-fixed self-reversible arcs" for each contributing cycle type.
  This involves tournaments on the fixed-point subset that are
  extendable to sigma-fixed tournaments on all n vertices.
""")

# How many permutations of each contributing type?
from math import comb as C
for n in range(3, 10):
    types = []
    # Generate all partitions of n with all odd parts and ≥2 parts equal to 1
    # Simple: k fixed points + rest in odd cycles ≥3
    for k in range(2, n+1):
        rem = n - k
        if rem == 0:
            types.append(((1,)*k, factorial(n) // (factorial(k) * 1)))
        elif rem >= 3 and rem % 2 == 1:
            # Single odd cycle of length rem
            n_perms = C(n, k) * factorial(k-1) * factorial(rem-1)  # wrong
            # Actually: choose k fixed points, then arrange remaining in one odd cycle
            n_perms = C(n, k) * factorial(rem-1)  # (rem-1)! ways to arrange cycle
            types.append(("1^%d + %d-cycle" % (k, rem), n_perms))
        # Could have multiple odd cycles, but let's start with single cycle
    if types:
        print("  n=%d: contributing types: %s" % (n, types))

print("\nDone. opus-2026-03-23-S215")
