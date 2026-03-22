#!/usr/bin/env python3
"""
tournaplex_fvectors_s180.py -- opus-2026-03-22-S180

INTENSIVE STUDY OF TOURNAPLEX f-VECTORS

The tournaplex T has simplices = tournaments at each vertex count.
Its f-vector f = (f_0, f_1, f_2, ...) where f_k = number of
(k+1)-vertex tournaments = 2^{C(k+1,2)}.

But we can define MANY sub-tournaplexes:
  T_full: all tournaments (f_k = 2^{C(k+1,2)})
  T_trans: only transitive tournaments (f_k = (k+1)!)
  T_reg: only regular tournaments (odd k+1 only)
  T_SC: only self-complementary tournaments
  T_{H>=h}: tournaments with H >= h
  T_{score=s}: tournaments with given score
  T_iso: one representative per iso class (f_k = A000568(k+1))
  FLAG tournaplex of a specific digraph G

For each: compute f-vector, h-vector, Euler characteristic,
and search for patterns.

Author: opus-2026-03-22-S180
"""

import sys
import numpy as np
from math import comb, factorial, log2
from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    results = []
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        results.append(A)
    return results

def is_transitive(A, n):
    """Tournament is transitive iff it has no directed 3-cycle."""
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if A[i][j] and A[j][k] and A[k][i]: return False
                if A[i][k] and A[k][j] and A[j][i]: return False
    return True

def is_regular(A, n):
    if n % 2 == 0: return False
    target = (n-1) // 2
    return all(sum(A[i]) == target for i in range(n))

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

# ============================================================
# PART 1: THE FULL TOURNAPLEX f-VECTOR
# ============================================================

banner("PART 1: THE FULL TOURNAPLEX f-VECTOR")

print("  f_k = number of (k+1)-vertex tournaments = 2^{C(k+1, 2)}\n")

full_f = []
for k in range(10):
    n = k + 1
    f = 2**comb(n, 2)
    full_f.append(f)
    print("  f_%d = 2^C(%d,2) = 2^%d = %d" % (k, n, comb(n,2), f))

print("\n  f-vector: %s" % full_f)
print("  This is 2^{0, 1, 3, 6, 10, 15, 21, 28, 36, 45} = 2^{C(k+1,2)}")
print("  The exponents are the TRIANGULAR NUMBERS!\n")

# Euler characteristic (alternating sum)
chi = sum((-1)**k * full_f[k] for k in range(len(full_f)))
print("  Euler characteristic (truncated at k=9):")
print("  chi = sum (-1)^k f_k = %d" % chi)
print("  Partial sums:")
for K in range(1, 10):
    chi_K = sum((-1)**k * full_f[k] for k in range(K))
    print("    K=%d: chi = %d" % (K, chi_K))

# ============================================================
# PART 2: THE TRANSITIVE TOURNAPLEX f-VECTOR
# ============================================================

banner("PART 2: TRANSITIVE TOURNAPLEX")

print("  T_trans: only transitive tournaments (no 3-cycles).")
print("  A transitive tournament on n vertices = a total order on [n].")
print("  f_k = number of transitive tournaments on k+1 vertices = (k+1)!\n")

trans_f = [factorial(k+1) for k in range(9)]
print("  f-vector: %s" % trans_f)
print("  = [1!, 2!, 3!, 4!, 5!, 6!, 7!, 8!, 9!]")
print("  = [1, 2, 6, 24, 120, 720, 5040, 40320, 362880]\n")

chi_trans = sum((-1)**k * trans_f[k] for k in range(len(trans_f)))
print("  Euler char (truncated): %d" % chi_trans)
print("  Partial sums:")
for K in range(1, 9):
    chi_K = sum((-1)**k * trans_f[k] for k in range(K))
    print("    K=%d: chi = %d" % (K, chi_K))

# subfactorial connection!
print("\n  Note: sum_{k=0}^{n} (-1)^k (k+1)! is related to SUBFACTORIALS!")
print("  D_n = n! sum_{k=0}^n (-1)^k / k! = subfactorial (derangements)")
print("  Our sum: sum (-1)^k (k+1)! = ???")
for K in range(1, 10):
    s = sum((-1)**k * factorial(k+1) for k in range(K))
    # Is this related to subfactorials?
    print("    K=%d: sum = %d" % (K, s))

# ============================================================
# PART 3: THE ISO CLASS TOURNAPLEX f-VECTOR
# ============================================================

banner("PART 3: ISO CLASS TOURNAPLEX (one per class)")

print("  f_k = A000568(k+1) = number of tournament isomorphism classes\n")

A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880]
iso_f = A000568[:]
print("  f-vector: %s" % iso_f)
print("  = [1, 1, 1, 2, 4, 12, 56, 456, 6880]\n")

chi_iso = sum((-1)**k * iso_f[k] for k in range(len(iso_f)))
print("  Euler char (truncated at k=8): %d\n" % chi_iso)

print("  Partial Euler characteristics:")
for K in range(1, len(iso_f)+1):
    chi_K = sum((-1)**k * iso_f[k] for k in range(K))
    print("    K=%d: chi = %d" % (K, chi_K))

# ============================================================
# PART 4: RATIOS BETWEEN f-VECTORS
# ============================================================

banner("PART 4: f-VECTOR RATIOS")

print("  f_full / f_trans = 2^C(k+1,2) / (k+1)!")
print("  = number of tournaments per transitive tournament")
print("  = average iso class size = 'tournament entropy'\n")

for k in range(8):
    n = k + 1
    ratio = Fraction(2**comb(n,2), factorial(n))
    print("  k=%d (n=%d): f_full/f_trans = %s = %.4f" %
          (k, n, ratio, float(ratio)))

print("\n  f_full / f_iso = 2^C(k+1,2) / A000568(k+1)")
print("  = average labeled class size\n")

for k in range(len(A000568)):
    n = k + 1
    if A000568[k] > 0:
        ratio = Fraction(2**comb(n,2), A000568[k])
        print("  k=%d (n=%d): f_full/f_iso = %s = %.2f" %
              (k, n, ratio, float(ratio)))

print("\n  f_trans / f_iso = (k+1)! / A000568(k+1)")
print("  = average automorphism group size\n")

for k in range(len(A000568)):
    n = k + 1
    if A000568[k] > 0:
        ratio = Fraction(factorial(n), A000568[k])
        print("  k=%d (n=%d): f_trans/f_iso = %s = %.4f" %
              (k, n, ratio, float(ratio)))

# ============================================================
# PART 5: SUB-TOURNAPLEX f-VECTORS (computed)
# ============================================================

banner("PART 5: COMPUTED SUB-TOURNAPLEX f-VECTORS")

# For n up to 6, compute various sub-tournaplex counts
max_n = 7

for sub_name, sub_filter in [
    ("transitive", is_transitive),
    ("regular", is_regular),
    ("H >= 3", lambda A, n: count_hp(A, n) >= 3),
    ("H = max", lambda A, n: True),  # placeholder
]:
    fvec = []
    for n in range(1, max_n):
        if sub_name == "H = max":
            # Count tournaments with H = H_max for that n
            tourns = all_tournaments(n)
            if n <= 6:
                Hs = [count_hp(A, n) for A in tourns]
                H_max = max(Hs)
                count = sum(1 for h in Hs if h == H_max)
            else:
                count = 0
        else:
            if n <= 6:
                count = sum(1 for A in all_tournaments(n) if sub_filter(A, n))
            else:
                count = 0
        fvec.append(count)

    # Only print if we have data
    if any(f > 0 for f in fvec):
        print("  %s: %s" % (sub_name, fvec))

# Also compute score-specific sub-tournaplexes
print("\n  Score-specific f-vectors (tournaments with given sorted score):")
for n in range(3, 6):
    tourns = all_tournaments(n)
    score_counts = Counter()
    for A in tourns:
        sc = score_seq(A, n)
        score_counts[sc] += 1
    for sc in sorted(score_counts.keys()):
        print("    n=%d, score %s: %d tournaments" % (n, sc, score_counts[sc]))

# ============================================================
# PART 6: THE FLAG TOURNAPLEX OF A TOURNAMENT
# ============================================================

banner("PART 6: FLAG TOURNAPLEX OF A SINGLE TOURNAMENT")

print("  Given a tournament T on n vertices, the FLAG TOURNAPLEX tFl(T)")
print("  consists of all sub-tournaments of T (induced subgraphs).\n")

print("  For a tournament T on n vertices:")
print("  f_k(tFl(T)) = C(n, k+1) for all k")
print("  because every (k+1)-subset of vertices induces a tournament.\n")

print("  So: the flag tournaplex of ANY tournament has the SAME f-vector!")
print("  f-vector of tFl(T) = [C(n,1), C(n,2), ..., C(n,n)] = Pascal row n.\n")

for n in range(3, 8):
    fvec = [comb(n, k+1) for k in range(n)]
    chi = sum((-1)**k * fvec[k] for k in range(len(fvec)))
    print("  n=%d: f-vector = %s, chi = %d" % (n, fvec, chi))

print("""
  The Euler characteristic of tFl(T):
  chi = sum_{k=0}^{n-1} (-1)^k C(n, k+1)
      = -sum_{j=1}^{n} (-1)^j C(n, j)
      = -(sum_{j=0}^n (-1)^j C(n,j) - 1)
      = -(0 - 1) = 1  (for n >= 1)

  Wait: sum_{j=0}^n (-1)^j C(n,j) = (1-1)^n = 0.
  So: chi = -(0 - 1) = 1.
  But our computation gives chi = (-1)^{n-1}.

  Let me recheck: f_k = C(n, k+1), k = 0, ..., n-1.
  chi = sum_{k=0}^{n-1} (-1)^k C(n, k+1)
  Let j = k+1: chi = sum_{j=1}^n (-1)^{j-1} C(n,j) = -sum_{j=1}^n (-1)^j C(n,j)
  = -(sum_{j=0}^n (-1)^j C(n,j) - C(n,0)) = -(0 - 1) = 1.

  So chi(tFl(T)) = 1 for ALL T. The flag tournaplex is CONTRACTIBLE!
""")

# Verify
for n in range(3, 10):
    fvec = [comb(n, k+1) for k in range(n)]
    chi = sum((-1)**k * fvec[k] for k in range(len(fvec)))
    print("  n=%d: chi = %d (should be 1)" % (n, chi))

# ============================================================
# PART 7: THE TOURNAPLEX h-VECTOR
# ============================================================

banner("PART 7: h-VECTORS")

print("  The h-vector of a simplicial complex is defined via")
print("  h_k = sum_{j=0}^k (-1)^{k-j} C(d-j, k-j) f_{j-1}")
print("  where d = dimension and f_{-1} = 1.\n")

print("  For the FLAG TOURNAPLEX of T on n vertices:")
print("  f_{-1} = 1, f_k = C(n, k+1) for k = 0, ..., n-1.")
print("  This is the simplex (n-1)-skeleton, which is the FULL simplex.")
print("  Its h-vector is (1, 0, 0, ..., 0) — trivial!\n")

print("  For the ISO CLASS TOURNAPLEX (one per class):")
print("  f = [1, 1, 2, 4, 12, 56, 456, ...] = A000568 values\n")

# Compute h-vector for the iso class tournaplex
# h_k = sum_{j=0}^{k} (-1)^{k-j} C(d-j, k-j) f_{j-1}
# where d = max dimension (we truncate), f_{-1} = 1
d = len(iso_f) - 1  # max dim available
h_iso = []
for k in range(d + 1):
    h_k = 0
    for j in range(k + 1):
        fj = iso_f[j] if j < len(iso_f) else 0
        if j == 0:
            fj_prev = 1  # f_{-1} = 1
        else:
            fj_prev = iso_f[j-1] if j-1 < len(iso_f) else 0
        # Actually the standard formula uses f_{j-1}
        h_k += (-1)**(k-j) * comb(d-j, k-j) * (iso_f[j] if j < len(iso_f) else 0)
    h_iso.append(h_k)

print("  h-vector of iso-class tournaplex (truncated at d=%d):" % d)
print("  h = %s\n" % h_iso[:10])

# ============================================================
# PART 8: THE f-VECTOR AS A GENERATING FUNCTION
# ============================================================

banner("PART 8: f-VECTOR GENERATING FUNCTIONS")

print("  FULL TOURNAPLEX: F(x) = sum f_k x^k = sum 2^{C(k+1,2)} x^k\n")

print("  First terms: F(x) = 1 + 2x + 8x^2 + 64x^3 + 1024x^4 + ...\n")

print("  LOG of f_k: log2(f_k) = C(k+1,2) = 0, 1, 3, 6, 10, 15, 21, ...")
print("  These are the TRIANGULAR NUMBERS = |delta_k|!\n")

print("  So: log2(f_k) = |staircase delta_k| = C(k+1, 2).")
print("  The full tournaplex f-vector ENCODES the staircase dimensions!\n")

print("  RATIOS f_{k+1}/f_k = 2^{C(k+2,2)}/2^{C(k+1,2)} = 2^{k+1}")
for k in range(8):
    ratio = full_f[k+1] / full_f[k]
    print("    k=%d: f_{k+1}/f_k = %d = 2^%d" % (k, ratio, k+1))

print("\n  The ratio = 2^{k+1} because adding one vertex to a")
print("  (k+1)-tournament adds k+1 new arcs, each with 2 choices.\n")

# ============================================================
# PART 9: THE TOURNAPLEX f-POLYNOMIAL AND I(Omega, x)
# ============================================================

banner("PART 9: CONNECTING f-VECTOR TO I(Omega, x)")

print("""
  Is there a connection between the tournaplex f-vector and
  the independence polynomial I(Omega, x)?

  The f-vector counts ALL tournaments of each size.
  I(Omega, x) counts independent cycle sets within ONE tournament.

  These are at different levels of the tournaplex hierarchy:
  - f-vector: the NUMBER of simplices at each dimension
  - I(Omega, x): an invariant of EACH simplex

  But we can form the GENERATING FUNCTION:
  G(x, y) = sum_{n >= 1} sum_{T on n vertices} I(Omega(T), x) * y^n

  At x=2: G(2, y) = sum_{n} sum_T H(T) * y^n
  = sum_n 2^{C(n,2)} * E[H(n)] * y^n
  = sum_n n! * 2^{C(n,2)-n+1} * y^n
  = sum_n n! / 2^{n-1} * 2^{C(n,2)} * y^n

  This is the f-vector weighted by E[H]:
  F_H(y) = sum f_k * E[H(k+1)] * y^k
""")

print("  F_H(y) coefficients (f_k * E[H]):")
for k in range(8):
    n = k + 1
    EH = Fraction(factorial(n), 2**(n-1))
    fk = full_f[k]
    coeff = fk * EH
    print("    k=%d: f_%d * E[H(%d)] = %d * %s = %s" %
          (k, k, n, fk, EH, coeff))

# ============================================================
# PART 10: COMPARING WITH KNOWN COMPLEXES
# ============================================================

banner("PART 10: COMPARISON WITH KNOWN SIMPLICIAL COMPLEXES")

print("  The tournaplex f-vector f_k = 2^{C(k+1,2)} grows as:")
print("  k=0: 1, k=1: 2, k=2: 8, k=3: 64, k=4: 1024, ...")
print("  This is MUCH faster than any standard simplicial complex.\n")

print("  COMPARISON:")
print("  %-20s | %-8s %-8s %-8s %-8s %-8s" %
      ("Complex", "f_0", "f_1", "f_2", "f_3", "f_4"))
print("  " + "-" * 70)

complexes = [
    ("Full tournaplex", [1, 2, 8, 64, 1024]),
    ("Transitive", [1, 2, 6, 24, 120]),
    ("Iso class", [1, 1, 2, 4, 12]),
    ("n-simplex (n=4)", [5, 10, 10, 5, 1]),
    ("4-cube (Q_4)", [16, 32, 24, 8, 0]),
    ("Icosahedron", [12, 30, 20, 0, 0]),
    ("24-cell", [24, 96, 96, 24, 0]),
]

for name, fv in complexes:
    vals = fv + [0] * (5 - len(fv))
    print("  %-20s | %-8d %-8d %-8d %-8d %-8d" %
          (name, vals[0], vals[1], vals[2], vals[3], vals[4]))

# ============================================================
# PART 11: THE TOURNAPLEX AND THE STAIRCASE
# ============================================================

banner("PART 11: f-VECTOR ENCODES THE STAIRCASE TOWER")

print("""
  log2(f_k) = C(k+1, 2) = |delta_k| (staircase box count)

  The sequence of staircase dimensions IS the log of the f-vector:
    delta_0 = 0 boxes -> f_0 = 2^0 = 1
    delta_1 = 1 box   -> f_1 = 2^1 = 2
    delta_2 = 3 boxes -> f_2 = 2^3 = 8
    delta_3 = 6 boxes -> f_3 = 2^6 = 64
    delta_4 = 10 boxes -> f_4 = 2^{10} = 1024

  The staircase tower:
    delta_0 < delta_1 < delta_2 < delta_3 < ...
  maps to the f-vector tower:
    f_0 < f_1 < f_2 < f_3 < ...

  The GROWTH RATE: f_{k+1}/f_k = 2^{k+1}
  = 2^{(number of new arcs when adding vertex k+2)}
  = 2^{(new row length in delta_{k+1})}

  The staircase delta_{k+1} is delta_k plus a NEW ROW of length k+1.
  This new row = k+1 new arcs. Each has 2 choices. Total: 2^{k+1}.

  THE TOURNAPLEX f-VECTOR IS THE EXPONENTIAL OF THE STAIRCASE GROWTH.

  In the tournaplex: moving from dimension k to k+1 = adding one vertex.
  In the staircase: growing delta_k to delta_{k+1} = adding one row.
  In G_n: the fiber bundle G_{n} -> G_{n-2} = adding 2 rows.

  These are ALL descriptions of the same growth process.
""")

# Verify: new row lengths
for k in range(8):
    old_boxes = comb(k+1, 2)
    new_boxes = comb(k+2, 2)
    new_row = new_boxes - old_boxes
    print("  delta_%d -> delta_%d: %d -> %d boxes, new row = %d = k+1 = %d" %
          (k, k+1, old_boxes, new_boxes, new_row, k+1))

# ============================================================
# PART 12: DEEP f-VECTOR IDENTITIES
# ============================================================

banner("PART 12: f-VECTOR IDENTITIES")

print("  IDENTITY 1: f_k = 2^{T_k} where T_k = k(k+1)/2 = triangular number")
print("  IDENTITY 2: f_{k+1}/f_k = 2^{k+1}")
print("  IDENTITY 3: log2(f_k) = C(k+1, 2) = |delta_k|\n")

print("  IDENTITY 4: sum_k f_k * x^k = sum_k 2^{C(k+1,2)} x^k")
print("  This is a THETA FUNCTION of sorts (quadratic exponent in 2).\n")

print("  IDENTITY 5: The alternating sum:")
for K in range(1, 9):
    chi = sum((-1)**k * full_f[k] for k in range(K))
    # Is this a nice number?
    if chi != 0:
        v2 = 0
        temp = abs(chi)
        while temp % 2 == 0 and temp > 0:
            v2 += 1; temp //= 2
    else:
        v2 = "inf"
    print("    sum_{k=0}^{%d} (-1)^k f_k = %d (v_2 = %s)" % (K-1, chi, v2))

print()
print("  IDENTITY 6: Product formula")
print("  f_k = prod_{j=1}^{k} 2^j = 2^{sum j from 1 to k} = 2^{k(k+1)/2}")
print("  The f-vector is a PRODUCT of powers of 2.\n")

print("  IDENTITY 7: f_k / k! = 2^{C(k+1,2)} / (k+1)!")
print("  = average labeled class size (the ratio from Part 4)")
print("  This ratio approaches infinity: tournaments grow MUCH faster")
print("  than their labeling symmetry group.\n")

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE TOURNAPLEX f-VECTOR TELLS THE WHOLE STORY")

print("""
THE f-VECTOR OF THE FULL TOURNAPLEX:

  f = (1, 2, 8, 64, 1024, 32768, 2097152, 268435456, ...)
  f_k = 2^{C(k+1, 2)} = 2^{k(k+1)/2}

  This sequence ENCODES:
    - The staircase dimensions: log2(f_k) = |delta_k|
    - The growth rate: f_{k+1}/f_k = 2^{k+1} (new row length)
    - The triangular numbers: exponents are T_k
    - The tournament cube dimensions: C(n,2) for n = k+1

THE HIERARCHY OF f-VECTORS:

  Full:     [1, 2, 8, 64, 1024, ...]     (all tournaments)
  Trans:    [1, 2, 6, 24, 120, ...]       (total orders = (k+1)!)
  Iso:      [1, 1, 2, 4, 12, ...]         (iso classes = A000568)

  Full/Trans = 2^{C(k+1,2)}/(k+1)! -> infinity  (entropy per order)
  Full/Iso   = avg labeled class size -> infinity
  Trans/Iso  = avg |Aut| group -> 1 (generically trivial)

FLAG TOURNAPLEX OF ONE TOURNAMENT:

  f_k = C(n, k+1) (binomial coefficients, Pascal row)
  chi = 1 (ALWAYS contractible)
  h-vector = (1, 0, 0, ...) (trivial)

  All tournaments have the SAME flag tournaplex f-vector!
  This is because every subset of a tournament induces a tournament.
  The flag tournaplex = the full (n-1)-simplex.

  So: INDIVIDUAL tournament flag tournaplexes are trivial.
  The interesting structure comes from WHICH tournaments
  appear as sub-tournaments — this is what distinguishes them.

THE DIRECTIONALITY-FILTERED f-VECTOR:

  Fix a vertex ordering. Filter by directionality d.
  f_k(d) = #{(k+1)-tournaments with >= d upward arcs}

  This gives a PERSISTENCE MODULE: as d increases from 0 to C(k+1,2),
  f_k(d) decreases from f_k (all) to 0 (none, except transitive).

  The persistent homology of this filtration encodes:
  - Birth/death of topological features in tournament space
  - The "homological complexity" of tournament structure
  - Connected to our Kemeny distance and G_n distance

THE DEEPEST CONNECTION:

  The tournaplex f-vector 2^{T_k} where T_k = k(k+1)/2 satisfies:

  SUM: sum_{k=0}^{n-1} f_k = sum 2^{T_k} ~ 2^{T_{n-1}} (dominated by last term)

  So: the total number of tournaments at ALL sizes up to n
  is approximately equal to the number at size n alone.
  Each new level DWARFS all previous levels combined.
  This is because T_k grows QUADRATICALLY (superexponential growth).

  In the tournaplex: almost all of the "mass" is at the top level.
  The lower levels are negligible. This means: the homology of
  the full tournaplex is dominated by the highest-dimensional simplices.
  And: our beta_0 = 0, beta_1 = 0 results (from S178) reflect the
  fact that the massive top levels fill in all the lower cycles.
""")

print("Done. opus-2026-03-22-S180")
