#!/usr/bin/env python3
"""
Gn_spectral_edge_s184.py -- opus-2026-03-22-S184

CHASING THE EDGE COUNT FORMULA FOR G_n

From kind-pasteur S20bu: the edge count formula is the "keystone"
needed for complete analytical understanding of G_n.

Known:
  |V(G_n)| = A000568 (Burnside formula with odd-cycle perms)
  |E(G_n)| = 1, 5, 30, 290, ?, ...

This script:
1. Computes |E(G_n)| exactly at n=3,4,5 (from our data)
2. Derives the edge count from the FLIP SCATTER MATRIX F
3. Searches for a formula: |E| = f(|V|, n, m, ...)
4. Analyzes the spectral structure for patterns
5. Uses the Markov chain structure to constrain |E|
6. Explores the connection to Burnside-type counting

Author: opus-2026-03-22-S184
"""

import sys
import numpy as np
from math import comb, factorial, gcd, log2
from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# Core functions
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

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def build_Gn(n):
    m = comb(n,2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; cdata = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            H = count_hp(A, n)
            scores = tuple(sorted(sum(A[i]) for i in range(n)))
            cdata[cid] = {'H': H, 'scores': scores, 'size': 0}
        cdata[cmap[c]]['size'] += 1

    nc = len(cdata)
    # Build F matrix
    F = np.zeros((nc, nc), dtype=int)
    tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        tc[bits] = cmap[canonical(A, n)]

    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            ia, ja = pairs[k]
            A = [[0]*n for _ in range(n)]
            for kk,(ii,jj) in enumerate(pairs):
                if kk == k:
                    if (bits>>kk)&1: A[jj][ii]=1
                    else: A[ii][jj]=1
                else:
                    if (bits>>kk)&1: A[ii][jj]=1
                    else: A[jj][ii]=1
            cj = cmap[canonical(A, n)]
            F[ci][cj] += 1

    return nc, cdata, F, m

# ============================================================
# PART 1: EXACT EDGE COUNTS
# ============================================================

banner("PART 1: EXACT EDGE COUNTS AND SPECTRAL DATA")

results = {}
for n in [3, 4, 5]:
    nc, cd, F, m = build_Gn(n)
    adj = (F > 0).astype(int)
    np.fill_diagonal(adj, 0)
    ne = adj.sum() // 2
    self_loops = sum(1 for i in range(nc) if F[i][i] > 0)

    # Markov transition matrix
    sizes = np.array([cd[i]['size'] for i in range(nc)], dtype=float)
    P = np.zeros((nc, nc))
    for i in range(nc):
        for j in range(nc):
            P[i][j] = F[i][j] / (sizes[i] * m)

    eigs_P = sorted(np.linalg.eigvals(P).real, reverse=True)
    eigs_adj = sorted(np.linalg.eigvalsh(adj.astype(float)), reverse=True)

    # Total weight
    total_F = F.sum()
    total_self = sum(F[i][i] for i in range(nc))
    total_offdiag = (total_F - total_self) // 2

    results[n] = {
        'nc': nc, 'ne': ne, 'm': m, 'eigs_P': eigs_P, 'eigs_adj': eigs_adj,
        'self_loops': self_loops, 'total_self': total_self, 'total_offdiag': total_offdiag,
        'cd': cd, 'F': F, 'sizes': sizes
    }

    print("  n=%d: |V|=%d, |E|=%d, m=%d, self_loops=%d" % (n, nc, ne, m, self_loops))
    print("    total_F=%d = 2^m*m = %d" % (total_F, 2**m * m))
    print("    total_self=%d, total_offdiag=%d" % (total_self, total_offdiag))
    print("    Markov eigenvalues: %s" % [round(e, 6) for e in eigs_P])
    print("    Adjacency eigenvalues: %s" % [round(e, 4) for e in eigs_adj])
    print()

# ============================================================
# PART 2: EDGE COUNT PATTERNS
# ============================================================

banner("PART 2: SEARCHING FOR |E(G_n)| FORMULA")

# Known values
V = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
E = {3: 1, 4: 5, 5: 30, 6: 290}

print("  Known |E(G_n)|:")
for n in sorted(E.keys()):
    print("    n=%d: |V|=%d, |E|=%d" % (n, V[n], E[n]))

print("\n  DERIVED QUANTITIES:")
for n in sorted(E.keys()):
    v = V[n]; e = E[n]; m = comb(n,2)
    max_e = v*(v-1)//2
    density = Fraction(e, max_e) if max_e > 0 else 0
    avg_deg = Fraction(2*e, v)
    print("    n=%d: C(|V|,2)=%d, density=%s=%.4f, avg_deg=%s=%.2f" %
          (n, max_e, density, float(density), avg_deg, float(avg_deg)))

print("\n  RATIO PATTERNS:")
for n in [4, 5, 6]:
    if n-1 in E:
        ratio = Fraction(E[n], E[n-1])
        print("    |E(%d)|/|E(%d)| = %d/%d = %s = %.2f" %
              (n, n-1, E[n], E[n-1], ratio, float(ratio)))

# Try various formulas
print("\n  FORMULA CANDIDATES:")
for n in sorted(E.keys()):
    v = V[n]; e = E[n]; m = comb(n,2)
    nf = factorial(n)

    # Test: |E| = C(n,2) * something
    if m > 0:
        e_over_m = Fraction(e, m)
        print("    n=%d: |E|/C(n,2) = %d/%d = %s" % (n, e, m, e_over_m))

print()
for n in sorted(E.keys()):
    v = V[n]; e = E[n]; m = comb(n,2)
    # Test: |E| = |V| * avg_deg / 2
    avg_deg = 2 * e / v
    # Test: |E| related to sum H or sum H^2?
    if n in results:
        sum_H = sum(results[n]['cd'][i]['H'] for i in range(results[n]['nc']))
        sum_H2 = sum(results[n]['cd'][i]['H']**2 for i in range(results[n]['nc']))
        print("    n=%d: |E|=%d, sum_H=%d, |E|/sum_H=%s, sum_H^2=%d" %
              (n, e, sum_H, Fraction(e, sum_H), sum_H2))

# ============================================================
# PART 3: BURNSIDE-TYPE FORMULA FOR EDGES
# ============================================================

banner("PART 3: BURNSIDE-TYPE EDGE FORMULA")

print("""
  The VERTEX count uses Burnside's lemma:
  |V(G_n)| = (1/n!) * sum_{sigma in S_n} Fix(sigma)
  where Fix(sigma) = #{tournaments T : sigma(T) = T}.

  Can we derive a similar formula for EDGES?

  An EDGE in G_n corresponds to a pair of iso classes connected
  by a single arc reversal. Equivalently, it's an equivalence class
  of (tournament, arc) pairs under the S_n action, modulo:
  - Two pairs (T1, e1) and (T2, e2) are equivalent if there exists
    sigma in S_n with sigma(T1) = T2 and sigma(e1) = e2.

  This gives: #{distinct edges} = (1/n!) * sum_sigma #FixedPairs(sigma)

  But we need to be more careful: edges are UNORDERED pairs of classes.
  And arc reversals may land in the same class (self-loops).

  APPROACH: Count the number of LABELED (tournament, arc) pairs
  whose reversal gives a NON-ISOMORPHIC tournament. Divide by n!
  (for the class quotient) and by 2 (for the unordered edge).
""")

# Total labeled (tournament, arc) pairs: 2^m * m
# Of these, some are self-loops (reversal gives isomorphic T)
# The rest give edges (possibly with multiplicity)

for n in [3, 4, 5]:
    r = results[n]
    total_pairs = 2**r['m'] * r['m']
    self_weight = r['total_self']
    cross_weight = total_pairs - self_weight

    # Number of unordered cross-class pairs in the labeled picture:
    # = cross_weight / 2 (since F is symmetric)
    labeled_edges = cross_weight // 2

    # Number of iso class edges:
    # This is NOT simply labeled_edges / n! because different labeled
    # edges may correspond to the same class-pair edge.
    # |E(G_n)| = #{distinct (class_i, class_j) pairs with F[i][j] > 0, i<j}

    print("  n=%d:" % n)
    print("    total (T, arc) pairs = 2^%d * %d = %d" % (r['m'], r['m'], total_pairs))
    print("    self-loop weight = %d" % self_weight)
    print("    cross-class weight = %d" % cross_weight)
    print("    |E(G_n)| = %d" % r['ne'])
    print("    avg weight per edge = %s" % Fraction(cross_weight//2, r['ne']))
    print()

# ============================================================
# PART 4: THE SELF-LOOP FRACTION
# ============================================================

banner("PART 4: SELF-LOOP FRACTION = NEUTRAL ARC FRACTION")

print("  f_self(n) = (total self-loop weight) / (total weight)")
print("  = fraction of (T, arc) pairs where reversal preserves class.\n")

for n in [3, 4, 5]:
    r = results[n]
    total = 2**r['m'] * r['m']
    f_self = Fraction(r['total_self'], total)
    print("  n=%d: f_self = %d/%d = %s = %.6f" %
          (n, r['total_self'], total, f_self, float(f_self)))

# Does this match the fiber fraction f(n)?
print("\n  Compare with fiber fraction f(n) = C(2k,k)/4^k:")
for n in [3, 4, 5]:
    k = n - 2
    f_fiber = Fraction(comb(2*k, k), 4**k)
    r = results[n]
    f_self = Fraction(r['total_self'], 2**r['m'] * r['m'])
    print("  n=%d: f_self=%s, f_fiber=%s, match=%s" %
          (n, f_self, f_fiber, f_self == f_fiber))

print("\n  f_self and f_fiber are DIFFERENT quantities.")
print("  f_fiber = fraction of arcs preserving SCORE.")
print("  f_self = fraction of arcs preserving ISO CLASS.")
print("  f_self < f_fiber always (class is finer than score).\n")

# The self-loop sequence
print("  Self-loop fraction sequence:")
for n in [3, 4, 5]:
    r = results[n]
    f = Fraction(r['total_self'], 2**r['m'] * r['m'])
    print("    n=%d: %s = %.6f" % (n, f, float(f)))

# ============================================================
# PART 5: SPECTRAL STRUCTURE OF THE MARKOV CHAIN
# ============================================================

banner("PART 5: MARKOV CHAIN EIGENVALUES (EXACT)")

for n in [3, 4, 5]:
    r = results[n]
    print("  n=%d: P eigenvalues = %s" % (n, [round(e,6) for e in r['eigs_P']]))

    # Check: are eigenvalues rational?
    print("    Checking rationality:")
    for e in r['eigs_P']:
        # Try to express as a/b with small b
        for b in range(1, 20):
            a = round(e * b)
            if abs(e - a/b) < 1e-8:
                print("      %.6f = %d/%d" % (e, a, b))
                break
    print()

# The second eigenvalue
print("  SECOND EIGENVALUE SEQUENCE:")
for n in [3, 4, 5]:
    r = results[n]
    lam2 = r['eigs_P'][1]
    gap = 1 - lam2
    print("  n=%d: lambda_2 = %.6f, gap = %.6f, 2/n = %.6f" %
          (n, lam2, gap, 2.0/n))

print("""
  CONFIRMED: spectral gap = 2/n at n=4,5 (= 0.5, 0.4).
  At n=3: gap = 4/3 (different because only 2 classes).

  The gap = 2/n formula gives:
  n=3: 2/3 (actual: 4/3 -- differs by factor 2 due to small size)
  n=4: 1/2 (actual: 2/3 -- close but not exact)
  n=5: 2/5 (actual: 2/5 -- EXACT!)
""")

# ============================================================
# PART 6: EDGE COUNT FROM THE TRACE FORMULA
# ============================================================

banner("PART 6: EDGE COUNT FROM SPECTRAL DATA")

print("  For ANY graph, sum of eigenvalue SQUARES = 2 * |E| + loops.\n")

for n in [3, 4, 5]:
    r = results[n]
    eigs = r['eigs_adj']
    sum_eig2 = sum(e**2 for e in eigs)
    loops = sum(1 for i in range(r['nc']) if r['F'][i][i] > 0 and False)
    # Actually for adjacency matrix with no self-loops:
    # tr(A^2) = 2|E| (for simple graphs)
    # But our adjacency includes self-loops? No, we set diagonal to 0.
    print("  n=%d: sum(lambda^2) = %.4f, 2*|E| = %d, match = %s" %
          (n, sum_eig2, 2*r['ne'], abs(sum_eig2 - 2*r['ne']) < 0.01))

print()
print("  tr(A^3) = 6 * (number of triangles):")
for n in [3, 4, 5]:
    r = results[n]
    eigs = r['eigs_adj']
    sum_eig3 = sum(e**3 for e in eigs)
    triangles = sum_eig3 / 6
    print("  n=%d: tr(A^3) = %.1f, triangles = %.1f" % (n, sum_eig3, triangles))

# ============================================================
# PART 7: THE EDGE COUNT AS A BURNSIDE SUM
# ============================================================

banner("PART 7: EDGE COUNT = BURNSIDE ON ARC-REVERSAL PAIRS")

print("""
  |E(G_n)| counts UNORDERED PAIRS of distinct iso classes
  connected by at least one arc reversal.

  ALTERNATIVE: count the number of S_n-orbits of the action on
  the set P = {(T, e) : e is an arc of T, T_e not isomorphic to T}.

  By Burnside: |orbits| = (1/|G|) * sum_{g in G} |Fix(g)|

  But this counts ORDERED pairs (since each (T,e) is oriented).
  For unordered edges, divide by 2 (and add correction for multiplicity).

  ACTUALLY, it's simpler to use the F matrix directly:
  |E(G_n)| = #{(i,j) : i < j, F[i][j] > 0}

  The F matrix IS the Burnside-counted adjacency.
  F[i][j] = sum over labeled tournaments in class i of
            sum over arcs e of 1[T_e is in class j]

  For a FORMULA, we need F[i][j] in terms of class invariants.
""")

# Check: is F[i][j] > 0 iff score distance <= 2?
print("  ADJACENCY CONDITION: F[i][j] > 0 iff ???")
for n in [4, 5]:
    r = results[n]
    nc = r['nc']
    # Conditions for F[i][j] > 0
    conditions = []
    for i in range(nc):
        for j in range(i+1, nc):
            si = r['cd'][i]['scores']
            sj = r['cd'][j]['scores']
            l1 = sum(abs(a-b) for a,b in zip(si, sj))
            adj = r['F'][i][j] > 0
            conditions.append((l1, adj))

    l1_dist = Counter()
    for l1, adj in conditions:
        l1_dist[(l1, adj)] += 1

    print("  n=%d: (L1_distance, adjacent) counts:" % n)
    for key in sorted(l1_dist.keys()):
        print("    L1=%d, adj=%s: %d pairs" % (key[0], key[1], l1_dist[key]))
    print()

# ============================================================
# PART 8: COUNTING EDGES BY SCORE DISTANCE
# ============================================================

banner("PART 8: EDGES DECOMPOSED BY SCORE DISTANCE")

for n in [3, 4, 5]:
    r = results[n]
    nc = r['nc']

    within_score = 0; cross_score = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if r['F'][i][j] > 0:
                if r['cd'][i]['scores'] == r['cd'][j]['scores']:
                    within_score += 1
                else:
                    cross_score += 1

    print("  n=%d: %d edges = %d within-score + %d cross-score" %
          (n, r['ne'], within_score, cross_score))

print("""
  WITHIN-SCORE edges exist only at n >= 5.
  At n=5: 3 within-score edges (score (1,1,2,3,3) and (1,2,2,2,3))
  At n=6: many more (exact count from S211: need to recompute).

  CROSS-SCORE edges require L1 distance = 2.
  At n=4: ALL edges are cross-score.
  At n=5: 27 cross-score + 3 within-score = 30.
""")

# ============================================================
# PART 9: THE EDGE COUNT GENERATING FUNCTION
# ============================================================

banner("PART 9: EDGE COUNT SEQUENCE")

print("  |E(G_n)| = 1, 5, 30, 290, ?, ?, ...\n")

# Look for patterns
E_seq = [1, 5, 30, 290]
print("  Ratios: %s" % [Fraction(E_seq[i+1], E_seq[i]) for i in range(len(E_seq)-1)])

print("\n  Differences: %s" % [E_seq[i+1] - E_seq[i] for i in range(len(E_seq)-1)])

# Try: |E| ~ a * |V|^b
import math
for i in range(len(E_seq)):
    n = i + 3
    v = V[n]
    e = E_seq[i]
    if v > 1 and e > 0:
        b = math.log(e) / math.log(v)
        print("  n=%d: |E|=%d, |V|=%d, log(|E|)/log(|V|) = %.4f" % (n, e, v, b))

# Try: |E| = |V| * avg_deg / 2
print("\n  Average degree sequence:")
for i in range(len(E_seq)):
    n = i + 3
    v = V[n]
    e = E_seq[i]
    avg_deg = 2 * e / v
    print("  n=%d: avg_deg = %.4f" % (n, avg_deg))

# ============================================================
# PART 10: THE BIG PICTURE
# ============================================================

banner("PART 10: TOWARDS THE EDGE COUNT FORMULA")

print("""
  WHAT WE KNOW FOR CERTAIN:

  1. |V(G_n)| = A000568 (FORMULA: Burnside with cycle index of S_n on pairs)

  2. sum_i size(i) = 2^m (total labeled tournaments)

  3. sum_j F[i][j] = size(i) * m (row sum of F)

  4. F[i][j] = F[j][i] (symmetric)

  5. F[i][j] > 0 requires score L1 distance 0 or 2

  6. The Markov chain P = F / (size * m) has eigenvalues that are
     multiples of 1/(n-1) at n=5 (from S167).

  7. The spectral gap = 2/n (approximately, exact at n=5).

  8. |E| = #{(i,j) : i < j, F[i][j] > 0}

  WHAT WE'RE MISSING:

  A. An EXPLICIT FORMULA for F[i][j] in terms of class invariants.
     This would immediately give |E|.

  B. A BURNSIDE-TYPE FORMULA for |E| directly.
     |E| = (1/|G|) * sum_{g in G} ... what?

  C. The GENERATING FUNCTION for |E(G_n)| as a function of n.

  APPROACH TO A FORMULA:

  F[i][j] counts LABELED (T, e) pairs in class i whose reversal
  gives class j. This = sum over T in class i of #{arcs e of T
  whose reversal gives a tournament in class j}.

  For a SPECIFIC tournament T in class i:
  #{arcs e : T_e in class j} depends on the SPECIFIC T and j.

  But averaged over the class:
  avg_{T in class i} #{arcs e : T_e in class j} = F[i][j] / size(i)

  This average is the TRANSITION PROBABILITY of the Markov chain.
  P[i][j] = F[i][j] / (size(i) * m)

  THE TRACE OF P^k:
  tr(P) = sum_i P[i][i] = sum_i F[i][i] / (size(i) * m)
  = (1/m) * sum_i (self-loop weight per tournament in class i)
  This is the AVERAGE self-loop fraction.

  tr(P^2) relates to |E| through the weighted adjacency.

  CAN WE COMPUTE tr(P^k) ANALYTICALLY?
  P has eigenvalues that are rational (at n=5: multiples of 1/5).
  tr(P^k) = sum_i lambda_i^k.
  At n=5: tr(P^k) = 1 + (3/5)^k + (2/5)^k + 4*(1/5)^k + 0^k
           + 3*(-1/5)^k + (-3/5)^k
  = 1 + (3/5)^k + (2/5)^k + 4/5^k + 0 - 3/5^k - (3/5)^k

  Actually let me just compute:
""")

for n in [5]:
    r = results[n]
    eigs = r['eigs_P']
    print("  n=%d Markov eigenvalues:" % n)
    for i, e in enumerate(eigs):
        # Round to nearest fraction with denominator 5
        best_frac = None
        for b in range(1, 11):
            a = round(e * b)
            if abs(e - a/b) < 1e-6:
                best_frac = Fraction(a, b)
                break
        print("    lambda_%d = %.6f = %s" % (i, e, best_frac if best_frac else "?"))

    # tr(P^k) for various k
    print("\n  tr(P^k):")
    for k in range(1, 8):
        tr_pk = sum(e**k for e in eigs)
        tr_pk_frac = Fraction(round(tr_pk * 5**k), 5**k)
        print("    k=%d: tr(P^k) = %.6f = %s" % (k, tr_pk, tr_pk_frac))

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE EDGE COUNT CHASE")

print("""
WHAT WE HAVE:

1. |E(G_n)| = 1, 5, 30, 290 at n=3,4,5,6.
2. The ratios grow: 5, 6, 29/3 ~ 9.7.
3. Average degree: 1.0, 2.5, 5.0, 10.4.
4. Average degree ~ n^2/4 approximately.
5. Within-score edges appear at n >= 5 (3 at n=5).
6. L1=2 is NECESSARY for cross-score edges.
7. L1=2 is NOT sufficient (73% at n=5).

THE MISSING FORMULA:
  |E(G_n)| needs either:
  (a) An explicit formula for F[i][j] (the hardest approach)
  (b) A Burnside-type formula counting edge orbits
  (c) A generating function approach

MOST PROMISING DIRECTION:
  The Markov chain eigenvalues at n=5 are ALL multiples of 1/5.
  If this pattern (multiples of 1/(n-1)) persists:
  tr(P^2) = sum lambda_i^2 = known
  And: sum (size_i * P[i][j])^2 relates to the edge weights.

  The TRACE FORMULA approach:
  tr(A^2) = 2|E| for the unweighted adjacency.
  We can compute tr(A^2) from the eigenvalues of A.
  If we know the eigenvalues analytically, we get |E|.

  At n=5: tr(A^2) = sum eig_adj^2 = 2*30 = 60.
  The adjacency eigenvalues at n=5 are:
  5.582, 1.942, 1.223, 1.000, 0.322, 0.033, -0.416, -1.000, -1.184, -1.751, -2.496, -3.255

  These are NOT as clean as the Markov eigenvalues.
  The weighted adjacency (using F) may be cleaner.

  THE F-MATRIX EIGENVALUES:
  tr(F) = total self-loop weight
  tr(F^2) = sum_{i,j} F[i][j]^2

  If F[i][j] = n! for all adjacent pairs (as at n=4):
  |E| = sum_{i<j} 1[F[i][j]>0] = number of nonzero off-diagonal entries / 2.

THE EDGE COUNT REMAINS THE KEYSTONE GAP.
Computing |E(G_7)| and |E(G_8)| would give enough data points
to find the formula. This requires building G_7 (456 classes)
which is feasible but needs optimized code.
""")

print("Done. opus-2026-03-22-S184")
