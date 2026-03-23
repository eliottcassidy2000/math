#!/usr/bin/env python3
"""
tiling_class_analysis_s204.py -- opus-2026-03-23-S204

THE META-GRAPH FROM THE TILING PERSPECTIVE

Focus: isomorphism classes and their tilings. NOT blue/black.

Each iso class C is a "cloud" of H(C)/|Aut(C)| tilings.
The tiling flip (bit complement) creates a PERFECT MATCHING
on the 2^m tilings. This matching pairs each tiling with exactly
one partner. The class structure emerges from how the matching
respects or crosses class boundaries.

KEY QUESTIONS:
1. How do tiling counts distribute across classes?
2. What is the "tiling profile" of G_n?
3. How does the flip matching organize the classes?
4. What sequences emerge from tiling counts?
5. How do tiling counts relate to other class invariants?

Author: opus-2026-03-23-S204
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

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

def aut_size(A, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: count += 1
    return count

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def build_tiling_data(n):
    """Build complete tiling-level data."""
    m = comb(n-1, 2)  # tiling bits
    pairs = []
    for i in range(n):
        for j in range(i+2, n):
            pairs.append((i,j))

    # Build class data from tilings
    cmap = {}; cdata = {}; creps = {}
    tiling_to_class = {}

    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for i in range(n-1): A[i][i+1] = 1
        for idx, (i,j) in enumerate(pairs):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            H = count_hp(A, n)
            aut = aut_size(A, n) if n <= 5 else 1
            cdata[cid] = {'H': H, 'aut': aut, 'scores': score_seq(A, n),
                          'tilings': [], 'tiling_count': 0}
            creps[cid] = [row[:] for row in A]
        cid = cmap[c]
        cdata[cid]['tilings'].append(bits)
        cdata[cid]['tiling_count'] += 1
        tiling_to_class[bits] = cid

    nc = len(cdata)
    return nc, cdata, creps, tiling_to_class, m

# ============================================================
# COMPUTE FOR n=3,4,5
# ============================================================

for n in [3, 4, 5]:
    banner("n = %d: TILING-CLASS STRUCTURE" % n)

    nc, cd, creps, tc, m = build_tiling_data(n)

    print("  Total tilings: 2^%d = %d" % (m, 2**m))
    print("  Iso classes: %d\n" % nc)

    # ===== 1. TILING COUNT DISTRIBUTION =====
    tiling_counts = sorted([cd[c]['tiling_count'] for c in range(nc)])
    print("  TILING COUNTS PER CLASS:")
    print("  Sorted: %s" % tiling_counts)
    print("  Sum: %d = 2^%d = %d" % (sum(tiling_counts), m, 2**m))
    print("  Min: %d, Max: %d, Mean: %.2f, Median: %d" %
          (min(tiling_counts), max(tiling_counts),
           sum(tiling_counts)/len(tiling_counts),
           tiling_counts[len(tiling_counts)//2]))

    # ===== 2. TILING COUNT = H/|Aut| VERIFICATION =====
    print("\n  VERIFICATION: tiling_count = H/|Aut|")
    all_match = True
    for cid in range(nc):
        d = cd[cid]
        predicted = d['H'] // d['aut']
        actual = d['tiling_count']
        if predicted != actual:
            print("    MISMATCH: class %d H=%d |Aut|=%d predicted=%d actual=%d" %
                  (cid, d['H'], d['aut'], predicted, actual))
            all_match = False
    print("  All match: %s" % all_match)

    # ===== 3. FLIP MATCHING STRUCTURE =====
    print("\n  FLIP MATCHING:")
    self_flip = 0  # tilings where flip partner is in same class
    cross_flip = 0

    flip_map = defaultdict(Counter)  # class -> {target_class: count}
    for bits in range(2**m):
        f = bits ^ ((1 << m) - 1)
        ci = tc[bits]; cj = tc[f]
        if bits < f:  # avoid double counting
            flip_map[ci][cj] += 1
            if ci != cj:
                flip_map[cj][ci] += 1
                cross_flip += 1
            else:
                self_flip += 1

    print("  Self-flip pairs (same class): %d" % self_flip)
    print("  Cross-flip pairs (different class): %d" % cross_flip)
    print("  Total pairs: %d = 2^%d / 2 = %d" %
          (self_flip + cross_flip, m, 2**(m-1)))

    # ===== 4. CLASS TABLE =====
    print("\n  CLASS TABLE:")
    print("  %-5s %-5s %-6s %-8s %-8s %-6s %-12s" %
          ("cls", "H", "|Aut|", "#tiling", "self_fl", "cross", "flip_targets"))

    for cid in range(nc):
        d = cd[cid]
        sf = sum(1 for bits in d['tilings']
                 if tc.get(bits ^ ((1<<m)-1), -1) == cid and bits < bits ^ ((1<<m)-1))
        targets = {}
        for bits in d['tilings']:
            f = bits ^ ((1 << m) - 1)
            cj = tc[f]
            if cj != cid:
                targets[cj] = targets.get(cj, 0) + 1

        total_cross = sum(targets.values())
        target_str = str(dict(sorted(targets.items())))[:40]

        print("  %-5d %-5d %-6d %-8d %-8d %-6d %s" %
              (cid, d['H'], d['aut'], d['tiling_count'], sf, total_cross, target_str))

    # ===== 5. TILING COUNT SEQUENCES =====
    print("\n  TILING COUNT SORTED BY H:")
    by_H = sorted(range(nc), key=lambda c: cd[c]['H'])
    h_tiling_seq = [(cd[c]['H'], cd[c]['tiling_count']) for c in by_H]
    print("  %s" % h_tiling_seq)

    # The sequence of tiling counts (sorted by H)
    just_tilings = [cd[c]['tiling_count'] for c in by_H]
    print("  Tiling count sequence: %s" % just_tilings)

    # ===== 6. DEGREE FROM TILING PERSPECTIVE =====
    print("\n  DEGREE = #{distinct target classes via flip}:")
    for cid in range(nc):
        targets = set()
        for bits in cd[cid]['tilings']:
            f = bits ^ ((1 << m) - 1)
            cj = tc[f]
            if cj != cid:
                targets.add(cj)
        deg = len(targets)
        print("    class %d (H=%d, %d tilings): flip-degree = %d, targets = %s" %
              (cid, cd[cid]['H'], cd[cid]['tiling_count'], deg,
               sorted([(t, cd[t]['H']) for t in targets], key=lambda x: x[1])))

    # ===== 7. THE "TILING WEIGHT" ON EACH EDGE =====
    print("\n  TILING WEIGHT ON EDGES:")
    edge_weights = defaultdict(int)
    for bits in range(2**m):
        f = bits ^ ((1 << m) - 1)
        if bits >= f: continue
        ci = tc[bits]; cj = tc[f]
        if ci != cj:
            edge = (min(ci,cj), max(ci,cj))
            edge_weights[edge] += 1

    for edge in sorted(edge_weights.keys()):
        ci, cj = edge
        w = edge_weights[edge]
        ti = cd[ci]['tiling_count']; tj = cd[cj]['tiling_count']
        print("    (%d,%d) H=(%d,%d) tilings=(%d,%d) weight=%d w/min=%s" %
              (ci, cj, cd[ci]['H'], cd[cj]['H'], ti, tj, w,
               Fraction(w, min(ti,tj))))

    # ===== 8. TILING DEGREE DISTRIBUTION =====
    print("\n  TILING-DEGREE vs H:")
    for cid in by_H:
        targets = set()
        for bits in cd[cid]['tilings']:
            f = bits ^ ((1 << m) - 1)
            cj = tc[f]
            if cj != cid: targets.add(cj)
        sf = cd[cid]['tiling_count'] - len([bits for bits in cd[cid]['tilings']
                                            if tc.get(bits ^ ((1<<m)-1)) != cid])
        # Actually count self-flip more carefully
        self_count = sum(1 for bits in cd[cid]['tilings']
                        if tc.get(bits ^ ((1<<m)-1)) == cid)
        cross_count = cd[cid]['tiling_count'] - self_count

        print("    H=%2d: %d tilings, %d self-flip, %d cross-flip -> %d targets" %
              (cd[cid]['H'], cd[cid]['tiling_count'], self_count, cross_count, len(targets)))

    print()

# ============================================================
# CROSS-n PATTERNS
# ============================================================

banner("CROSS-n PATTERNS IN TILING COUNTS")

print("""
  TILING COUNT SEQUENCES (sorted by H):

  n=3: [1, 1]                                    (H = 1, 3)
  n=4: [1, 1, 1, 5]                              (H = 1, 3, 3, 5)
  n=5: [1, 1, 1, 1, 5, 5, 9, 9, 11, 13, 5, 3]   (H sorted)

  PATTERNS IN TILING COUNTS:
  1. The TRANSITIVE class always has 1 tiling.
  2. Classes with |Aut|=3 have H/3 tilings (always integer since 3|H).
  3. Classes with |Aut|=5 have H/5 tilings (regular at n=5: 15/5=3).
  4. Generic classes (|Aut|=1) have #tilings = H.

  THE TILING COUNT SEQUENCE IS:
  - Sorted by H: a non-decreasing sequence (mostly)
  - But JUMPS occur when |Aut| changes

  THE SELF-FLIP FRACTION:
  = fraction of tilings whose flip partner is in the same class
  = measure of "class stability" under flip

  Self-flip fraction per class:
  - Classes with 1 tiling: self-flip iff flip(tiling) is in same class
    This requires the tiling AND its complement to be in the same class
    = the class is "blueself" or "blackself"
  - Classes with many tilings: self-flip is common (many tilings, some pair up)
""")

# ============================================================
# THE TOTAL TILING WEIGHT = 2^{m-1}
# ============================================================

banner("THE TOTAL TILING WEIGHT")

print("  Total flip pairs = 2^m / 2 = 2^{m-1}.\n")

for n in [3, 4, 5]:
    m = comb(n-1, 2)
    total = 2**(m-1)

    nc, cd, creps, tc, m_val = build_tiling_data(n)
    self_total = 0; cross_total = 0
    for bits in range(2**m):
        f = bits ^ ((1 << m) - 1)
        if bits < f:
            ci = tc[bits]; cj = tc[f]
            if ci == cj: self_total += 1
            else: cross_total += 1

    print("  n=%d: 2^{m-1} = %d = %d self + %d cross" %
          (n, total, self_total, cross_total))
    print("    Self fraction: %s" % Fraction(self_total, total))
    print("    This = 1 - (2*|E_tiling_flip|) / 2^m")
    print()

# ============================================================
# THE TILING COUNT DISTRIBUTION
# ============================================================

banner("THE TILING COUNT DISTRIBUTION")

for n in [3, 4, 5]:
    nc, cd, creps, tc, m = build_tiling_data(n)
    counts = sorted([cd[c]['tiling_count'] for c in range(nc)])
    count_dist = Counter(counts)

    print("  n=%d: tiling count distribution:" % n)
    for tc_val in sorted(count_dist.keys()):
        freq = count_dist[tc_val]
        print("    %d tilings: %d classes (%.1f%%)" %
              (tc_val, freq, 100*freq/nc))

    # What determines the tiling count?
    # tiling_count = H / |Aut|
    # So: classes with the same H/|Aut| ratio have the same tiling count
    print("    Distinct tiling counts: %d" % len(count_dist))

    # The SORTED tiling count sequence
    print("    Sorted: %s" % counts)

    # Is the tiling count sequence related to a known sequence?
    print("    Sum: %d = 2^%d" % (sum(counts), m))
    print()

# ============================================================
# THE FLIP GRAPH ON TILINGS
# ============================================================

banner("THE FLIP GRAPH ON TILINGS")

print("""
  The FLIP GRAPH on tilings is a PERFECT MATCHING on 2^m nodes.
  Every tiling pairs with exactly one other.

  The INDUCED GRAPH on iso classes: weight = #{tiling pairs connecting them}.
  This IS the G_n graph with weights from the tiling model.

  KEY STRUCTURAL FACT:
  Each class C has tiling_count(C) = H(C)/|Aut(C)| tilings.
  Of these, some flip to the same class (self-flip), others cross.
  self_flip(C) = #{tilings in C whose flip is also in C}
  cross_flip(C) = tiling_count(C) - self_flip(C)

  DEGREE: degree(C) = #{distinct classes reached by cross-flips from C}
  WEIGHT: weight(C,D) = #{tilings in C whose flip is in D}

  The weight is ALWAYS <= min(tiling_count(C), tiling_count(D)).
  From S202: for small-tiling classes (1-2 tilings), weight = min exactly.
""")

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE TILING PERSPECTIVE ON G_n")

print("""
  THE META-GRAPH FROM THE TILING PERSPECTIVE:

  LEVEL 1: THE TILING SPACE
    2^m tilings, each a binary string of length m = C(n-1,2).
    A PERFECT MATCHING: every tiling paired with its bit-complement.
    The matching partitions 2^m tilings into 2^{m-1} pairs.

  LEVEL 2: THE ISO CLASS PARTITION
    The 2^m tilings partition into nc = A000568(n) iso classes.
    Class C has |C| = H(C)/|Aut(C)| tilings.
    Sum of all |C| = 2^m.

  LEVEL 3: THE FLIP MATCHING ON CLASSES
    Each tiling pair either:
    (a) Both in same class -> SELF-FLIP (contributes to class stability)
    (b) In different classes -> CROSS-FLIP (contributes to G_n edge)

    Self-flip count: sum_C floor(|C|_self/2)
    Cross-flip count: 2^{m-1} - self_flip

  LEVEL 4: THE WEIGHTED CLASS GRAPH
    Vertices = iso classes
    Edge (C,D) exists iff some tiling in C flips to D
    Weight = #{such tiling pairs}
    Total weight = cross-flip count

  THE TILING COUNT IS THE "MASS" OF EACH NODE:
    Heavy nodes (large H/|Aut|) have more tilings -> more flip partners
    -> more connections -> higher degree.
    Light nodes (H=1 or large |Aut|) have few tilings -> degree 1.

  THE PERFECT MATCHING CONSTRAINT:
    Every tiling has exactly one partner. This means:
    sum over all edges of weight(C,D) + sum_C self_flip(C) = 2^{m-1}
    This is a GLOBAL CONSERVATION LAW on the graph.

  THE TILING COUNT DETERMINES EVERYTHING:
    Given the tiling counts {|C| = H(C)/|Aut(C)|} and the flip matching,
    the entire G_n structure is determined. The tiling counts are
    themselves determined by H and |Aut|, which are class invariants.

  OPEN QUESTION:
    Can we predict the DEGREE of each class from its tiling count alone?
    From S195: deg(transitive) = C(n-1,2) = m.
    The transitive has 1 tiling and degree = m (the maximum!).
    So: small tiling count does NOT mean small degree.
    In fact: the TRANSITIVE (1 tiling) has the HIGHEST degree
    among classes with |Aut|=1 (since it has m non-neutral arcs,
    each leading to a distinct class).

  THE PARADOX:
    More tilings -> more total flip partners -> more WEIGHT.
    But NOT necessarily more DISTINCT targets.
    A class with 100 tilings might flip 50 pairs to class A
    and 50 to class B -> degree only 2, despite high weight.
    A class with 1 tiling flips that 1 tiling to some class -> degree 1.
    UNLESS the tiling is the transitive (|Aut|=1 forces all arcs distinct)
    then the ARC-REVERSAL degree (not flip-degree) = C(n-1,2).

  RESOLUTION: The FLIP DEGREE (from the tiling model) and the
  ARC-REVERSAL DEGREE (from G_n) are DIFFERENT quantities.
  - Flip degree of transitive = 1 (one flip partner, one target class)
  - Arc-reversal degree of transitive = C(n-1,2) (many arc reversals)
  The tiling-flip is a GLOBAL operation (complement ALL non-base arcs).
  Arc reversal is a LOCAL operation (flip ONE arc).
""")

print("Done. opus-2026-03-23-S204")
