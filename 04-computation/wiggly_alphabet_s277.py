#!/usr/bin/env python3
"""
wiggly_alphabet_s277.py -- opus-2026-03-24-S277

THE WIGGLY ALPHABET

Each wiggly class = a letter. m = C(n,2) letters.
A wiggly path = a word over this alphabet.
S_n permutes the letters (by relabeling vertices → relabeling cells).

ALTERNATIVE MERGINGS for the wiggly meta-graph:

1. STANDARD MERGING (complement = Z_2):
   Merge tiling with its complement. Same as current G_n/Z_2.
   The alphabet is UNCHANGED by complement — complement flips ALL bits
   but doesn't relabel cells.

2. ALPHABET MERGING (S_n on letters):
   Two wiggly paths are equivalent iff one is obtained from the other
   by permuting the letters. This gives "wiggly words up to S_n."
   This IS the existing iso class structure.

3. LETTER-CLASS MERGING:
   Group letters by their STAIRCASE RANGE.
   Range-1 letters (adjacent pairs): n-1 of them
   Range-2 letters: n-2 of them
   ...
   Range-(n-1) letter (apex): 1 of them
   Merge all letters at the same range into one "super-letter."
   This gives an alphabet of size n-1 (one per range level).

4. POSITION-FREE WORDS:
   A wiggly path of length k = a sequence of k letters.
   If we forget the ORDER of letters and just count:
   how many times each letter appears → a "wiggly multiset."
   The multiset determines which cells were flipped (with multiplicity).

Author: opus-2026-03-24-S277
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
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

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

banner("THE WIGGLY ALPHABET AT n=%d" % n)

print("  Alphabet size: m = C(%d,2) = %d letters\n" % (n, m))

# Label letters by their staircase range
for k, (i,j) in enumerate(pairs):
    r = j - i
    print("  Letter %d: pair (%d,%d), range %d" % (k, i, j, r))

# Group by range
range_groups = defaultdict(list)
for k, (i,j) in enumerate(pairs):
    range_groups[j-i].append(k)

print("\n  RANGE GROUPS:")
for r in sorted(range_groups.keys()):
    letters = range_groups[r]
    print("  Range %d: %d letters (%s)" % (r, len(letters), letters))

# Build iso classes
cmap = {}; sizes = {}; tc = {}; reps = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    sizes[cmap[c]] += 1; tc[bits] = cmap[c]
nc = len(cmap)

info = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                 'SC': cmap[comp_c] == cid, 'size': sizes[cid]}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

H_vec = [info[merged[mi][0]]['H'] for mi in range(nm)]

# ============================================================
banner("WIGGLY WORDS: PATHS IN THE META-GRAPH AS WORDS OVER THE ALPHABET")
# ============================================================

# A wiggly path of length k from tiling T:
# (T, letter_1, letter_2, ..., letter_k) → T' = T flipped at letters 1..k
# The final tiling T' = T ^ (1<<letter_1) ^ (1<<letter_2) ^ ...
# If any letter appears TWICE, it cancels (XOR).
# So effectively, a path of length k with distinct letters gives
# T' = T ^ (letters as bit mask).
# A path with repeated letters gives T' = T ^ (XOR of letters).

# The KEY: the set of reachable tilings from T by k wiggly flips
# = all tilings at Hamming distance ≤ k from T in Q_m.

# But the PATH matters for the meta-graph traversal: at each step,
# the intermediate tiling must be a tournament (always true in Q_m).

print("""
  WIGGLY WORDS:

  A word of length k = a sequence of k letters from the alphabet.
  Each word traces a path in Q_m (the m-cube of all tilings).

  PROPERTIES:
  - The same letter twice CANCELS (XOR): ...X...X... returns to before X.
  - A word of length k with all distinct letters = Hamming distance k in Q_m.
  - The ISO CLASS at each step determines the meta-graph path.

  WORD EQUIVALENCES:
  - S_n acts on the alphabet: σ relabels letters.
    Words w and σ(w) give ISOMORPHIC paths.
  - Complement (Z_2) flips all bits: equivalent to prepending "flip all."
    This doesn't change the meta-graph path at all (merging).

  The WIGGLY WORD LANGUAGE:
  L_n = {words w : w spells a path in the meta-graph from class C to class D}
  This language is REGULAR (determined by the finite meta-graph automaton).
""")

# ============================================================
banner("RANGE-REDUCED ALPHABET")
# ============================================================

# The user suggests: think of the alphabet in terms of RANGES.
# Instead of m individual letters, use n-1 "super-letters" (one per range).

# At n=5: 4 super-letters:
# Range 1: 4 letters (adjacent pairs)
# Range 2: 3 letters (skip-1)
# Range 3: 2 letters (skip-2)
# Range 4: 1 letter (apex)

# A range-reduced word = sequence of ranges, forgetting which specific
# letter within the range.

# The RANGE-REDUCED META-GRAPH: two merged classes are "range-X connected"
# if some tiling in one can reach the other by flipping a range-X letter.

# From S_n transitivity: ALL ranges see ALL edges (every range connects
# every pair of classes). So the range-reduced meta-graph is the same
# for every range — the full meta-graph.

# BUT: the MULTIPLICITIES differ by range.

print("  RANGE-REDUCED ALPHABET: %d super-letters (one per range)\n" % (n-1))

# Weight per range per edge
W_range = [np.zeros((nm, nm), dtype=int) for _ in range(n)]
for bits in range(2**m):
    mi = mid_map[tc[bits]]
    for k in range(m):
        r = pairs[k][1] - pairs[k][0]
        mj = mid_map[tc[bits ^ (1 << k)]]
        W_range[r][mi][mj] += 1

print("  EDGE WEIGHTS BY RANGE:")
print("  Edge                       ", end="")
for r in range(1, n):
    print("| r=%d     " % r, end="")
print("| total")
print("  " + "-"*70)

# Build adjacency
A_adj = np.zeros((nm, nm), dtype=int)
for mi in range(nm):
    for mj in range(nm):
        if mi != mj and sum(W_range[r][mi][mj] for r in range(1,n)) > 0:
            A_adj[mi][mj] = 1

for mi in range(nm):
    for mj in range(mi+1, nm):
        if not A_adj[mi][mj]: continue
        Hi = H_vec[mi]; Hj = H_vec[mj]
        print("  {%d(H=%d),%d(H=%d)}" % (mi, Hi, mj, Hj), end="")
        total = 0
        for r in range(1, n):
            w = W_range[r][mi][mj] + W_range[r][mj][mi]
            total += w
            print("| %-8d" % w, end="")
        print("| %d" % total)

# ============================================================
banner("SELF-LOOP WEIGHTS BY RANGE")
# ============================================================

print("  Node (H)              ", end="")
for r in range(1, n):
    print("| r=%d     " % r, end="")
print("| total")
print("  " + "-"*60)

for mi in range(nm):
    sl = sum(W_range[r][mi][mi] for r in range(1, n))
    if sl == 0: continue
    H = H_vec[mi]
    print("  %d (H=%d)               " % (mi, H), end="")
    for r in range(1, n):
        print("| %-8d" % W_range[r][mi][mi], end="")
    print("| %d" % sl)

# ============================================================
banner("THE RANGE RATIO: A UNIVERSAL PATTERN?")
# ============================================================

# At n=5, from S274: weights by range follow 4:3:2:1 ratio
# = number of letters at each range.
# Is this universal? Let me check the per-edge ratios.

print("  For each edge, the weight ratio by range:")
print("  (Normalized so range-1 = 1.000)\n")

for mi in range(nm):
    for mj in range(mi+1, nm):
        if not A_adj[mi][mj]: continue
        weights = [W_range[r][mi][mj] + W_range[r][mj][mi] for r in range(1, n)]
        if weights[0] == 0: continue
        ratios = [w / weights[0] for w in weights]
        expected = [(n-1-r+1) / (n-1) for r in range(n-1)]
        # Actually expected = (n-r) / (n-1) normalized to range-1
        expected_norm = [(n-1-r) / (n-1-0) for r in range(n-1)]
        print("  {%d,%d}: ratios=%s, expected=%s, match=%s" %
              (mi, mj, [round(r, 3) for r in ratios],
               [round(e, 3) for e in expected_norm],
               all(abs(ratios[i] - expected_norm[i]) < 0.01 for i in range(n-1))))

banner("THE WIGGLY WORD AUTOMATON")

print("""
  THE META-GRAPH AS A FINITE AUTOMATON:

  States: merged iso classes (V_merged nodes)
  Alphabet: wiggly classes {A, B, C, ...} (m letters)
  Transitions: δ(state, letter) = set of states reachable by flipping letter

  Wait — from a single merged class, flipping letter X doesn't give a
  SINGLE target class. Different labeled tilings in the source class
  may land in DIFFERENT target classes when flipping the same letter.

  So: δ(C, X) = {D : ∃ T ∈ C with T^X ∈ D}

  This is a NONDETERMINISTIC automaton (NFA).

  But from the local=global theorem: δ(C, X) = set of ALL neighbors of C.
  So for every letter X: δ(C, X) = {D : {C,D} is a meta-graph edge} ∪ {C if neutral}.

  This means: THE NFA IS THE SAME FOR EVERY LETTER.
  Every letter gives the same transition function!
  The wiggly word language is LETTER-INDEPENDENT.

  IMPLICATION: The meta-graph path depends on the SEQUENCE of classes
  visited, not on WHICH specific letter was flipped at each step.
  The alphabet labels carry NO information about the meta-graph dynamics.

  The information carried by the alphabet labels is WITHIN-CLASS:
  which specific tiling within the class is selected. This is the
  fiber structure (permutohedron at the transitive class, etc.).

  THE RANGE RATIO:
  Although every range sees every edge, the WEIGHT differs.
  Range-r letters contribute weight proportional to (n-r) = #{arcs at range r}.
  This is because S_n acts transitively on pairs, and there are (n-r)
  pairs at distance r in the staircase.

  So: the weight of edge {C,D} at range r = (n-r) × (base weight of {C,D}).
  The base weight depends only on the edge, not the range.
  This is a PERFECT FACTORIZATION: weight(edge, range) = w(edge) × (n-r).
""")

# Verify the perfect factorization
print("  VERIFYING: weight(edge, range) = w(edge) × (n-r)\n")

factorization_holds = True
for mi in range(nm):
    for mj in range(mi+1, nm):
        if not A_adj[mi][mj]: continue
        weights = [W_range[r][mi][mj] + W_range[r][mj][mi] for r in range(1, n)]
        # If factorization holds: w_r / (n-r) should be constant
        base_weights = [w / (n-r) for r, w in enumerate(weights, 1) if w > 0]
        if base_weights and max(base_weights) - min(base_weights) > 0.01:
            factorization_holds = False
            print("  FAIL at {%d,%d}: base_weights=%s" %
                  (mi, mj, [round(b, 2) for b in base_weights]))

if factorization_holds:
    print("  ✓ PERFECT FACTORIZATION VERIFIED for all edges at n=%d" % n)
    print("    weight(edge, range r) = base_weight(edge) × (n - r)")
    print("    where (n-r) = #{pairs at distance r in the staircase}")

# Compute base weights
print("\n  BASE WEIGHTS (range-independent edge weights):")
for mi in range(nm):
    for mj in range(mi+1, nm):
        if not A_adj[mi][mj]: continue
        w_r1 = W_range[1][mi][mj] + W_range[1][mj][mi]
        base = w_r1 / (n - 1)
        Hi = H_vec[mi]; Hj = H_vec[mj]
        print("    {%d(H=%d),%d(H=%d)}: base = %.0f" % (mi, Hi, mj, Hj, base))

print("\nDone. opus-2026-03-24-S277")
