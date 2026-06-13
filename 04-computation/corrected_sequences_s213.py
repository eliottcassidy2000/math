#!/usr/bin/env python3
"""
corrected_sequences_s213.py -- opus-2026-03-23-S213

THE CORRECTED META-GRAPH SEQUENCES

After MISTAKE-029 (formula E=(T-D)/2 is wrong), we now have CORRECT
values for all quantities through n=7, and partially through n=8.

This script compiles all known exact sequences and analyzes patterns.

Author: opus-2026-03-23-S213
"""

import sys
from math import comb, factorial, gcd
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
banner("ALL KNOWN EXACT META-GRAPH SEQUENCES")
# ============================================================

# V = iso classes (A000568)
V = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536}

# T = total arc-orbits (Burnside on (tournament, arc) pairs)
T = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}

# E = edges in simple undirected meta-graph G_n
E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161}

# G = gap = T - 2E = SL + excess
G = {n: T[n] - 2*E[n] for n in E}

# SL = self-loop arc-orbits (CORRECT, from orbit computation)
SL = {3:2, 4:6, 5:16, 6:58, 7:324}

# excess = cross-orbit excess beyond 2 per edge
excess = {3:0, 4:0, 5:12, 6:66, 7:416}

# D = degeneracy = sum C(mult, 2)
D = {3:1, 4:6, 5:28, 6:122, 7:632}

# S_2 = second moment = T + 2D (CORRECT values)
S2 = {n: T[n] + 2*D[n] for n in D}

print("  THE CORRECTED META-GRAPH SEQUENCES")
print()
print("  n  | V     | T      | E     | G=T-2E | SL   | excess")
print("  ---|-------|--------|-------|--------|------|--------")
for n in range(3, 8):
    print("  %d  | %5d | %6d | %5d | %5d  | %4d | %5d" %
          (n, V[n], T[n], E[n], G[n], SL[n], excess[n]))
print("  8  | %5d | %6d | %5d | %5d  |  ??? |   ???" %
      (V[8], T[8], E[8], G[8]))

print()
print("  n  | D     | S_2    | S_2/T  | D/V   | G/T")
print("  ---|-------|--------|--------|-------|------")
for n in range(3, 8):
    print("  %d  | %5d | %6d | %.4f | %.3f | %.4f" %
          (n, D[n], S2[n], S2[n]/T[n], D[n]/V[n], G[n]/T[n]))

# ============================================================
banner("CORRECTED S_2 vs WRONG S_2")
# ============================================================

S2_wrong = {n: 3*T[n] - 4*E[n] for n in E}

print("  n  | S_2_correct | S_2_wrong | difference")
print("  ---|-------------|-----------|----------")
for n in range(3, 8):
    s2c = S2[n]
    s2w = S2_wrong[n]
    print("  %d  | %10d  | %9d | %d" % (n, s2c, s2w, s2w - s2c))

print("""
  At n=3,4,5: the wrong formula gives the same answer (coincidence).
  At n=6: off by 4 (952 vs 948).
  At n=7: off by 216 (10392 vs 10176).
  The discrepancy grows rapidly.
""")

# ============================================================
banner("SL SEQUENCE ANALYSIS")
# ============================================================

sl_seq = [SL[n] for n in range(3, 8)]
print("  SL sequence: %s" % sl_seq)

print("\n  SL/V (avg self-loops per class):")
for n in range(3, 8):
    print("  n=%d: SL/V = %d/%d = %s" % (n, SL[n], V[n], Fraction(SL[n], V[n])))

print("\n  SL/m (fraction of arcs that self-loop, per tournament):")
for n in range(3, 8):
    m = comb(n, 2)
    # For |Aut|=1: SL per tournament = SL(C) = arcs targeting same class.
    # Average across all classes: SL / V (but weighted by what?)
    # Total self-loop arcs across all tournaments = sum_C size(C) * arcs_to_self(C)
    # = sum_C F[C][C]. And SL(C) = F[C][C]/size(C) for |Aut|=1.
    # So total_self_loop_arcs = sum_C size(C) * SL(C) approximately.
    pass

print("\n  Consecutive ratios SL(n+1)/SL(n):")
for i in range(len(sl_seq)-1):
    r = sl_seq[i+1] / sl_seq[i]
    print("  SL(%d)/SL(%d) = %d/%d = %.4f" % (i+4, i+3, sl_seq[i+1], sl_seq[i], r))

print("\n  SL factorizations:")
for n in range(3, 8):
    s = SL[n]; temp = s; factors = []
    for p in range(2, min(temp+1, 10000)):
        if p*p > temp and temp > 1: factors.append(temp); break
        while temp % p == 0: factors.append(p); temp //= p
    print("  n=%d: SL=%d = %s" % (n, s, " * ".join(str(f) for f in factors)))

# SL / (n-1)!:
print("\n  SL / (n-1)!:")
for n in range(3, 8):
    print("  n=%d: SL/(n-1)! = %d/%d = %s" %
          (n, SL[n], factorial(n-1), Fraction(SL[n], factorial(n-1))))

# ============================================================
banner("EXCESS SEQUENCE ANALYSIS")
# ============================================================

ex_seq = [excess[n] for n in range(3, 8)]
print("  excess sequence: %s" % ex_seq)

print("\n  excess factorizations:")
for n in range(3, 8):
    s = excess[n]
    if s == 0:
        print("  n=%d: excess=0" % n)
        continue
    temp = s; factors = []
    for p in range(2, min(temp+1, 10000)):
        if p*p > temp and temp > 1: factors.append(temp); break
        while temp % p == 0: factors.append(p); temp //= p
    print("  n=%d: excess=%d = %s" % (n, s, " * ".join(str(f) for f in factors)))

print("\n  excess/E (fraction of edges with multi-edges):")
for n in range(3, 8):
    if excess[n] == 0:
        print("  n=%d: excess/E = 0" % n)
    else:
        print("  n=%d: excess/E = %d/%d = %s" %
              (n, excess[n], E[n], Fraction(excess[n], E[n])))

# ============================================================
banner("THE GAP SEQUENCE G(n) = SL + excess")
# ============================================================

g_seq = [G[n] for n in range(3, 9)]
print("  G = %s" % g_seq)

print("\n  G/(n-1)! ratios:")
for n in range(3, 9):
    r = Fraction(G[n], factorial(n-1))
    print("  n=%d: G/(n-1)! = %s = %.6f" % (n, r, float(r)))

print("\n  G/(n-1)! - 1:")
for n in range(3, 9):
    r = Fraction(G[n], factorial(n-1)) - 1
    print("  n=%d: %s = %.6f" % (n, r, float(r)))

# G/n! gives:
print("\n  G/n!:")
for n in range(3, 9):
    r = Fraction(G[n], factorial(n))
    print("  n=%d: %s = %.6f" % (n, r, float(r)))

# ============================================================
banner("BURNSIDE INTERPRETATION OF SL")
# ============================================================

print("""
  SL = total self-loop arc-orbits
     = #{(class C, arc-orbit O of C) : reversing O gives back class C}

  Burnside form:
  SL = (1/n!) * sum_{sigma in S_n} #{(T, e) : sigma(T)=T, sigma(e)=e, T_e ~ T}

  The condition T_e ~ T means: reversing arc e in T gives a tournament
  isomorphic to T itself. Call this a "SELF-REVERSIBLE arc."

  For sigma = identity: every (T, e) where T_e ~ T contributes.
  #{self-reversible arcs across all T} = (1/n!) * [main term] + corrections.

  COUNTING SELF-REVERSIBLE ARCS:
  For a tournament T, arc e=(a->b) is self-reversible iff:
  exists tau in S_n with tau(T_e) = T.
  Equivalently: tau reverses the arc e back and maps everything else correctly.
  So tau maps (b,a) to some arc consistent with T.

  This is a FIXED-POINT condition: T_e = tau^{-1}(T), so
  tau is a "near-automorphism" that fixes T except at edge {a,b}.

  SIMPLER: e is self-reversible iff T and T_e are isomorphic.
  This means the UNLABELED tournament doesn't change when we flip one arc.
  The set of self-reversible arcs gives the SELF-ADJACENCY of the class
  in the meta-graph.
""")

# What fraction of arcs are self-reversible?
# At n=3: transitive has 2 self-reversible arcs (of 3). 3-cycle has 0.
# Weighted avg: (6*2 + 2*0) / (6*3 + 2*3) = 12/24 = 1/2.
# SL=2 = sum over classes of (self-reversible orbits).

# ============================================================
banner("THE KEY OBSERVATION: SL AND LEVEL EDGES")
# ============================================================

print("""
  Self-reversible arcs are arcs connecting vertices at the "same level"
  in the tournament structure. Reversing such an arc doesn't change
  the iso class because the reversal is "internal" — it permutes
  equivalent arcs without changing the overall shape.

  At n=3: in the transitive tournament 0->1->2 (with 0->2):
  - Reverse (0,1): get 1->0, still transitive (top/bottom vertices swap roles).
    This is a SELF-LOOP: H doesn't change. The arc (0,1) connects vertices
    of score 2 and 1 — their roles are "equivalent" under some relabeling.
  - Reverse (0,2): get 3-cycle. NOT self-reversible.
  - Reverse (1,2): get transitive (bottom/middle swap). SELF-LOOP.

  So the 2 self-reversible arcs are (0,1) and (1,2) — the "consecutive" arcs.
  The 1 non-self-reversible arc is (0,2) — the "skip" arc.

  PATTERN: Self-reversible arcs are those whose reversal only
  "shuffles" equivalent substructure.

  AT LARGE n: As tournaments get more "generic," fewer arcs are
  self-reversible (most reversals change the iso class).
  SL/T -> 0 confirms this.
""")

# Fraction of total arc-orbits that are self-loops
print("  SL/T (fraction of orbits that are self-loops):")
for n in range(3, 8):
    print("  n=%d: SL/T = %d/%d = %s = %.4f" %
          (n, SL[n], T[n], Fraction(SL[n], T[n]), SL[n]/T[n]))

# ============================================================
banner("FINAL CORRECTED PICTURE")
# ============================================================

print("""
  THE META-GRAPH G_n HAS:

  V(n) = iso classes (A000568)
  T(n) = total arc-orbits = V * m at large n (Burnside, EXACT)
  E(n) = undirected edges in simple graph (COMPUTED, no closed form)
  G(n) = T - 2E = self-loops + multi-edge excess
  SL(n) = self-loop orbits
  excess(n) = directed multi-edge excess

  KNOWN SEQUENCES:
  V:      2,    4,     12,     56,    456,    6880,   191536
  T:      4,   16,     88,    704,   8912,  188288,  6847200
  E:      1,    5,     30,    290,   4086,   91161
  G:      2,    6,     28,    124,    740,    5966,    85698
  SL:     2,    6,     16,     58,    324
  excess: 0,    0,     12,     66,    416
  D:      1,    6,     28,    122,    632
  S_2:    6,   28,    144,    948,  10176

  KEY RELATIONS:
  T = 2E + G                    (definition)
  G = SL + excess               (definition)
  D = (S_2 - T) / 2             (exact)
  S_2 = T + 2D                  (exact)
  E ≠ (T - D) / 2               (WRONG, see MISTAKE-029)
  E = (T - SL - excess) / 2     (correct, requires orbit data)
  E ~ T/2 as n → ∞              (asymptotic, since G/T → 0)

  THE OPEN PROBLEM:
  Can G(n), SL(n), or excess(n) be computed from Burnside theory
  without building the full meta-graph?
""")

# Check: do E values need correction for n=9?
# E_known at n=9 was E=3380751 (from edge_formula_final_s188).
# But that might have been computed from the wrong formula!
# Need to check the SOURCE of E(9).

# Actually: E was computed at n≤8 via nauty (gn_n8_nauty_s216).
# E at n=9 was claimed in some scripts. Let me check if it's trustworthy.

print("  STATUS OF E(9):")
print("  The value E(9) = 3380751 appears in several scripts.")
print("  If computed from T and D via wrong formula: INVALID.")
print("  If computed from direct adjacency: VALID.")
print("  Need to verify source. G(9) = T-2E = 6847200-6761502 = 85698.")
print("  (If E=3380751, G=6847200-6761502=85698.)")
print("  Actually 2*3380751 = 6761502. T-6761502 = 85698. Consistent with G sequence.")
print("  But is this E correct? Unclear without verification.")

print("\nDone. opus-2026-03-23-S213")
