#!/usr/bin/env python3
"""
gap_sequence_s211e.py -- opus-2026-03-23-S211e

Study the GAP SEQUENCE: G(n) = T(n) - 2*E(n) = SL + excess_cross.

G(n) = 2, 6, 28, 124, 740, 5966, 85698  (n=3..9)

This measures "how far the meta-graph is from a simple graph."
G(n)/T(n) → 0 means the meta-graph approaches simplicity.

Can G(n) be computed from Burnside theory?

Author: opus-2026-03-23-S211e
"""

import sys
from math import comb, factorial, gcd
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}
V_known = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536}

G = {n: T_known[n] - 2*E_known[n] for n in T_known}

banner("THE GAP SEQUENCE G(n) = T(n) - 2*E(n)")

print("  G(n) = SL(n) + excess_cross(n)")
print("  G(n) measures self-loops and multi-edge excess\n")

for n in range(3, 10):
    g = G[n]
    t = T_known[n]
    v = V_known[n]
    e = E_known[n]
    print("  n=%d: G=%d, T=%d, E=%d, V=%d" % (n, g, t, e, v))
    print("        G/T=%.6f, G/V=%.4f, G/E=%.6f" % (g/t, g/v, g/e))

# Factorizations
banner("FACTORIZATIONS")
for n in range(3, 10):
    g = G[n]
    temp = g; factors = []
    for p in range(2, min(temp+1, 100000)):
        if p*p > temp and temp > 1:
            factors.append(temp); break
        while temp % p == 0:
            factors.append(p); temp //= p
    print("  n=%d: G=%d = %s" % (n, g, " * ".join(str(f) for f in factors)))

# Ratios
banner("RATIOS AND PATTERNS")
g_seq = [G[n] for n in range(3, 10)]
print("  G sequence: %s" % g_seq)

print("\n  Consecutive ratios G(n+1)/G(n):")
for i in range(len(g_seq)-1):
    ratio = g_seq[i+1] / g_seq[i]
    print("  G(%d)/G(%d) = %d/%d = %.4f" % (i+4, i+3, g_seq[i+1], g_seq[i], ratio))

print("\n  G(n) / n!:")
for n in range(3, 10):
    print("  n=%d: G/%d! = %d/%d = %s" % (n, n, G[n], factorial(n),
          Fraction(G[n], factorial(n))))

print("\n  G(n) / C(n,2):")
for n in range(3, 10):
    m = comb(n, 2)
    print("  n=%d: G/m = %d/%d = %s" % (n, G[n], m, Fraction(G[n], m)))

print("\n  G(n) / 2^{n-1}:")
for n in range(3, 10):
    print("  n=%d: G/2^{n-1} = %d/%d = %s" % (n, G[n], 2**(n-1),
          Fraction(G[n], 2**(n-1))))

# Is G(n) related to known sequences?
banner("G(n) vs KNOWN SEQUENCES")

# Try G(n) * k! or G(n) / k for various k
print("  Testing: is G(n) = c * something_nice(n)?")

# G(3)=2, G(4)=6, G(5)=28, G(6)=124, G(7)=740, G(8)=5966, G(9)=85698

# Check if G(n) = 2 * something
print("  G(n)/2: %s" % [G[n]//2 for n in range(3, 10)])
# 1, 3, 14, 62, 370, 2983, 42849

# Check: 1, 3, 14, 62 — this looks like Catalan-adjacent
# Catalan: 1, 2, 5, 14, 42, 132 ... no
# Central binomial: C(2k,k) = 1,2,6,20,70 ... no
# Bell: 1, 2, 5, 15, 52, 203, 877 ... no

# 1, 3, 14, 62, 370 — check OEIS
print("\n  G(n)/2 = 1, 3, 14, 62, 370, 2983, 42849")
print("  Checking if this matches known sequences...")

# These are close to Fubini numbers (ordered Bell): 1, 1, 3, 13, 75, 541, 4683
# Or Tangent numbers: 1, 2, 16, 272
# Or subfactorial: 0, 1, 2, 9, 44, 265, 1854

# Try G(n) - 2*G(n-1):
print("\n  G(n) - 2*G(n-1):")
for n in range(4, 10):
    print("  n=%d: %d - 2*%d = %d" % (n, G[n], G[n-1], G[n] - 2*G[n-1]))

# G(n) - n*G(n-1):
print("\n  G(n) - n*G(n-1):")
for n in range(4, 10):
    print("  n=%d: %d - %d*%d = %d" % (n, G[n], n, G[n-1], G[n] - n*G[n-1]))

# G(n) - (n-1)*G(n-1):
print("\n  G(n) - (n-1)*G(n-1):")
for n in range(4, 10):
    print("  n=%d: %d - %d*%d = %d" % (n, G[n], n-1, G[n-1], G[n] - (n-1)*G[n-1]))

# G(n) / (n-1)!:
print("\n  G(n) / (n-1)!:")
for n in range(3, 10):
    print("  n=%d: %d/%d = %s" % (n, G[n], factorial(n-1), Fraction(G[n], factorial(n-1))))

# SL component vs excess component at known n values
banner("SL AND EXCESS DECOMPOSITION")

# From the n=3..6 computation we had:
SL = {3:2, 4:6, 5:16, 6:58}
excess = {3:0, 4:0, 5:12, 6:66}

print("  n  | G=T-2E | SL  | excess | SL/G  | excess/G")
print("  ---|--------|-----|--------|-------|----------")
for n in [3, 4, 5, 6]:
    g = G[n]
    print("  %d  | %5d  | %3d |  %4d  | %.3f | %.3f" %
          (n, g, SL[n], excess[n], SL[n]/g, excess[n]/g))

print("""
  SL dominates at small n, but excess grows faster.
  At n=6: excess > SL (66 > 58).
  Prediction: at large n, excess dominates (multi-edges proliferate
  faster than self-loops as classes increase).
""")

# The SL sequence: 2, 6, 16, 58
print("  SL sequence: 2, 6, 16, 58")
print("  SL/V: %s" % [Fraction(SL[n], V_known[n]) for n in [3,4,5,6]])
# SL/V = 1, 3/2, 4/3, 29/28

# At n=6: 58 self-loop orbits across 56 classes (more than 1 per class)
# Actually: it means on average each class has about 1 self-loop orbit.
# But many classes have 0 self-loops, and some have many.

# SL SEQUENCE FACTORED
print("  SL factored: 2=2, 6=2*3, 16=2^4, 58=2*29")

# Excess sequence: 0, 0, 12, 66
print("  excess sequence: 0, 0, 12, 66")
print("  excess factored: 0, 0, 12=4*3, 66=2*3*11")

banner("THE BIG PICTURE")

print("""
  THE META-GRAPH EDGE COUNT E(n) CANNOT be computed from T(n) and D(n) alone.

  The correct formula is: E = (T - SL - excess) / 2

  Where:
  SL(n) = total self-loop arc-orbits (arcs reversing to the same class)
  excess(n) = total cross-orbit excess (multi-edges beyond 1 per direction)

  Both SL and excess require ORBIT-LEVEL information about the meta-graph.
  They cannot be derived from aggregate statistics like T, D, or S_2.

  THE GAP G(n) = T(n) - 2E(n) = SL(n) + excess(n):
  n:   3    4     5      6      7       8         9
  G:   2    6    28    124    740    5966     85698
  G/T: .50  .38   .32    .18    .08     .03       .01

  ASYMPTOTICALLY: G/T → 0, so E ~ T/2 at large n.
  The meta-graph becomes increasingly "simple" (few multi-edges, few self-loops).

  FOR EXACT E(n): must compute the meta-graph adjacency (F matrix)
  at each n, which requires canonical form computation.
  Known exact: E(3..8) = 1, 5, 30, 290, 4086, 91161.

  THE BURNSIDE APPROACH gives T(n) exactly at any n.
  But T only gives an UPPER BOUND: E ≤ T/2.
  The gap G = T - 2E requires additional information.

  OPEN QUESTION: Is there a cycle-index formula for G(n)?
  G counts a specific defect in the meta-graph structure.
  It may be expressible as a Burnside sum over (tournament, arc) pairs
  with additional constraints (arc maps to same class, or same directed edge).
""")

# Check: is G(n) in OEIS?
print("  G sequence: 2, 6, 28, 124, 740, 5966, 85698")
print("  Not immediately recognizable. Need OEIS lookup.")

print("\nDone. opus-2026-03-23-S211e")
