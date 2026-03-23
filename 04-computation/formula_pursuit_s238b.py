#!/usr/bin/env python3
"""
formula_pursuit_s238b.py -- opus-2026-03-23-S238b

PURSUIT OF FORMULAS for the merged meta-graph.

Target: express E_merged, #BLUE, #BLACK, #GREEN, #triangles as functions of n.

Approach: collect sequences, search for patterns, test conjectures.

Known exact sequences:
  V_merged:  2, 3, 10, 34, 272     (n=3..7, NOT in OEIS)
  E_merged:  1, 3, 21, 143, 2123   (n=3..7)
  BLUE:      1, 1, 12, 13, 174     (n=3..7, SC-SC edges)
  BLACK:     0, 2, 8, 45, 550      (n=3..7, SC-NS edges)
  GREEN:     0, 0, 1, 85, 1399     (n=3..7, NS-NS edges)
  Tri:       ?, ?, 12, 139, ?      (n=3..7, triangles in merged graph)
  Width:     1, 1, 2, 3, 8, 25     (n=3..8, max level size)
  AMBER:     0, 0, 3, 25, 393      (n=3..7, score-preserving edges)
  LEVEL:     0, 0, 1, 5, 71        (n=3..7, dH=0 edges)

Author: opus-2026-03-23-S238b
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
banner("ALL KNOWN MERGED META-GRAPH SEQUENCES")
# ============================================================

# A000568 values (iso classes)
V = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}
SC = {3:1, 4:2, 5:8, 6:12, 7:88, 8:456}  # SC classes (unmerged)
# V_merged = (V + SC) / 2
VM = {n: (V[n] + SC[n]) // 2 for n in V}

# E values (unmerged)
E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161}
# Merged edges
EM = {3:1, 4:3, 5:21, 6:143, 7:2123, 8:45550}

# BLUE/BLACK/GREEN (merged)
BLUE = {3:1, 4:1, 5:12, 6:13, 7:174, 8:319}
BLACK = {3:0, 4:2, 5:8, 6:45, 7:550, 8:1894}
GREEN = {3:0, 4:0, 5:1, 6:85, 7:1399, 8:43337}

# Other
AMBER = {5:3, 6:25, 7:393}
LEVEL = {5:1, 6:5, 7:71}
WIDTH = {3:1, 4:1, 5:2, 6:3, 7:8, 8:25}
TRI = {5:12, 6:139}

print("  n  | VM  | EM    | BLUE | BLACK | GREEN | AMBER | LEVEL | WIDTH | 2E/VM")
print("  ---|-----|-------|------|-------|-------|-------|-------|-------|------")
for n in range(3, 9):
    vm = VM.get(n, '?')
    em = EM.get(n, '?')
    bl = BLUE.get(n, '?')
    bk = BLACK.get(n, '?')
    gr = GREEN.get(n, '?')
    am = AMBER.get(n, '?')
    lv = LEVEL.get(n, '?')
    wd = WIDTH.get(n, '?')
    avg_d = "%.1f" % (2*em/vm) if isinstance(em, int) and isinstance(vm, int) else '?'
    print("  %d  | %s | %s | %s | %s | %s | %s | %s | %s | %s" %
          (n, vm, em, bl, bk, gr, am, lv, wd, avg_d))

# ============================================================
banner("SEQUENCE ANALYSIS: V_merged")
# ============================================================

vm_seq = [VM[n] for n in range(3, 8)]
print("  V_merged: %s" % vm_seq)
print("  Ratios: %s" % [round(vm_seq[i+1]/vm_seq[i], 2) for i in range(len(vm_seq)-1)])
print("  V_merged / V: %s" % [Fraction(VM[n], V[n]) for n in range(3, 8)])
print("  V_merged / SC: %s" % [Fraction(VM[n], SC[n]) for n in range(3, 8)])

# V_merged = (A000568 + SC) / 2
# SC is A000171 (self-complementary tournaments)
# Check OEIS for these

# ============================================================
banner("SEQUENCE ANALYSIS: E_merged")
# ============================================================

em_seq = [EM[n] for n in range(3, 9)]
print("  E_merged: %s" % em_seq)
print("  Ratios: %s" % [round(em_seq[i+1]/em_seq[i], 2) for i in range(len(em_seq)-1)])
print("  E_merged / VM: %s" % [round(EM[n]/VM[n], 2) for n in range(3, 8)])
print("  E_merged / E: %s" % [Fraction(EM[n], E[n]) for n in range(3, 8)])

# E_merged vs E: ratio approaches 1/2 as complements merge
print("  E_merged/E approaching 1/2:")
for n in range(3, 9):
    if n in EM and n in E:
        print("    n=%d: %d/%d = %.4f" % (n, EM[n], E[n], EM[n]/E[n]))

# ============================================================
banner("SEQUENCE ANALYSIS: BLUE (SC-SC edges)")
# ============================================================

blue_seq = [BLUE[n] for n in range(3, 9)]
print("  BLUE: %s" % blue_seq)
print("  BLUE / SC^2: %s" % [round(BLUE[n] / (SC[n]**2/2), 4) if SC[n] > 0 else '?' for n in range(3, 8)])
print("  BLUE / C(SC_merged, 2): ", end="")
SC_merged = {n: SC[n] for n in SC}  # In merged graph, SC classes stay as-is
for n in range(3, 8):
    sc_m = SC_merged[n]
    max_blue = sc_m * (sc_m - 1) // 2
    print("%s " % Fraction(BLUE[n], max_blue) if max_blue > 0 else "N/A ", end="")
print()

# Blue density = BLUE / C(SC_merged, 2) tells us how complete the SC backbone is
print("  BLUE density (fraction of possible SC-SC edges):")
for n in range(3, 9):
    if n not in SC: continue
    sc = SC[n]
    max_e = sc * (sc-1) // 2
    if max_e > 0 and n in BLUE:
        print("    n=%d: %d / %d = %.4f" % (n, BLUE[n], max_e, BLUE[n]/max_e))

# ============================================================
banner("SEQUENCE ANALYSIS: GREEN (NS-NS edges)")
# ============================================================

green_seq = [GREEN[n] for n in range(3, 9)]
print("  GREEN: %s" % green_seq)

NS = {n: VM[n] - SC[n] for n in VM}
print("  NS_merged: %s" % [NS[n] for n in range(3, 8)])
print("  GREEN density (fraction of possible NS-NS edges):")
for n in range(3, 9):
    if n not in NS or n not in GREEN: continue
    ns = NS[n]
    max_e = ns * (ns-1) // 2
    if max_e > 0:
        print("    n=%d: %d / %d = %.4f" % (n, GREEN[n], max_e, GREEN[n]/max_e))

# ============================================================
banner("LOOKING FOR SIMPLE FORMULAS")
# ============================================================

# Test: E_merged ~ V_merged * m / 2
m_vals = {n: comb(n, 2) for n in range(3, 9)}
print("  E_merged vs VM * m / 2:")
for n in range(3, 8):
    predicted = VM[n] * m_vals[n] / 2
    actual = EM[n]
    print("    n=%d: predicted=%.0f, actual=%d, ratio=%.4f" %
          (n, predicted, actual, actual/predicted if predicted > 0 else 0))

# Test: E_merged ~ VM * (VM-1) * p for some edge probability p
print("\n  E_merged / C(VM, 2) (edge density):")
for n in range(3, 8):
    max_e = VM[n] * (VM[n]-1) // 2
    if max_e > 0:
        density = EM[n] / max_e
        print("    n=%d: %d / %d = %.4f" % (n, EM[n], max_e, density))

# Test: BLUE ~ SC * (SC-1) / 6 (roughly 1/3 of max)
print("\n  BLUE / SC ~ constant?")
for n in range(3, 9):
    if n in BLUE and n in SC:
        print("    n=%d: BLUE/SC = %d/%d = %.2f" % (n, BLUE[n], SC[n], BLUE[n]/SC[n]))

# Test: BLACK ~ SC * NS_merged
print("\n  BLACK / (SC * NS_merged):")
for n in range(3, 8):
    if n in BLACK:
        product = SC[n] * NS[n]
        if product > 0:
            print("    n=%d: %d / %d = %.4f" % (n, BLACK[n], product, BLACK[n]/product))

# ============================================================
banner("RATIO SEQUENCES (looking for convergence)")
# ============================================================

print("  GREEN / EM (fraction that is NS-NS):")
for n in range(3, 9):
    if n in GREEN and n in EM:
        print("    n=%d: %.4f" % (n, GREEN[n]/EM[n]))

print("\n  BLUE / VM (avg blue degree / 2):")
for n in range(3, 9):
    if n in BLUE and n in VM:
        print("    n=%d: %.4f" % (n, BLUE[n]/VM[n]))

print("\n  Avg degree = 2*EM/VM:")
for n in range(3, 9):
    if n in EM and n in VM:
        print("    n=%d: %.2f (m=%d)" % (n, 2*EM[n]/VM[n], m_vals[n]))

print("\n  Avg degree / m:")
for n in range(3, 9):
    if n in EM and n in VM:
        print("    n=%d: %.4f" % (n, 2*EM[n]/(VM[n]*m_vals[n])))

# ============================================================
banner("THE KEY OPEN SEQUENCES")
# ============================================================

print("""
  SEQUENCES NOT IN OEIS (candidates for submission):

  1. V_merged (n=3..8): 2, 3, 10, 34, 272, 3668
     = (A000568 + A000171) / 2

  2. E_merged (n=3..8): 1, 3, 21, 143, 2123, 45550
     = edges in merged arc-reversal graph on tournament iso classes

  3. BLUE (n=3..8): 1, 1, 12, 13, 174, 319
     = SC-SC edges in merged meta-graph

  4. BLACK (n=3..8): 0, 2, 8, 45, 550, 1894
     = SC-NS edges in merged meta-graph

  5. GREEN (n=3..8): 0, 0, 1, 85, 1399, 43337
     = NS-NS edges in merged meta-graph

  6. AMBER (n=5..7): 3, 25, 393
     = score-preserving edges in merged meta-graph

  7. LEVEL (n=5..7): 1, 5, 71
     = same-H edges (dH=0) in merged meta-graph

  8. Triangles (n=5,6): 12, 139
     = triangles in merged meta-graph

  9. SL (n=3..7): 2, 6, 16, 58, 324
     = self-loop arc-orbits

  10. G = T - 2E (n=3..8): 2, 6, 28, 124, 740, 5966
      = gap sequence (self-loops + multi-edge excess)

  TOTAL: 10+ new integer sequences arising from the merged meta-graph.
""")

# ============================================================
banner("FACTORIZATION PATTERNS")
# ============================================================

sequences = {
    'V_merged': [2, 3, 10, 34, 272, 3668],
    'E_merged': [1, 3, 21, 143, 2123, 45550],
    'BLUE': [1, 1, 12, 13, 174, 319],
    'BLACK': [0, 2, 8, 45, 550, 1894],
    'GREEN': [0, 0, 1, 85, 1399, 43337],
    'LEVEL': [0, 0, 1, 5, 71],
}

for name, seq in sequences.items():
    print("  %s:" % name)
    for i, val in enumerate(seq):
        n = i + 3
        if val == 0:
            print("    n=%d: 0" % n)
            continue
        temp = val; factors = []
        for p in range(2, min(abs(temp)+1, 10000)):
            if p*p > abs(temp) and abs(temp) > 1:
                factors.append(abs(temp)); break
            while temp % p == 0:
                factors.append(p); temp //= p
        print("    n=%d: %d = %s" % (n, val, " * ".join(str(f) for f in factors) if factors else str(val)))
    print()

print("Done. opus-2026-03-23-S238b")
