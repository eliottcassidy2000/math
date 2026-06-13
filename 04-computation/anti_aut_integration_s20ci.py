#!/usr/bin/env python3
"""
anti_aut_integration_s20ci.py -- kind-pasteur-2026-03-22-S20ci

INTEGRATING OPUS S214: Anti-automorphism Burnside formula.

5 NEW EXACT SEQUENCES (from opus S214):

SC_n (SC iso classes): 2, 2, 8, 12, 88, 176, 2752, 8784, 279968
                       n=3, 4, 5, 6, 7,  8,   9,  10,   11

V_merged = (A000568 + SC)/2:
  n=3: (2+2)/2 = 2
  n=4: (4+2)/2 = 3
  n=5: (12+8)/2 = 10
  n=6: (56+12)/2 = 34
  n=7: (456+88)/2 = 272
  n=8: (6880+176)/2 = 3528

Author: kind-pasteur-2026-03-22-S20ci
"""
import sys
from math import comb, log2
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("  ANTI-AUTOMORPHISM INTEGRATION: COMPLETE PICTURE")
print("=" * 70)

# Exact data
A000568 = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}
SC = {3:2, 4:2, 5:8, 6:12, 7:88, 8:176, 9:2752, 10:8784, 11:279968}
E_orig = {3:1, 4:5, 5:30, 6:290, 7:4086}
collapsed = {3:0, 4:0, 5:0, 6:5}
twin = {3:0, 4:2, 5:9, 6:142}

# ================================================================
# 1. THE COMPLETE TABLE
# ================================================================
print(f"\n{'='*70}")
print(f"  1. THE COMPLETE TABLE")
print(f"{'='*70}\n")

print(f"  {'n':>3s} {'A000568':>8s} {'SC':>6s} {'NS':>6s} {'NS/2':>6s} {'V_merg':>7s} {'E_orig':>7s} {'E_merg':>7s} {'Coll':>5s} {'Twin':>5s}")
for n in range(3, 11):
    V = A000568.get(n, '?')
    sc = SC.get(n, '?')
    if isinstance(V, int) and isinstance(sc, int):
        ns = V - sc
        ns_pairs = ns // 2
        v_merged = sc + ns_pairs
    else:
        ns = '?'
        ns_pairs = '?'
        v_merged = '?'

    e = E_orig.get(n, '?')
    c = collapsed.get(n, '?')
    t = twin.get(n, '?')
    if isinstance(e, int) and isinstance(c, int) and isinstance(t, int):
        e_merged = e - c - t
    else:
        e_merged = '?'

    print(f"  {n:>3d} {str(V):>8s} {str(sc):>6s} {str(ns):>6s} {str(ns_pairs):>6s} {str(v_merged):>7s} {str(e):>7s} {str(e_merged):>7s} {str(c):>5s} {str(t):>5s}")

# ================================================================
# 2. SC FRACTION ANALYSIS
# ================================================================
print(f"\n{'='*70}")
print(f"  2. SC FRACTION ANALYSIS")
print(f"{'='*70}\n")

print(f"  {'n':>3s} {'SC/V':>8s} {'SC frac':>9s} {'Odd n?':>7s}")
for n in range(3, 11):
    V = A000568.get(n)
    sc = SC.get(n)
    if V and sc:
        frac = sc / V
        odd = n % 2 == 1
        print(f"  {n:>3d} {sc:>4d}/{V:<4d} {frac:>8.4f} {'ODD' if odd else 'even':>7s}")

print(f"""
  SC FRACTION PATTERN:
    n=3: 100%   (odd, all SC)
    n=4: 50%    (even)
    n=5: 66.7%  (odd)
    n=6: 21.4%  (even)
    n=7: 19.3%  (odd) -- DRAMATIC DROP from n=5's 66.7%!
    n=8: 2.6%   (even)
    n=9: 1.4%   (odd)
    n=10: 0.09% (even)

  The SC fraction -> 0 at ALL parities.
  Odd n has higher SC fraction than even n, but both decrease.
  At n=10: only 0.09% of iso classes are SC!
""")

# ================================================================
# 3. MERGED VERTEX SEQUENCE
# ================================================================
print(f"{'='*70}")
print(f"  3. MERGED VERTEX SEQUENCE")
print(f"{'='*70}\n")

print(f"  V_merged = (A000568 + SC) / 2:")
for n in range(3, 11):
    V = A000568.get(n)
    sc = SC.get(n)
    if V and sc:
        vm = (V + sc) // 2
        ratio_to_V = vm / V
        print(f"  n={n}: V_merged = ({V} + {sc})/2 = {vm} ({100*ratio_to_V:.1f}% of A000568)")

print(f"""
  V_merged / A000568 -> 50% as SC fraction -> 0.
  At large n: V_merged ~ A000568 / 2 (every class pairs with complement).
""")

# ================================================================
# 4. CORRECTED n=7 PREDICTIONS
# ================================================================
print(f"{'='*70}")
print(f"  4. CORRECTED n=7 PREDICTIONS (SC=88, not 240)")
print(f"{'='*70}\n")

n = 7
V = 456
sc = 88
ns = V - sc  # = 368
ns_pairs = ns // 2  # = 184
v_merged = sc + ns_pairs  # = 272

print(f"  n=7: V={V}, SC={sc}, NS={ns}, NS_pairs={ns_pairs}")
print(f"  V_merged = {v_merged}")
print(f"  SC fraction: {sc/V:.3f} = {100*sc/V:.1f}%")
print()

# My earlier estimate was SC~240, giving V_merged~348.
# Actual SC=88, V_merged=272. Off by ~28%.
print(f"  CORRECTION: My estimate was SC~240, actual is SC=88.")
print(f"  V_merged estimate was 348, actual is 272.")
print(f"  The SC fraction at odd n drops MUCH faster than I predicted.")

# Updated descent ratio
print(f"\n  DESCENT RATIOS (corrected):")
v_merged_seq = {3: 2, 4: 3, 5: 10, 6: 34, 7: 272}
for n in [5, 6, 7]:
    if n-2 in v_merged_seq:
        ratio = v_merged_seq[n] / v_merged_seq[n-2]
        a_ratio = A000568[n] / A000568[n-2]
        print(f"  n={n}: V_merged ratio = {v_merged_seq[n]}/{v_merged_seq[n-2]} = {ratio:.1f}, A000568 ratio = {a_ratio:.1f}")

# ================================================================
# 5. THE BLUE FRACTION PREDICTION (corrected)
# ================================================================
print(f"\n{'='*70}")
print(f"  5. BLUE FRACTION (corrected for SC=88)")
print(f"{'='*70}\n")

# At n=7: 88 SC classes, 368 NS classes.
# NS-NS connections are BLUE. SC-NS connections are BLACK.
# SC-SC connections are BLUE.

# The NS fraction = 368/456 = 80.7%.
# If random connections: blue fraction ~ NS^2 + SC^2 = 0.807^2 + 0.193^2 = 0.689
# This is a LOWER BOUND because NS classes have higher degree.

ns_frac = ns / V
sc_frac = sc / V
blue_rand = ns_frac**2 + sc_frac**2
print(f"  NS fraction: {ns_frac:.3f}")
print(f"  SC fraction: {sc_frac:.3f}")
print(f"  Random blue fraction: ns^2 + sc^2 = {blue_rand:.3f}")
print(f"  Previous data: n=6 blue=69%, n=5 blue=47%")
print(f"  Predicted n=7 blue: ~{100*blue_rand:.0f}%-80%")

# ================================================================
# 6. THE GRAND ANALYTICAL PICTURE (UPDATED)
# ================================================================
print(f"\n{'='*70}")
print(f"  6. THE GRAND ANALYTICAL PICTURE")
print(f"{'='*70}\n")

print(f"""  WITH THE ANTI-AUTOMORPHISM FORMULA, WE NOW HAVE:

  EXACT FORMULAS (computable at any n from cycle type enumeration):
    1. A000568(n) = tournament iso class count [Burnside, odd-cycles]
    2. SC(n) = SC iso class count [anti-aut Burnside, opus S214]
    3. V_merged(n) = (A000568 + SC) / 2 [from 1 and 2]
    4. f(n) = (1/2)_{{n-2}} / (n-2)! [fiber fraction]
    5. Width(n) = C(n-2, (n-2)//2) [central binomial]
    6. Fix(sigma) for tournaments [only odd cycles contribute]
    7. Fix_anti(sigma) for anti-automorphisms [opus S214]
    8. Tilings * |Aut| = H [orbit-stabilizer]
    9. Weight symmetry W[i,j] = W[j,i]
   10. Down edges = 0 always (DAG property)

  NEAR-FORMULAS (asymptotically exact):
   11. E ~ V*m*(1-f)/2 [edge approximation, 95%+ at n>=6]
   12. Blue fraction -> 1 [NS-NS dominates]
   13. T_n/(2*E_n) -> 1 [transition orbits approach 2*edges]

  THE n -> n-2 DESCENT:
   14. G_n/Z_2 -> G_(n-2)/Z_2 via PoS projection
   15. Descent ratio ~ A000568(n)/A000568(n-2)

  PHASE TRANSITIONS AT n=6:
   16. H-convexity FAILS
   17. Degree-H correlation flips sign
   18. Assortativity flips
   19. Alpha_2 onset
   20. Morse secondary peak

  THIS IS A NEAR-COMPLETE ANALYTICAL DESCRIPTION OF G_n FOR ALL n.
  The only remaining unknowns are:
    - Exact edge count formula (have near-formula)
    - Spectral structure
    - Collapsed/twin edge formulas
""")
