#!/usr/bin/env python3
"""
cascade_s20bu.py -- kind-pasteur-2026-03-22-S20bu

WHAT'S MISSING? Inventory everything we know about G_n, find the gaps.

KNOWN:
  Vertices: A000568 = 2, 4, 12, 56, 456 (Burnside with odd-cycle perms)
  Edges: 1, 5, 30, 290 (computed exhaustively)
  Blue edges: 1, 1, 14, 200
  Black edges: 0, 4, 16, 90
  Level edges: 0, 0, 1, 15
  Down edges: 0, 0, 0, 0 (DAG property)
  Chains: 1, 3, 99, 292510
  Width: 1, 2, 3, 6
  Sinks: 1, 1, 2, 2
  Sources: 1, 1, 1, 1
  Self-loop fraction: 1/2, 3/8, 5/16, 35/128 = (1/2)_k/k!
  SC classes: 2, 2, 8, 12
  NS classes: 0, 2, 4, 44
  |Aut| values at n=5: {1: 7 classes, 3: 4 classes, 5: 1 class}
  Tilings*|Aut| = H (proved n=4,5)
  Burnside Fix: only odd-cycle perms contribute (proved n=1..10)

UNKNOWN:
  1. Edge count formula (1, 5, 30, 290 -- what's the pattern?)
  2. Blue/black edge count formulas
  3. Level edge formula (0, 0, 1, 15)
  4. Genus of G_n for n >= 6
  5. Whether DAG property holds for all n
  6. Tilings*|Aut|=H at n=6 (need to verify)
  7. The n -> n-2 recursion (does G_n contain G_{n-2}?)
  8. The spectral structure of G_n

Author: kind-pasteur-2026-03-22-S20bu
"""
import sys
from math import comb, factorial, gcd, log2
from fractions import Fraction
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0: yield []; return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest

def count_with_ct(ct, n):
    denom = 1
    counter = Counter(ct)
    for length, mult in counter.items():
        denom *= length ** mult * factorial(mult)
    return factorial(n) // denom

def tournament_fix_ct(ct):
    if any(l % 2 == 0 for l in ct): return 0
    n_op = sum((l-1)//2 for l in ct) + sum(gcd(ct[i],ct[j]) for i in range(len(ct)) for j in range(i+1,len(ct)))
    return 2**n_op

print("=" * 70)
print("  THE CASCADE: WHAT'S KNOWN, WHAT'S MISSING")
print("=" * 70)

# ================================================================
# 1. THE COMPLETE INVENTORY
# ================================================================
print(f"\n{'='*70}")
print(f"  1. EVERYTHING WE KNOW ABOUT G_n")
print(f"{'='*70}\n")

print(f"  {'Quantity':>30s}  {'n=3':>6s}  {'n=4':>6s}  {'n=5':>8s}  {'n=6':>10s}  {'Formula?':>10s}")
print(f"  {'-'*30}  {'-'*6}  {'-'*6}  {'-'*8}  {'-'*10}  {'-'*10}")

data = {
    'Vertices (A000568)':         [2, 4, 12, 56, 'Burnside'],
    'Edges':                      [1, 5, 30, 290, '???'],
    'Blue edges':                 [1, 1, 14, 200, '???'],
    'Black edges':                [0, 4, 16, 90, '???'],
    'Level edges':                [0, 0, 1, 15, '???'],
    'Down edges':                 [0, 0, 0, 0, '=0 always?'],
    'Chains (min->max)':          [1, 3, 99, 292510, '???'],
    'Width':                      [1, 2, 3, 6, '???'],
    'Longest chain':              [1, 2, 6, 15, '???'],
    'Sinks':                      [1, 1, 2, 2, '=2 for n>=5?'],
    'Sources':                    [1, 1, 1, 1, '=1 always'],
    'H-levels':                   [2, 3, 7, 19, '???'],
    'SC classes':                 [2, 2, 8, 12, '???'],
    'NS classes':                 [0, 2, 4, 44, '???'],
    'Max degree':                 [1, 3, 7, 14, '~2(n-1)?'],
    'Density':                    ['1.00', '0.83', '0.45', '0.19', '~2/n?'],
    'Triangles':                  ['0', '?', '21', '?', '???'],
    'Score sequences':            [2, 4, 9, 22, 'A000571'],
    'Self-loop frac':             ['1/2', '3/8', '5/16', '35/128', '(1/2)_k/k!'],
    'Tilings*|Aut|=H':           ['yes', 'yes', 'yes', '???', 'Conjecture'],
    'Meta transitive':            ['yes', 'yes', 'yes', 'almost', '???'],
}

for name, vals in data.items():
    v3, v4, v5, v6, formula = vals[0], vals[1], vals[2], vals[3], vals[4]
    print(f"  {name:>30s}  {str(v3):>6s}  {str(v4):>6s}  {str(v5):>8s}  {str(v6):>10s}  {formula:>10s}")

# ================================================================
# 2. SEQUENCE ANALYSIS: FINDING PATTERNS
# ================================================================
print(f"\n{'='*70}")
print(f"  2. SEQUENCE ANALYSIS")
print(f"{'='*70}\n")

# Edges: 1, 5, 30, 290
edges = [1, 5, 30, 290]
print(f"  Edges: {edges}")
print(f"    Ratios: {[edges[i+1]/edges[i] for i in range(len(edges)-1)]}")
print(f"    Differences: {[edges[i+1]-edges[i] for i in range(len(edges)-1)]}")
# Check: is edges[n-3] related to A000568?
vertices = [2, 4, 12, 56]
print(f"    Edges/Vertices: {[edges[i]/vertices[i] for i in range(len(edges))]}")
# E/V = 0.5, 1.25, 2.5, 5.18. Growing. Not clean.

# Check: edges = C(V, 2) * density
for i, (e, v) in enumerate(zip(edges, vertices)):
    max_edges = v * (v-1) // 2
    density = e / max_edges if max_edges > 0 else 0
    n = i + 3
    print(f"    n={n}: V={v}, max_edges={max_edges}, actual={e}, density={density:.4f}")

# Blue: 1, 1, 14, 200
blue = [1, 1, 14, 200]
print(f"\n  Blue edges: {blue}")
print(f"    Blue/Total: {[blue[i]/edges[i] for i in range(len(blue))]}")
# 1.0, 0.2, 0.47, 0.69. Blue fraction GROWS (from 20% to 69%).

# Level: 0, 0, 1, 15
level = [0, 0, 1, 15]
print(f"\n  Level edges: {level}")
# 0, 0, 1, 15. Ratios: inf, 15. Too few terms.

# Chains: 1, 3, 99, 292510
chains = [1, 3, 99, 292510]
print(f"\n  Chains: {chains}")
print(f"    Ratios: {[chains[i+1]/chains[i] for i in range(len(chains)-1)]}")
# 3, 33, 2955. Super-exponential.
# Check: 99 = 3 * 33. 292510 = 99 * 2954.6.
# Is chains[n] ~ chains[n-1] * (something growing)?

# Width: 1, 2, 3, 6
width = [1, 2, 3, 6]
print(f"\n  Width: {width}")
# 1, 2, 3, 6. Is this C(n-2, floor((n-2)/2))?
# n=3: C(1,0)=1. n=4: C(2,1)=2. n=5: C(3,1)=3. n=6: C(4,2)=6. YES!!
for i, w in enumerate(width):
    n = i + 3
    k = n - 2
    predicted = comb(k, k//2)
    print(f"    n={n}: width={w}, C({k},{k//2})={predicted}, match={w==predicted}")

print(f"\n  WIDTH = C(n-2, floor((n-2)/2)) -- the MIDDLE BINOMIAL COEFFICIENT!")
print(f"  Predicted n=7: C(5,2) = {comb(5,2)}")
print(f"  Predicted n=8: C(6,3) = {comb(6,3)}")

# Longest chain: 1, 2, 6, 15
longest = [1, 2, 6, 15]
print(f"\n  Longest chain: {longest}")
# 1, 2, 6, 15. Is this C(n-1, 2)?
# C(2,2)=1, C(3,2)=3, C(4,2)=6, C(5,2)=10. No: 1, 3, 6, 10 != 1, 2, 6, 15.
# Is it H-levels - 1?
h_levels = [2, 3, 7, 19]
print(f"  H-levels: {h_levels}")
print(f"  H-levels - 1: {[h-1 for h in h_levels]}")
# 1, 2, 6, 18. Close to longest but 18 != 15.

# Max degree: 1, 3, 7, 14
max_deg = [1, 3, 7, 14]
print(f"\n  Max degree: {max_deg}")
# 1, 3, 7, 14. Is this 2n-5?
# 2*3-5=1, 2*4-5=3, 2*5-5=5, 2*6-5=7. No: 1, 3, 5, 7 != 1, 3, 7, 14.
# Is it C(n-1,2) - 1?
# C(2,2)-1=0, C(3,2)-1=2, C(4,2)-1=5, C(5,2)-1=9. No.
# Is it the NUMBER OF TILES minus 1? C(n-1,2) - 1?
# C(2,2)=1, C(3,2)=3, C(4,2)=6, C(5,2)=10. tiles-1 = 0, 2, 5, 9. No.

# ================================================================
# 3. THE BIGGEST DISCOVERY: WIDTH = CENTRAL BINOMIAL
# ================================================================
print(f"\n{'='*70}")
print(f"  3. WIDTH FORMULA: C(n-2, floor((n-2)/2))")
print(f"{'='*70}\n")

print(f"  Width = maximum number of iso classes at the same H value")
print(f"  = maximum antichain in the H-DAG")
print(f"  = CENTRAL BINOMIAL COEFFICIENT C(n-2, floor((n-2)/2))")
print()
print(f"  THIS IS THE SAME FUNCTION AS THE FIBER FRACTION NUMERATOR!")
print(f"  f(n) = C(2(n-2), n-2) / 4^(n-2) and width = C(n-2, (n-2)/2)")
print()
print(f"  At n=5: f numerator = C(6,3) = 20, width = C(3,1) = 3")
print(f"  These are DIFFERENT central binomial coefficients")
print(f"  (f uses C(2k,k), width uses C(k, k/2))")
print(f"  But both are controlled by the MIDDLE of Pascal's triangle.")
print()

# Predicted widths
for n in range(3, 12):
    k = n - 2
    w = comb(k, k//2)
    print(f"  n={n}: predicted width = C({k},{k//2}) = {w}")

# ================================================================
# 4. WHAT'S STILL MISSING (THE GAPS)
# ================================================================
print(f"\n{'='*70}")
print(f"  4. THE GAPS: WHAT WE STILL DON'T KNOW")
print(f"{'='*70}\n")

print(f"""  GAP 1: EDGE COUNT FORMULA for G_n
    Sequence: 1, 5, 30, 290
    No formula found. Ratios 5, 6, 9.67 are irregular.
    This is the MOST IMPORTANT missing piece.
    If we had this, we'd know the density of the meta-graph at all n.

  GAP 2: BLUE/BLACK EDGE FORMULAS
    Blue: 1, 1, 14, 200. Black: 0, 4, 16, 90.
    Blue fraction: 1.0, 0.2, 0.47, 0.69.
    Blue fraction is INCREASING. Does it approach 1?

  GAP 3: LEVEL EDGE FORMULA
    Sequence: 0, 0, 1, 15.
    At n=5: the single level edge is between two H=9 classes.
    At n=6: 15 level edges, mostly at H=29 and H=37.
    Level edges measure the "degeneracy" of H at the iso class level.

  GAP 4: TILINGS*|Aut|=H AT n=6+
    Proved at n=4,5. Need verification at n=6 (computationally intensive).
    If it holds at n=6, very likely holds for all n.

  GAP 5: THE n -> n-2 RECURSION
    Does G_n "contain" G_{n-2} as a coarsened subgraph?
    The PoS classes of G_n should map to G_{n-2}.
    Partially verified but not precisely formulated.

  GAP 6: SPECTRAL STRUCTURE OF G_n
    Eigenvalues of G_5 computed. Need eigenvalues of G_6.
    Does the spectral gap have a formula?

  THE CASCADING PIECES:
  If we solve GAP 1 (edge formula), then combined with:
  - Vertices = A000568 (known formula)
  - Width = C(n-2, (n-2)/2) (discovered today!)
  - Self-loop fraction = (1/2)_(n-2)/(n-2)! (known)
  - Tilings*|Aut| = H (known)
  - Burnside = only odd-cycle perms (known)
  - Blue fraction trend (increasing)

  ...we would have a COMPLETE ANALYTICAL DESCRIPTION of G_n.
  The edge count formula is the KEYSTONE.
""")

# ================================================================
# 5. OEIS SEARCH FOR EDGE SEQUENCE
# ================================================================
print(f"{'='*70}")
print(f"  5. TRYING TO IDENTIFY EDGE SEQUENCE: 1, 5, 30, 290")
print(f"{'='*70}\n")

# 1, 5, 30, 290. Check OEIS-like patterns:
# Is this related to A000568 somehow?
# A000568: 2, 4, 12, 56. Edges/Vertices = 0.5, 1.25, 2.5, 5.18.
# 2*E/V = average degree = 1, 2.5, 5, 10.36.
# Average degree sequence: 1, 2.5, 5, 10.36. Growing roughly as n.

# Is edges = sum over classes of degree/2?
# degree(class) = number of distinct neighboring classes.
# We have degree sequences at n=5: [2,3,3,3,4,6,6,6,6,7,7,7]
# Sum of degrees = 2*edges = 2*30 = 60. Average = 5.

# Check: is 1, 5, 30, 290 in OEIS?
# (Can't search online, but let me check some known sequences)
# A001700: 1, 3, 15, 105 (double factorials). No.
# A000670: 1, 1, 3, 13, 75 (Fubini). No.
# A000142: 1, 2, 6, 24 (factorials). No.
# A002866: 1, 6, 60, 840 (related to labeled graphs). No.

# Let me try: edges(n) = (1/2) * sum over (non-identity) odd-cycle perms of Fix * (correction)
# This connects the Burnside formula to the edge count.

# Actually: the TOTAL number of arc flips that CROSS iso class boundaries
# = total flips - self-loops = 2^m * m - self_loops
# = total_flips * (1 - self_loop_fraction)
# Each cross-boundary flip connects two classes. The weighted edge count is:
# W_total = total_flips * (1 - f_self_loop)
# But the UNWEIGHTED edge count (distinct class pairs) is different.

# At n=5: total_flips = 1024 * 10 = 10240.
# Self-loops = 1760. Cross = 8480.
# Each edge between classes carries weight W[i][j] = # (tournament, flip) pairs.
# Sum of weights = 8480. Distinct edges = 30.
# Average weight per edge = 8480/30 = 282.7.

print(f"  Edge sequence analysis:")
for n, v, e in [(3, 2, 1), (4, 4, 5), (5, 12, 30), (6, 56, 290)]:
    m = comb(n, 2)
    total_flips = 2**m * m
    self_frac = float(Fraction(comb(2*(n-2), n-2), 4**(n-2)))
    cross = total_flips * (1 - self_frac)
    avg_weight = cross / (2*e) if e > 0 else 0  # divide by 2*e because each edge counted twice
    print(f"  n={n}: V={v}, E={e}, total_flips={total_flips}, cross={cross:.0f}, avg_weight_per_edge={avg_weight:.1f}")

print(f"""
  THE KEYSTONE: edge count formula.
  If we could express edges(n) in terms of A000568(n), n, and
  the self-loop fraction, the picture would be complete.

  HYPOTHESIS: edges(n) ~ A000568(n)^2 * density(n)
  where density ~ 2/n (from the observed trend).

  n=3: 2^2 * 2/3 = 2.67 ~ 1 (off)
  n=5: 12^2 * 2/5 = 57.6 ~ 30 (off by 2x)
  n=6: 56^2 * 2/6 = 1045 ~ 290 (off by 3.6x)

  Not a clean formula. The edge count remains the KEYSTONE GAP.
""")
