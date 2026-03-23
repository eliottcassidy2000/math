#!/usr/bin/env python3
"""
blue_black_degree_s183.py -- opus-2026-03-22-S183

WHAT DETERMINES BLUE DEGREE AND BLACK DEGREE IN G_n?

From existing data (S211):
  n=5: SC classes avg blue_deg=3.0, avg black_deg=2.0
       Non-SC classes avg blue_deg=1.0, avg black_deg=4.0
  n=6: SC classes avg blue_deg=2.2, avg black_deg=7.5
       Non-SC classes avg blue_deg=8.5, avg black_deg=2.0

The REVERSAL between n=5 and n=6 is striking.
At n=5 (odd): SC classes have MORE blue edges.
At n=6 (even): Non-SC classes have MORE blue edges.

HYPOTHESIS: This is controlled by the GS (grid-symmetric) tiling count.
Blue edges = GS flip pairs. Black edges = non-GS flip pairs.
GS tilings count = 2^k where k = (m + floor((n-1)/2))/2.

For n=7: We need a FAST computation strategy.
n=7 has 2^21 = 2M tournaments, 456 iso classes.
Full enumeration is feasible (minutes, not hours).

Author: opus-2026-03-22-S183
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
import time
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
# PART 1: ANALYZE EXISTING DATA
# ============================================================

banner("PART 1: PATTERNS IN BLUE/BLACK DEGREE (from S211 data)")

# Extracted from the S211 output
data = {
    5: {
        'SC_blue': [2, 2, 2, 3, 3, 4, 4, 4],
        'SC_black': [0, 0, 2, 2, 2, 2, 4, 4],
        'NSC_blue': [1, 1, 1, 1],
        'NSC_black': [2, 2, 6, 6],
        'n_SC': 8, 'n_NSC': 4, 'total_blue': 14, 'total_black': 16,
    },
    6: {
        'SC_blue': [1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3],
        'SC_black': [2, 2, 4, 8, 8, 8, 8, 10, 10, 10, 10, 10],
        'NSC_blue': [3,3,3,3,4,4,5,5,6,6,6,6,7,7,7,7,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,12,13,13],
        'NSC_black': [0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,6,6],
        'n_SC': 12, 'n_NSC': 44, 'total_blue': 200, 'total_black': 90,
    },
}

for n in [5, 6]:
    d = data[n]
    print("  n=%d:" % n)
    print("    SC classes: %d" % d['n_SC'])
    print("      avg blue_deg = %.2f, avg black_deg = %.2f" %
          (np.mean(d['SC_blue']), np.mean(d['SC_black'])))
    print("      blue_deg range: [%d, %d]" % (min(d['SC_blue']), max(d['SC_blue'])))
    print("      black_deg range: [%d, %d]" % (min(d['SC_black']), max(d['SC_black'])))
    print("    Non-SC classes: %d" % d['n_NSC'])
    print("      avg blue_deg = %.2f, avg black_deg = %.2f" %
          (np.mean(d['NSC_blue']), np.mean(d['NSC_black'])))
    print("    Total: %d blue, %d black edges" % (d['total_blue'], d['total_black']))
    print("    Blue fraction: %.1f%%" % (100 * d['total_blue'] / (d['total_blue'] + d['total_black'])))
    print()

print("""
  KEY OBSERVATION: The blue fraction INCREASES with n.
    n=3: 100% blue (1/1)
    n=4: 20% blue (1/5)
    n=5: 47% blue (14/30)
    n=6: 69% blue (200/290)

  At even n: blue fraction is LOW (n=4: 20%).
  At odd n: blue fraction is MODERATE (n=3: 100%, n=5: 47%).
  But n=6 has HIGH blue fraction (69%)!

  THIS SUGGESTS: blue fraction grows with n, overriding the
  even/odd parity effect at larger n.
""")

# ============================================================
# PART 2: WHAT DETERMINES BLUE DEGREE?
# ============================================================

banner("PART 2: WHAT DETERMINES BLUE DEGREE?")

print("""
  A blue edge between classes i and j exists iff:
  There exist GS (grid-symmetric) tilings T_a in class i and T_b in class j
  such that T_b = flip(T_a).

  So: blue_degree(i) = number of distinct classes j (j != i) connected
  to i by at least one GS flip pair.

  FACTORS CONTROLLING BLUE DEGREE:

  1. GS TILING COUNT: How many GS tilings does class i contain?
     More GS tilings = more chances for blue connections.
     #GS tilings in class i = #GS_tilings_total * (fraction in class i)

  2. GS FLIP DIVERSITY: How many DIFFERENT classes does flipping
     the GS tilings of class i produce?
     If all GS flips land in the same class: blue_deg = 1.
     If they spread across many classes: blue_deg is large.

  3. SCORE COMPATIBILITY: A GS tiling has a palindromic score
     (score(v) + score(n-1-v) = n-1). GS flips preserve this
     palindromic structure. So blue edges connect classes whose
     scores are compatible with palindromic structure.

  4. SC STATUS: GS tilings produce SC (self-complementary) tournaments.
     So blue edges connect SC classes to SC classes.
     Non-SC classes can only have blue edges if they share GS tilings
     with SC classes... wait, non-SC classes have NO GS tilings!
     But they CAN be connected by blue edges if a GS tiling in an
     SC class flips to a tiling in the non-SC class.

  ACTUALLY: kind-pasteur proved GS tilings are CLOSED under flip.
  So flip(GS) = GS. Both endpoints of a blue edge must contain
  GS tilings, hence both must be SC classes... unless?

  Wait: from S211 data, non-SC classes at n=5 have blue_deg > 0
  (value 1 for all 4 non-SC classes). How?

  The answer: the FULL adjacency graph includes BOTH arc reversals
  (from our G_n computation) and the blue/black classification from
  the TILING MODEL. These are DIFFERENT notions of "blue" vs "black"!

  In the ARC REVERSAL model:
  - An edge is "blue" if it involves a GS arc reversal
  - An edge is "black" if it involves a non-GS arc reversal

  In the TILING FLIP model:
  - A blue line connects GS tilings (bit complement of GS = GS)
  - A black line connects non-GS tilings

  The S211 script uses the TILING model, not the arc reversal model.
  Let me verify.
""")

# ============================================================
# PART 3: FORMULA FOR BLUE DEGREE
# ============================================================

banner("PART 3: BLUE DEGREE FORMULA CANDIDATES")

print("  From the n=5 data:")
print("  SC classes (H, blue_deg, black_deg):")
n5_sc = [(1,2,4), (3,2,2), (9,4,2), (9,3,4), (11,4,2), (13,4,2), (15,3,0), (15,2,0)]
for H, bd, kd in n5_sc:
    total = bd + kd
    print("    H=%2d: b_deg=%d, k_deg=%d, total=%d, b_frac=%.2f" %
          (H, bd, kd, total, bd/(bd+kd) if bd+kd > 0 else 0))

print("\n  Checking: is blue_deg related to H?")
Hs = [h for h,_,_ in n5_sc]
bds = [b for _,b,_ in n5_sc]
corr = np.corrcoef(Hs, bds)[0,1]
print("    corr(H, blue_deg) = %.4f" % corr)

print("\n  Checking: is blue_deg related to c3 (3-cycle count)?")
# From S211 data: c3 for SC classes at n=5
c3s = [0, 1, 3, 3, 4, 4, 4, 5]
corr_c3 = np.corrcoef(c3s, bds)[0,1]
print("    corr(c3, blue_deg) = %.4f" % corr_c3)

print("\n  For n=6 SC classes:")
n6_sc = [(1,2,8), (5,2,10), (9,1,2), (17,2,10), (17,3,8), (29,3,10),
         (29,2,10), (33,3,10), (37,3,8), (41,2,8), (45,2,2), (45,1,4)]
for H, bd, kd in n6_sc:
    print("    H=%2d: b_deg=%d, k_deg=%d, total=%d" % (H, bd, kd, bd+kd))

Hs6 = [h for h,_,_ in n6_sc]
bds6 = [b for _,b,_ in n6_sc]
corr6 = np.corrcoef(Hs6, bds6)[0,1]
print("    corr(H, blue_deg) for n=6 SC = %.4f" % corr6)

# ============================================================
# PART 4: GS TILING COUNT PER CLASS
# ============================================================

banner("PART 4: GS TILING COUNT PER CLASS")

print("  The number of GS tilings in each iso class determines")
print("  how many BLUE flip pairs originate from that class.\n")

print("  Total GS tilings: 2^k where k depends on n.")
for n in [3, 4, 5, 6, 7, 8]:
    m = comb(n-1, 2)
    # k = number of free bits in GS tiling
    # From S164: k = (m + floor((n-1)/2)) / 2
    # Actually from S176: GS count = 2^k where k = diag + pairs
    # floor((n-1)/2) diagonal boxes, (m - floor((n-1)/2))/2 off-diagonal pairs
    # k = floor((n-1)/2) + (m - floor((n-1)/2)) // 2
    diag = (n-1) // 2
    pairs = (m - diag) // 2
    k = diag + pairs
    gs_count = 2**k
    print("  n=%d: m=%d, diag=%d, pairs=%d, k=%d, 2^k=%d GS tilings" %
          (n, m, diag, pairs, k, gs_count))

print()
print("  At n=7: 512 GS tilings among 2^15=32768 total tilings.")
print("  At n=8: 4096 GS tilings among 2^21=2097152 total tilings.\n")

print("  The GS fraction = 2^k / 2^m = 2^{k-m}:")
for n in [3, 4, 5, 6, 7, 8]:
    m = comb(n-1, 2)
    diag = (n-1) // 2
    pairs = (m - diag) // 2
    k = diag + pairs
    frac = Fraction(2**k, 2**m)
    print("    n=%d: GS fraction = 2^%d / 2^%d = 2^{%d} = %s" %
          (n, k, m, k-m, frac))

# ============================================================
# PART 5: PREDICTIONS FOR n=7
# ============================================================

banner("PART 5: PREDICTIONS FOR n=7")

print("""
  At n=7:
    456 iso classes (A000568)
    88 SC classes, 184 NSC pairs (from S25h data)
    512 GS tilings = 2^9

  PREDICTION: Blue edges at n=7
    At n=5: 14 blue edges out of 30 total = 47%
    At n=6: 200 blue edges out of 290 total = 69%
    Trend: blue fraction increases with n.
    Prediction: n=7 blue fraction ~ 75-85%?

  PREDICTION: Blue degree
    At n=5, max blue_deg(SC) = 4
    At n=6, max blue_deg(SC) = 3
    The blue degree of SC classes DECREASES while total degree increases.
    This is because the number of SC classes grows (8 -> 12 -> 88)
    but the GS tilings only grow as 2^k.

  GS tilings per SC class = 512/88 ~ 5.8 at n=7.
  At n=5: 16/8 = 2 GS tilings per SC class.
  At n=6: 64/12 ~ 5.3 GS tilings per SC class.

  So the average GS tiling count per SC class is GROWING.
  Each GS tiling connects to its flip partner.
  The flip partner lands in SOME class.
  If the distribution is uniform: blue_deg ~ #distinct classes hit.

  With 5.8 GS tilings per class and ~88 target classes:
  If flips are "random": blue_deg ~ 5.8 * (1 - (87/88)^5.8) ~ 5.8 * 0.064 ~ 0.37
  Wait, that gives very low blue degree. Let me reconsider.

  Actually: each of the 5.8 GS tilings has ONE flip partner.
  That partner lands in some SC class (since GS->GS under flip).
  If the partner classes are all DIFFERENT: blue_deg = 5.8.
  If some partners are the SAME class (self-flip): blue_deg < 5.8.

  At odd n (7): NO self-flip (blueself requires even n from THM-023).
  So all 5.8 flip partners go to OTHER classes.
  If they're all distinct: blue_deg ~ 5.8 ~ 6.
  If some coincide: blue_deg < 6.

  PREDICTION: avg blue_deg(SC) at n=7 ~ 5-6.
""")

# ============================================================
# PART 6: WHAT DETERMINES BLACK DEGREE?
# ============================================================

banner("PART 6: WHAT DETERMINES BLACK DEGREE?")

print("""
  A black edge between classes i and j exists iff:
  There exist NON-GS tilings T_a in class i and T_b in class j
  such that T_b = flip(T_a).

  Non-GS tilings: total - GS = 2^m - 2^k.
  At n=7: 32768 - 512 = 32256 non-GS tilings.

  Each non-GS tiling has a flip partner (also non-GS, since GS is
  closed under flip). So: 32256/2 = 16128 non-GS flip pairs.

  These 16128 pairs define the BLACK edges of G_n.
  Each pair connects some class i to some class j (possibly i=j
  for blackself classes).

  The black_degree(i) = number of distinct classes j connected to i
  by at least one non-GS flip pair.

  WHAT DETERMINES BLACK DEGREE:
  1. Tiling count: class i has #tilings = H/|Aut| tilings.
     Of these, ~H/|Aut| * (1 - 2^{k-m}) are non-GS.
  2. Flip diversity: how many different classes do the non-GS flips hit?
  3. Score compatibility: non-GS flips can change scores in ways
     that GS flips cannot. This allows connections to MORE classes.

  KEY DIFFERENCE GS vs non-GS:
  GS tilings have palindromic scores -> connect to palindromic classes.
  Non-GS tilings have arbitrary scores -> connect to any compatible class.
  So: black edges have MORE connectivity options than blue edges.

  At n=5: avg black_deg(SC) = 2.0 < avg blue_deg(SC) = 3.0
  At n=6: avg black_deg(SC) = 7.5 >> avg blue_deg(SC) = 2.2

  The REVERSAL happens because at n=6, there are MANY more non-GS
  tilings (1024-64=960) than GS tilings (64), so black edges dominate.
  But at n=5, the ratio is more balanced (64-16=48 non-GS vs 16 GS).

  PREDICTION FOR n=7:
  Non-GS tilings: 32256. GS tilings: 512. Ratio: 63:1.
  Black edges should VASTLY outnumber blue edges in total weight.
  But in the UNWEIGHTED degree count, blue edges may still be
  significant because GS flip partners are very SPREAD across classes.
""")

# ============================================================
# PART 7: THE DEGREE SUM FORMULA
# ============================================================

banner("PART 7: DEGREE SUM CONSTRAINTS")

print("  The HANDSHAKING LEMMA constrains the total degree.\n")

for n_val, data_n in [(5, data[5]), (6, data[6])]:
    all_blue = data_n['SC_blue'] + data_n['NSC_blue']
    all_black = data_n['SC_black'] + data_n['NSC_black']
    sum_blue = sum(all_blue)
    sum_black = sum(all_black)
    print("  n=%d:" % n_val)
    print("    sum blue_deg = %d = 2 * blue_edges = %d" % (sum_blue, 2*data_n['total_blue']))
    print("    sum black_deg = %d = 2 * black_edges = %d" % (sum_black, 2*data_n['total_black']))
    print("    Verified: %s, %s" % (sum_blue == 2*data_n['total_blue'],
                                     sum_black == 2*data_n['total_black']))
    print()

# ============================================================
# PART 8: FORMULAS
# ============================================================

banner("PART 8: EMERGING FORMULAS")

print("""
  BLUE EDGE COUNT FORMULA:
  At n=3: 1 blue edge, 0 black. GS tilings: 2.
  At n=4: 1 blue edge, 4 black. GS tilings: 4.
  At n=5: 14 blue, 16 black. GS tilings: 16.
  At n=6: 200 blue, 90 black. GS tilings: 64.

  Blue edges = GS_tiling_pairs / 2 - blueself_pairs
  = 2^{k-1} - (# blueself GS tiling pairs)

  At odd n: 0 blueself (THM-023). So blue_edges = 2^{k-1}.
  At n=3: 2^0 = 1. MATCH.
  At n=5: 2^3 = 8. Does NOT match 14.

  Hmm, 14 != 8. The blue edge count is NOT simply 2^{k-1}.
  The issue: multiple GS flip pairs can connect the SAME pair of classes.
  Each such connection counts as ONE edge but MULTIPLE GS flip pairs.

  Total GS flip pairs = 2^{k-1} = 8 at n=5.
  But there are 14 blue EDGES (distinct class pairs).
  That means some blue edges come from MULTIPLE GS pairs?
  No: 8 GS flip pairs can give AT MOST 8 distinct edges.
  14 > 8 is IMPOSSIBLE.

  RESOLUTION: The S211 script counts blue edges differently.
  It's counting ARC REVERSALS, not tiling flips.
  An arc reversal edge is "blue" if the reversed arc is a GS arc.
  This is a DIFFERENT definition from tiling flip!

  Let me reconcile:
  In the ARC REVERSAL model:
  - C(n,2) arcs total
  - Each arc reversal changes one arc
  - The arc is "blue" if it's a GS-symmetric arc (lies on the
    diagonal of the staircase, i.e., arc (v, n-1-v))
  - Otherwise "black"

  In the TILING FLIP model:
  - Flip = reverse ALL non-base arcs simultaneously
  - Blue = both tilings are GS
  - Black = at least one is non-GS

  These are COMPLETELY DIFFERENT! The S211 script uses arc reversals.
  The S181 script uses tiling flips.

  CLARIFICATION: What the user described IS the tiling flip model.
  The S211 data is the arc reversal model.
  These give different blue/black counts!
""")

print("  TILING FLIP blue edges (from S181 data):")
print("  n=3: 1 blue, 0 black (from 1 GS flip pair)")
print("  n=4: 2 blue, 2 black (from 4 GS flip pairs, 2 same-class)")
print("  n=5: 8 blue, ? black (from 8 GS flip pairs)")
print("  These are DIFFERENT from S211's arc-reversal blue/black counts.\n")

print("  AT ODD n: TILING FLIP blue edges = 2^{k-1}")
print("    n=3: 2^0 = 1. MATCH with S181.")
print("    n=5: 2^3 = 8. MATCH with S181.")
print("    n=7: 2^8 = 256 predicted tiling-blue edges.\n")

print("  AT EVEN n: TILING FLIP blue edges = 2^{k-1} - #blueself_pairs")
print("    n=4: 2^1 - 0 = 2. From S181: 2 cross-class blue + 1 self = 3 total pairs.")
print("    n=6: 2^5 - ? = ?. Need to compute.")

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: BLUE AND BLACK DEGREE")

print("""
CRITICAL DISTINCTION: TWO DIFFERENT DEFINITIONS OF "BLUE"

1. TILING FLIP BLUE: flip = bit complement of ALL non-base arcs.
   Blue iff both tilings are GS. This gives the user's picture:
   each tiling has ONE partner, connected by one blue or black line.

2. ARC REVERSAL BLUE: reverse ONE arc at a time.
   Blue iff the reversed arc is on the GS diagonal.
   This gives the G_n arc-reversal graph from S211.

These give DIFFERENT edge counts and degree sequences!

FOR THE TILING FLIP MODEL (user's picture):
  - At odd n: 2^{k-1} blue edges, 2^{m-1} - 2^{k-1} black edges.
    n=3: 1 blue, 0 black
    n=5: 8 blue, 24 black
    n=7: 256 blue, 16128 black (PREDICTED)
  - Blue fraction in TILING model: 2^{k-1} / 2^{m-1} = 2^{k-m}
    = the GS fraction! This DECREASES with n.
    n=3: 1/1 = 100%
    n=5: 8/32 = 25%
    n=7: 256/16384 = 1.6%

  So in the TILING model, blue edges become RARE at large n.

FOR THE ARC REVERSAL MODEL (S211):
  - Blue fraction INCREASES with n (47% at n=5, 69% at n=6).
  - This is because many arc reversals preserve the SC structure.

FOR n=7 TILING MODEL:
  Total tilings: 2^15 = 32768
  GS tilings: 2^9 = 512
  Tiling flip pairs: 16384 total
  Blue flip pairs: 256
  Black flip pairs: 16128
  Blue fraction: 1.6%

  Blue edges (distinct class pairs): <= 256
  SC classes: 88
  Blue degree per SC class: ~256/88 ~ 2.9 on average

  Black edges: much larger, connecting the 456 classes via
  the 16128 non-GS flip pairs.

THE TILING-LEVEL BLUE DEGREE at n=7:
  avg ~ 2.9 for SC classes (from 256 blue pairs among 88 SC classes)
  Each SC class has ~5.8 GS tilings, each with 1 blue partner.
  If partners spread evenly: blue_deg ~ 5.8.
  If concentrated: blue_deg ~ 1-3.

THE ARC-REVERSAL BLUE DEGREE at n=7:
  Will be MUCH higher because there are C(7,2)=21 arcs,
  and many arc reversals count as "blue" in the S211 sense.
""")

print("Done. opus-2026-03-22-S183")
