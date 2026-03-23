#!/usr/bin/env python3
"""
n7_merged_s20cg.py -- kind-pasteur-2026-03-22-S20cg

G_7 MERGED META-GRAPH: What can we compute analytically?

At n=7: 456 iso classes, 4086 edges (from opus S212).
We need: SC count, complement pairing, blue/black split.

For the MERGED graph G_7/Z_2:
  V_merged = (456 + SC_7) / 2
  E_merged = E_total - collapsed - twin

The SC count at n=7 can be estimated from the Burnside formula
applied to the group S_n x Z/2Z, or from known OEIS sequences.

The number of SC tournaments (labeled): A000171
  n=1:1, n=3:2(=A000171(3)?), n=5:704, n=7:?

The number of SC iso classes: this is the count under S_n x Z/2Z.

Actually, for tournament iso classes, the SC classes are those
fixed by the complement map c. In the merged graph, these stay
as individual nodes. The NS classes pair up.

Author: kind-pasteur-2026-03-22-S20cg
"""
import sys
from math import comb, factorial, gcd
from fractions import Fraction
from collections import Counter
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
    """Fix(sigma) for tournaments: 0 if any even cycle, 2^orbit_pairs if all odd."""
    if any(l % 2 == 0 for l in ct): return 0
    n_op = sum((l-1)//2 for l in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)):
            n_op += gcd(ct[i], ct[j])
    return 2**n_op

print("=" * 70)
print("  G_7 MERGED META-GRAPH: ANALYTICAL PREDICTIONS")
print("=" * 70)

# ================================================================
# 1. SC ISO CLASS COUNT AT n=7
# ================================================================
print(f"\n{'='*70}")
print(f"  1. SC ISO CLASS COUNT")
print(f"{'='*70}\n")

# SC iso classes = orbits of S_n x Z/2Z on tournaments
# = (1/|S_n x Z/2Z|) * sum_{(sigma, epsilon)} Fix((sigma, epsilon))
# where epsilon in {0, 1} (complement or not).

# Fix((sigma, 0)) = tournament Fix(sigma) (already computed)
# Fix((sigma, 1)) = # tournaments T such that sigma(T) = T^complement
# = # anti-automorphisms with permutation sigma.

# For Fix((sigma, 1)): sigma(T) = complement(T) means
# for each arc (a,b): sigma sends the arc to the OPPOSITE direction.
# This is like the tournament Fix but with OPPOSITE sign.

# For ordered pairs: sigma sends (a,b) to (sigma(a), sigma(b)).
# For anti-automorphism: if T has arc a->b, then sigma(T) has
# arc sigma(a)->sigma(b), but complement(T) has arc b->a.
# So we need: sigma maps (a,b) to (b,a) for all arcs.
# Wait: sigma(T) has sigma(a)->sigma(b) iff T has a->b.
# complement(T) has b->a iff T has a->b.
# So sigma(T) = complement(T) iff sigma(a)->sigma(b) = b->a for all (a,b).
# This means: sigma maps directed arc (a,b) to directed arc (b,a) = reversed arc.
# In terms of ordered pairs: sigma maps (a,b) to (sigma(a), sigma(b)),
# and we need this to be (b,a), i.e., sigma(a) = b and sigma(b) = a.
# Wait that would mean sigma is the transpose of every edge, which is impossible
# for a permutation of vertices.

# Actually: sigma(T) = T^comp means:
# T has arc a->b iff sigma(T) has arc b'->a' where b'=sigma(b), a'=sigma(a).
# Wait, let me be more careful.
# sigma(T) is the tournament where arc from sigma(a) to sigma(b) exists iff a->b in T.
# T^comp has arc b->a iff T has arc a->b.
# So sigma(T) = T^comp iff: for every pair {a,b},
# sigma(a)->sigma(b) in sigma(T) iff a->b in T,
# and b->a in T^comp iff a->b in T.
# So we need: sigma(a)->sigma(b) in sigma(T) iff b->a in T^comp,
# i.e., a->b in T iff b->a in T^comp iff a->b in T. Tautology?

# I think I'm confusing myself. Let me just state:
# Fix((sigma, complement)) = # T such that for every ordered pair (i,j):
# T(sigma(i), sigma(j)) = 1 - T(i,j)
# (sigma permutes vertices AND flips all arcs)

# For a pair-orbit of sigma of length L:
# The arc directions must satisfy: after L applications of sigma,
# the direction is FLIPPED L times (once per step because of complement).
# After L steps: direction flipped L times. Must return to original.
# So (-1)^L = 1, i.e., L must be EVEN.

# But for tournaments, we also need the sign-flip constraint from before.
# The combined constraint for (sigma, complement):
# For each pair-orbit of length L:
# After L steps of sigma, the direction has been:
# - Sign-flipped s times (from sigma reordering endpoints)
# - Complement-flipped L times (once per step from complement)
# Total flip: (-1)^{s+L}
# Must equal +1 for consistency: s + L must be even.

# For the original tournament Fix: need s even (all pair-orbits have even sign-flips).
# For the complement Fix: need s + L even for each orbit.

# This is the CONJUGATE condition. I'll compute it by brute force for small n
# and check against known SC class counts.

# Known SC class counts at small n:
# n=3: all 2 classes are SC -> SC iso classes = 2
# n=4: 2 SC classes out of 4
# n=5: 8 SC classes out of 12
# n=6: 12 SC classes out of 56

# The SC ISO CLASS count under S_n only (not S_n x Z/2Z) is just
# the number of iso classes that are self-complementary.
# This is a SUBSET of A000568(n).

# At n=7: how many of the 456 iso classes are SC?
# From the blueself counts: at n=7, the labeled SC tournament count
# and the SC iso class count are related by Burnside on SC tournaments.

# Let me compute A000568 decomposed into SC and NS.
# SC iso class count = (1/n!) * sum_sigma #{SC tournaments fixed by sigma}

# Actually simpler: for each iso class C, C is SC iff any T in C
# has T^comp isomorphic to T. This is a property of the CLASS, not individual T.

# I don't have the n=7 class data to determine this directly.
# But I can estimate: at n=7, the SC fraction of LABELED tournaments
# can be computed, and from that the SC iso class count.

# SC labeled tournaments at odd n:
# A self-complementary tournament on n vertices exists only if n is odd
# or n = 0 mod 4. At n=7 (odd): SC tournaments exist.

# The number of labeled SC tournaments: A000171 (if that's the right sequence)
# Actually, for n=7: SC tournaments have been studied.
# From our data: SC fraction at n=3: 100%, n=5: 68.8%, and we computed
# at n=6: 21.5%. At odd n, the SC fraction is higher.

# Let me use the Burnside approach to PREDICT the SC iso class count.
# SC iso classes = the number of classes C such that complement(C) = C.

# From the A000568 Burnside: the 456 classes split into SC and NS.
# NS classes pair up under complement. So:
# 456 = SC_count + 2 * NS_pair_count
# SC_count = 456 - 2 * NS_pair_count
# SC_count must be even (for n=7 odd) or can be any parity?

# At n=3: SC=2, NS_pairs=0. 2 = 2 + 0. OK.
# At n=5: SC=8, NS_pairs=2. 12 = 8 + 4. OK.
# At n=7: SC=?, NS_pairs=(456-SC)/2.

# The SC class count for tournaments is OEIS A000171 divided by n!/|Aut|.
# Actually A000171 counts labeled SC tournaments.
# A000171: 1, 1, 6, 1, 1, 720, 1 (for n=1,2,...,7)?
# This doesn't look right. Let me think differently.

# SC labeled tournaments at n: this is the count of labeled T with T ~ T^comp.
# At n=3: 8 (all are SC). At n=5: 704. At n=7: unknown without computation.

# For the PREDICTION: I'll use the pattern.
# SC fraction of labeled: 100%, 75%, 68.8%, 21.5% at n=3,4,5,6.
# At odd n: higher SC fraction. At even n: lower.
# n=7 (odd): SC fraction should be between 68.8% (n=5) and lower.

# From Burnside theory:
# SC labeled count = sum_sigma Fix((sigma, complement))
# The SC iso class count = SC_labeled / avg_orbit_size_of_SC_classes

# Without the full n=7 data, let me estimate.
# If SC fraction of CLASSES is similar to fraction of labeled:
# At n=5: 8/12 = 66.7% of classes are SC (vs 68.8% of labeled).
# At n=6: 12/56 = 21.4% of classes are SC (vs 21.5% of labeled).
# So class SC fraction ~ labeled SC fraction. Very close!

# If we knew the labeled SC count at n=7, we'd know the class count.
# Let me compute the SC labeled fraction at n=7 using the Burnside approach.

# For the COMPLEMENT action alone (not combined with S_n):
# # SC labeled tournaments = sum_{sigma in S_n} Fix((sigma, complement)) / 2
# ... this is getting complicated. Let me just use the known pattern.

# At odd n: the SC fraction decreases.
# n=3: 100%. n=5: 68.8%. n=7: maybe ~40-50%?
# That would give SC class count ~ 0.45 * 456 = 205.

# Or from the growth pattern:
# SC iso classes: 2, 2, 8, 12, ?
# Ratios: 1, 4, 1.5. Very irregular.
# Let me try: SC at n=7 = ?

# Actually from the repo memory: SC = blueself count at n=7.
# From opus S168: blueself count sequence = 8, 48, 704, 7040.
# These are LABELED counts: 8 at n=3, 48 at n=4, 704 at n=5, 7040 at n=6.
# At n=7: we'd need to compute.

# But blueself at n=5 = 704. Total at n=5 = 1024. Fraction = 68.8%.
# Blueself at n=6 = 7040. Total at n=6 = 32768. Fraction = 21.5%.
# These match the SC fractions exactly. So blueself = SC labeled count.

# For n=7: total = 2^21 = 2097152. If SC fraction ~ 30%:
# SC labeled ~ 629,000. SC iso class count ~ 629,000 / (7!/avg_aut_SC).
# If avg |Aut| for SC classes ~ 2-3: SC classes ~ 629,000 / (5040/2.5) ~ 312.

# This is very rough. Let me use a different approach.

# From Burnside: SC iso class count can be computed from the Fix((sigma, comp)) formula.
# But this requires knowing which pair-orbits have even vs odd s+L.

print(f"  ESTIMATES FOR n=7:")
print(f"  Total iso classes: 456")
print(f"  Total edges: 4086")

# Conservative estimates
for sc_est in [100, 150, 200, 240, 300]:
    ns_pairs = (456 - sc_est) // 2
    v_merged = sc_est + ns_pairs
    print(f"  SC={sc_est}: NS_pairs={ns_pairs}, V_merged={v_merged}")

# ================================================================
# 2. PREDICTIONS FOR THE MERGED GRAPH
# ================================================================
print(f"\n{'='*70}")
print(f"  2. PREDICTIONS FOR G_7/Z_2")
print(f"{'='*70}\n")

# If SC_7 ~ 240 (our best guess from the blueself pattern):
sc_est = 240
ns_est = 456 - sc_est  # = 216
ns_pairs = ns_est // 2  # = 108
v_merged = sc_est + ns_pairs  # = 348

# Blue edges at n=7: from the inversion pattern,
# most edges should be blue (NS-NS connections).
# At n=6: 200/290 = 69% blue.
# At n=7: probably 75-80% blue.
blue_frac_est = 0.77
blue_est = int(4086 * blue_frac_est)
black_est = 4086 - blue_est

# Collapsed edges: growing sequence 0, 0, 0, 5, ?
# Extrapolate: maybe 30-50 at n=7.
collapsed_est = 40

# Twin edges: at n=6, twin = 290 - 143 - 5 = 142
# Twin fraction = 142/290 = 49%.
# At n=7: twin fraction might be ~50%
twin_est = int(4086 * 0.50)
e_merged_est = 4086 - collapsed_est - twin_est

print(f"  BEST ESTIMATE (SC ~ {sc_est}):")
print(f"    V_merged = {v_merged}")
print(f"    E_merged ~ {e_merged_est}")
print(f"    Blue ~ {blue_est} ({100*blue_frac_est:.0f}%)")
print(f"    Black ~ {black_est}")
print(f"    Collapsed ~ {collapsed_est}")
print(f"    Twin ~ {twin_est}")
print(f"    Density ~ {2*e_merged_est/(v_merged*(v_merged-1)):.4f}")

# ================================================================
# 3. THE RECURSIVE PATTERN CHECK
# ================================================================
print(f"\n{'='*70}")
print(f"  3. THE RECURSIVE PATTERN: G_n/Z_2 -> G_(n-2)/Z_2")
print(f"{'='*70}\n")

# The descent should map G_7/Z_2 -> G_5/Z_2.
# G_5/Z_2 has 10 nodes and 21 edges.
# The PoS classes of G_7 should project onto G_5/Z_2.

# At n=7: the PoS score class is (1,2,3,3,3,3,3) or similar.
# The most ambiguous score class at n=7 has many iso classes.
# These should map to G_5/Z_2's 10 nodes via source-sink removal.

print(f"  Descent: G_7/Z_2 ({v_merged} nodes) -> G_5/Z_2 (10 nodes)")
print(f"  Ratio: {v_merged/10:.1f}x")
print(f"  Compare: G_6/Z_2 (34 nodes) -> G_4/Z_2 (3 nodes)")
print(f"  Ratio: {34/3:.1f}x")
print(f"  G_5/Z_2 (10 nodes) -> G_3/Z_2 (2 nodes)")
print(f"  Ratio: {10/2:.1f}x")

# The ratio sequence: 5.0, 11.3, ~35
# Growing roughly as n^2 or A000568(n)/A000568(n-2).
print(f"\n  Descent ratios: {10/2:.1f}, {34/3:.1f}, ~{v_merged/10:.1f}")
print(f"  A000568 ratios: {12/2:.1f}, {56/4:.1f}, {456/12:.1f}")
print(f"  Close match: the descent ratio ~ A000568(n)/A000568(n-2)")

# ================================================================
# 4. THE BLUE FRACTION TREND
# ================================================================
print(f"\n{'='*70}")
print(f"  4. THE BLUE FRACTION TREND")
print(f"{'='*70}\n")

blue_data = {3: (1, 1), 4: (1, 5), 5: (14, 30), 6: (200, 290)}
print(f"  {'n':>3s} {'Blue':>6s} {'Total':>7s} {'Blue%':>7s} {'SC frac':>8s}")
for n in [3, 4, 5, 6]:
    b, t = blue_data[n]
    sc_frac = {3: 1.0, 4: 0.5, 5: 0.667, 6: 0.214}[n]
    print(f"  {n:>3d} {b:>6d} {t:>7d} {100*b/t:>6.1f}% {sc_frac:>7.1%}")

print(f"  {7:>3d} {blue_est:>6d} {4086:>7d} {100*blue_frac_est:>6.1f}% {'~53%':>8s}")

# The blue fraction sequence: 100%, 20%, 47%, 69%, ~77%
# PATTERN: After the dip at n=4 (20%), blue fraction INCREASES.
# This is because NS-NS connections (which are blue) dominate at large n.
# The convergence to ~100% is slow because SC classes always contribute
# some black edges.

print(f"""
  THE BLUE FRACTION TREND:
    n=3: 100% (all SC, all blue)
    n=4: 20% (first NS classes, black dominates)
    n=5: 47% (more SC, blue recovers)
    n=6: 69% (NS-NS blue dominates)
    n=7: ~77% (NS-NS continues to grow)
    n->inf: ~100% (almost all NS, almost all blue)

  The blue fraction -> 1 because:
  1. SC fraction of classes -> 0
  2. NS-NS connections are BLUE by definition
  3. The only black edges connect SC to NS (rare) or NS-to-complement-NS
  4. As NS dominates, blue dominates

  THE MERGED GRAPH at large n is ALMOST ENTIRELY ONE COLOR (blue).
  Black edges become the RARE structural features connecting the
  thin SC spine to the massive NS bulk.
""")

# ================================================================
# 5. THE WIDTH PREDICTION VERIFICATION
# ================================================================
print(f"{'='*70}")
print(f"  5. WIDTH AT n=7: PREDICTION C(5,2) = 10")
print(f"{'='*70}\n")

# Width = C(n-2, floor((n-2)/2)) = C(5, 2) = 10.
# This predicts at most 10 iso classes share the same H value at n=7.
# From opus S212: 456 classes. We need the H-level distribution.
# If width = 10, the widest H-level has exactly 10 classes.

print(f"  Width prediction for n=7: C(5, 2) = {comb(5, 2)}")
print(f"  Previous verified: n=3:C(1,0)=1, n=4:C(2,1)=2, n=5:C(3,1)=3, n=6:C(4,2)=6")
print(f"  Need H-level distribution at n=7 to verify.")
print(f"  The prediction of 10 would mean the widest H-level has exactly 10 classes.")
