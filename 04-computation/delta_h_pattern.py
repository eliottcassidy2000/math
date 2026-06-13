#!/usr/bin/env python3
"""
delta_h_pattern.py — opus-2026-03-14-S71f

Pattern analysis of ΔH gaps under single arc flip.

n=4: ΔH ∈ {-4,-2,0,2,4} — NO gaps
n=5: ΔH ∈ {-12,-8,-6,-4,-2,0,2,4,6,8,12} — gap at ±10
n=6: ΔH ∈ {-32,...,-2,0,2,...,32} — gaps at ±28, ±30

Question: What is the pattern of ΔH gaps? Is it related to H-spectrum gaps?

Also: ΔH values are always even. What is the complete set of achievable ΔH/2?
"""

from collections import Counter

# ΔH/2 values from previous analysis:
# n=4: {0, 1, 2}
# n=5: {0, 1, 2, 3, 4, 6} — skips 5
# n=6: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16} — skips 14, 15

# Pattern of |ΔH/2| max values:
# n=4: max = 2, skipped = none
# n=5: max = 6, skipped = {5}
# n=6: max = 16, skipped = {14, 15}

# The maximum ΔH/2 for arc flip: each arc is in n-2 triples,
# flipping can change at most n-2 3-cycles, plus 5-cycles etc.

# At n=5: max |ΔH/2| = 6 = 3 (from t₃) + 3 (from d₅)
# At n=6: max |ΔH/2| = 16: how?

# Let's understand the n=6 maximum: ΔH = ±32
# ΔH/2 = 16 = Δα₁ (change in total directed odd cycle count)
# This means flipping one arc can create 16 new directed odd cycles!
# At n=6, cycles can be of length 3 or 5.
# The arc (i,j) appears in:
#   - n-2 = 4 triples (3-cycle candidates)
#   - C(n-2, 3) = 4 five-element sets containing both i,j (5-cycle candidates)
#     Wait: 5-element sets containing BOTH i and j from {0,...,5}:
#     choose 3 from {0,...,5}\{i,j} = C(4,3) = 4

# For 3-cycles: each of the 4 triples {i,j,k} can gain/lose 1 directed 3-cycle
# So max |Δt₃| = 4

# For 5-cycles: each of the 4 five-element sets can gain/lose multiple directed 5-cycles
# At n=5, one 5-element set can have 0,1,2,3 directed 5-cycles
# So max |Δd₅| per set = 3, times 4 sets = max |Δd₅| = 12

# Total: max |Δα₁| ≤ 4 + 12 = 16. This matches!

# The gap at ±28 and ±30:
# |ΔH/2| = 14 or 15 would need Δα₁ = 14 or 15
# With max Δt₃ = 4 and max Δd₅ = 12: Δα₁ ≤ 16
# But Δα₁ = 14 or 15 requires Δt₃ + Δd₅ = 14 or 15
# This needs extreme coordination: e.g., Δt₃=4, Δd₅=10 or 11

print("=" * 60)
print("ΔH Gap Analysis: Theoretical Bounds")
print("=" * 60)

for n in range(3, 9):
    max_dt3 = n - 2  # triples containing the flipped arc
    # 5-element sets containing both endpoints: C(n-2, 3)
    from math import comb
    five_sets = comb(n-2, 3) if n >= 5 else 0
    max_dd5 = 3 * five_sets  # each set can gain/lose up to 3 directed 5-cycles

    # 7-element sets containing both endpoints: C(n-2, 5)
    seven_sets = comb(n-2, 5) if n >= 7 else 0
    # Directed 7-cycles per set: up to 6!/7 = ~102 (too many to enumerate)
    # But the actual max change is bounded

    max_da1 = max_dt3 + max_dd5  # ignoring 7-cycles
    print(f"  n={n}: max|Δt₃|={max_dt3}, five_sets={five_sets}, max|Δd₅|≤{max_dd5}, max|Δα₁|≤{max_da1}")

# Now: why does |Δα₁|=5 skip at n=5?
print(f"\n{'='*60}")
print("Why |Δα₁|=5 skips at n=5")
print(f"{'='*60}")
print("""
At n=5:
  max|Δt₃| = 3 (3 triples contain the flipped arc)
  max|Δd₅| = 3 (1 five-element set, at most 3 directed 5-cycles change)
  max|Δα₁| = 6

  |Δα₁| = 5 requires |Δt₃| + |Δd₅| = 5 with SAME SIGN
  Possibilities: (3,2), (2,3)

  But from the joint distribution (computed in delta_h_gap.py):
    Only (Δt₃, Δd₅) pairs where signs match give high |Δα₁|:
    (-3,-3)→6, (-3,-1)→4, (-2,-2)→4, (-1,-2)→3, (1,2)→3, (2,2)→4, (3,1)→4, (3,3)→6

  Missing: (3,2) or (2,3) or (-3,-2) or (-2,-3)

  WHY? Because Δt₃=3 (all 3 triples gain a cycle) constrains the
  tournament structure so strongly that Δd₅ can only be 1 or 3, not 2.

  Similarly, Δd₅=3 means the 5-cycle gains 3 directed cycles,
  which constrains Δt₃ to be exactly 3 (forced by the same structure).

  The (3,2) and (2,3) combinations require conflicting structural constraints.
""")

# This is a beautiful LOCAL parity constraint!
# It means the 3-cycle and 5-cycle changes are CORRELATED.

# Does this correlation persist at n=6?
print(f"{'='*60}")
print("At n=6: gaps at |Δα₁| = 14, 15")
print(f"{'='*60}")
print("""
At n=6:
  max|Δt₃| = 4
  max|Δd₅| = 12 (4 five-sets × 3 max change each)
  max|Δα₁| = 16

  |ΔH/2| takes values {0,1,...,13,16} — skipping 14, 15.

  |Δα₁| = 14 would need (4,10) or (3,11) or (2,12) — extreme near-max 5-cycle changes
  |Δα₁| = 15 would need (3,12) or (4,11)
  |Δα₁| = 16 achieved: needs (4,12) — ALL triples and ALL five-sets max out

  So the gap is at |Δα₁| = max - 1 and max - 2.
  At n=5: gap at |Δα₁| = max - 1 (= 5)
  At n=6: gaps at |Δα₁| = max - 2 and max - 1 (= 14, 15)

  This suggests: the EXTREME simultaneous change can happen (max),
  but slightly-less-than-max requires contradictory constraints.
""")

# Predict gaps for n=7:
print("Prediction for n=7:")
n = 7
max_dt3 = n - 2  # 5
five_sets = comb(n-2, 3)  # 10
max_dd5 = 3 * five_sets  # 30
seven_sets = comb(n-2, 5)  # 1
# 7-cycle contributions: each set can have many directed 7-cycles
# At n=7, the full set {0,...,6} has many directed 7-cycles
# For now, ignore 7-cycles for the prediction

print(f"  max|Δt₃| = {max_dt3}")
print(f"  max|Δd₅| ≤ {max_dd5}")
print(f"  max|Δα₁| ≤ {max_dt3 + max_dd5}")
print(f"  Expected gaps near {max_dt3 + max_dd5 - 1} and {max_dt3 + max_dd5 - 2}")
print()
print("Note: 7-cycles add more, so actual max |Δα₁| is higher.")
print("The gap pattern should still appear near the max.")
