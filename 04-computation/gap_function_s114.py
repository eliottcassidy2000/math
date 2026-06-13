#!/usr/bin/env python3
"""
gap_function_s114.py — opus-2026-03-21-S114

THE GAP FUNCTION AND FORBIDDEN VALUES: WHY 7 AND 21?

Forbidden H values at n=7: 7, 21, 63, 107, 119, 149, 161-187.
The first three form a geometric sequence: 7, 3×7, 9×7 = 7×3^k.

CONNECTION TO THE OCR/NESTING:
The OCR residual comes from the POS, which nests tournaments recursively.
Forbidden values are structural constraints on I(Omega,2).
Can we explain the forbidden values from the nesting structure?

KEY INSIGHT FROM THM-029/079:
  H=7 impossible because i₁=3, i₂=0 forces c₅≥1 → α₁≥4 → H≥9
  H=21 impossible because I(Omega,2)=21 requires specific Omega structures
  that are unrealizable as conflict graphs of actual tournaments.

The connection: forbidden values are OCF OBSTRUCTIONS.
The independence polynomial I(Omega,2) cannot equal certain values
because tournament structure constrains Omega.

Author: opus-2026-03-21-S114
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import comb, factorial
import random
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

print("=" * 72)
print("  THE GAP FUNCTION AND FORBIDDEN VALUES")
print("  opus-2026-03-21-S114")
print("=" * 72)

# ========================================================================
# PART 1: Complete gap function at n=3,...,6 (exhaustive)
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: COMPLETE GAP FUNCTION")
print("=" * 72)

for n in range(3, 8):
    if n >= 7:
        # Use known data from the exhaustive computation
        # At n=7: gaps = [7, 21, 63, 107, 119, 149, 161-187]
        print(f"\n  n=7 (from prior exhaustive computation):")
        print(f"    Achievable: 77 values in [1, 189]")
        print(f"    Forbidden (low): 7, 21, 63")
        print(f"    Forbidden (high): 107, 119, 149, 161, 163, 165, 167, 169, 173, 177, 179, 181, 183, 185, 187")
        continue

    print(f"\n--- n = {n} ---")
    H_vals = set()
    for A in all_tournaments(n):
        H_vals.add(count_hp(A, n))

    H_sorted = sorted(H_vals)
    H_max = max(H_sorted)
    all_odd = set(range(1, H_max + 1, 2))
    gaps = sorted(all_odd - H_vals)

    print(f"  Achievable: {H_sorted}")
    print(f"  Gaps (forbidden odd): {gaps}")
    print(f"  H_max = {H_max}")

# ========================================================================
# PART 2: The 7, 21, 63 pattern — is it 7 × 3^k?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: THE 7 × 3^k PATTERN")
print("=" * 72)

forbidden_low = [7, 21, 63]
print(f"  First three structural gaps at n=7: {forbidden_low}")
print(f"  Ratios: {forbidden_low[1]}/{forbidden_low[0]} = {forbidden_low[1]/forbidden_low[0]:.1f}, "
      f"{forbidden_low[2]}/{forbidden_low[1]} = {forbidden_low[2]/forbidden_low[1]:.1f}")
print(f"  Pattern: 7 × 3^k for k=0,1,2")
print(f"  Predicted next: 7 × 3³ = 189 = H_max(n=7). Is 189 achievable? YES (Paley).")

print(f"""
  OBSERVATION: 7 × 3^k for k=0,1,2 gives 7, 21, 63.
  The NEXT value would be 189, which IS achievable (it's the maximum H at n=7).

  So the geometric sequence 7 × 3^k terminates at k=2 because 7 × 3³ = 189
  is the PALEY MAXIMUM.

  At n=5: forbidden = [7]. Gap sequence: just 7 = 7 × 3⁰.
  At n=6: gaps = [7, 21, 35, 39]. Is 35 = 7 × 5? Yes!
  And 39 = 3 × 13.

  Let me check n=6 gaps more carefully.
""")

# n=6 exhaustive
H6 = set()
for A in all_tournaments(6):
    H6.add(count_hp(A, 6))

H6_sorted = sorted(H6)
all_odd_6 = set(range(1, max(H6_sorted) + 1, 2))
# At n=6, H can be even too!
all_vals_6 = set(range(1, max(H6_sorted) + 1))
gaps_6 = sorted(all_vals_6 - H6)
# But actually H is always ODD at odd n. At n=6 (even n), H can be even.
print(f"  n=6: achievable H = {H6_sorted}")
print(f"  Gaps: {gaps_6}")
print(f"  Note: at even n, H can be even!")

# Wait — is H always odd? Rédei's theorem says H is always odd.
# So gaps should be among odd numbers only.
gaps_odd_6 = sorted(set(range(1, max(H6_sorted) + 1, 2)) - H6)
gaps_even_6 = sorted(set(range(2, max(H6_sorted) + 1, 2)) - H6)
print(f"  Odd gaps: {gaps_odd_6}")
print(f"  Even gaps: {gaps_even_6}")

# ========================================================================
# PART 3: The OCF decomposition of forbidden values
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: WHY 7 AND 21 ARE FORBIDDEN — THE OCF VIEW")
print("=" * 72)
print(f"""
  From OCF: H = I(Omega, 2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...

  H=7: Requires 1 + 2α₁ + 4α₂ + ... = 7.
    Only solution: α₁=3, α₂=0, α₃=0, ... (since 1+6=7, no room for 4α₂)
    This means: exactly 3 odd cycles, ALL pairwise sharing a vertex.
    THM-029 proves: 3 pairwise-sharing triples on 5+ vertices ALWAYS force
    a 5-cycle → α₁≥4 → H≥9.

  H=21: Requires 1 + 2α₁ + 4α₂ + ... = 21 → α₁ + 2α₂ + 4α₃ + ... = 10.
    Possible decompositions:
      (10, 0, 0): 10 cycles, all pairwise conflicting. 10 is VERY large.
      (8, 1, 0): 8 cycles, 1 disjoint pair.
      (6, 2, 0): 6 cycles, 2 disjoint pairs.
      (4, 3, 0): 4 cycles, 3 disjoint pairs.
      (2, 4, 0): 2 cycles, 4 disjoint pairs? Impossible (max α₂ with α₁=2 is 1).
    Plus decompositions with α₃ > 0.
    THM-079 eliminates all of these via tournament structure constraints.

  H=63: Requires α₁ + 2α₂ + 4α₃ + ... = 31.
    63 = 1 + 2×31. Many decompositions possible.
    Is 63 provably forbidden? Or just absent at n≤7?
""")

# Check if 63 is achievable at n=8
print(f"  Checking H=63 at n=8 (sampled)...")
rng = random.Random(42)
found_63 = False
for _ in range(100000):
    A = [[0]*8 for _ in range(8)]
    for i in range(8):
        for j in range(i+1, 8):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    H = count_hp(A, 8)
    if H == 63:
        found_63 = True
        break

print(f"  H=63 at n=8: {'FOUND' if found_63 else 'NOT FOUND in 100K samples'}")

# Check H=7 and H=21 at n=8
found_7 = False
found_21 = False
rng = random.Random(42)
for _ in range(100000):
    A = [[0]*8 for _ in range(8)]
    for i in range(8):
        for j in range(i+1, 8):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    H = count_hp(A, 8)
    if H == 7:
        found_7 = True
    if H == 21:
        found_21 = True

print(f"  H=7 at n=8: {'FOUND' if found_7 else 'NOT FOUND in 100K'}")
print(f"  H=21 at n=8: {'FOUND' if found_21 else 'NOT FOUND in 100K'}")

# ========================================================================
# PART 4: The gap function and the nesting
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: GAP FUNCTION AND THE RECURSIVE NESTING")
print("=" * 72)
print(f"""
  From S110: at the POS, H = 14 - H_inner + 4×sbs (at n=5).
  The achievable H values in the POS are {{11, 13, 15}}.

  The GAP at H=7 is OUTSIDE the POS — it can't come from the nesting formula.
  H=7 would require: 14 - H_inner + 4×sbs = 7
    If sbs=0: H_inner = 7. Is H_inner=7 achievable at n=3? NO! H(n=3) ∈ {{1, 3}}.
    If sbs=1: H_inner = 7 - 4 = 3. H_inner=3 IS achievable. But sbs=1 forces
    inner = 3-cycle, and H_inner = 3 for 3-cycle.
    So H = 14 - 3 + 4 = 15. Not 7!

  H=7 is forbidden because there's NO WAY to get it from the nesting:
  - The nesting formula at the POS gives H ∈ {{11, 13, 15}}
  - Non-POS score classes give different H ranges but still skip 7

  The gap at H=7 is a GLOBAL constraint, not a POS-specific one.
  It holds across ALL score classes, ALL nesting depths.

  From the OCF: H = 1 + 2α₁ + 4α₂ + ...
  H=7 requires α₁=3, α₂=0. But tournament structure forces α₁≥4
  whenever 3 pairwise-conflicting cycles exist.

  THIS IS A CONSTRAINT ON THE CONFLICT GRAPH Omega, not on the nesting.
  The nesting determines WHICH tournaments appear in the POS.
  The gap function determines WHICH H values are achievable GLOBALLY.

  The intersection: forbidden values that would fall in the POS range.
  At n=5: POS range [11,15], forbidden=7. No intersection.
  At n=7: POS range is wider. Does any forbidden value fall in the POS range?
""")

# ========================================================================
# PART 5: The forbidden values as OCF obstructions
# ========================================================================
print("=" * 72)
print("  PART 5: FORBIDDEN VALUES AS OCF OBSTRUCTIONS")
print("=" * 72)
print(f"""
  Each forbidden H value corresponds to an UNREALIZABLE OCF decomposition.

  H = 1 + 2α₁ + 4α₂ + 8α₃ + 16α₄ + ...

  The achievable (α₁, α₂, α₃, ...) tuples are constrained by:
  1. α_k is the count of independent k-sets in Omega(T)
  2. Omega(T) must be the conflict graph of SOME tournament
  3. The conflict graph structure restricts which (α₁,α₂,...) are possible

  For H=7: only (3,0,...) works. But 3 pairwise-conflicting odd cycles
           always create more cycles → α₁≥4. BLOCKED.

  For H=21: only (10,0), (8,1), (6,2), (4,3) work.
            Each requires specific Omega topology that can't be realized.
            THM-079's 464-line proof eliminates all cases.

  For H=63: 63 = 1 + 62 → α₁+2α₂+... = 31.
            Many decompositions. Is there a short obstruction?

  THE KEY: 7 = 1 + 2×3, and α₁=3 with α₂=0 is the UNIQUE minimal
  obstruction. For 21 = 1 + 2×10, the obstructions are more complex.

  PATTERN: The forbidden values 7, 21 are of the form 1 + 2k where
  k is a value that CANNOT be achieved as α₁+2α₂+4α₃+... with all the
  tournament structure constraints.

  7: k=3 forbidden because α₁=3,α₂=0 impossible
  21: k=10 forbidden because all decompositions blocked
  63: k=31 forbidden because... ?
""")

# ========================================================================
# PART 6: The connection between gaps and the OCR minimum at n=7
# ========================================================================
print("=" * 72)
print("  PART 6: GAPS AND THE OCR MINIMUM AT n=7")
print("=" * 72)
print(f"""
  At n=7: 18 forbidden values out of 95 odd values (18.9%).
  At n=5: 1 forbidden value out of 8 odd values (12.5%).
  At n=6: {len(gaps_odd_6)} forbidden values out of {max(H6_sorted)//2 + 1} odd values.

  n=5: 1/8 = 12.5% gaps
  n=6: {len(gaps_odd_6)}/{max(H6_sorted)//2+1} = {len(gaps_odd_6)/(max(H6_sorted)//2+1)*100:.1f}% gaps
  n=7: 18/95 = 18.9% gaps

  The gap DENSITY increases from n=5 to n=7!
  This is related to the OCR dip: more gaps → more structural constraints
  → more "forced" H values → less freedom → lower OCR?

  ACTUALLY: the relationship is more subtle. The OCR measures within-class
  variance, while gaps measure global achievability. But they're connected
  through the OCF:

  The POS contributes to the OCR residual when its score class allows
  multiple H values. The gap function says which H values CAN'T appear.
  If the POS range includes a gap, there are fewer achievable H values,
  potentially REDUCING the OCR residual.

  At n=7: the POS has wide H ranges (e.g., 16 distinct values in one class).
  The gaps remove some of these, but 77 achievable values remain.
  The OCR residual is high because the achievable values are SPREAD OUT.

  THE CONNECTION: the gap density and the OCR residual are both measures
  of how CONSTRAINED the H landscape is by tournament structure.
  High gap density → heavily constrained → the structure "wants" specific values.
  High OCR residual → scores don't capture these constraints.

  Both peak near n=7 because n=7 is where the tournament structure
  (conflict graphs, cycle forcing, matching constraints) is maximally
  complex relative to the number of possible values.
""")

# ========================================================================
# PART 7: What makes 7 special? The 7-ness of the forbidden value
# ========================================================================
print("=" * 72)
print("  PART 7: WHAT MAKES 7 SPECIAL?")
print("=" * 72)
print(f"""
  7 = 1 + 2 × 3. The number 3 is special because:
  1. 3 is the MINIMUM number of pairwise-conflicting odd cycles that
     creates a forcing argument.
  2. With 1 or 2 cycles: no constraint (1+2=3 ✓, 1+4=5 ✓).
  3. With 3 cycles, all pairwise conflicting: the overlap structure
     FORCES additional cycles (THM-029).

  21 = 1 + 2 × 10 = 3 × 7. The number 10 is special because:
  1. 10 is the minimum α₁+2α₂ value that requires COMPLEX Omega topology.
  2. All four viable decompositions (10,0),(8,1),(6,2),(4,3) are blocked.
  3. The 464-line proof of THM-079 shows each decomposition leads to
     a tournament structure contradiction.

  63 = 1 + 2 × 31 = 9 × 7. Is 31 an obstruction like 3 and 10?
  31 in binary: 11111. It's 2⁵-1.
  The OCF formula: 31 = α₁ + 2α₂ + 4α₃ + 8α₄ + 16α₅.
  With α₅≥1: 16+... ≤ 31, so α₅=1, remaining 15 = α₁+2α₂+4α₃+8α₄.
  This gives MANY decompositions. Harder to block.

  THE PATTERN 7, 21, 63:
  7 = 2³ - 1
  21 = ? Not 2^k - 1.
  63 = 2⁶ - 1

  Actually: 7=7, 21=3×7, 63=9×7. So the pattern is 7 × {{1, 3, 9}} = 7 × 3^k.

  WHY 7 × 3^k?
  - H=7: needs 3 all-conflicting cycles → forces extra cycle → H≥9
  - H=21=3×7: needs I(component,2)=7 in disconnected Omega → impossible
    (same argument as H=7 applied to a component)
  - H=63=9×7: needs I(component,2)=7 or 21 in some factorization?
    63 = 7×9 or 3×21 or 1×63.
    If Omega has components with I=7: impossible (THM-029).
    If Omega has components with I=21: impossible (THM-079).
    If Omega connected with I=63: more complex.

  THE RECURSIVE STRUCTURE:
  H=7 is forbidden because I=7 is impossible for ANY Omega.
  H=21=3×7 is forbidden because BOTH factors (3,7) create problems:
    factor 7 → impossible component; factor 3 → forces more structure → contradiction.
  H=63=9×7: factor 7 → impossible; factor 9 → needs 4 pairwise-conflicting cycles?

  THE MULTIPLICATION BY 3 IS THE OCF PRODUCT:
  I(G₁∪G₂, 2) = I(G₁, 2) × I(G₂, 2).
  If Omega has 2 components, H = I₁ × I₂.
  H=21 = 3×7: one component has I=3 (one cycle), other has I=7 (impossible).
  H=63 = 9×7: one component has I=9, other has I=7 (impossible).
  Also: H=63 = 3×21: one has I=3, other has I=21 (impossible).

  SO: 7 IS THE ATOMIC FORBIDDEN VALUE.
  21 = 3 × 7 inherits the impossibility via the product structure.
  63 = 9 × 7 also inherits it. As does 63 = 3 × 21.

  ANY MULTIPLE OF 7 THAT FACTORS AS k × 7 WHERE k IS ACHIEVABLE
  BY A SINGLE OMEGA COMPONENT IS FORBIDDEN (via disconnected Omega).
  But for CONNECTED Omega, 63 might still be achievable...

  From the n=7 data: H=63 is absent. Is it provably forbidden?
""")

# Check which multiples of 7 are forbidden
print(f"  Multiples of 7 up to 189:")
achievable_n7 = set([1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,35,37,39,41,43,
                     45,47,49,51,53,55,57,59,61,65,67,69,71,73,75,77,79,81,83,
                     85,87,89,91,93,95,97,99,101,103,105,109,111,113,115,117,
                     121,123,125,127,129,131,133,135,137,139,141,143,145,147,
                     151,153,155,157,159,171,175,189])

for k in range(1, 28):
    val = 7 * k
    if val > 189:
        break
    if val % 2 == 0:
        continue
    status = "ACHIEVABLE" if val in achievable_n7 else "FORBIDDEN"
    print(f"  7 × {k} = {val}: {status}")

# ========================================================================
# PART 8: The gap function as an OCR constraint
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 8: THE GAP FUNCTION CONSTRAINS THE OCR")
print("=" * 72)
print(f"""
  SYNTHESIS: How forbidden values connect to the OCR.

  1. The OCR residual = Var(H | scores) / Var(H).
  2. Var(H | scores) depends on how SPREAD OUT H is within each score class.
  3. Forbidden values CREATE GAPS in the H distribution within score classes.
  4. These gaps REDUCE the possible spread → HELP the OCR.

  EXAMPLE at n=5:
  POS class (1,2,2,2,3) has H ∈ {{11, 13, 15}}.
  The gap at H=7 is BELOW this range, so it doesn't directly help.
  But it means no tournament with this score has H=7, limiting ambiguity.

  AT n=7: more gaps exist in the MIDDLE of the H range.
  H=63 is forbidden in the range where many score classes have H values.
  This means some expected H values are "missing," constraining the spread.

  THE PARADOX: More gaps (n=7) correspond to LOWER OCR (n=7 minimum).
  This seems contradictory — gaps should reduce spread.

  RESOLUTION: The gaps at n=7 are STRUCTURAL (forced by Omega topology),
  while the OCR residual comes from STATISTICAL variation within score
  classes. The gaps remove specific values but don't reduce the overall
  VARIANCE, because the remaining values are still spread out.

  In fact: the gaps may INCREASE variance by creating a multimodal
  distribution within score classes (values cluster around the gap,
  with nothing in between).

  THE PICTURE:
  - Few gaps (n=3,4): H is tightly constrained → OCR = 1
  - More gaps (n=5-7): gaps create structure but don't reduce variance → OCR dips
  - Many gaps near H_max (n=8+): gaps compress the high-H tail → OCR recovers

  The gap density near the MODE of the H distribution matters most.
  At n=7: gaps at 7, 21, 63 are in the LEFT TAIL (few tournaments there).
  At larger n: gaps would need to be in the BULK to significantly affect OCR.
""")

print("\nDone. opus-2026-03-21-S114")
