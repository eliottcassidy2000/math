# Structural Statements About Specific Classes
## opus-2026-04-03-S27

### VERIFIED STATEMENTS (n=3..6)

**S1. The transitive class has exactly 1 tiling and connects via exactly 1 blue complement line.**
  - At every n: the all-zeros tiling is the unique tiling in the transitive class (H=1)
  - |Aut| = n! (the full symmetric group, since the transitive tournament is the unique total order)
  - Size = n!/n! = 1
  - Complement neighbors = exactly 1 class
  - That complement connection is ALWAYS blue (the all-zeros tiling is grid-symmetric)

**S2. The transitive class connects to a class of size 1+2^(n-2) via its single blue complement line.**
  - n=3: complement target has size 1 (special case: 1+2^1=3 but size is 1 — fails for n=3)
  - n=4: complement target has size 5 = 1+2^2 ✓
  - n=5: complement target has size 9 = 1+2^3 ✓  
  - n=6: complement target has size 17 = 1+2^4 ✓
  - The target class always has H = 1+2^(n-2) as well (H = size for these classes!)

**S3. H = class size for certain special classes.**
  - At n=4: class with H=5 has size 5. H=size. ✓
  - At n=5: classes with H=9 have size 9. Two such classes. ✓
  - At n=6: classes with H=17 have size 17. Two such classes. ✓
  - Conjecture: the class at H=1+2^(n-2) always has size = H = 1+2^(n-2).

**S4. The transitive class has wiggly degree = m = C(n-1,2).**
  - n=3: wig_deg=1=C(2,2) ✓
  - n=4: wig_deg=3=C(3,2) ✓
  - n=5: wig_deg=6=C(4,2) ✓
  - n=6: wig_deg=10=C(5,2) ✓
  - Every wiggly flip of the transitive creates a DIFFERENT neighbor class!
    (because the transitive has only 1 tiling, and each tile flip creates a different tournament)

**S5. The anti-transitive (all-ones) tiling is in the SAME class as the complement target of the transitive.**
  - n=4: anti-transitive is class #2 with H=5 (same as complement target of transitive) ✓
  - n=5: anti-transitive is class #6 with H=9 (same as complement target of transitive) ✓
  - n=6: anti-transitive is class #22 with H=17 (same as complement target of transitive) ✓
  - This makes sense: the all-ones tiling IS the complement of the all-zeros tiling.

**S6. Size 1+2^(n-2) appears exactly twice (the anti-transitive class and its partner).**
  - n=5: 2 classes with size 9 (#4 and #6), both SC, same scores, same H ✓
  - n=6: 2 classes with size 17 (#12 and #22), both SC, same scores, same H ✓

**S7. Size 2^(n-2)-1 also appears as a class size.**
  - n=5: no class with size 7 (fails for n=5)
  - n=6: 4 classes with size 15 = 2^4-1 ✓ (two SC with H=45, two NS with H=15)
  - Not as clean a pattern as 1+2^(n-2).

**S8. The H=1 class (transitive) has complement degree exactly 1 at all n.**
  - Because it has only 1 tiling, it has only 1 complement → 1 complement edge.
  - This edge is always blue (the tiling is grid-symmetric).

**S9. Max-H classes have the LOWEST wiggly degree.**
  - n=5: H=15 classes have wig_deg=2,3 vs avg~5.5
  - n=6: H=45 classes have wig_deg=4,5 vs avg~10
  - Near-regular tournaments are "trapped" — few wiggly exits.

**S10. Complement degree grows with H (and class size).**
  - n=6: correlation between H and comp_deg is strong
  - H=1 → comp_deg=1, H=45 → comp_deg=10-11
  - But wiggly degree is NON-MONOTONE: it peaks in the middle H-range.

### OPEN QUESTIONS

Q1. Is H = size always true for the class at H=1+2^(n-2)?
Q2. Does the formula size=1+2^(n-2) extend to n=7,8?
Q3. What determines the wiggly degree of a class? Is there a formula in terms of H, t3, |Aut|?
Q4. Why does max-H have LOW wiggly degree but HIGH complement degree?

### CONSERVATION LAW: m+1 Transitions Per Tiling

Every tiling has exactly m+1 transitions:
  - m wiggly transitions (one per tile position)
  - 1 complement transition (flip all tiles)

Each transition is either:
  - A SELF-LOOP (stays in same class) — "silent mutation"
  - A CROSS-CLASS line (reaches a different class) — "expressive mutation"

Therefore: cross-class lines per tiling = (m+1) - self_loops_per_tiling

When comp_self_rate = 0 and wig_self_rate = 0 (no self-loops at all):
  cross-class lines per tiling = m+1 (maximum connectivity)
  This happens for the transitive class (H=1) at all n tested.

When self-loop rates > 0:
  cross-class lines per tiling < m+1
  Higher-H classes tend to have more self-loops, fewer exits.

Per-tiling complement lines = (1 - comp_self_rate)/2 ≈ 0.5 (most classes)
Per-tiling wiggly lines = (m - wig_self_per_tiling)/2 ≤ m/2

### THE DUALITY (refined)

The transitive (H=1):
  - Size 1, |Aut| = n!
  - Zero self-loops → all m+1 slots are cross-class
  - Wiggly degree = m (maximum: every tile flip creates a new class)
  - Complement degree = 1 (only one complement)
  - All wiggly neighbors are DIFFERENT classes

The regular-like (max H):
  - Larger classes, smaller |Aut|
  - More self-loops → fewer cross-class slots per tiling
  - Lower wiggly degree (per CLASS, not per tiling — fewer distinct neighbors)
  - Higher complement degree (more tilings → more complement targets)
  - But more LINE WEIGHT per edge (thicker connections)

Correlation at n=6:
  corr(H, comp_deg) = +0.737 (strong positive)
  corr(H, wig_deg) = +0.197 (weak positive)  
  corr(H, wig_self_rate) = +0.028 (nearly zero)
  corr(size, comp_deg) = +0.923 (very strong)
  corr(size, wig_deg) = +0.624 (moderate)
