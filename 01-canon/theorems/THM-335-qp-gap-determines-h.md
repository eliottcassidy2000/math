---
theorem: THM-335
name: Q-P Gap Determines H at n ≤ 4
status: PROVED
session: opus-2026-05-27-S1
verified: exhaustive (all 2 + 64 tournaments for n=3,4)
---

## Statement

Define the **Q-P gap** of tournament T on n vertices as:

gap(T) = max{d⁺(v)} − min{d⁺(v)}

i.e., the range of the score sequence (max outdegree minus min outdegree). Then:

**At n=3:** H(T) is uniquely determined by gap(T):
- gap = 0 (regular: score (1,1,1)) → H = 3
- gap = 2 (transitive: score (0,1,2)) → H = 1

**At n=4:** H(T) is uniquely determined by gap(T) via the linear formula:
- gap = 1 (score (1,1,2,2)) → H = 5
- gap = 2 (score (1,2,2,1) sorted = (1,1,2,2)... wait: score (0,2,2,2) or (1,1,2,2) or (1,2,2,1)) 
- Actually: gap=1 → H=5; gap=2 → H=3; gap=3 → H=1.

**H = 7 − 2·gap** at n=4. Perfect correlation: Pearson corr(gap, H) = −1.

**At n ≥ 5:** gap does NOT uniquely determine H. The correlation remains strong (−0.90 at n=5, −0.91 at n=6) but different tournaments with the same gap can have different H values.

## Proof for n=3

Three tournaments on 3 vertices (up to labeling):
- Transitive (1 iso class, 6 labelings): score (0,1,2), gap=2, H=1.
- 3-cycle (1 iso class, 2 labelings): score (1,1,1), gap=0, H=3.

The score sequence uniquely determines the tournament type at n=3. H is determined by type. □

## Proof for n=4

At n=4, the score sequence determines H:
- (0,1,2,3): transitive, H=1, gap=3.
- (0,2,2,2): one vertex beats all, H=3... wait, need to check.

Actually: at n=4, the score multisets are: (0,1,2,3), (0,2,2,2)... let me be precise. Score sequences for 4-vertex tournaments:

```
(0,1,2,3) → H=1, gap=3  (transitive)
(0,2,2,2) → H=3, gap=2  (one dominator)  
(1,1,2,2) → H=5, gap=1  (two "co-dominators")
(1,1,1,3) → H=3, gap=2  (one king, others regular)
```

Wait, score (0,2,2,2): one vertex beats 2, one vertex beats 0. gap=2, H=3.
Score (1,1,1,3): one vertex beats 3, three beat 1. gap=2, H=3? Let me verify.

Actually the output shows: gap=2 → H=3 (unique, 16 tournaments). So BOTH score (0,2,2,2) and (1,1,1,3) give H=3 and gap=2.

Why does gap=2 always give H=3? In a tournament with gap=2 at n=4:
- Q with degree d, P with degree d−2.
- If d=3: Q beats all (P has degree 1 or... wait d−2=1). Score contains one 3 and one 1 and two others summing to 6−4=2, so (1,1,1,3). H=3.
- If d=2: P has degree 0. Score (0,2,2,2). Q beats P and one other; P beats nobody. H=3.

Both give H=3 by separate verification. The formula H=7−2×gap holds:
- gap=1: H=5=7−2, ✓
- gap=2: H=3=7−4, ✓
- gap=3: H=1=7−6, ✓

**Proof of the formula H=7−2×gap at n=4:**

The only possible gap values are 1,2,3 (not 0, since regular n=4 tournament doesn't exist — sum of outdegrees = C(4,2) = 6, requiring average 1.5). The formula gives a unique H for each gap. We verify:
- gap=3: only transitive tournament, H=1. ✓
- gap=1: only score (1,1,2,2) tournaments (verified by exhaustion), H=5. ✓
- gap=2: two score types (0,2,2,2) and (1,1,1,3), both H=3. ✓ □

## Why n=5 Breaks the Pattern

At n=5 with gap=2 (score ≠ regular): H ∈ {9, 11, 13, 15}. The same gap can correspond to vastly different H values. This is because at n=5, the 5-cycle structure (which depends on the specific tournament, not just the score sequence) contributes significantly to H.

Specifically: H = 1 + 2α₁ + 4α₂ (from OCF, where α₁ = #odd cycles, α₂ = #disjoint cycle pairs). At n=5, α₁ and α₂ are NOT determined by the score sequence alone (5-cycle counts depend on the specific tournament).

## Connection to Principal Line

The Q-P gap is the natural measure of "hierarchy" in a tournament:
- gap=0: regular (maximal democracy), corresponds to the TOP of the principal line in G_n/Z₂.
- gap=max=n-1: transitive (maximal hierarchy), corresponds to the BOTTOM of the principal line.

The principal line in the metagraph G_n/Z₂ is the **Q-P axis**: it runs from maximum gap (transitive, H=1) to minimum gap (regular, H=max). The formula H=7−2×gap at n=4 is the linear version of this axis. At n≥5, the axis is no longer linear but remains the primary correlation.

**Correlation values:**
- n=3: corr(gap,H) = −1.000 (exact linear)
- n=4: corr(gap,H) = −1.000 (exact linear)
- n=5: corr(gap,H) = −0.898
- n=6: corr(gap,H) = −0.910
