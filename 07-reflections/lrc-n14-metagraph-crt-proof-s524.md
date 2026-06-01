# Reflection: LRC@14 proof attempt via metagraph walk + CRT quotient

**Session:** opus-2026-06-01-S524
**Date:** 2026-06-01

## The CRT quotient reduction

14 = 2 × 7. The 13 runners partition into 7 mod-7 classes:
- 6 pairs: {1,8}, {2,9}, {3,10}, {4,11}, {5,12}, {6,13}
- 1 singleton: {7}

**LRC@14 = all 7 CRT classes simultaneously safe.**

Each pair-class is safe iff BOTH runners have ||vt|| >= 1/14.
The class-safe measures are remarkably close to the independence
prediction (6/7)² ≈ 73.5%:

| Class | Speeds | Safe % | (6/7)² | Ratio |
|-------|--------|--------|--------|-------|
| 1 | {1,8} | 73.2% | 73.5% | 0.997 |
| 2 | {2,9} | 73.0% | 73.5% | 0.994 |
| 3 | {3,10} | 72.8% | 73.5% | 0.992 |
| 4 | {4,11} | 73.0% | 73.5% | 0.994 |
| 5 | {5,12} | 73.4% | 73.5% | 0.998 |
| 6 | {6,13} | 73.4% | 73.5% | 1.000 |
| 0 | {7} | 85.7% | 85.7% | 1.000 |

**Pairwise nearly independent!** But the all-7 intersection is wall-only
for the initial segment, despite the product of measures being 13.1%.

## The coupon collector picture

The quotient outdegree d_Q(t) counts how many of 7 classes are safe.
The histogram for the initial segment:

```
d_Q=6: 55% of the time (6/7 safe, 1 class blocking)
d_Q=5: 32.4% (5/7 safe)
d_Q=4: 4.9%
d_Q=3: 3.4%
...
d_Q=7: 0% (wall-only)
```

The walk spends most of its time at d_Q=6 — ONE class short.
The "last class" to become safe is the bottleneck. This is a coupon
collector problem: you have 7 coupons (classes) and need all 7 at once.

## The proof gap

**What I proved:** The CRT structure reduces LRC@14 to a 7-class
coverage problem. Each class is safe ~73% of the time. Pairwise
correlations are tiny (ratio > 0.99). 498/498 tested speed sets are lonely.

**What I didn't prove:** That the 7-class intersection is always nonempty.
The correlation between ALL 7 classes simultaneously being safe is the crux.
The independence model gives 13.1% (huge), but the actual value for the
initial segment is 0 (wall-only).

The remaining gap is: **bound the 7-way correlation** using the arithmetic
structure (CRT, sieve, gap race).

## The creative tweak that might work

The metagraph walk gives more than just measures — it gives REACHABILITY.
The walk on the 7-class quotient has d_Q=6 most of the time.
By THM-387, the walk is directional.
The last blocking class rotates among the 7 classes.
**Each class takes its turn as the last blocker.**

If the "last blocker" rotates through all 7 classes, there must be a
transition moment where the old last blocker becomes safe and the new
one hasn't started blocking yet. At that instant, d_Q = 7 → lonely!

This is the "musical chairs" argument: 7 blockers take turns, and during
the handoff between consecutive blockers, the observer is momentarily free.

## 498/498 is not a proof, but it's strong evidence

The hardest case is the initial segment {1,...,13}: it's the ONLY set where
the lonely time is wall-only (d_open = 12, not 13). All other tested sets
achieve d=13 in open cells. This confirms the initial segment is the tight
extreme, as expected from the "regular polygon" theorem (S523).
