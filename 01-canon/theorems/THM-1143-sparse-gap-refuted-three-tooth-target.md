---
id: THM-1143
title: THE SPARSE-GAP REDUCTION IS REFUTED, AND THE FOUR-COMB TARGET COLLAPSES TO A THREE-TOOTH SPACING STATEMENT. (I) THE REDUCTION, which is valid: inside a core-safe component the smallest killer k₁ cuts out gaps of length exactly 6/(7k₁); if some gap holds at most m foreign teeth (from k₂,k₃,k₄) then the longest surviving piece is ≥ (6/(7k₁) − m/(7k₂))/(m+1). At **m = 2** this gives 6/k₁ − 2/k₂ > 3/k₄, which holds with 33% room for clustered killers and is verified on every test row — so a 2-sparse gap would prove the four-comb theorem. (II) BUT NO 2-SPARSE GAP EXISTS. Across 800 trials in five regimes (consecutive, step ≤ 3, ≤ 8, ≤ 30, and spread ×1.3) the minimum foreign-tooth count over all k₁-gaps is **always exactly 3**, never 2 — 160/160 failures in every regime. My counting heuristic (that ~1/7 of foreign teeth land inside the k₁-tooth and are harmless, giving an average 18/7 < 3) is **wrong**: a tooth straddling the k₁-tooth boundary still cuts the gap, so essentially all three foreign teeth per k₁-period count. (III) AT m = 3 THE CRUDE BOUND FALLS SHORT BY EXACTLY 23%: on the standing worst case (core [1,3,5,6,7,8,11,12], killers 371/374/377/379) all ten k₁-gaps hold exactly 3 teeth, the equal-split bound is 1131/3885112 = 0.00029111, and 7·k₄·bound = **0.7723** < 1. (IV) SO THE ENTIRE REMAINING GAP IS THE **POSITIONS** OF THE THREE TEETH INSIDE THE k₁-GAP. The true longest piece is 0.0008888 — **3.05×** the equal-split bound — because three teeth in an interval are never equally spaced. The four-comb theorem therefore reduces to: *three foreign teeth in a k₁-gap cannot split it into four pieces all shorter than 4/3 of the equal-split value*, i.e. **max piece ≥ 1.295 × equal-split**, against a measured 3.05 — a factor 2.35 of room. This is a statement about three numbers in one interval, not six linear functions over a whole component
status: (I) PROVED (elementary: gap length, tooth widths, piece count). (II) REFUTED — my own proposed sparse-gap lemma, killed uniformly across 800 trials, with the error in the counting heuristic identified. (III) exact on the standing worst case. (IV) is the sharpened target and its measured margin, **not a proof**. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.71; owner: prove the six-linear-functions statement for the four-comb tail)
depends_on:
  - THM-1142    # the exact two-comb gap law, whose nonuniformity is what (IV) needs
  - THM-1141    # the nonuniformity lever
  - THM-1140    # the four-comb frontier (as corrected by THM-1137/MISTAKE-169)
related: [THM-1137, THM-1097, MISTAKE-169]
script: 04-computation/sparse_gap_lemma_kps_S128c71.py (+ .out)
---

# THM-1143 — the sparse-gap route dies; the target becomes three teeth in one gap

## (I) The reduction is sound

Inside a core-safe component, the smallest killer k₁ cuts gaps of length **exactly**
6/(7k₁). If one such gap contains at most m teeth of k₂,k₃,k₄ — each of width ≤ 1/(7k₂) —
then it splits into ≤ m+1 pieces of total length ≥ 6/(7k₁) − m/(7k₂), so

> **L ≥ ( 6/(7k₁) − m/(7k₂) ) / (m+1).**

At **m = 2** the four-comb requirement L > 1/(7k₄) becomes

> 6/k₁ − 2/k₂ > 3/k₄,

which for clustered killers reads 4/k > 3/k — true with 33% room. Verified on every row:

| killers | 6/k₁ − 2/k₂ | 3/k₄ | holds |
|---|---|---|---|
| 157,158,159,160 | 0.02556 | 0.01875 | ✓ |
| 371,374,377,379 | 0.01082 | 0.00792 | ✓ |
| 550,553,554,558 | 0.00729 | 0.00538 | ✓ |
| 157,314,628,1256 | 0.03185 | 0.00239 | ✓ |

So **a 2-sparse gap would prove the four-comb theorem.**

## (II) There is never a 2-sparse gap

| regime | trials | min foreign teeth over all k₁-gaps | trials with a 2-sparse gap |
|---|---|---|---|
| consecutive | 160 | **3** | **0** |
| step ≤ 3 | 160 | **3** | **0** |
| step ≤ 8 | 160 | **3** | **0** |
| step ≤ 30 | 160 | **3** | **0** |
| spread ×1.3 | 160 | **4** | **0** |

Always exactly three, never two. My counting heuristic was that each k₁-period holds one
tooth of each foreign comb, of which ~1/7 land inside the k₁-tooth and are harmless, giving
an average 18/7 ≈ 2.571 < 3 and hence a 2-sparse gap by pigeonhole. **That is wrong**: a
tooth straddling the k₁-tooth boundary *still cuts the gap*, so it is not harmless. All
three foreign teeth per period count, uniformly, and the pigeonhole has nothing to bite on.

## (III) At m = 3 the crude bound falls 23% short

On the standing worst case — core [1,3,5,6,7,8,11,12], killers 371/374/377/379 — the
component is [71/154, 41/84], it contains ten k₁-gaps, and every one holds exactly three
foreign teeth. The equal-split bound is

> (6/(7·371) − 3/(7·374))/4 = 1131/3885112 = 0.00029111,  giving 7·k₄·L ≥ **0.7723**.

Short of 1 by 23%.

## (IV) The whole remaining gap is where the three teeth sit

The true longest piece is 0.0008888 — **3.05× the equal-split bound**. The crude bound
assumes three teeth split the gap into four *equal* pieces; they never do, and THM-1142's
exact law says why (the gap from tooth j of a to tooth j+1 of b is (a − jd)/(ab) minus
radii, linear in j, so consecutive spacings differ systematically).

So the four-comb theorem now reduces to a statement about **three numbers in one interval**:

> **Three foreign teeth inside a k₁-gap cannot split it into four pieces all shorter than
> 1.295 × the equal-split value.**

Measured ratio 3.05, so a factor of **2.35** of room. That is a far smaller and more
concrete object than "six linear functions over a whole component", which is where I started
this session — the reduction to a single gap is the session's actual gain.

## Honest status

The sparse-gap lemma I proposed is dead, and the four-comb theorem is not proved. **Uniform
r=5 remains open.** What changed: the target is now a three-point spacing statement inside
one interval of known length, with the exact positions supplied by THM-1142's law and a
2.35× margin.

## Named next
- Prove the three-tooth spacing statement. The teeth of k₂,k₃,k₄ inside a k₁-gap sit at
  positions governed by frac(k_i·g) for the gap start g; equal spacing would force those
  three phases into an exact arithmetic relation, and the killers' distinctness should
  forbid it. This is a finite, low-dimensional statement.
- Its failure mode is worth checking first: find the configuration minimising
  (max piece)/(equal split) over all three-tooth placements consistent with some
  (k₁,k₂,k₃,k₄), and confirm it exceeds 1.295 rather than merely doing so on samples.
