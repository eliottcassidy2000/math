---
source: opus-2026-07-09-S198
status: FORMALIZED the disjoint-block burden bound (LEM-016(i) dominant-gap case, kernel-pure):
  burden_ge_of_dominant_gap -- a 7-set split at a gap dominating both spans (g > D/2) has burden >= 14,
  so B <= 13 => no dominant gap, reducing LEM-016(i) to its sliver. Via three_block_card (three sum-blocks
  in disjoint ranges) + restrictedSum_card_ge' (within) + Mathlib Cauchy-Davenport for Z (cross). AND
  noticed a genuine convergence: kps-S126 independently hit the SAME k>=5 threshold on the E3 axis
  (dist<=deficit false for k<=4, holds k>=5) that I hit on the burden axis (MISTAKE-133, n>=5). Both
  Freiman axes fail on the same small sets.
tags:
  - lrc14
  - freiman
  - disjoint-block
  - lem-016
  - cauchy-davenport
  - convergence
  - k5-threshold
---

# Both Freiman axes fail at k ≤ 4; and the disjoint-block bound, formalized

**opus-2026-07-09-S198.** Owner: work the still-open burden axis. I formalized the dominant-gap case of
LEM-016(i), and — pulling kps-S126 mid-session — found that the two Freiman axes of the THM-681 ladder
share a threshold that is not a coincidence.

## The disjoint-block bound (LEM-016(i), g > D/2 case), formalized

`LRCFreiman.burden_ge_of_dominant_gap` (`LRCBurdenGap.lean`, kernel-pure): ANY `k`-set `s = L ⊔ R` split
at a gap dominating both sides' spans (`span(L) + span(R) < gap`, equivalently `g > D/2`) has descent
burden `≥ 3k − 7` (`k = 7` corollary `≥ 14`; and — a bonus of proving it for general `k` — `k = 13` gives
`≥ 32`, exactly monad's THM-682 "core `B ≥ 32`" bound on the dominant-gap branch). The engine is
`three_block_card`: the three sum-blocks `L +̂ L`, `L + R` (Minkowski), `R +̂ R` occupy disjoint integer
ranges — the ranges' separation follows from the gap condition by `omega` — so the restricted sumset of
`L ∪ R` contains all three *disjointly*, and

  `|A +̂ A| ≥ |L +̂ L| + |L + R| + |R +̂ R| ≥ (2|L|−3) + (|L|+|R|−1) + (2|R|−3) = 3·7 − 7 = 14 > 12`.

The within-block sizes are `restrictedSum_card_ge'` (my S188 Freiman lower bound, made unconditional); the
cross block `|L + R| ≥ |L| + |R| − 1` is **Mathlib's Cauchy–Davenport for ℤ**
(`cauchy_davenport_add_of_linearOrder_isCancelAdd` — the linearly-ordered-cancellative form, which ℤ
satisfies). Finding that lemma already in Mathlib is what made this feasible in one session: the cross-block
bound is the one genuinely additive-combinatorial ingredient, and it was a free import.

**Contrapositive: `B ≤ 13 ⟹` no gap exceeds `D/2`.** So the Lean proof reduces LEM-016(i) to exactly the
sliver mac-mini flagged on paper (`D > 60`, all gaps `≤ D/2`), and no further — an honest Lean mirror of
the paper's honest gap.

## The convergence: both Freiman axes fail at k ≤ 4, hold at k ≥ 5

kps-S126 (`LRCSchurPeel.lean`) materialized the E3-axis stability ladder and flagged HYP-5852: their
capstone `dist ≤ deficit` is **FALSE for k ≤ 4** (counterexample `{1,4,5}`) and **holds for k ≥ 5** — and
they noted this is the "SAME threshold as opus `ap_of_min_burden`," whose `n ≥ 5` hypothesis I was forced
to add after `{0,1,3,4}` broke the naive general claim (MISTAKE-133).

This is not a coincidence of two off-by-one bookkeeping choices. Both axes measure Freiman rigidity, and
both lose it in the same place:

- **Burden axis** (translation-invariant, restricted sumset): at `k ≤ 4` the minimal restricted sumset is
  achieved by bi-arithmetic sets (`{0,1,3,4}`) that are not APs. The even/odd gap tie needs a 5th element
  (`a₄`) to break.
- **E3 axis** (translation-broken, Schur incidences): at `k ≤ 4` the deficit-distance capstone fails on
  `{1,4,5}`; the peeling recursion needs enough elements for the reflection-symmetry base to bite.

`k = 5` is where a 7-set-style Freiman argument first has enough room for its second-order structure — the
smallest size at which "few sums / high incidence ⟹ dilated interval" becomes rigid rather than
accidentally satisfiable. Two independent formalizations, on two complementary measures, drew the boundary
at the same integer. When the mathematics repeats a constant across frameworks, that constant is telling
you something structural — here, that Freiman rigidity for these small sets is a `k ≥ 5` phenomenon on
*both* axes, and the `k ≤ 4` failures are the same "too small to be rigid" phenomenon wearing two costumes.
The LRC application lives at `k = 7`, safely past the boundary on both axes.

## Ledger

- FORMALIZED (kernel-pure) `burden_ge_of_dominant_gap` (LEM-016(i), `g > D/2` case) + `three_block_card` +
  the block-range helpers (`LRCBurdenGap.lean`). `B ≤ 13 ⟹` no dominant gap; reduces LEM-016(i) to its
  sliver. Cross-block bound = Mathlib Cauchy–Davenport for ℤ.
- NOTICED the cross-axis `k ≥ 5` convergence (burden: MISTAKE-133 / `ap_of_min_burden`; E3: kps HYP-5852 /
  `LRCSchurPeel`). Both Freiman axes fail on the same small sets.
- Canon: LEM-016 now records the `g > D/2` Lean formalization.
- Files: `LRCBurdenGap.lean` (+root). Remaining on the axis: the `D > 60` all-gaps-`≤ D/2` sliver (open on
  paper too) and `AP ⟹ non-covering` (klein-S211). → LEM-016, THM-676(iv), kps-S126/THM-681 ladder,
  MISTAKE-133, opus-S195/S196/S197.
