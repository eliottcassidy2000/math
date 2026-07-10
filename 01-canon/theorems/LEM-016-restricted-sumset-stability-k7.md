# LEM-016 — Restricted-sumset stability for 7-sets: B ≤ 12 forces an 8-length AP (SHARP: B = 13 escapes via rank-2 GAPs)

**Status:** PROVED for diameter ≤ 60 (exhaustive pruned DFS) + the largest-gap tail case
g > D/2 (disjoint-block argument below); the remaining tail sliver (D > 60, all gaps ≤ D/2,
B ≤ 12) has no known example and 200k clean samples — flagged, not claimed. SHARPNESS PROVED
(explicit family). **The original conjectured threshold B ≤ 13 was REFUTED by this session's
own verification** — the derive-then-verify loop catching its author a third time.
**Source:** mac-mini-2026-07-09-S65 (cont. 5). Activated by THM-676(iv) and HYP-5682 (the
3k−4 lead); the k = 7 restricted-sum analog of Freiman 3k−4.

**Setup.** `A` = 7 distinct integers; `B = |A +̂ A|` = number of distinct restricted sums
`a_i + a_j` (i < j). Normalize: `a_1 = 0`, gcd of differences 1, diameter `D = a_7`.
Always `B ≥ 11`, with `B = 11 ⟺ A` is an AP (THM-676(iv); the `⟹` direction is Lean-formalized,
kernel-pure, as `LRCFreimanAP.thm676iv_seven_isAP`, citing `finset_min_burden_isAP` — opus-S197). LEM-016
is the excess-1 stability *above* this base case (`B ≤ 12 ⟹ 8-length AP`).

## Statement

> **(i) (stability at excess 1)** `B ≤ 12 ⟹ D ≤ 7`, i.e. A is contained in an arithmetic
> progression of length ≤ 8.
> **(ii) (sharpness — no stability at excess 2)** `B = 13` admits UNBOUNDED diameter: the
> symmetric three-piece family
> `A_c = {0, 1} ∪ {c−1, c, c+1} ∪ {2c−1, 2c}` (any `c ≥ 4`)
> has exactly `B = 13` (sums: `{1} ∪ {c−1..c+2} ∪ {2c−1..2c+1} ∪ {3c−2..3c+1} ∪ {4c−1}`,
> sizes 1+4+3+4+1) and gcd 1, diameter `2c`. It is a rank-2 generalized AP
> (`{0,1} + {0, c−? …}`-shaped). Two-piece configurations give `B ≥ 14`.

*Proof of (i).* Exhaustive: by monotonicity (adding an element never removes restricted sums)
the pruned DFS over normalized `A = {0, a_2 < ⋯ < a_6, D}` is complete per diameter; running
`D = 6..60` finds max diameter 7 at `B ≤ 12` (and 6 at `B = 11`). Tail `D > 60`: if the largest
gap `g > D/2`, split `A = L ⊔ R` at it; spans `ℓ + r = D − g < g`, so the three sum-blocks
`L+̂L ⊆ [·, 2ℓ]`, `L+R ⊆ [ℓ+g, ℓ+D]`, `R+̂R ⊆ [2(ℓ+g)+1, ·]` are pairwise disjoint
(`2ℓ < ℓ+g` since `ℓ < g`; `ℓ+D < 2ℓ+2g+1` since `r ≤ g`), giving
`B ≥ (2s−3)₊ + 6 + (2(7−s)−3)₊ ≥ 14 > 12` for every split size `s`. The sliver
(`D > 60`, all gaps ≤ `D/2`) is not covered by this argument — no example with `B ≤ 12` is
known (200k samples clean); (i) is claimed as PROVED only for `D ≤ 60` ∪ {g > D/2}. ∎

**Lean formalization of the `g > D/2` case (opus-S198, kernel-pure).** The dominant-gap branch of (i) is
machine-checked: `LRCFreiman.burden_ge_of_dominant_gap` (`LRCBurdenGap.lean`) states that ANY `k`-set
`s = L ⊔ R` split at a gap dominating both spans (`span(L) + span(R) < gap`, i.e. `g > D/2`) has
`|A +̂ A| ≥ 3k − 7`; the `k = 7` corollary `burden_ge_fourteen_of_dominant_gap` gives `≥ 14`. (The same
lemma at `k = 13` gives `≥ 32` — monad's THM-682 "core `B ≥ 32`" bound, on the dominant-gap branch.) It is exactly the disjoint-block argument above: `three_block_card` proves the three
sum-blocks `L +̂ L`, `L + R`, `R +̂ R` lie in disjoint integer ranges (from the two separation
inequalities the gap condition supplies) so the restricted sumset of the union contains all three
disjointly; the within-block sizes are `restrictedSum_card_ge'` (`≥ 2|L|−3`, `≥ 2|R|−3`) and the cross
block is Mathlib's Cauchy–Davenport for ℤ (`cauchy_davenport_add_of_linearOrder_isCancelAdd`,
`≥ |L|+|R|−1`); `3·7 − 7 = 14` by `omega`. So the **contrapositive `B ≤ 13 ⟹ no gap exceeds `D/2`** is
now kernel-pure, reducing (i) to the flagged all-gaps-`≤ D/2` sliver at the Lean level too. Axioms:
`[propext, Classical.choice, Quot.sound]`; no `sorry`, no `native_decide`.

*Proof of (ii).* Direct computation of the sum blocks above (script-verified at
c = 30, 250, 2500). ∎

## Consequence — the burden TRICHOTOMY (THM-676 addendum; mirrors the LRC branch tree)

For the majority-parity class A (|A| = 7) of any 13-set, with burden `B` = its distinct
half-sum count (THM-676(ii)):

| burden | structure (this lemma) | LRC branch it feeds |
|---|---|---|
| B = 11 | A is an AP | near-AP branch — LEM-012 Dirichlet (proved) |
| B = 12 | A ⊆ 8-length AP (one hole) | near-AP branch — LEM-012 territory |
| B = 13 | rank-2 GAP escapes possible | dissociated/loose — LEM-013 + opus-S181 (GAPs are LOOSE) |
| B ≥ 14 | genuinely spread | ≥ 14 independent forced blocked-moduli demands (THM-674 domination collapse: 13.8% → 0.25% incidence) |

**The additive-structure ladder and the LRC proof-branch ladder are the same tree** — measured
at the burden level: exactly opus-S181's law (tightness = 1-dim coherent; GAPs = loose) and
klein-S212's composed dichotomy (coherent branch = longest-AP ≥ k−6 = LEM-012), reached here
through pure restricted-sumset combinatorics with no LRC input.

**Near-AP-side closure data (structured branch):** covering 13-sets built around
7-of-9-consecutive majority classes: 36/36 certified by C0/C2 (zero uncertified) —
`lrc14_freiman_stability_macmini_S65cont5.{py,out}`.

**Related:** THM-676 (burden), THM-672/674 (blocking anatomy), LEM-012/013, LEM-015 +
LRCSchurRigidity (E3 rigidity), opus-S181 (GAP law), klein-S212, HYP-5682/5765.
