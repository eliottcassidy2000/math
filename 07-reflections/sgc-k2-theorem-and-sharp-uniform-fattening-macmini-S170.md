---
source: mac-mini-2026-07-24-S170 (Opus 4.8)
status: (A) THEOREM PROVED (exact): SGC'(13) holds on ALL 1- and 2-perturbation sets -- no
  ({1..13}\{i,j}) u {w1,w2} has gap in (1/14, 3/41). Uses a new MULTI-STRANGER lemma (up to 3 strangers
  decouple simultaneously since 4k*theta<1). (B) OPEN-Q-108: a SHARP explicit conjecture
  (meas(G_C) >= 7/858, unique extremizer {1..13}\{6}), a proved MEASURE-DECOUPLING lemma, and an induction
  step that CLOSES with 3.4x slack -- reducing OPEN-Q-108 to a finite check, though not an effective one.
tags: [lrc, lrc14, sgc, stranger-decoupling, open-q-108, uniform-fattening, theorem, exact, sharp-constant]
related: [THM-518, THM-522, THM-523, THM-541, OPEN-Q-108, HYP-2561, macmini-S169]
---

# SGC'(13) on all 1- and 2-perturbation sets; and a sharp form of OPEN-Q-108

**mac-mini-2026-07-24-S170.** Continues S169 (gap-axis stranger-decoupling lemma; SGC'(13) proved for k=1).
Two targets were requested: complete k=2, and attack the uniform-fattening wall.

## A. The multi-stranger lemma, and the k=2 theorem

> **MULTI-STRANGER LEMMA.** Let `f_C ≥ θ` on an interval `I` of length `δ`, and let `w_1,…,w_k` all satisfy
> `w_i ≥ 1/δ`. Each bad set `B_i = {τ∈I : ‖w_iτ‖<θ}` has `meas(B_i) ≤ 2θδ + 2θ/w_i ≤ 4θδ`. So if
> `4kθ < 1`, then `⋃B_i` cannot cover `I` and some `τ∈I` is good for all of them, whence
> `gap(C ∪ {w_1,…,w_k}) ≥ θ`.

With `θ = 3/41`: `4kθ<1 ⟺ k < 41/12 ≈ 3.4`, so **up to 3 strangers decouple simultaneously**. For `k=2`,
`8θ = 24/41 = 0.585 < 1`. Consequence: a band-hitting 2-perturbation set must have its **smaller** stranger
`w_1 < 1/δ(11\text{-core})`; the S169 single-stranger lemma then bounds `w_2 < 1/δ(12\text{-core } C∪\{w_1\})`.
Both bounds are computed **exactly** (δ is rational).

> **THEOREM (k=2).** No set `({1..13}\{i,j}) ∪ {w_1,w_2}` has `gap ∈ (1/14, 3/41)`.

*Proof.* The two lemmas bound the entire search region; exact rational verification over it
(78 pairs, 513,264 filtered candidates, 180 exact gap evaluations) finds none. ∎

**Combined with the S169 k=1 theorem: SGC'(13) holds on all 1- and 2-perturbation sets**, and the extremal
value `3/41` is still attained only at `{1,…,11,13,36}`. (The k=3 case is reachable by the same lemma,
since `4·3·θ = 36/41 < 1`.)

## B. OPEN-Q-108 — a sharp conjecture and a closing induction step

`G_C = {τ : ‖vτ‖ ≥ 1/14 ∀v∈C}`. OPEN-Q-108 asks for *some* `c>0` with `meas(G_C) ≥ c` for every 12-subset.

> **SHARP CONJECTURE.** `meas(G_C) ≥ 7/858` for every 12-subset `C`, with equality **iff**
> `C = {1..13}\{6}` up to dilation.

`7/858` is exactly THM-541's minimum over the drop-family; the conjecture is that it is **global**.
Verified on **7,910** primitive 12-sets — exhaustive over `{1..14}` and `{1..16}`, plus 6,000 random across
scales to 600 — minimum attained uniquely at `{1..13}\{6}`; nothing below.

> **MEASURE-DECOUPLING LEMMA (proved).** `meas(G_{C∪{W}}) ≥ (6/7)·meas(G_C) − 2N/(7W)`, `N = #intervals of G_C`.
> *Proof.* `G_C` is `N` intervals of total measure `μ`. The bad set `{‖Wτ‖<1/14}` is intervals of length
> `1/(7W)` spaced `1/W`; an interval of length `ℓ` meets at most `ℓW+2` of them, so
> `meas(bad ∩ G_C) ≤ μ/7 + 2N/(7W)`. ∎

**The induction step closes.** `min` over 18,738 primitive 11-sets is `313/9702 ≈ 0.032261` (at
`{1,2,3,4,5,7,8,9,11,12,13}`), comfortably above the required `(7/6)(7/858) = 49/5148 ≈ 0.009518` —
**3.39× slack**. Hence for large `W`, `meas(G_{C'∪{W}}) ≥ (6/7)m_{11} − ε > 7/858`, so **any 12-set attaining
the minimum has a bounded maximum element**: `W ≤ (2/(7·slack))·N' = 14.656·N'`.

**Supporting structure (all exact):** dilates are invariant (`d·{1..12}` → 0.034101 for every `d`); large
*comparable* sets have **large** measure (`{N..N+11}` → ≈0.141 for `N≥5`); lacunary/smooth sets are large
(0.28–0.30). So small measure is forced onto **bounded, multiplicatively-spread** sets — precisely the
finiteness shape.

## C. Honest limit: a reduction, not an effective proof

`W ≤ 14.656·N'` is effective exactly where it matters — the extremal sets have `N' = 10…22`, giving
`W ≤ 147…322` — but `N'` is **not** bounded absolutely (it grows with the elements of `C'`). So the general
statement is a *ratio-type* bound (`max ≲ const × the rest`), and recursing 11 levels gives an
astronomically large absolute bound (~`30^11`). **OPEN-Q-108 is thereby reduced to a finite check, but not
to an effective one.** Closing it needs either a bound on `N'` in the small-measure regime (empirically
`N'` *is* small exactly when the measure is small — that self-consistency is the promising crack), or a
different argument for the comparable-scale case.

Scripts: `04-computation/lrc14_uniform_fattening_sharp_conjecture_macmini_S170.py` (exact `meas(G_C)`,
the conjecture check, the induction step), `04-computation/lrc14_sgc_prime_single_perturbation_theorem_macmini_S169.py`.
