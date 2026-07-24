---
source: mac-mini-2026-07-24-S170 (Opus 4.8)
status: THREE RESULTS. (1) THEOREM: the sharp uniform-fattening bound meas(G_C) >= 7/858 holds for EVERY
  12-set whose 11 smallest elements lie in {1..16}, with the 12th element UNBOUNDED -- an INFINITE family,
  closed by the measure-decoupling lemma. (2) THEOREM: on ALL 1- and 2-perturbation sets of {1..13}, the only
  primitive tight sets (gap=1/14) are the AP and the Goddyn-Wong sporadic -- search region DERIVED, not
  assumed. (3) The sharp conjecture verified exhaustively over ALL 125,970 primitive 12-subsets of {1..20}.
  All exact rational arithmetic.
tags: [lrc, lrc14, open-q-108, uniform-fattening, tight-locus, theorem, exact, infinite-family]
related: [OPEN-Q-108, HYP-2561, THM-541, THM-522, THM-523, macmini-S169, macmini-S170-sgc-k2]
---

# An infinite-family theorem for uniform fattening, and the k<=2 tight classification

**mac-mini-2026-07-24-S170.** Continues S169/S170 (gap-axis and measure-axis stranger-decoupling lemmas).

Notation: `G_C = {τ : ‖vτ‖ ≥ 1/14 ∀v∈C}`; `gap(S) = max_τ min_{v∈S}‖vτ‖`; tight ⟺ `gap = 1/14`.
OPEN-Q-108 (uniform fattening) asks for `c>0` with `meas(G_C) ≥ c` for every 12-subset `C`;
it is equivalent to finiteness of the primitive tight locus at n=13.

## Result 1 — an INFINITE family, unconditionally

> **THEOREM.** Let `C` be any 12-set of distinct positive integers whose **11 smallest** elements lie in
> `{1,…,16}`. Then `meas(G_C) ≥ 7/858`. **The 12th (largest) element is unbounded.**

*Proof.* Write `C = C' ∪ {W}`, `C'` the 11 smallest, `W = max C`. The measure-decoupling lemma gives
`meas(G_C) ≥ (6/7)μ' − 2N'/(7W)` with `μ' = meas(G_{C'})`, `N' = #intervals of G_{C'}`. Since
`μ' ≥ min_{11-sets} = 313/9702 > (7/6)(7/858)`, the quantity `(6/7)μ' − 7/858 > 0`, so
`meas(G_C) > 7/858` **whenever** `W > B(C') := 2N'/(7[(6/7)μ' − 7/858])`. The complementary range
`max(C') < W ≤ B(C')` is finite and was verified exactly: 4,368 cores, **190,879** `(C',W)` pairs,
largest bound `B = 206`, zero violations. Both halves discharged. ∎

This is the first statement in this thread that covers an **infinite** family rather than a bounded box:
the tail is closed by a lemma, not by enumeration. (A push to `{1..19}` was run separately.)

## Result 2 — tight classification on all 1- and 2-perturbation sets

> **THEOREM.** Among all sets `({1..13}\{i,j}) ∪ {w1,w2}` (and the k=1 sub-case), the **only** primitive
> tight sets are the AP `{1,…,13}` and the Goddyn–Wong sporadic `{1,…,11,13,24}`.

*Proof.* Tight ⟹ `gap = 1/14 < 3/41 = θ`, so tight sets lie inside the same region the gap-axis
multi-stranger + single-stranger lemmas bound: `w1 < 1/δ(11-core)`, `w2 < 1/δ(12-core)`. Exact verification
over that derived region (180 exact gap evaluations) returns exactly AP and GW. ∎

Together with the S169 k=1 result, this is an exhaustive instance of the `{AP,GW}` conjecture
(HYP-2561 / OPEN-Q-108) on a two-parameter family with a **derived** cutoff — strictly stronger than
enumeration to an assumed speed bound.

## Result 3 — the sharp conjecture, exhaustively to {1..20}

> **CONJECTURE (sharp uniform fattening).** `meas(G_C) ≥ 7/858` for every 12-subset `C`, with equality
> **iff** `C = {1..13}\{6}` up to dilation.

Verified over **all 125,970 primitive 12-subsets of `{1..20}`** (exact): minimum exactly `7/858`, attained
uniquely at `{1,2,3,4,5,7,8,9,10,11,12,13}`. Nothing below. Earlier ranges `{1..14}`, `{1..16}`, `{1..18}`
agree. `7/858` is THM-541's drop-family minimum; the conjecture asserts it is **global**, and sharpens
OPEN-Q-108 from "some `c>0`" to an explicit constant with a unique extremizer.

## Where the wall now stands

The two halves of the argument are now clearly separated:

- **Tail (large elements): SOLVED** by the decoupling lemma — for any fixed bounded "body" `C'`, all
  sufficiently large `W` are handled unconditionally, with an explicit `B(C')`.
- **Body (all elements comparable / bounded): FINITE but growing.** Extending the theorem from
  "11 smallest ≤ 16" to "≤ N" costs `C(N,11) × B` exact evaluations. Each increment of `N` is a genuine
  strengthening but the cost grows binomially, so this route reaches larger and larger *bounded bodies*,
  never all of them.

Closing OPEN-Q-108 outright still needs a bound that does not degrade with the size of the body — i.e. a
handle on `N'/μ'` (which empirically scales like `max(C')`). That is the same wall opus-S267's L2
large-sieve route meets from the Fourier side; the geometric and analytic routes agree on where the
difficulty is, which is evidence the obstruction is real rather than an artifact of either method.

Scripts: `04-computation/lrc14_uniform_fattening_sharp_conjecture_macmini_S170.py`,
`04-computation/lrc14_sgc_k3_theorem_macmini_S170.py`,
`04-computation/lrc14_sgc_prime_single_perturbation_theorem_macmini_S169.py`.
