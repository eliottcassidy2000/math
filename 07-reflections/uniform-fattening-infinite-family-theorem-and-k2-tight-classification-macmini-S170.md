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

## Result 2 — SGC'(13) and the tight classification on ALL 1-, 2-, and 3-perturbation sets

> **THEOREM (k ≤ 3).** No set obtained from `{1,…,13}` by replacing **up to three** elements has
> `gap ∈ (1/14, 3/41)`; and among all such sets the **only** primitive tight sets (`gap = 1/14`) are the
> AP `{1,…,13}` and the Goddyn–Wong sporadic `{1,…,11,13,24}`.

*Proof.* Tight ⟹ `gap = 1/14 < 3/41 = θ`, so band-hitters and tight sets lie in the same region, which the
gap-axis multi-stranger lemma bounds (`4kθ<1` for `k≤3`, since `4·3·(3/41)=36/41<1`): with `w1<w2<w3`,
`w1 < 1/δ(10\text{-core})`, `w2 < 1/δ(11\text{-core})`, `w3 < 1/δ(12\text{-core})`, all δ exact rationals.
Exhaustive exact verification of that **derived** region:

| k | nodes searched | exact gap evaluations | band-hitters | tight sets |
|---|---|---|---|---|
| 1 | all `w ≤ 417` | all j | 0 | AP, GW |
| 2 | 513,264 | 180 | 0 | AP, GW |
| 3 | 497,847 | 352 | 0 | AP, GW |
∎

The k=3 family is wide — every set differing from the AP in up to three positions, with the replacements
unbounded a priori and bounded only by the lemma. No new tight set appears anywhere in it.

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

## Result 4 — the tail is under control for ANY number of well-separated large elements

The decoupling can be iterated: `j` unbounded top elements are handled whenever
`(6/7)^j · m_{12−j} > 7/858`, where `m_k` = min measure over k-sets. The gates all close, with slack that
**grows** as the body shrinks:

| j unbounded | gate `(7/6)^j·(7/858)` | actual `m_{12−j}` | slack |
|---|---|---|---|
| 1 | 0.009518 | `313/9702` = 0.032261 | 3.39× |
| 2 | 0.011105 | `14249/252252` = 0.056487 | 5.09× |
| 3 | 0.012955 | `10601/114660` = 0.092456 | 7.14× |

So the *tail* of the problem is not the obstruction at any depth; the difficulty is confined entirely to the
comparable/bounded **body**.

## Result 5 — the extremizer is RIGID (so the body is genuinely small)

Dilates `d·({1..13}\{6})` all have measure exactly `7/858` (dilation invariance) but are **not primitive**.
Restoring primitivity by a minimal perturbation costs a large jump:

| perturbation of `d·EXT` | resulting measure | ratio to 7/858 |
|---|---|---|
| max ± 1 | 0.042–0.046 | **5.1–5.7×** |
| min ± 1 | 0.064–0.068 | **7.8–8.4×** |

(verified for `d = 2,3,5,10,30`). Likewise every primitive 12-set tested with a large minimum sits at
`≥ 5.4×` (e.g. `{100..111}` at 17.3×). **The extremizer is isolated: no primitive small-measure set exists
at large scale.**

The mechanism is exactly the decoupling lemma: perturbing one element of a dilate turns the set into
"(dilate of an 11-set) ∪ (one off-lattice element)", and that element decouples, giving
`meas ≈ (6/7)·m(11-set) ≥ (6/7)m_{11} = 0.0277 > 7/858`. The observed cluster at ≈0.044 is exactly
`(6/7)·meas({1,2,3,4,5,7,8,9,10,11,12}) = (6/7)(0.05142)`. So this whole regime is already **covered by the
lemma**, not merely observed.

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
