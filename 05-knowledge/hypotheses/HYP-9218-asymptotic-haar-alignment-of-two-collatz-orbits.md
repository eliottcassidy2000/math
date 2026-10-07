---
id: HYP-9218
title: "Refuted full-past Haar-depth formulation; open coarse-conditioning replacement for coupled Collatz overlap depths"
status: >
  REFUTED as a full-joint-past asymptotic distribution claim; a precisely specified coarse-conditioning replacement remains OPEN.
  Correction 2026-10-07: depth is measurable in the full past, so its conditional tail at d=2 differs from 1/2 by exactly 1/2.
  NUMERICAL (exact integer orbits, 3000 pairs, 20,000-bit sources, 1.66M alignments):
  unconditional excess eps(l) = E[2^(1-delta) | L_s = l] - 2/3 = +0.043, -0.103, +0.029, -0.094 (l = 1..4, parity-alternating),
  |eps(l)| <= 0.024 for 5 <= l <= 8, within about 1 sigma (0.006) for 9 <= l <= 40;
  P(delta >= d)/2^(1-d) = 1.000 +- 0.001 for L_s >= 8;
  joint law of aligned exponents chi^2 = 66 on 48 dof;
  CONDITIONAL numerical test: consecutive fresh alignments of one pair (both at L >= 8, n = 638,702) are consistent with Haar x Haar
  (chi^2 = 32.6 on 35 dof; corr(2^(1-delta_j), 2^(1-delta_(j+1))) = -0.0004 +- 0.0013; P(delta' >= d | delta >= 3)/2^(1-d) = 1.000, 0.999, 1.000, 0.993, 0.981, 0.978).
  The structure it rests on is PROVED (THM-4564).
source: mac-mini-2026-10-07-oaimath3, 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
related:
  - 01-canon/theorems/THM-4564-two-orbits-read-one-tape-alignment-coupling-of-collatz-exponent-streams.md
  - 05-knowledge/hypotheses/HYP-9217-mersenne-debt-recurrence-law-t-minus-half.md (the target)
  - 05-knowledge/hypotheses/HYP-9214-reset-two-debt-resolves-with-probability-tending-to-one.md (per-visit success)
---

# HYP-9218 — full-past formulation refuted; coarse-conditioning direction open

**Correction, 2026-10-07.** The alignment depth is a function of the full
joint past. Thus `P(delta>=2 | joint past)` is an indicator, with error
exactly `1/2` from the geometric target. Its error cannot decay to zero
along unbounded alignment gaps. The proposed full-past statement below
is refuted in its intended asymptotic sense. The numerical measurements
average over coarser bins and remain observations of those bins.
Such gaps are nonvacuous: raw y=3 mod4 and v=4k give an alignment at
(s,t)=(0,2k), gap4k+1, and deterministic depth2+v2(k), on a cylinder of
conditional probability1/2. This also excludes a uniform geometric law
at every initial alignment; any replacement needs explicit sampling.
See the [proof and controls](../results/collatz_overlap_kernel_integration_20261007.md).

**Original statement, retained as the failed formulation.** In the Haar
model, there were proposed constants `C` and `λ < 1` such that for every
exact alignment `(s, t)` and every `d ≥ 1`:

    | P(δ ≥ d | joint past, L_s = l) − 2^(1−d) | ≤ C λ^l .

**Strongest surviving research direction.** Choose a coarser conditioning
field that retains the desired debt/return obligation but does not already
determine the depth, and formulate an averaged depth-tail estimate there.
THM-4564 covers exact alignments; THM-4565 extends the kernel to all
overlapping windows, including unequal starts. A geometric depth mixture
makes the normalized pair independent within the stated mixing class.
The transition preserves a Haar depth input only under its explicit input
law, and its equality branch can exceed the clock depth. Cross-time
conditioning, an invariance principle, recurrence, and a usable per-visit
merge chance are additional obligations; none is supplied by a one-pair
moment or by the refuted formulation.

**Evidence.**
* See the status. Detectable deviations in these tables are concentrated at `l ≤ 8`, with a parity pattern; no exact cutoff is proved.
* Even gaps phase-lock: merges need `L` even (S19, Theorem D), and the sibling relation `x = 4^k y + (4^k − 1)/3` sits there.
* At `L_s ≥ 8`, the sampled distributions of fresh alignments (n = 1.06M) and continuations (n = 0.53M) agree with the tested Haar tails to about 0.1%.
* Successive fresh depths are numerically consistent with independence: `χ² = 32.6` on 35 dof (n = 638,702), correlation `−0.0004 ± 0.0013` (script `twoorbit_conditional_haar_20261007.py`).

**What a replacement needs.** A specified distribution of alignment or
overlap events, a smaller conditioning sigma-field, and a quantitative
bound on the resulting depth tails. Equidistribution of low debt digits
given an archimedean size is a possible related formulation, not a proved
equivalence. The exact normalized debt retains a forcing term beyond the
comparison perpetuity. The full depth profile, or the diagonal joint
probabilities that recover it, distinguishes laws which agree on covariance
and other single-moment summaries.

## Concurrent update: THM-4569 supplies count recurrence separately

The Terras-clock odd-count difference is a simple random walk on a
predictable disagreement clock. On nonmerge paths its count is recurrent;
this conclusion does not use the refuted full-past hypothesis. A valid
coarse-conditioning replacement could help with the clock rate or the
translation debt, but recurrence of a count alone is not archimedean box
recurrence or almost-sure merging. The inverse-clock transfer is recorded
in the current session synthesis. The rate remains a separate target.

## Concurrent update (mac-mini-2026-10-07-oaimath3): box recurrence and almost-sure merging are PROVED without any depth hypothesis (THM-4581)

* **What THM-4581 adds.** The archimedean control that the paragraph above says is missing is now supplied. In the Terras clock, the weight `|e·3^(−max(k,0))|^θ s^|k|`, with `s = 2^θ(1 − √(1 − (3/4)^θ))`, is a martingale on flips: a move toward `k = 0` is exactly the ×3/2 branch. Runs at level `h` cost a fresh fair coin per continuation beyond `v_2(3^h − 1) − 1` steps.
* **Consequences.**
  * `E|e|^θ` contracts by `ρ = 1 − √(1 − (3/4)^θ)` (0.634 at `θ = 1/2`) along returns, plus a constant.
  * Every return state reaches `(0, 0)` with positive probability, so the pair merges almost surely (THM-4581 (3)).
  * `c T^(−1/2) ≤ P(no merge by T) ≤ C T^(−1/2) (log T)^2` (THM-4581 (4)).
* **What is left for any coarse-conditioning replacement of this hypothesis.** Only fine structure: an invariance principle for `L_t` with variance 4 per odd step, the exact disagreement density, and the constant of the `T^(−1/2)` law. These are no longer needed for recurrence, merging, HYP-9213, HYP-9214 or HYP-9220.
* **The refutation above is accepted.** The full-past formulation was this session's error. See the MISTAKES entry by codex-tiling.
