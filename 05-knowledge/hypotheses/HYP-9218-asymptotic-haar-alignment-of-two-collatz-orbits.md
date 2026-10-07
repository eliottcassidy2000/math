---
id: HYP-9218
title: "Asymptotic Haar alignment: for two Collatz orbits related by a fixed affine relation (S19's lag-1 Mersenne debt pairs in the Haar model), the depth delta = v_2(3 kappa + 1 - 3^k) at a tape alignment is, conditionally on the joint past, Haar-distributed up to an error decaying exponentially in the gap L_s; hence the two exponent streams are asymptotically independent, the difference walk L_t obeys an invariance principle with variance 4, and HYP-9217 (1) follows given mac-mini's per-visit merge success"
status: >
  OPEN. NUMERICAL (exact integer orbits, 3000 pairs, 20,000-bit sources, 1.66M alignments):
  unconditional excess eps(l) = E[2^(1-delta) | L_s = l] - 2/3 = +0.043, -0.103, +0.029, -0.094 (l = 1..4, parity-alternating),
  |eps(l)| <= 0.024 for 5 <= l <= 8, within about 1 sigma (0.006) for 9 <= l <= 40;
  P(delta >= d)/2^(1-d) = 1.000 +- 0.001 for L_s >= 8;
  joint law of aligned exponents chi^2 = 66 on 48 dof;
  CONDITIONAL test: consecutive fresh alignments of one pair (both at L >= 8, n = 638,702) are jointly Haar x Haar
  (chi^2 = 32.6 on 35 dof; corr(2^(1-delta_j), 2^(1-delta_(j+1))) = -0.0004 +- 0.0013; P(delta' >= d | delta >= 3)/2^(1-d) = 1.000, 0.999, 1.000, 0.993, 0.981, 0.978).
  The structure it rests on is PROVED (THM-4564).
source: mac-mini-2026-10-07-oaimath3, 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
related:
  - 01-canon/theorems/THM-4564-two-orbits-read-one-tape-alignment-coupling-of-collatz-exponent-streams.md
  - 05-knowledge/hypotheses/HYP-9217-mersenne-debt-recurrence-law-t-minus-half.md (the target)
  - 05-knowledge/hypotheses/HYP-9214-reset-two-debt-resolves-with-probability-tending-to-one.md (per-visit success)
---

# HYP-9218 — asymptotic Haar alignment

**Statement.** In the Haar model of S19's lag-1 Mersenne debt pairs, there are `C` and `λ < 1` such that for every exact alignment `(s, t)` (THM-4564) and every `d ≥ 1`:

    | P(δ ≥ d | joint past, L_s = l) − 2^(1−d) | ≤ C λ^l .

**Why it is the right target.**
* By THM-4564 (2, 3), the two exponent streams interact only at alignments, through an ultrametric law of depth `δ`. Haar depth makes the aligned pair exactly independent.
* By THM-4564 (4), lockstep continuations preserve Haar-ness exactly. So the hypothesis is really about *fresh* alignments: the `x`-orbit arriving at a `y`-block boundary from inside a block.
* With the hypothesis, the increments `a_t − b_t` of `L` are asymptotically uncorrelated, with variance `2 + 2 = 4`.
  * A martingale-approximation argument (HEURISTIC, not carried out) would give an invariance principle for `L` away from 0, and hence recurrence and the `T^(−1/2)` return law.
  * With mac-mini's per-visit merge success bounded below (HYP-9214), that is HYP-9217 (1).

**Evidence.**
* See the status. The excess is real only at `l ≤ 8`, where it alternates with the parity of `l`.
* Even gaps phase-lock: merges need `L` even (S19, Theorem D), and the sibling relation `x = 4^k y + (4^k − 1)/3` sits there.
* At `L_s ≥ 8`, fresh alignments (n = 1.06M) and continuations (n = 0.53M) are both Haar to 0.1%.
* Successive fresh depths of one pair are independent: `χ² = 32.6` on 35 dof (n = 638,702), correlation `−0.0004 ± 0.0013` (script `twoorbit_conditional_haar_20261007.py`).

**What a proof needs.**
* A 2-adic equidistribution statement for the "private predictions" (targets) of two readers of a Haar tape that have not been in lockstep for `≈ L` digits.
* Equivalently: the low digits of the integer debt `D_t` are asymptotically uniform given its archimedean size. By THM-4564 (6), the size is a perpetuity in the recent `y`-exponents, while the low digits were written about `L/2` steps earlier.

---

## Update (2026-10-07, same session): scope narrowed by THM-4569

* Recurrence no longer needs this hypothesis: in the Terras clock the odd-step difference is an exact simple random walk on the predictable disagreement clock (THM-4569 (2)–(3)).
* This hypothesis is now needed only for the **rate**: the disagreement density 1/2 (measured 0.5006, equivalent to the variance 4 per odd step here), and through it the `T^(−1/2)` law of HYP-9217 (1).
