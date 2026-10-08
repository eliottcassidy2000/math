# Audit C — THM-4603, HYP-9242, HYP-9240 remark, results note §§0–3, 5, 6

Independent adversarial audit, 2026-10-07, for session mac-mini-2026-10-07-twoanchor (continuation). Saved by the main session from the auditor's final message; the harness blocked subagent report writes. Scripts in this directory, each with a matching `.out`.

## Verdicts and corrections (all applied)

**1. Drift lemma — CONFIRMED** (2,748 exact tests: k in [−9, 9], boundary precision, integral cycle points, the halving point 0).
* The lemma holds for any T-periodic points c, c′, including 0.
* `v2 ≥ jΛ` suffices for the Terras-step conclusions; `>` is needed only for an exact last U-letter.
* The session's own test covered only k ∈ [0, 6].

**2. Ladder completeness — CONFIRMED on 2,948 merges** (10 named states plus 1,500 random admissible states). Scope correction:
* The state must be admissible (3^max(0,−k) e ∈ Z), and each orbit must have an odd value before the merge.
* Counterexamples otherwise:
  * the reset state (1,1) merges at s = 1 with no odd value before s;
  * inadmissible starts give half-ladders with odd b.
* For non-integer 2-adic orbits the evenness of the letter difference follows from admissibility mod 3.
* The ladder 4z+1 and the predecessor formula are classical folklore.
* The universal-state frequencies replicate (i = 1 share 0.936; source share 0.431).

**3. Transitions.**
* All 532 table lines re-derive. Seven transitions were checked on actual integers.
* *Classification.* It is conditional on eventual periodicity of the rational limit orbit, an open 3x+d question. For c′ = −1 this contains Collatz on the negative integers. All 550 computed cases are periodic.
* *Realisability.* The universal-state lists were incomplete, and 4 of the 8 absorptions (source (2,4), (2,6), (2,2,4), (2,3,3)) start with letter 2, which is not realisable after a two-run. The halving point 0 is not realisable either.
* *SHIFT label.* 224 of the 280 SHIFT outcomes are different cycles of equal odd density, not the same cycle out of phase.
* *Heavy tail (overclaim).* Positive-drift runs explain unbounded debt, not the s^(−1/2) tail. From the universal state, survival to s = 3000 is 0.316 overall and 0.312 among paths without a source ones-run ≥ 14.

**4. HYP-9242.**
* *Data.* Mersenne data are byte-identical to an independent recomputation (K ≤ 12800, 135 orphans).
* *Partner criterion.* It is exact for every K: aligned orbits that first reach 1 together agree from 8 on. 173 partner and 155 non-partner pairs were checked directly. There are no sporadic coincidences.
* *Exponent (overclaim of "about 1/2").*
  * Mersenne fit 0.59, 95% CI [0.42, 0.75].
  * General sources, 2000 per row, B ≤ 3200: about 0.42 for K = 2, 3, and 0.54–0.63 for K = 9–33, excluding 1/2.
  * The effective exponent grows with the number of children, and small-L estimates are biased downward (K = 2 measures 0.42, against a theoretical 1/2).
* *"10× enrichment among long orbits" is a selection artifact.* Partners share σ_T − K, so non-orphans inherit their length from a smaller exponent.
  * A null model with the class structure and random-orbit lengths reproduces the terciles.
  * The genuine signal is that orphans' post-run lengths are +1.0 SD long against random orbits of the same size (median +0.94, 86% positive, n = 83). Non-orphans sit at −0.01.
* *Leading digits.* The frac(K log2 3) quarters hold for K in [600, 6400]; 1.6–2.2% for [1000, 6400].
* *Confirmed:* nearest partner D = 1 for 92%; the shells; ones-runs (z = −0.38); 3-adic rescue 20/53 against a population rate of 0.403 (depth 3 ⟺ K ≡ 5 mod 6); the 33 doubly uncovered exponents; the general-source table.
* *Constants.* The √(K+B)-scaled constants are not constant (3.0 → 2.5 for K = 9; 2.7 → 1.7 for K = 17).

**5. Partial result — CONFIRMED at sketch level.**
* Say "t uniform among the odd L-bit integers in the residual class", not "Haar".
* Inputs: THM-4581 (4)'s constant is at most linear in the debt; P(J ≥ m) = 4^(1−m); the forced prefix is (0,1) or (0,0,1); the D = 1 law does not depend on K.
* P(not absorbed by L−K−2)·√L = 7.1 → 14.3 for L = 100 → 1600, which is pre-asymptotic.

**6. HYP-9240 remark.**
* CONFIRMED: N_D is a function of D and the first s_0 binary digits of 3^(−D); N_D + 1 < (3/2)^D strictly; observed landings are far smaller.
* The "zeroless-power" sentence is an analogy, not a reduction (mild overclaim).

**7. Typing and prior art.** Parts (1) and (2) are elementary bookkeeping on classical facts: the 4z+1 ladder, Böhm–Sontacchi cycle points, and 2-adic continuity of parity vectors.

## Scripts
* `c1_drift_lemma.py`
* `c2_ladder_completeness.py`, `c2b_universal_freq.py`
* `c3_transitions.py`, `c3_heavy_tail.py`
* `c4_mersenne_orbits.py`, `c4_orphan_stats.py`, `c4_partner_merge_direct.py`, `c4_sporadic.py`, `c4_rescue.py`, `c4_ones_runs.py`, `c4_length_selection.py`, `c4_general_sources.py`, `c4_general_exponent.py`, `c4_fit_session_general.py`
* `c5_partial_result.py`
* `c6_hyp9240_remark.py`
