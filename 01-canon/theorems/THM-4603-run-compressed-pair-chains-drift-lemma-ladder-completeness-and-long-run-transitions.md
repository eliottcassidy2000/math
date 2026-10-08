---
id: THM-4603
title: "Run-compressed pair chains: drift lemma, ladder completeness, and long-run transitions conditional on an authenticated periodic limit pair"
status: >
  PROVED: (1) and (2) (elementary). CONDITIONAL: (3) the classification requires a supplied or independently proved eventual-periodicity witness for the rational limit pair; bounded denominators alone do not provide one. OPEN: unconditional eventual periodicity for every such limit pair. FINITE-EXACT: (1) on 720 exact
  integer tests over 36 word pairs; (2) on all 1651 merges observed from the universal state in 3000 residual sources (child ladders i = 1, 2, 3:
  870, 50, 18; source ladders i = 1, 2, 3, 4: 672, 33, 7, 1); (3) the exact limit-transition table for 5 states x 55 cycles (U-words with sum <= 8,
  length <= 4, and the halving point 0), with the two universal-state transitions re-checked on 40 + 40 actual sources.
session: mac-mini-2026-10-07-twoanchor (continuation)
source: 05-knowledge/results/runcompress_orphans_cayley_20261007.md
scripts:
  - 04-computation/experiments/runcompress_20261007/ladder_and_drift.py (ALL CHECKS PASSED)
  - 04-computation/experiments/runcompress_20261007/transitions.py (+ transitions.out), drift_check_actual.py (ALL CHECKS PASSED)
related:
  - THM-4600 (run transparency: the case w = w'); THM-4601 (the residual branch; (iv) the ladder compiler, made complete here)
  - THM-4581 (the pair chain); THM-4594 (class-decided certificates)
---

# THM-4603 — run-compressed pair chains

**Scope correction, 2026-10-07:** part (3)'s original bounded-denominator
argument did not prove eventual repetition. Its general classification is
conditional on an authenticated finite lasso. Also, a zero-drift periodic
pair can occupy distinct cycles of equal odd density; it need not be the
same cycle out of phase. The finite transition values remain unchanged.
See [the focused scope controls](../../05-knowledge/results/collatz_run_classification_scope_20261007e.md).
Parts (1) and (2) are not modified by this correction.

## Setting

* T is the Terras map on `Q ∩ Z_(2)`, and U is the odd-to-odd map.
* A pair-chain state is the relation `u = A(v) = 3^k v + e`, with Terras times aligned and debt k. The debt changes by `par(u) − par(v)` per step.
* For a U-word w, let `c_w = B_w/(2^(Σw) − 3^|w|)` be its rational cycle point. This is the point at which w is read; its odd density is `|w|/Σw`.

## Statements

**(1) Drift lemma.**
* Let w, w′ be U-words, `c = c_w`, `c′ = c_(w′)`, and `Λ = lcm(Σw, Σw′)`.
* Suppose odd u, v satisfy `u − c = 3^k (v − c′)` and `v2(v − c′) > jΛ`. Then:
  * for the first jΛ Terras steps, u follows `w^(jΛ/Σw)` and v follows `w′^(jΛ/Σw′)`;
  * `T^(jΛ)(u) − c = 3^(k′) (T^(jΛ)(v) − c′)` with `k′ = k + jΛ(|w|/Σw − |w′|/Σw′)`.
* Both deviations lose one bit of 2-adic precision per Terras step, so the two runs end together.
* THM-4600 is the case w = w′ (zero drift).
* The case of the source on +1 (density 1/2) and the child near the halving point 0 (density 0) gives drift +1/2. This is exactly the collapse of the reset children D ≤ 2 in THM-4601 (ii)/(v), whose debt grows like J.

**(2) Ladder completeness.**
* Suppose two orbits u, v in any pair-chain state first coincide at Terras time s with zero debt.
* Let Z be their first common odd value, and z_u, z_v their last odd values before s.
* Then `z_u ≠ z_v`, and both are U-preimages of Z. Hence for some i ≥ 1, either:
  * `z_v = 4^i z_u + (4^i−1)/3` (child ladder), with letters c and c+2i at z_u and z_v; or
  * the same relation with u and v exchanged (source ladder).
* Writing u's and v's words before z_u, z_v as heads `h_u`, `h_v` (Terras lengths `Σh_u`, `Σh_v`):
  * `|h_v| = |h_u| + k_0`, where k_0 is the debt at the start;
  * the equal-time condition is `Σh_u = Σh_v + 2i` (child ladder) or `Σh_v = Σh_u + 2i` (source ladder);
  * the merge time is `s = max(Σh_u, Σh_v) + 1`, independent of the final letter c.
* For a state anchored at a point c* the head identity reads `F_(h_v)(c*) = 4^i F_(h_u)(c*) + (4^i−1)/3`, or its mirror for source ladders.
* So the ladder compiler of THM-4601 (iv) (anchor +1, debt 3), with every index i and both directions, is complete: its pairs are exactly the absorbing post-run patterns. The same holds for the compiler at −1 of `collatz_reset2_rules_20261007.md` once every i and the source ladders are included.

**(3) Long-run transitions, conditional on a periodic limit witness.**
* Let a state A be given, and let one orbit run near the periodic orbit of a cycle point c′ for L Terras steps (2-adically within `2^(−L)`).
* The other orbit is then near the rational point `A(c′)` (child run) or `A^(−1)(c′)` (source run). For s < L the chain follows the exact limit pair.
* Additionally suppose the rational limit pair has a supplied or independently
  proved finite eventual-periodicity witness. Then, after its transient τ
  independent of L, the marked limit transition has one of the following forms:
  * **absorption at a fixed time**: a periodic head;
  * **re-anchoring**: both limits on the same cycle in phase, so the state at the run end is a fixed anchored state, depending only on L mod the period;
  * **noncoincident zero drift (SHIFT)**: bounded periodic debt, either on the
    same cycle out of phase or on distinct cycles with equal odd density;
  * **nonzero drift**: the limit cycles have different odd densities. The debt
    at the run end is `k_τ + ρ(L − τ) + O(1)`, with ρ their density difference.
* Without such a witness, the outcome is **UNRESOLVED**. A finite search bound
  or a fixed denominator does not supply the missing witness.
* From the universal residual state (3, 1−27) of THM-4601, among the transitions for the 55 cycles (exact table in transitions.out):
  * **child ones-run**: re-anchors at −1 with debt −5 after 12 steps, via the limit pair (−53, −1);
  * **source ones-run**: drift +1/2 after an 8-step transient (debt 6 at step 8), via the child limit 25/27 on the trivial cycle;
  * **child runs of letter 3 or 4**: re-anchor with debts 4 and 9;
  * **positive drift**: child runs of letters 5, 7, 8 and (1,6), (3,5), between +0.27 and +0.47 per step;
  * **negative drift (debt paid down)**: source runs of letters 5, 7 and of (1,5), (3,4), (3,5), between −0.18 and −0.43 per step;
  * **absorption at a fixed time**: source runs of (2,4), (2,6), (1,2,5), (1,5,2), (2,2,4), (2,3,3), and child runs of (1,1,4), (1,1,4,2) (times 17–37).

## Proofs

**(1).** Repeat THM-4600 (2) along both cycles:
* `F_w(u) − c = (3^|w|/2^(Σw))(u − c)`, and likewise for v.
* Over Λ steps the factors are `3^(|w|Λ/Σw)/2^Λ` and `3^(|w′|Λ/Σw′)/2^Λ`. Their ratio is `3^(drift·Λ)`.
* The precision condition is THM-4600 (2) iterated j times.

**(2).** Strictly before the coincidence the orbits differ.
* Just before a coincidence at an odd value w, one orbit is at 2w and the other at the odd value (2w−1)/3.
* At an even value, the common halving chain leads to Z.
* Either way, z_u and z_v are distinct odd U-preimages of Z.
* The odd U-preimages of Z are `(2^a Z − 1)/3` for a in an arithmetic progression with difference 2. Hence `z_(j+1) = 4z_j + 1`, and `z′ = 4^i z + (4^i−1)/3`.
* The equal-time and debt conditions give the head relations. The anchored form follows from the slope identities of THM-4601 (iv).

**(3), repaired.** Read the supplied eventual-periodicity witness and verify
every exact rational transition. The difference of odd-step counts over its
period gives the debt increment, and division by the period gives ρ.
Classify coincidence, phase and density separately. The finite-L shadowing
then transfers this authenticated limit calculation to the guarded long run.
Bounded denominators do not imply bounded numerators or repetition. In fact
`A(v)=v+n+1` and the fixed child limit `v=-1` give the limit pair `(n,-1)`
for any positive odd n; an unconditional eventual-periodicity conclusion
would settle the open no-divergence question for arbitrary integer orbits.
The original universal inference is therefore not supplied by this proof. ∎

For the zero-drift boundary, the actual universal state `(3,-26)` with
child cycle point `7/55` reaches the limit pair `(2/55,7/55)` at time24,
with debt8. The two cycles have periods12 and6, both odd density1/3,
and are disjoint. This is a SHIFT entry in the existing exact table, not
an example of two phases of one cycle.

## Reading

* The pair chain is a debt walk whose long-run increments are differences of odd densities.
  * Re-anchoring keeps the debt bounded.
  * Drift moves it linearly in the run length.
* The heavy tail of THM-4601's post-run problem therefore has an identified source: post-run runs with positive drift. Above all, a long source ones-run is a new Mersenne-like excursion, during which the debt grows by L/2.
* A complete variable-depth rule needs a mechanism that pays such debts. Candidates are negative-drift runs and re-anchoring.
* Ladder completeness turns "find a certificate" into "find a ladder pair", at every anchor.
