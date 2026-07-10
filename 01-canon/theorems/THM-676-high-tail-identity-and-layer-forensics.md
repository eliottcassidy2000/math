---
id: THM-676
title: The high-tail identity and the moment-layer forensics — for odd D, B_D = LM − Σ_{c>D} C(c−1,D)·hist[c] EXACTLY (partial-binomial identity), so the entire Bonferroni error is the HIGH-COVERAGE mass (c ≥ 6 for the quintic — the discrete apex-7 object: points killed by ≥6 classes at once); layer forensics show the generic branch survives by CROSS-LAYER cancellation (layers are coupled cumulants; per-layer deviations ±0.1–0.25 with total deficit only 0.02–0.08) while the coherent branch has cancel-ratio → 1.0 at d ≥ 3 (deviations 2.3/5.4/8.5 — pure constructive coherence); the one-sided margins LM/penalty = 1.6–81× on generic vs 0.006 on the near-dilation
status: MIXED. PROVED: the identity (f_D(c) = Σ_{d≤D}(−1)^d C(c,d) = (−1)^D C(c−1,D) for c ≥ 1 — standard partial alternating binomial sum — hence B_D = hist[0] − Σ_{c>D} C(c−1,D) hist[c] for odd D; verified exactly on every (S,q) tested); corollaries: the odd-depth ladder's monotonicity/exactness in closed form, and the cheapest certificate form (hist[0] itself = LM, consumed by kps-S117's mreach_ge_of_blocked_lt). MEASURED: the forensics and margins tables. HONESTLY OPEN: the a-priori generic-branch bound — the identity RELOCATES the two-sided moment control (L_2, L_4 lower bounds = the forced-independence direction) into the high-tail form but does not remove it; cross-layer cancellation is systematic (coupled layers + alternating signs) but proving it is the same interference problem. The generic a-priori bound is THE last mathematical item of the modular route.
source: klein-2026-07-09-S213 (HYP-5770; owner-directed "interference-aware generic bound / finish LRC(14)")
depends_on:
  - THM-671   # the B5 certificate (its error now in closed form)
  - THM-675   # the composed dichotomy; coherence law
  - kps-S117 LRCLedgerConsumer  # mreach_ge_of_blocked_lt — the Lean consumer of hist-counts
related:
  - THM-604 (iid values), THM-673 (aggregated skeleton), OPEN-Q-108 (apex-7 — the high-coverage mass is its discrete analog)
---

# THM-676 — the high-tail identity and the layer forensics

## 1. The identity (proved)

For c ≥ 1 and any D: Σ_{d≤D} (−1)^d C(c,d) = (−1)^D C(c−1, D). Hence for odd D,
summing over p ≠ 0:

> **B_D(S,q) = LM(q) − Σ_{c>D} C(c−1, D)·hist_q[c].**

The Bonferroni error is EXACTLY the high-coverage mass: points p killed by more
than D classes simultaneously. For D = 5 the penalty starts at c = 6 — the
discrete cousin of the apex-7 threshold (≥ 6–7 danger sets overlapping at once).
Corollaries: the odd ladder B5 ≤ B7 ≤ … ≤ B13 = LM with differences = re-weighted
high tails; and the terminal certificate is hist[0] itself (exact LM — one
histogram pass, native_decide-shaped, consumed by `mreach_ge_of_blocked_lt`).

## 2. Layer forensics (measured; the binomial-histogram identity L_d = Σ_p C(C(p), d))

Generic covering instances (V = 120–280): per-layer deviations from iid are
±0.05–0.25 with cancel ratios 0.13–0.81, and the B5-deficit decomposition shows
CROSS-LAYER cancellation (e.g. +0.155 at d=2 against −0.126 at d=3): the layers
are coupled cumulants of one overlap structure, and the alternating signs damp
their common movement. Total deficits 0.02–0.08 < 0.1221 ⟹ avg-B5 > 0.
The near-dilation (V=260): cancel ratio 0.97–1.00 at d ≥ 2 with deviations
+0.30/+2.26/+5.42/+8.49 — fully coherent, no cancellation, the catastrophe in
closed view. Practical branch statistic: the support-3 relation count at
height ≤ 3 (generic: 52–132; @91-comb: 104; near-dilation: 672).

## 3. The one-sided margins (measured)

avg LM/q vs avg penalty/q over q ∈ (V, 2V]: 81.4× (adv-worst-120), 1.6×
(adv-200), 4.4× (adv-280), 7.7× (@91), 0.006× (near-dilation). The generic
inequality LM > penalty holds with room; the object is comfortably true.

## 4. Honest scope — what remains for LRC(14)

The identity does NOT remove the two-sided moment control: lower-bounding LM
through any truncation still needs L_2, L_4 from below (forced independence —
THM-604's hard direction) or equivalently an upper bound on the high tail that
itself engages the same overlaps. Cross-layer cancellation is systematic but
proving it IS the interference problem. **The a-priori generic-branch bound
remains the single open mathematical item of the modular route**, now in its
sharpest form: *show the coverage histogram of a low-coherence covering set
stays sub-binomial at c = 0 and c ≥ 6 simultaneously, for some q ∈ (V, 2V]*.
Realistic paths to a COMPLETE LRC(14) proof from tonight's map:
 (i) a genuinely new idea for the generic histogram bound (the item above);
 (ii) brute-scale: native_decide enumeration (kps-S115 pattern) pushed as far
     as compute allows + the coherent branch by rigidity + (i) only for the
     asymptotic tail;
 (iii) the τ-line composed route (LEM-012 + boxeph/death-star's composed
     realization + the mid-band via boxeph-S3's comb-discrepancy D_m, measured
     ≤ 4% vs needed 5.7% — a 1.4× gap, the closest analytic margin anywhere).
All three are live; none is closed tonight. The certificate infrastructure
(dispatch + ledger consumer + histogram) is complete and Lean-ready regardless.

## Files

`04-computation/lrc14_moment_layer_forensics_klein_S213.py` (+ `.out`),
`05-knowledge/results/lrc14_hightail_penalty_klein_S213.out` (identity check
12/12 instances-wide + margins).
