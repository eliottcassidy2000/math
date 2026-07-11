---
id: THM-693
title: The finite-V two-scale witness — for S_V = P ∪ {V−e : e ∈ E} with test data (q* ∈ [8,13], a, missed class c) satisfying (p·a) % q* ≠ 0 for every small speed and (c − e·a) % q* ≠ 0 for every co-offset, the EXPLICIT multiplier w = a·V + ((c − a·V) mod q*) is STRICTLY LIVE at modulus q*·V whenever V > 2184 and V > 168·e per co-offset — the residues are r_p·V + p·Δ (small) and s_e·V − e·Δ (cluster), both strictly inside (q*V/14, 13q*V/14); composed with kps's chain: an explicit rational witness time w/(q*V) for every large-V two-scale family. FORMALIZED SORRY-FREE (LRCTwoScaleWitness.lean, kernel-pure, green on first compile) — the test-point program made constructive at finite V: no limits, no transfer, no measure theory
status: PROVED AND FORMALIZED (the Lean file IS the proof: small_speed_band + cluster_speed_band + twoScale_strictlyLive + twoScale_strictWitness, all [propext, Classical.choice, Quot.sound]; demo family {1,…,5} ∪ {V−e : e ∈ {0,1,2,3,5,8,13,21}} at V = 10000 witnessed at (q, w) = (130000, 10001) through the theorem with decide-checked hypotheses; Python cross-check all 13 residues strictly in-band with the tight cases exactly at s_e ∈ {1, q*−1} as the proof predicts).
source: klein-2026-07-10-S244 (HYP-5945; the named finite-V composition from S243)
depends_on:
  - THM-690/691/692 (the test-point structure: the missed class c is pigeonhole_missed_class's output)
  - kps LRCStrictRuler (StrictlyLive, strictWitness_of_strictlyLive)
related:
  - LRCMeasureTransfer (S242 — the certificate route; THIS theorem is the direct route, no interval needed)
  - LRCTestPointCore (S243 — supplies the hypotheses: residue_in_range, pigeonhole_missed_class)
---

# THM-693 — the finite-V two-scale witness

## Statement

Let S_V = P ∪ {V − e : e ∈ E} (13 speeds), q* ∈ [8,13], a ∈ [1, q*−1], and
c a residue with:
- (p·a) % q* ≠ 0 for every p ∈ P (the q*-test-point P-side; automatic from
  coprimality when q* ∉ P — THM-691(A)/LRCTestPointCore.qstar_p_nonzero);
- (c − e·a) % q* ≠ 0 for every e ∈ E (c = a missed class of E mod q* —
  supplied by pigeonhole_missed_class whenever |E| < q*);
- V > 2184 and V > 168·e for every e ∈ E.

Then with Δ = (c − a·V) mod q* ∈ [0, q*) and **w = a·V + Δ**, the pair
(q, p) = (q*·V, w) is a strictly-live ruler for S_V:

> (p·w) mod q*V = r_p·V + p·Δ,   r_p = (p·a) mod q* ∈ [1, q*−1];
> ((V−e)·w) mod q*V = s_e·V − e·Δ,   s_e = (c − e·a) mod q* ∈ [1, q*−1];

both strictly inside (q*V/14, 13·q*V/14): the small side because
14·r_p ≥ q* + 1 beats the p·Δ ≤ 156 perturbation once V > 2184; the cluster
side because 14·s_e·V ≥ (q*+1)·V beats 14·e·Δ ≤ 168·e once V > 168·e.
Hence (kps's strictWitness_of_strictlyLive) S_V has a strict witness at the
explicit rational time w/(q*·V). ∎ (The Lean file is the complete proof.)

## What this closes

The witness is DIRECT: no limit measure, no C/V transfer, no certificate
interval — the test-point structure (P-side nonvanishing + one missed
class) hands the finite family an explicit modulus and multiplier. This is
the constructive finite-V form of THM-690/691/692: those theorems' μ∞ > 0
said witnesses exist in the limit; this one writes the witness down. The
thresholds are absolute constants (V > 2184) and per-co-offset linear
(V > 168·e); below them, the bounded-window machinery (THM-686, banks)
applies. For the wall on two-scale slices the chain is now fully
constructive end to end: pigeonhole → missed class → w = aV + Δ →
StrictlyLive → StrictWitness → lonely.

## Formalization

`04-computation/lean/TournamentH7/TournamentH7/LRCTwoScaleWitness.lean`
(kernel-pure, root-wired, green on first compile): `small_speed_band`,
`cluster_speed_band` (the two residue computations — Int.ediv_add_emod
decompositions identified via Int.add_mul_emod_self_left, band bounds by
nlinarith), `twoScale_strictlyLive` (THE THEOREM), `twoScale_strictWitness`
(the kps composition), `demoTwoScale_strictlyLive/strictWitness` (the
(130000, 10001) demo through the theorem). Python cross-check:
`05-knowledge/results/lrc14_twoscale_witness_demo_klein_S244.out`.
