---
id: HYP-9802
title: "Literal span-3 growth certificates: (i) every cheap primitive cycle of period at most 6 letters (growth below 1.72 per letter) contains a rank-one factor of at most 9 letters, so a targeted suffix dictionary of 2,527 strings cuts them; (ii) letter-level cone certificates cannot cap the potential by H on raw generators at any rate, because one letter can split a path fragment, so they need the bounded-growth lemma H(W l) <= K H(W) (observed K = 39); (iii) the remaining slow bottlenecks have periods of 7-12 letters and need rank-adaptive memory"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-frame, wave 39 of chanpak_waves_20261113. UPDATE (klein-2026-11-15-ribbon, THM-5243 /
  HYP-9803): the lemma of (ii) is exactly Z_l(W) <= (K - 2) H(W) for the bridged two-covers (PROVED decomposition); it holds with
  K = 39 for every word with H(W) <= 100 (closed set of 100,203 states) and at all lengths <= 8 (P105 confirmed at length 8); the constant is sharp at
  unboundedly large H (sealing by three empty letters); the annihilator LP and the coordinate-cap cone both fail to certify it
  (PROVED obstructions), so a proof must use fragment ages.
  (i) FINITE-EXACT: 3,343 cheap cycles of period <= 6, and 1,303 cycles of period <= 7 with growth below 1.3 (each
  has a rank-one factor of <= 6 letters). Ranks are exact (nonnegative integer columns, primitive-direction test);
  growth is floating-point.
  (ii) The obstruction is PROVED in outline and FINITE-EXACT (the LP with the cap is infeasible for every tested rate
  >= 1.2 on a 6,262-state automaton). The lemma is CONJECTURE: the one-vertex ratio maximum is 3, 5, 11, 23, 39, 39 at
  word lengths 2-7, and the injection proof is only a sketch.
  (iii) is the concurrent finding of klein-2026-11-15-ocf (LITERAL_CONES.md: cheap cyclic code words needing up to 12
  code units of memory; rank as a potential), restated at the letter level.
  The literal certificate itself is NOT built. HYP-9742's literal span-3 case remains open.
session: mac-mini-2026-11-15-frame
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 39)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w39_lit_setup.py, w39_automaton.py, w39_automaton_extra.py,
    w39_cheap_cycles.py, w39_cone_lp.py, w39_cone_lp_growth_only.py, w39_verify_literal_cert.py, w39_receipts.txt
related:
  - HYP-9742, THM-5212 (E), chanpak_mod8_synthesis_span3_20261115/SPAN3_LITERAL.md,
    chanpak_apex_ocf_literal_cones_20261115/LITERAL_CONES.md
---

# HYP-9802 - what a literal span-3 growth certificate needs

**Setting.**
- Literal span-3 words: a first vertex, then letters ℓ ⊆ {1, 2, 3}.
- The 119-state profile transfer A_ℓ (THM-4880). After the first vertex the minimal representation has dimension 61
  (62 with the empty-prefix state).
- Bottleneck-free means every factor has rank ≥ 2 there.

**Target.** H ≥ c·ρⁿ for all bottleneck-free words, with ρ > 5^(1/3) (HYP-9742). THM-5212 (E) proves ρ = 1.7699 for
words over the 13 floor units.

**The cone method.** Use states given by a deterministic suffix automaton, with memory = the longest suffix in a
dictionary closed under deleting the last letter. The cones are images of the nonnegative orthant under the remembered
suffix, and potentials are linear functionals. The LP needs:
- growth ρ along every edge on every raw generator;
- a final comparison with H.

## (i) Short cheap cycles are cut within nine letters (FINITE-EXACT)

| period ≤ | growth below | cheap primitive cycles | shortest rank-one factor | factor dictionary |
|---|---|---:|---|---:|
| 4 | 1.72 | 104 | ≤ 9 letters | 411 |
| 5 | 1.72 | 690 | ≤ 9 letters | 1,463 |
| 6 | 1.72 | 3,343 | ≤ 9 letters | 2,527 |
| 7 | 1.30 | 1,303 | ≤ 6 letters | 164 |

Every cheap cycle in these ranges is a bottleneck, and its shortest bottleneck factor is short. The 13 floor units
are among them: their powers become rank one at 2–5 repetitions.

Memory must be measured in letters, for every phase. The concurrent code-unit census measures powers from code-unit
phases. That is why it reports bottleneck lengths of up to 12 code units for periods up to 12 letters.

## (ii) The readout obstruction and the bounded-growth lemma

**Obstruction (PROVED in outline).** Suppose the final comparison is imposed as φ_s ≤ H on the raw generators of each
cone, as in the 13-unit certificates. Then no rate ρ > 0 is feasible at the letter level:
1. A single-fragment raw configuration with a young tail is split into two fragments by any letter containing that
   gap. Such a letter gives readout 0, so φ ≤ 0 is forced on the predecessor generator.
2. After one letter, every single-fragment configuration has a young tail.
3. Multi-fragment generators have readout 0.
4. So φ ≤ 0 on every generator, which contradicts φ ≥ 1 on initial states.

Floor units escape this because a whole unit always re-merges fragments.

The LP confirms it: the 6,262-state automaton is infeasible at ρ = 1.2, 1.4, 1.55, 1.65 and 1.705.

**Fix (the lemma, CONJECTURE).** For span-3 ribbon words, H(W ℓ) ≤ K·H(W) for an absolute K.
- *Observed.* The maximum ratio over all words of length 1–7 is 1, 3, 5, 11, 23, 39, 39. So K = 39 is the
  conjectured value.
- *Proof sketch.* Remove the new vertex from a Hamiltonian path of T + v.
  - Either a Hamiltonian path of T remains, at most w + 2 = 5 insertion slots per path,
  - or two paths remain whose junction vertices lie among the last 2w = 6 vertices. These are mapped with bounded
    multiplicity to Hamiltonian paths of T minus a bounded young zone, which is at most H(T) by monotonicity. The
    merging step needs care and is not written out.
- *Why it suffices.* The potential is a combination of Hankel coordinates (O X)_i = H(W w_i) for continuation words
  w_i of bounded length. So the lemma gives |Φ| ≤ C·H on actual states, and no readout cap is needed in the LP.

## (iii) What remains

- *Slow bottlenecks of longer period.* Cheap cycles of period 7–12 letters can need much longer memory before their
  first rank-one factor (concurrent LITERAL_CONES.md: up to 12 code units).
- *The growth-only LP.* On the 6,262-state automaton with the period-≤ 6 dictionary it did not terminate in 25 minutes
  at ρ = 1.2, so it is inconclusive.
- *The promising route.* klein's rank potential: states (short tag, current rank), a bounded loss per rank drop.
  Combined with the lemma of (ii), it gives a letter-level certificate scheme with no readout cap and finitely many
  rank drops.

**Predictions.**
- *P105.* The one-vertex ratio maximum stays 39 for all word lengths ≥ 6.
- *P106.* Every cheap primitive cycle of period ≤ 8 letters with growth below 1.72 has a rank-one factor of at most 12
  letters.
