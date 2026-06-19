---
id: HYP-2640
title: LRC(14) correction rank scaling - Freiman relation rank versus signed visible rank
status: CLAIMED; computation pending
source: codex-2026-06-19-S28
depends_on:
  - HYP-2639
  - HYP-2638
  - HYP-2637
  - HYP-2635
  - HYP-2606
related:
  - HYP-2607
  - HYP-2604
  - THM-531
  - OPEN-Q-108
---

# HYP-2640 - LRC(14) Correction Rank Scaling

## Claim Being Tested

The correction term

```text
Corr_y(E) = L_y(E) - L_y^inf(k)
```

should scale with relation structure, but not with a naive scalar rank.

There are two different rank notions that must not be conflated:

1. The Fourier offset-relation lattice rank for the one-parameter orbit
   `sum n_i e_i=0`.  For integer offsets this rank is essentially fixed
   (`k-2` after the pinned zero coordinate); its covolume and short vectors
   vary, but the rank alone cannot explain the correction.
2. The Freiman / summand relation rank: the rank of pair-sum and bounded
   weighted-summand relations among the finite set elements.  This rank varies
   from `0` on dissociated rows to `k-2` on AP rows, and is the natural
   "pocket dimension" complement `rank = k-1-d`.

The test is whether the positive correction is bounded by a typed rank packet:

```text
total Freiman relation rank
+ visible fold rank
+ hidden balanced shell rank
+ signed cancellation / parity type
```

HYP-2639 predicts that total rank alone is too coarse: AP and shifted AP can
share scalar energy, but differ in observer-visible fold payload and LRC
hardness.  HYP-2640 will measure this in the sector correction itself.

## Planned Computation

Add a scout that:

1. computes exact `L_y`, `p0=meas(S7)`, and independent baselines;
2. computes Freiman relation rank from pair-sum equality vectors;
3. splits that rank into visible-fold and hidden balanced-shell spans;
4. optionally adds bounded weighted relation rank for small coefficients;
5. groups exact bounded banks by these rank packets and records max correction.

The desired theorem shape is:

```text
high total rank + high visible signed rank -> near-AP finite table;
high total rank + low visible signed rank -> large signed cancellation / safe;
low total rank -> independent-limit or peel.
```

## Assumption Challenge

The challenged assumption is that relation-lattice rank is a single scalar.
The Fourier rank, Freiman rank, visible shell rank, and signed reciprocal rank
are different quotients.  The useful tournament vertices are therefore proof
obligations:

```text
signed_visible_rank
> total_freiman_rank
> hidden_balanced_rank
> correction_per_rank
> covolume_short_vector_data
> raw_additive_energy
> raw_runner_vertices
```
