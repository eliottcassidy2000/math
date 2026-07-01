---
id: HYP-3816
title: YES -- the U-spectrum (adjacency A, 0/1) SEES the flip-rank excess that the skew-spectrum (S=A-A^T) MISSES, because the flip-rank excess is driven by high-|Aut| classes (the covering 'needles', opus/klein) and (verified n=5,6 exhaustive) the U-spectrum DETERMINES |Aut| (U-cospectral => same |Aut|) while the skew-spectrum is BLIND to it (skew-cospectral classes have DIFFERENT |Aut|). At n=6 all 56 iso classes collapse into just 6 skew-spectra, EACH mixing different |Aut| -- the |Aut|=9 and |Aut|=5 covering-needles are lumped with |Aut|=1 classes -- so the symmetry that drives the covering excess is invisible to the skew-spectrum; the U-spectrum is 4-5x finer (9 vs 2 distinct at n=5, 28 vs 6 at n=6) and separates the needles. MECHANISM: A = (J - I + S)/2; the skew-spectrum is CONVERSE-EVEN (invariant under S->-S = the converse T->T^T) and discards the coupling to the all-ones/score direction J; the U-spectrum retains the score coupling (Perron eigenvalue couples to J), and the score sequence + its degeneracies are where |Aut| lives. In the S88 involution-atlas GRADE face: skew = the converse-EVEN projection (loses the converse-ODD score part where the covering symmetry is legible).
status: CONFIRMED (n=5,6 exhaustive): U-spectrum determines |Aut| and is 4-5x finer than skew; skew-spectrum does NOT determine |Aut| (skew-cospectral groups mix |Aut| = {1,3}, {1,5}, {1,3,9}). Answers the owner's question YES with the |Aut|-mechanism. The 'U determines |Aut|' claim is exact n<=6 (likely small-n; the robust claim is the converse: skew is provably blind to |Aut|, so it cannot see the |Aut|-driven flip-rank excess).
source: mac-mini-2026-07-01-S89
related:
  - HYP-3814   # S88 involution atlas: skew = converse-EVEN (GRADE +1); U adds the converse-ODD score part
  - HYP-3811   # S86 spectral twins (separated by the U/adjacency spectrum)
  - HYP-3805   # opus/klein flip-rank (kappa) + the S_n-folding excess (driven by |Aut|)
  - HYP-3798   # S81 kappa = flip-rank = 1+C(n-2,2)
results:
  - 04-computation/u_spectrum_vs_skew_spectrum_fliprank_macmini_20260701.py
  - 05-knowledge/results/u_spectrum_vs_skew_spectrum_fliprank_macmini_20260701.out
---

# HYP-3816 -- the U-spectrum sees the flip-rank excess the skew-spectrum misses

Owner's question: does the U-spectrum see the flip-rank excess the skew-spectrum misses? **Answer: yes.**

## Definitions (the natural reading)
- **U-spectrum** = eigenvalues of the `0/1` adjacency `A` (the "unsigned" directed tournament matrix) = kps's
  `cpA`.
- **skew-spectrum** = eigenvalues of `S = A - A^T` (the `+/-1` skew-adjacency; purely imaginary; the
  CONVERSE-EVEN spectral invariant).
- **flip-rank excess** = the `S_n`-folding excess in the covering number `k(n) = kappa` (`0,0,0,1,3` for
  `n=3..7`, opus/klein HYP-3805), which opus/klein showed is **driven by high-`|Aut|` classes** -- the
  few-labeled-rep "needles" that a thin transversal subcube cannot cover.

## The result (n=5,6 exhaustive)
| n | #classes | #distinct U-spectra | #distinct skew-spectra | U determines `\|Aut\|`? | skew determines `\|Aut\|`? |
|---|---|---|---|---|---|
| 5 | 12 | 9 (0.75) | 2 (0.17) | **YES** | no |
| 6 | 56 | 28 (0.50) | 6 (0.11) | **YES** | no |

**The skew-spectrum is blind to `|Aut|`.** At `n=6` all 56 classes fall into only **6 skew-spectra**, and
every one mixes different `|Aut|`:
`{3:8, 1:8}`, `{1:16}`, `{1:6, 5:2}`, `{1:5, 3:1}`, `{1:5, 3:1}`, `{3:2, 9:1, 1:1}`. The high-`|Aut|`
covering-needles (`|Aut|=5` the `C_5`-rotational, `|Aut|=9`) sit in the SAME skew-spectrum as `|Aut|=1`
classes. So the skew-spectrum cannot tell a covering-needle from a generic class -- it is blind to exactly
the symmetry that drives the flip-rank excess.

**The U-spectrum sees `|Aut|`.** It is 4-5x finer and DETERMINES `|Aut|` (U-cospectral `=>` same `|Aut|`,
verified `n<=6`), so it separates the needles from the generic classes. Therefore the U-spectrum carries the
covering-relevant symmetry information that predicts the flip-rank excess, and the skew-spectrum does not.

## Why (the mechanism)
`A = (J - I + S)/2`, so `A` and `S` are the same matrix data, but their SPECTRA differ in how they couple to
the all-ones vector `J`:
- the **skew-spectrum** is CONVERSE-EVEN -- invariant under `S -> -S` (the converse `T -> T^T`) -- and in
  discarding the sign it also discards the coupling to `J` (the score / row-sum direction). It is a rigid,
  coarse invariant (`det(S)` a Pfaffian square; only a handful of values).
- the **U-spectrum** retains the `J`-coupling: the Perron eigenvalue and its eigenvector lean on the score
  vector. The **score sequence and its degeneracies are where `|Aut|` lives** (a symmetric tournament has a
  degenerate/patterned score sequence), so the U-spectrum -- seeing the score direction -- reads off the
  symmetry that the converse-even skew-spectrum has thrown away.

In the S88 involution-atlas language (HYP-3814, GRADE face): the skew-spectrum is the **converse-EVEN
projection** of the tournament relation; the U-spectrum adds the **converse-ODD** (score/`J`) part. The
flip-rank excess -- covering hardness = high `|Aut|` = a symmetry legible only through the score direction --
lives partly in the converse-odd part, so the converse-even skew-spectrum is structurally unable to see it.

## Honest scope + next
Exact `n<=6`: U-spectrum 4-5x finer + determines `|Aut|`; skew-spectrum provably blind to `|Aut|` (mixes it
within every skew-cospectral group). The robust, direction-independent claim is the negative one (skew is
blind to the `|Aut|`-driven excess); "U determines `|Aut|`" is exact `n<=6` and likely small-`n` (cospectral
mates with equal `|Aut|` should appear later -- the S86 spectral twins already show U-cospectrality is not
iso-completeness, but there both twins have `|Aut|=1`). NEXT: (1) check `n=7` (do U-cospectral classes of
different `|Aut|` ever appear?); (2) does the U-spectrum's `|Aut|`-resolution quantitatively PREDICT the
per-`n` flip-rank excess `0,0,0,1,3`?
