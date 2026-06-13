# Good-Cut Height and Projection Polarity

**Session:** opus-2026-05-29-S14

This session fused two threads that looked adjacent but not yet coupled:

- kind-pasteur's projection-defect profile between `G_n/Z_2` and `E_n`;
- the S13 good-cut bucket coordinate `g(τ)`.

The first clean fact is almost tautological once seen, but it was not yet
measured:

> For single-tile flips through n=6, every line with `|Delta g| > 0`
> changes the merged tournament class.

Equivalently, even-only projection defects live in the `g`-neutral layer. This
is exactly what HYP-1764 would predict if `g` descends to `G_n/Z_2`: if the
quotient class is unchanged, `g` cannot move. The computation turns HYP-1764
from a static class-purity statement into a dynamic finite-difference test.

## Range Parity

The second fact is more interesting. Single-tile projection polarity tracks the
parity of the tile range:

- even tile ranges are even-graph biased;
- odd tile ranges are tournament-class biased.

At n=6:

- range 2: defect `-0.0742`;
- range 3: defect `+0.1615`;
- range 4: defect `-0.1094`;
- range 5: defect `+0.1523`.

The fundamental cycle of a tile with range `r` has length `r+1`, so this is
also a parity statement about the cycle-space generator. Odd fundamental cycles
lean toward the even-graph quotient; even fundamental cycles lean toward the
tournament quotient.

That is a surprisingly crisp bridge between the cut-space/cycle-space split
and the merged tournament/even graph projection split.

## Engineering Reading

For Tournament TDA, `g` is not just one scalar feature. The derivative of `g`
under named perturbations acts like a cheap certificate:

- `Delta g != 0` means the perturbation crossed a merged tournament boundary
  in all computed cases;
- `Delta g = 0` is where even-only defects can hide;
- tile range parity separates local cycle-space motion from tournament-class
  motion.

So a practical feature extractor should include not only `g(τ)`, but also
bucket-transition histograms stratified by tile range parity.
