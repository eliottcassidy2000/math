---
id: HYP-3539
title: The n=6 axis-aligned face-compression overhead (7 free tiles where information allows 6) is the cost of resolving the R-odd (complement-signed) coordinate
status: REFUTED (klein-2026-06-29-S1) — covering only the complement-MERGED classes also needs a 7-face at n=6, so the overhead lives in the R-even bulk, not the R-odd coordinate
source: klein-2026-06-29-S1
depends_on:
  - THM-584
related:
  - HYP-3538
  - HYP-2685
results:
  - 04-computation/r-eigenspace-merged-compression.py
  - 05-knowledge/results/r-eigenspace-merged-compression.out
  - 04-computation/iso-class-compression-faces.py
---

# HYP-3539 — Is the n=6 compression overhead the R-odd coordinate? (REFUTED)

## Background

The tiling cube has `m = C(n-1,2)` free tiles. The minimum-dimension **axis-aligned face** (fix some
tiles, vary the rest) that still hits every iso class is:

| n | info bound ⌈log2 A000568⌉ | min covering face | tight? |
|---|----|----|----|
| 4 | 2 | 2 | yes |
| 5 | 4 | 4 | yes |
| 6 | 6 | **7** | **no (+1)** |

(Exhaustive over all `C(m,k)` faces and all fixings; `iso-class-compression-faces.py`.) So at n=6 the
clean information-tight compression first fails by one bit.

## Claim tested

The +1 overhead is the cost of distinguishing the two members of a non-self-complementary (NS)
complement-pair — i.e. of resolving the **R-odd / `eps=-1`** coordinate (THM-584). Prediction: if you
only need to cover the complement-**merged** classes (where `+` and `−` collapse, so the R-odd bit is
irrelevant), the min face should drop to the information bound (tight).

## Result: REFUTED

`r-eigenspace-merged-compression.py` computed the min covering face for the merged class set:

| n | ALL classes: info / min face | MERGED classes: info / min face |
|---|----|----|
| 4 | 2 / 2 | 2 / 2 |
| 5 | 4 / 4 | 4 / 4 |
| 6 | 6 / **7** | 6 / **7** |

At n=6 the **merged** classes (34 of them, info bound 6) **also require a 7-face**. Collapsing the
R-odd coordinate does not remove the overhead. Therefore the overhead is a packing/covering deficiency
in the **R-even bulk** (the merged metagraph itself), not the cost of the signed R-odd coordinate.

## Why the failure is informative

Read against HYP-3538 this is the *consistent* outcome, not a surprise: the bulk is the R-even side,
and the compression hardness turned out to live in the bulk. The R-odd eigenspace is small (dim 22 of
56 at n=6, dim 2 of 12 at n=5) and cheap to carry; the expensive object is packing the R-even classes
into axis-aligned faces. The overhead is "Brouwer-flavored" (covering geometry), not
"Borsuk–Ulam-flavored" (signed obstruction).

## Open residue

Whether a general **GF(2)-affine** slice (a coset of a linear subcode, not just coordinate-fixing)
removes the n=6 overhead is still open — faces are a special case of affine slices, so affine ≤ face;
a structured search (not the weak randomized one already tried, which can't sample the rare covers) is
the next step. Logged as a lead in INVESTIGATION-BACKLOG.
