---
id: THM-755
title: The CAPPED-ENVELOPE theorem -- disc_v <= 4 r_P |G'_P|/(pi v) + 2|G'_P|^2 in the moderate regime (each Fourier coefficient obeys BOTH the trivial cap |c_m| <= |G'| and the jump envelope r/(pi m); split at the crossover), hence klein-S306's (H) HOLDS for every v > v* = r_P/(pi |G'_P|) -- unconditional, no alignment hypothesis; the (H)-band shrinks to the per-core finite interval (maxP, v*], nonempty only for tight cores (band edges: deep-well core 112 vs isolated 156; loose core 52 vs 182; spread core 644 vs 4771)
status: PROVED (six-line envelope split, below) + VERIFIED-EXACT (zero violations across four cores x six v-values each, exact Bernoulli disc vs the bound; band edges fire as computed)
source: opus-2026-07-14-S289 (owner prompt: prove the (H) inequality with the perspective frame)
depends_on: []
related:
  - klein-S306 / THM-753 (the one-step peel PROVED modulo (H): this theorem discharges (H) above v* per core; the skeleton's residual is now the finite band)
  - THM-731/732 (the certificate this feeds; the crude bound is the l* = 1 endpoint of the same split)
  - THM-752 (the fine-comb witness -- the band's other closer from below), kps exact-Q disc (per-body closer inside the band)
---

# THM-755 -- the capped-envelope theorem and the (H)-band

## Statement

Let G'_P (measure |G'|, r components) be the 1/14-good set of a core P, and for v >= 1 let
disc_v = sum_{l != 0} |c_{lv}|^2 (the THM-731 grid discrepancy).  Then for every integer
l* >= 1:

>  disc_v <= 2 l* |G'|^2 + 2 r^2 / (pi^2 v^2 l*),

and with l* = ceil(r/(pi v |G'|)) (the envelope crossover):

>  **disc_v <= 4 r |G'| / (pi v) + 2 |G'|^2**   (moderate regime r >= pi v |G'|),
>  disc_v <= r^2/(3 v^2)                        (l* = 1: the THM-732 crude bound recovered).

COROLLARY ((H) above the band edge): since (H) asks disc_v < 6|G'|^2, and
4 r |G'|/(pi v) + 2|G'|^2 < 6|G'|^2  <=>  v > r/(pi |G'|) =: v*,

>  **(H) holds for every v > v* = r_P/(pi |G'_P|)** -- unconditionally (no alignment, no
>  isolation).  The open part of (H) shrinks to the per-core FINITE band v in (maxP, v*].

## Proof

(i) |c_m| <= |G'| (triangle inequality on the integral).  (ii) integrating by parts,
c_m = (2 pi i m)^{-1} sum_p sigma_p e(-m x_p) over the 2r jumps: |c_m| <= r/(pi m).
(iii) split the sum at l*: the first 2l* terms take cap (i), the tail takes (ii) with
sum_{l > l*} l^{-2} <= 1/l*.  (iv) substitute the crossover l*; each half equals
2 r |G'|/(pi v) (up to the +2|G'|^2 rounding of ceil).  QED.

The perspective reading: the trivial cap is the ORIGIN's view (a set of measure |G'| cannot
correlate more than |G'| with any clock); the jump envelope is the SPOKE view (2r boundary
spokes, each 1/(2 pi m)); the theorem is the pair of perspectives spliced at their crossover
-- neither alone reaches past the crude bound, together they eat the moderate band.

## Verified band edges (exact Bernoulli disc, zero violations)

| core | \|G'\| | r | v* = r/(pi\|G'\|) | old isolated threshold 13 maxP | band reduction |
|---|---|---|---|---|---|
| deep-well {1..12} | 0.03410 | 12 | **112** | 156 | 1.4x |
| residue {1..11,13} | 0.01216 | 4 | **105** | 169 | 1.6x |
| loose {2..12,14} | 0.08512 | 14 | **52** | 182 | 3.5x |
| spread (klein-stall core) | 0.13743 | 278 | **644** | 4771 | 7.4x |

Inside the residual band the true disc is 10-100x below 6|G'|^2 everywhere measured; the
closers there are kps's exact-Q disc (per body), THM-752's fine-comb witness, and klein's
iterated peel -- all already operational.  The band is nonempty only for TIGHT cores, whose
near-AP rigidity (THM-724/730) is exactly the proved structure.

## Position

klein-S306's THM-753 skeleton ("LRC(14) covering = one inequality atop the proved finite
skeleton") now needs (H) only on the finite per-core band (maxP, v*] -- an exact-computable
edge.  The one inequality has become: finitely many v per core, each decidable.

## Files

05-knowledge/results/lrc14_capped_envelope_thm755_opus_S289.out (exact verification, four
cores, band-edge calibration).
