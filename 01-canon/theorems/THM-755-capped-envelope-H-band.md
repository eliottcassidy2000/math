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

## Lean (opus-S291/S292) -- COMPLETE, sorry-free, kernel-pure

LRCClosedBudget.lean now contains the FULL spectral THM-755, all audited
[propext, Classical.choice, Quot.sound]:
- `tail_inv_sq_le_sub` -- the telescoping Ioc tail (strengthened 1/m - 1/N form);
- `capped_envelope_kernel` -- the abstract splice: origin cap |c l| <= A + spoke envelope
  l|c l| <= B give Sum (c l)^2 <= m A^2 + B^2/m;
- `fourierCoeff` -- the coefficient of an n-interval family (interval integrals of
  exp(-2 pi m x I));
- `fourierCoeff_norm_le_measure` -- the ORIGIN CAP instantiated (norm-of-integral <= measure);
- `fourierCoeff_norm_le_env` -- the SPOKE ENVELOPE instantiated (closed-form antiderivative
  via `integral_exp_mul_complex`: each interval <= 2/(2 pi m); total n/(pi m));
- `spectral_thm755` -- the composition: Sum_{l=1}^{N} ||c(l v)||^2 <= m G^2 + (n/(pi v))^2/m.

THE POISSON BRIDGE IS NOW FORMALIZED TOO (opus-S293, sorry-free kernel-pure): B2R (the
periodized second Bernoulli), B2R_add_int / B2R_neg, both Gauss sums, raabe_shift_one /
raabe_shift_int (cyclic reindexing), raabe_base (fundamental-cell algebra), **raabe_B2** (the
multiplication formula Sum_{i<v} B2({y+i/v}) = (1/v) B2({vy}) -- the finite Poisson atom), and
**grid_deficit** (integral-free E-linearity: h = C + Sum w_r B2({. - y_r}) has grid deficit
(1/v^2) Sum w_r B2({v y_r})).  THE DECOMPOSITION CHAIN IS NOW ALSO FORMALIZED (opus-S294): the exact form derived and
tent-verified -- A(tau) = |G|^2 + (1/2) Sum_{p,q} sigma_p sigma_q B2({tau + x_q - x_p}) --
with, in Lean (sorry-free, kernel-pure): grid_deficit GENERALIZED to arbitrary Fintype
indices; acorrModel (the B2-model) and discB (THM-732's exact Bernoulli form) defined;
**discB_eq_grid_deficit PROVED** ((1/v) Sum_i model(i/v) - C = discB -- the geometric
grid-mean identity for the model); and two machine-checked rational instances of the
model-equals-autocorrelation identity (tent interior: [0,1/3] at tau = 1/4 gives 1/12; wrap
regime: tau = 5/6 gives 1/6).  THE LAST BRICK IS LAID (opus-S295): **pair_overlap_B2 PROVED in Lean, sorry-free,
kernel-pure** -- the single-pair overlap identity in its normalized one-parameter form:
max 0 (min l1 (l2-t)) + max 0 (min l1 (1-t+l2) - (1-t)) = l1 l2 + (1/2)(B2({t}) - B2({t-l2})
- B2({t+l1}) + B2({t+l1-l2})), for l1, l2 in [0,1], t in [0,1).  Six-leaf case cascade over
the fract/overlap thresholds; every leaf closes by ring after floor-computation resolution.
(The S295 pre-verification caught a reversed orientation convention before the formalization
-- 8392 exact configurations, 0 mismatches on the corrected form; the double sum over ordered
pairs symmetrizes either orientation, which is why the diagonal instances could not detect
it.)  The (H)-edge chain is now machine-checked END TO END at the identity level: kernel ->
Fourier envelopes -> spectral -> Raabe -> grid deficit -> discB -> single-pair overlap.  THE ASSEMBLY IS DONE (opus-S296): acorr_eq_model PROVED (Acorr -- the geometric
autocorrelation as the sum of pair overlaps -- equals acorrModel over the jump index
Fin n + Fin n, via B2R_fract_add de-fracting, the prodComm orientation swap, and the four-block
Fintype.sum expansion), and the CAPSTONE geometric_disc_eq_discB PROVED: for any interval
family, the v-grid mean of the geometric autocorrelation minus |G|^2 equals THM-732's exact
Bernoulli jump-pair discrepancy.  **THE GEOMETRIC THM-755 IS COMPLETE IN LEAN**: every
statement from the raw overlap geometry to the spectral band edge is machine-checked,
sorry-free, kernel-pure.

## Files

05-knowledge/results/lrc14_capped_envelope_thm755_opus_S289.out (exact verification, four
cores, band-edge calibration).
