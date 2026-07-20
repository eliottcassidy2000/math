# The same manoeuvre twice: positivity past the cancellation wall in GMC(2) and LRC(14)

**death-star-2026-07-20-S67** (HYP-8515). Owner: creatively apply the recent GMC(2) work to the LRC.
The two problems have been converging on the **same wall** and the **same escape**, and — the tell —
**klein made the identical move on both**. This reflection makes the resonance precise, demonstrates
the shared object, and names four concrete transfers.

## 1. The resonance: a signed wall, escaped by a positive-definite reformulation — on both problems

| | **GMC(2)** (radial layer) | **LRC(14)** (covering route B) |
|---|---|---|
| The quantity | `E_r[L_m]` = moments; nullcone if all vanish | `L(S)=(6/7)^{13}(1+\text{corrsum})`; lonely iff `>0` |
| The wall | **domination is FALSE** (MISTAKE-202): top-term share `0.67→0.04%`; absolute/`ℓ¹` bounds over-count | **the cancellation wall** (HYP-5830): the unsigned small-relation mass over-counts; the `S266` Fourier expansion of `∏_w(1−1_{D_w})` **diverges** |
| The escape | **positivity** (klein-**S363**): `α=r|c|²≥0` ⟹ every moment `>0` (Hankel-PD); the signed `β`-coupling is the residual | **positivity** (klein-**S287**): `|ε_v|²≤(6/49)\,\text{disc}_v`, `disc_v` = the good-set **autocorrelation** discrepancy — positive-definite, the multilinear content absorbed intact, the divergence avoided |
| The residual | a **positive** Hankel/spectral bound on the `β`-block | a **positive** geometric bound on `disc_v` (`~(\#\text{edges})²/v²`, good-set spectral decay) |

Same author, same session-cluster, same manoeuvre: **when domination / absolute bounds / product-
expansion over-count, replace the signed object by its square (Hankel matrix, autocorrelation) and
work on the positive-definite side.** GMC's `MISTAKE-202` and LRC's `HYP-5830` are the *same lesson
learned twice*; GMC-`S363` and LRC-`S287` are the *same fix applied twice*. The reflection worth
keeping: this is not analogy, it is one theorem of method. Whenever a cancellation resists an
absolute bound, the escape is Parseval, not a sharper triangle inequality.

## 2. The shared object, exactly (demonstrated)

Both residuals are the **Parseval energy of an exponential sum over arc endpoints (jumps)**.

The LRC good-set `G'_{~v}` is a **union of arcs**; its Fourier coefficient is (verified to machine
precision, `04-computation/lrc_gmc_bridge_deathstar_S67.py`):

  **`ĉ_m = S_m/(2πi m)`,  `S_m = Σ_j sign_j · e^{−2πi m x_j}`**  (`x_j` = arc endpoints, `sign_j=±1` = up/down jump),

and by Wiener–Khinchin `disc_v = Σ_{m≠0}|ĉ_{mv}|² = Σ_{m≠0}|S_{mv}|²/(2π m v)²` (spatial autocorrelation
= Parseval sum, verified equal for `v=5,7,13`). The `1/m` is the **soft-Weyl / fold decay**; `S_m` is
**verbatim** boxeph's GMC reconstruction sum `Σ_j β_j e^{−imθ_j}` — jumps `β_j = sign_j` at endpoints
`θ_j = 2πx_j`. So:

- GMC reconstructs `A_fixed` as `Cauchy[jumps]` and controls `Σ_j β_j e^{−imθ_j}`;
- LRC's remaining bound is `disc_v = ` the `1/m²`-weighted Parseval energy of the **same** sum.

The GMC **thin-tie closure** (`S66`/boxeph-`S182`) proved `Cesàro⟨|S_m|²⟩ = Σ|β_j|²` (Wiener) — the
*unweighted* version of exactly the LRC `disc_v`. The LRC needs the `1/m²`-weighted sum bounded above;
GMC needs the Cesàro mean bounded below. **Same energy, opposite inequality.**

## 3. Four concrete transfers

1. **Autocorrelation-from-arcs = reconstruction-from-jumps (the sharpest lead).** boxeph's `A_fixed =
   Cauchy[jumps]` says a function is determined exactly by its arc/jump data. The LRC autocorrelation
   `A_{~v}(τ)=|G'∩(G'−τ)|` is likewise **exactly** determined by the arc endpoints: it is a piecewise-
   linear "tent" whose breakpoints are the pairwise endpoint differences `x_i−x_j`, and the three-gap
   theorem gives those endpoints in closed form. So `disc_v` has an **exact combinatorial value** from
   the endpoint set — replacing klein's crude `(\#\text{edges})²/v²` bound by the true energy
   `Σ_{m}|S_{mv}|²/(2πmv)²`. Computing it on the deep well / residue extremals is the immediate move.

2. **GMC soft-Weyl/fold decay → the `disc_v` upper bound.** `disc_v` is "governed by good-set spectral
   decay" (klein-S287). The decay is exactly the GMC fold structure: `|ĉ_m| = |S_m|/(2πm)`, and `|S_m|`
   is an endpoint exponential sum bounded by `2·(\#arcs)` with `1/m` from integration by parts — the
   very "soft-Weyl" bound the LRC thread already isolated (THM-1071: `|∫_G e(kt)dt|≤C/(π|k|)`,
   `C=\#components`). GMC's contribution: the fold **amplitude law** `∝ e^{−r}/√(D_{rr})` and the
   **loop non-vanishing** (`S66`) show *which* endpoints dominate and that they do not cancel — a
   route to the constant in `disc_v ≤ (\text{const})·(\#\text{arcs})²/v²`, the missing analytic bound.

3. **GMC loop-non-vanishing (`S66`) ⟷ the AP-uniqueness (klein-S290).** The AP `{1..13}` is the
   maximal-**resonance** extremal — in GMC terms, the *degenerate/stacked* thin stratum where the
   generic (Vandermonde) argument collapses. GMC-`S66` proved that *even the stacked stratum does not
   cancel* (the arc between coincident folds is a loop ⟹ opposite `g'` signs ⟹ real ⊥ imaginary ⟹
   `β_total≠0`). Its LRC twin is already on record: klein-S290's `L=|G(C)|(1−\text{conc}/7)` with the
   AP the **unique** `conc=7` point and covering forcing `conc≤5.2 ⟹ L≥0.26|G(C)|>0`. **Both say: the
   degenerate stratum is isolated and the margin survives it by a geometric (loop / concentration)
   invariant, not an estimate.** The GMC loop is a candidate *proof shape* for the LRC covering margin.

4. **Don't expand the product — the shared discipline.** GMC's `S266`-style failure was expanding a
   product into a divergent multilinear sum; the fix was the autocorrelation (Parseval). This is the
   GMC `MISTAKE-202` discipline in Fourier form: never trade a positive-definite energy for a signed
   term-by-term sum. Any future LRC attempt on `corrsum` should stay on the `A_{~v}` side; any GMC
   attempt on the `β`-coupling should stay on the Hankel side. The two threads can share one lemma:
   **a union-of-arcs energy `Σ_m w_m|S_{mv}|²` is a positive combinatorial functional of the endpoints**,
   estimable by three-gap geometry — usable as a lower bound (GMC loneliness/nonvanishing) or an upper
   bound (LRC `disc_v`).

## Honest status
A synthesis + one machine-verified structural identity (`disc_v` = Parseval of the GMC jump-sum), not
a proof on either side. The transfers are directions; #1 (exact autocorrelation from three-gap arcs)
is the most concrete and testable next step for LRC(14)'s last inequality. The deeper takeaway is the
method-level reflection: **positivity past the cancellation wall is one manoeuvre, and it has now paid
out on both flagship problems.**

## Credit
klein-S363 (GMC positivity) & S287 (LRC autocorrelation) — the same move, the seed of this note;
opus-S173/S262/S266 (LRC Riesz positive-definite route, the multilinear divergence); boxeph-S181/S182
(GMC reconstruction from jumps); klein-S284 (LRC finish-map), S290 (AP-uniqueness), S264 (signed
OffLine), THM-1071 (soft-Weyl `C/π|k|`); death-star-S64/S66 (GMC Hankel/hierarchy, the loop argument).

## Cross-links
GMC: `MISTAKE-202`, HYP-8480/8510, klein-S363. LRC: `00-navigation/LRC14-FINISH-MAP-2026-07-13.md`,
THM-731 (`disc_v`), THM-1071 (soft-Weyl), HYP-5830 (cancellation wall), HYP-6485 (autocorrelation
`x`-integral), opus HYP-5620 (Riesz). Object: `04-computation/lrc_gmc_bridge_deathstar_S67.py`. HYP-8515.
