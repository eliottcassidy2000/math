---
id: HYP-2840
title: The resonance bound p0<=p0_decorr -- multi-far residual via VITALI COVERING (geometric, bypasses the divergent Fourier/lattice envelope)
status: STRATEGY + foundational ATOM FORMALIZED (LRCMarginalUniform, sorry-free); single-far DONE (THM-546), multi-far = the residual where Fourier diverges; Vitali covering = the geometric route, tight step = quasi-independence (open)
source: kind-pasteur-2026-06-22-S31
related:
  - HYP-2839   # hp0cap cores + the resonance residual isolated
  - THM-546    # single-far comb bound (Fourier/BV), PROVED
  - HYP-2643   # the absolute multi-D lattice envelope DIVERGES
  - HYP-2757   # the commensurable curve atlas (resonant far pairs)
  - THM-534    # decorrelated closed form (moment dual)
---

# HYP-2840 -- the resonance bound via Vitali covering

## The target (the one analytic residual of hp0cap)
`p0(E) <= p0_decorr(E) + o(1)`, where `p0_decorr = sum_t P_t^{(r)} p_t(B)` (THM-534) is the
decorrelated coverage (far phases independent-uniform). With `p0_decorr <= Q(k-1) < cap` (finite +
rational, machine-checkable), this closes hp0cap. So the residual is purely: **the far phases
decorrelate** (the true joint coverage <= the product/IE value).

## The split (where the spectral approach lives and dies)
- **SINGLE-FAR: DONE (THM-546).** Isolate one far element `w`: `Delta_w = p0(E'u{w}) - p0(E') -
  (1/7)p1(E')` obeys `|Delta_w| <= kappa V(E')/(pi^2 w)`, elementary Fourier/BV. Converges.
- **MULTI-FAR: the residual.** Iterating fails when far elements CLUSTER (resonance `w_i/w_j` near
  rational): the **absolute multi-D Fourier/lattice envelope DIVERGES** (HYP-2643) -- the spectral
  sum over the resonance lattice blows up. This is exactly where a NON-spectral tool is needed.

## VITALI COVERING (the owner's hint): bound the multi-far coverage GEOMETRICALLY
Write the coverage set as a union of "resonance boxes" indexed by sector->speed assignments:
`coverSet = U_{sigma injective} Box_sigma`,  `Box_sigma = cap_j {x: frac(sigma(j) x) in I_j}`
(each `Box_sigma` a finite union of intervals). The decorrelated value is the inclusion-exclusion of
the `Box_sigma` measures. **Vitali's covering theorem** (mathlib `Vitali.exists_disjoint_subfamily_
covering_enlargement`, 1-D `Mathlib/.../Covering/OneDim`) extracts a DISJOINT subfamily covering the
union, giving `meas(coverSet) <= C * sum_{disjoint} |Box|` WITHOUT summing the divergent Fourier
lattice -- a purely geometric (real-space) bound, immune to the HYP-2643 divergence.

### What Vitali delivers, and the honest tightness gap
- **Promise:** bypasses the spectral divergence -- the multi-far bound becomes a real-space
  disjointification, always convergent (intervals on `[0,1)`).
- **Limitation (verified by counting):** the crude constant is LOSSY for the tight margin. The
  first-order union bound `sum_sigma meas(Box_sigma) = (#assignments)*(box measure)` ALREADY
  overshoots: k=8 gives ~0.171 < cap 0.381 (ok), but k=9 gives ~0.514 > cap 0.494 (FAILS). So
  subadditivity/crude-Vitali is too lossy; the TIGHT bound is the inclusion-exclusion (decorrelated
  value, which SUBTRACTS the box overlaps). Vitali covering must therefore be refined to capture the
  overlap structure = **quasi-independence** (bounded box overlap).
- **The refinement = the actual decorrelation.** In 1-D the resonance boxes 2-color into disjoint
  classes by max-overlap m (interval graphs are perfect, chi = m); `meas(coverSet) <= m * (disjoint
  class sum)`. For GENERIC (Weyl-equidistributed) far elements the boxes are near-disjoint (m small,
  bound ~ tight); for RESONANT (commensurable) far elements the boxes align (m large) BUT the count
  drops proportionally (the curve reduction HYP-2757: the far pair sweeps a 1-D geodesic, |R|<=2
  coverage). So: **Vitali covering + the max-overlap/quasi-independence refinement** is the geometric
  route, with the resonant case handed to the finite curve atlas (HYP-2757).

## Vitali CONVERGENCE (a separate, cleaner sub-use)
For the single-far LIMIT `p0(E u {w}) -> p0(E) + p1(E)/7` as `w->infty`: the coverage indicators are
bounded by 1 (uniformly integrable on the finite measure `[0,1)`), and the cross term decorrelates by
equidistribution of `frac(w x)` (Weyl). NB the indicators converge only WEAKLY (oscillate), so this
is equidistribution, NOT mathlib's `tendstoInMeasure_iff_tendsto_Lp_finite` (strong) -- Vitali
CONVERGENCE does not apply to the raw indicator. Vitali COVERING (above) is the relevant Vitali.

## Formalizable pieces (sorry-free targets)
1. **Marginal uniformity** `meas{x in [0,1): frac(w x) in I} = |I|` (the base case; w-fold cover or
   `AddCircle.measurePreserving_zsmul`). The atom every box measure is built from. -- NEXT.
2. **Union bound** `meas(coverSet) <= sum_sigma meas(Box_sigma)` (subadditivity) -- the loose Vitali
   skeleton; sorry-free once Box_sigma is defined.
3. The max-overlap/quasi-independence refinement (the tight step) -- the genuine residual.

## CONVERGENCE with mac-mini's L_y route (S28) -- the decorrelation IS the far-element half
mac-mini (concurrent) found the CLEANER hp0cap reframing via THM-534: `p0(E) <= L_y(E)` is PROVED
(moment-LP dual, EXACT Bonferroni -- NO decorrelation needed for this step), so hp0cap reduces to the
SCALAR extremality **"consec maximizes L_y(E)"**. This is strictly cleaner than the measure-decorrelation
route for the `p0 <= (something)` step. BUT the remaining scalar piece "consec maximizes L_y" SPLITS:
- **bounded-spread**: finite exact check (consec argmax over span<=2k), VERIFIED.
- **far elements drop L_y**: mac-mini measured "ANY far element drops L_y by >=0.044 (decorrelation
  safe)". THIS half IS the decorrelation -- and the Vitali/marginal-uniformity work here is its
  rigorous backbone: far phases equidistribute (marginal uniformity, `LRCMarginalUniform`, PROVED) =>
  the far-element moments S_r decorrelate => L_y drops. So the two routes COMBINE: L_y reduces hp0cap
  to consec-extremality; its far-element half = this decorrelation (Vitali covering for the multi-far
  geometric bound where the Fourier envelope diverges). Defer the `p0<=L_y` step to THM-534 (cleaner);
  keep the Vitali/decorrelation for the far-element L_y drop.

## Status
PIVOT: the L_y route (THM-534, mac-mini S28) supersedes the measure-decorrelation for the `p0<=cap`
backbone (cleaner, exact). The resonance/decorrelation content survives as the FAR-ELEMENT half of
"consec maximizes L_y" (far elements decorrelate-safe). Formalized: the marginal-uniformity atom
(`LRCMarginalUniform`, sorry-free) -- the backbone of the far-element L_y drop. Open: the scalar
consec-extremality (mac-mini's L_y route) + the quasi-independence refinement (this Vitali route).
-> HYP-2839, THM-534 (L_y route), THM-546, HYP-2643, HYP-2757.
