# GMC(2) proven unconditionally, kernel-pure — the Omega-wiring closed

*kind-pasteur-2026-07-22-S128c153. Owner: finish the Omega-wiring + remaining parts. It closed —
`GMC2DvdKOmegaWiring.gmc2_unconditional` is on `origin/main`, `#print axioms = [propext,
Classical.choice, Quot.sound]`, no `sorry`. This ends the multi-agent GMC(2) formalization.*

## The theorem

```lean
theorem gmc2_unconditional (P Q : MvPolynomial (Fin 2) ℂ)
    (hnull : ∀ m : ℕ, 1 ≤ m → GMC2.E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, GMC2.E (Q * P ^ m) = 0
```

— GMC(2): if every central power `E(Pᵐ)` vanishes, then `E(Q·Pᵐ)` vanishes eventually. It is
`gmc2_of_crux` (boxeph) applied to the now-proved `singlePolyCrux_holds`.

## The Omega-wiring (`singlePolyCrux_holds`) — what closed the last hypothesis

`SinglePolyCrux`: for `Φ = Phi R M`, under the DvdK vanishing, the small-root packet product is `c·t` in
the splitting field. The obstacle was a **single field `Ω` holding both** the splitting field (for the
packet) and death-star's `LaurentSeries ℂ` divisibility with `Pω` split. `Ω = LaurentSeries ℂ` fails (the
small roots need ramification `x^{1/M}`). The resolution:

**`Ω := AlgebraicClosure (LaurentSeries ℂ)`.** The non-synthesizable `Algebra (RatFunc ℂ)(LaurentSeries ℂ)`
is provided locally via death-star's `rfToL.toAlgebra`; then `AlgebraicClosure.instAlgebra` lifts it to
`Ω`, `IsScalarTower.of_algebraMap_eq rfl` aligns the tower, and `ψ := IsAlgClosed.lift` embeds the
splitting field. The instance diamond that blocked everyone **resolves by `rfl`** once the base algebra is
handed over as a ring hom — the key finding.

Then the three landed pieces compose:
- **`Pω := (smallRootFactor R M).map(ofPowerSeries).map(algebraMap _ Ω)`**, dividing `Φ.map` via
  death-star's `smallRootFactor_map_dvd_phiVieta_map` + `map_dvd_map'` (goal-side `map_map` to avoid the
  ambiguous nested-map rewrite);
- **`Pω.coeff 0 = algebraMap v`, `v = −C(r₀)·X`**, from mac-mini's `smallRootFactor_coeff0_of_vanish`
  (`= −X·r₀`), transported by `ofPowerSeries_comp_C` and `rfToL_comp_algebraMap` (both `ofPS(algMap r₀)`
  and `rfToL(C r₀)` collapse to `algebraMap ℂ (LaurentSeries ℂ) r₀`);
- **boxeph's `hS_of_dvd_value`** then gives `∏_{β∈S} β = algebraMap((-1)ᴹ·v) = algebraMap(C cc·X)`,
  `cc = (-1)^{M+1} r₀ ≠ 0` (since `r₀ = R.coeff 0 ≠ 0`).

The `_x0` witness is a root of `Φ` in its splitting field (`IsSplittingField.splits` +
`Splits.roots_ne_zero` + `mem_rootSet`).

**Honest note on process.** I twice mis-scoped this: first as "trivial composition" (retracted — the `Ω`
crux is real), then feared it intractable (wrong — the instance path resolved on the first build). The
lesson: the viability was only knowable by *attempting the build*; the instance diamond that looked
fatal in prose was a `rfl`.

## Credit — this was the whole fleet

mac-mini (Weierstrass preparation, `hderiv_final`, the `−X·r₀` value, `hconst`), death-star (the
`(LaurentSeries F)⟦t⟧` frame, `xCoeff0`, the transpose `φ`, `rfToL`/`PhiCoincide`, the x-degree toolkit,
`hconst`), boxeph (the frame bridge, `smallRootFactor_dvd_PhiPoly`, `hS_of_dvd_value`, `gmc2_of_crux`),
codex (Check A, nodal Lagrange). My through-line: the F=D_m leg, the h-side (disk/annulus insight),
`logDeriv_map`, the `hderiv` assembly backbone (`hderiv_of_transpose_glue`, which `hderiv_final` is built
on), and this Omega-wiring.

## Cross-links
`GMC2DvdKOmegaWiring` (this) · `GMC2DvdKUnivariateReduction.gmc2_of_crux` · `GMC2FrameBridgeAssembly` ·
`GMC2DvdKPhiCoincide` · `GMC2DvdKTransposeAssembly` (hderiv_final, value) · `GMC2DvdKHderivAssembly`
(my backbone) · `GMC2DvdKFrameHSide` / `GMC2DvdKFrameExtraction` (my legs) · THM-1550 · HYP-9020.
