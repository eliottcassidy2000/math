# GMC(2) endgame: the cohesive assembly map — every link from concrete glue to GMC(2)

*kind-pasteur-2026-07-22-S128c153. Owner: work the (c) finish + final assembly wiring; pull at many times;
integrate all agent work into one cohesive picture. This is that picture — the verified file→theorem→
author→hypothesis chain, with the exact residual glue named.*

## The full chain (all kernel-pure `[propext, Classical.choice, Quot.sound]` unless marked ⟳=in progress)

```
GMC(2)
 ⟸ gmc2_of_crux (SinglePolyCrux)                     boxeph  GMC2DvdKUnivariateReduction   DONE
 ⟸ SinglePolyCrux  (∏_{S} β = c·t in splitting field)
     ⟸ boxeph frame bridge: smallRootFactor_dvd_PhiPoly + false_of_frame_data / hS_of_*   boxeph  DONE
     ⟸ Π = c·t  (Weierstrass frame:  smallRootFactor.coeff 0 = −t·r₀)
         ⟸ smallRootFactor_coeff0_eq_of_derivative_vanishes'   death-star/kps  GMC2DvdKUnitOrigin   DONE
             ⟸ hderiv:  derivativeFun (unitCoeff0 R M) = 0
                 ⟸ **hderiv_of_transpose_glue**            kps  GMC2DvdKHderivAssembly   DONE (this session)
                     internal: ha  = kps xCoeff0_logDeriv_map_ofPowerSeries   (GMC2DvdKFrameHSide)   DONE
                     internal: hF1 = kps xCoeff0_xM_div_PhiFrame_eq_one_of_vanish (GMC2DvdKFrameExtraction) DONE
                     internal: hderiv_of_frame (master identity + logDeriv split) death-star GMC2DvdKHderiv DONE
                     ── residual hypotheses (the ONLY open glue) ──────────────────────────
                     hfact : PhiFrame Rl M = Pfr · phi Wu     ⟸ phi_Phi = PhiFrame   death-star  ⟳ (generators phi_X/phi_C_X/phi_C_C done; polynomial-image wiring left)
                     hPu   : IsUnit Pfr,  Pfr := phi(smallRootFactor)                death-star  ⟳
                     hc    : xCoeff0(logDeriv Pfr) = 0  ⟸ xCoeff0_logDeriv_eq_zero_of_monic  death-star DONE, needs APPLYING to Pfr (coeff0=single M 1, xdeg<M)
                     hg    : IsUnit (xCoeff0 (phi Wu))                                (g(0)=1)   open, small
                     hbridge: xCoeff0 (phi Wu) = unitCoeff0 R M                        death-star  ⟳
                     hvanish: ∀ m≥1, (Rl^m).coeff(M·m) = 0  ⟸ R→Rl coeff transport + boxeph generatingFunction_eq_one (D_m=0)   open
```

## What this session contributed to the assembly

- **`hderiv_of_transpose_glue`** (GMC2DvdKHderivAssembly, kernel-pure): the integration backbone. It
  discharges the two *analytic* inputs of `hderiv_of_frame` internally — my h-side `ha` and my `F=1`
  `hF1` — using that `phi Wu = map(ofPowerSeries)(tau Wu)` is definitionally the disk form. So `hderiv`
  (hence GMC(2)) is now reduced to **purely algebraic transpose glue**: `hfact, hPu, hc-application, hg,
  hbridge` (death-star's `phi`↔Weierstrass connectors) and `hvanish` (the R→Rl moment transport).
- The h-side `ha` itself (`logDeriv_map` + disk `[x⁰]`=ring-hom) and `hF1` (geometric inverse ⇒ moment
  series ⇒ `F=1`) were the prior two sessions' deliverables; they now plug in with no glue.

## The residual is NOT analytic — it is transpose bookkeeping

Every remaining item is a **ring-hom-image identity** (`phi` applied to a specific Weierstrass object) or a
**coefficient-preservation** fact (`(Rl^m).coeff = (R^m).coeff` under the polynomial→LaurentSeries embedding).
No residues, no Puiseux, no new analysis. death-star owns the `phi`-connectors (`phi_Phi`, `Pfr`, `hbridge`);
`hvanish` is the R→Rl transport composing with boxeph's `generatingFunction_eq_one`. Once `phi_Phi = PhiFrame`
lands, `hfact` is `phi(Φ)=phi(P·h)=phi(P)·phi(h)` and the whole chain closes to **unconditional GMC(2)**.

## Cross-links
`GMC2DvdKHderivAssembly` (this session) · `GMC2DvdKFrameHSide` / `GMC2DvdKFrameExtraction` (kps ha/hF1) ·
`GMC2DvdKHderiv` (death-star assembly) · `GMC2DvdKTranspose` (death-star phi) · `GMC2DvdKFrameDegree` (c) ·
`GMC2DvdKUnitOrigin` (Π=c·t from hderiv) · `GMC2FrameBridge*` (boxeph) · `GMC2DvdKUnivariateReduction`
(gmc2_of_crux) · HYP-9016.
