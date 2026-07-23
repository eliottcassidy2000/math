import TournamentH7.GMC2DvdKOmegaWiring

/-!
# Independent axiom audit of the completed GMC(2) chain (mac-mini S167)

Run with:  `lake env lean AxiomCheckGMC2MacMiniS167.lean`
(or build the module: `lake build TournamentH7.GMC2DvdKOmegaWiring`).

Every line below must report exactly `[propext, Classical.choice, Quot.sound]` —
Mathlib's standard axioms — with **no `sorryAx`** and no project-local axiom.
Anything else falsifies the "GMC(2) proven unconditionally" claim.

The chain being certified, bottom to top:

* `heightWitnessSupplier_holds` — the normalized height package (proved, not assumed);
* `nc2_of_dvdK1` / `gmc2_of_dvdK1` — endpoints needing only the DvdK input;
* `singlePolyCrux_holds` — kind-pasteur's Ω-wiring (S128c153): the last analytic input,
  built over `Ω = AlgebraicClosure (LaurentSeries ℂ)` from death-star's divisibility
  connector, mac-mini's `smallRootFactor_coeff0_of_vanish` value, and boxeph's
  `hS_of_dvd_value` packet bridge;
* `dvdK1_of_crux` — the one-variable Duistermaat–van den Kallen statement is now
  **derived**, so the external citation is removed;
* `gmc2_unconditional` — the capstone.

This file is an audit only: it proves nothing and adds no dependency.
-/

-- the height package (formerly a visible premise)
#print axioms GMC2NC2.heightWitnessSupplier_holds

-- height-discharged endpoints from the DvdK input
#print axioms GMC2NC2.nc2_of_dvdK1
#print axioms GMC2NC2.gmc2_of_dvdK1

-- the Ω-wiring: the last remaining input, now proved
#print axioms GMC2DvdKOmegaWiring.singlePolyCrux_holds

-- DvdK1 is *derived* from the crux (citation removed, not assumed)
#print axioms GMC2DvdKUnivariateReduction.dvdK1_of_crux

-- load-bearing upstream inputs
#print axioms GMC2DvdKTransposeAssembly.smallRootFactor_coeff0_of_vanish
#print axioms GMC2DvdKPhiCoincide.smallRootFactor_map_dvd_phiVieta_map
#print axioms GMC2FrameBridgeAssembly.hS_of_dvd_value

-- the capstone
#print axioms GMC2DvdKOmegaWiring.gmc2_unconditional
