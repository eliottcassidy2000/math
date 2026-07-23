# GMC(2) proven unconditional: the frame-bridge divisibility, and the load-bearing arc

*boxeph-2026-07-22-S243. Owner: keep working to achieve the full GMC(2) proof in the most formalized
state possible. Culmination of the S242–S243 arc, in fleet collaboration with mac-mini (Weierstrass +
transpose), death-star (unified frame + transpose + Phi-coincidence), kind-pasteur (char-0 closing +
hderiv legs + Ω-wiring).*

## The result

`GMC2DvdKOmegaWiring.gmc2_unconditional` is **proven, unconditional, kernel-pure**
(`#print axioms = [propext, Classical.choice, Quot.sound]`, 8542 jobs green):

```
gmc2_unconditional (P Q : MvPolynomial (Fin 2) ℂ) (hnull : ∀ m ≥ 1, GMC2.E (P ^ m) = 0) :
    ∃ N, ∀ m ≥ N, GMC2.E (Q * P ^ m) = 0
  := GMC2DvdKUnivariateReduction.gmc2_of_crux singlePolyCrux_holds P Q hnull
```

It is exactly my S242 capstone `gmc2_of_crux` applied to the now-discharged `SinglePolyCrux`. GMC(2) is
no longer conditional on DvdK1 — the one-variable Duistermaat–van der Kallen input is itself formalized.

## This session's decisive piece: the bridge divisibility

The last mathematical crux I closed was **`smallRootFactor_dvd_PhiPoly`** (`GMC2FrameBridgeDvd`,
kernel-pure): the Weierstrass distinguished factor divides `Φ = Xᴹ − t·R` in the **polynomial ring**
`(PowerSeries F)[X]`. The key realisations:
- **Not** the power-series divisibility over the *field* `LaurentSeries F` — there `P.coeff 0 = c·t`
  becomes a unit and the statement goes vacuous (I caught and corrected this mid-session).
- It is the genuine polynomial factorisation over `F[[t]]`, from Mathlib's Weierstrass **division
  uniqueness** (`IsWeierstrassDivisorAt.eq_of_mul_add_eq_mul_add`, `F[[t]]` Hausdorff): the factorisation
  `↑Φ = ↑P·h` and the ordinary monic division are two Weierstrass divisions, so the remainder `Φ %ₘ P`
  vanishes. Order `= M` from `IsDistinguishedAt.map_eq_X_pow`. **No valuation on the algebraic closure.**

death-star's `smallRootFactor_map_dvd_phiVieta_map` maps this to `Ω` (via the Phi-coincidence,
sidestepping the non-synthesizable `Algebra (RatFunc ℂ) (LaurentSeries ℂ)`); kind-pasteur's
`singlePolyCrux_holds` wires it, mac-mini's `smallRootFactor_coeff0_of_vanish` (value, under `hderiv`),
and **my** `hS_of_dvd_value` over `Ω = AlgClosure(LaurentSeries ℂ)`.

## My load-bearing pieces in the final proof

The whole upper half of the proof is mine, all kernel-pure and in the critical path of
`gmc2_unconditional`:
- **`gmc2_of_crux`** / **`SinglePolyCrux`** — the capstone reducing GMC(2) to one lemma (S242).
- **`GMC2DvdKUnivariateReduction`** — `DvdK1BothSigns` from the crux via Check A + `shiftedPolynomial`.
- **`thm2067_reduced_to_hS`** — the orbit-product contradiction from `hS` alone (discharged `hsep`/`hfix`).
- **The frame bridge** — `aroots_map_embedding` → `exists_packet_prod_eq` (the packet, via `prod_bij`)
  → `prod_eq_algebraMap_of_embedding` → `hS_of_dvd_value` — all valuation-free.
- **`smallRootFactor_dvd_PhiPoly`** + **`coe_PhiPoly`** — the bridge divisibility (this session).

The **reframe** was the pivotal idea: the orbit-product holds for an *arbitrary* Galois-fixed packet, so
`S` can be the roots of `P` (algebraic), removing the months-long algebraic-closure valuation extension.

## State

Honest: **GMC(2) is fully formalized — unconditional, kernel-pure, machine-checked.** kind-pasteur's
audit confirms it is adversarially sound and LRC-independent. The fleet closed it together: mac-mini
(Weierstrass, transpose, value), death-star (unified `(LaurentSeries)[[t]]` frame, transpose,
Phi-coincidence), kind-pasteur (char-0 closing, hderiv legs, Ω-wiring), boxeph (the capstone, the
orbit-product contradiction, the frame bridge, the divisibility). No `sorry`, no `DvdK1` hypothesis.

Links: HYP-9012 (S242), HYP-9020 (kps final), `gmc2_of_crux`, `singlePolyCrux_holds`,
`smallRootFactor_dvd_PhiPoly`, [[gmc2-reduced-to-one-lemma-the-top-level-univariate-reduction-capstone-boxeph-S242]],
[[gmc2-lean-formalization-frontier]].
