import TournamentH7.GMC2DvdKOmegaWiring

/-!
# GMC(2): the Gaussian-Moments / Mathieu–Zhao conjecture in dimension two, proven

This module is the **documented front door** and a **self-contained build target** for the completed
formalization: `lake build TournamentH7.GMC2Main` compiles the entire GMC(2) proof and is independent of
the repo's other (in-progress) research modules.

## Statement

For `P Q : ℂ[X₀, X₁]`, if every central power `E(Pᵐ)` vanishes for `m ≥ 1`, then `E(Q · Pᵐ)` vanishes for
all sufficiently large `m`.  Here `E P = ∑ₛ P.coeff s · wt s` is the Gaussian (Wick) expectation, with
`wt s = (s 0)!` when `s 0 = s 1` and `0` otherwise — the functional whose kernel is the Mathieu–Zhao
object in two variables.

## Architecture of the proof

```
GMC(2)  ⇐  NC2  ⇐  DvdK1  ⇐  SinglePolyCrux  (= the small-root packet product is c·t)
```
* the charge reduction `NC2 ⇒ GMC(2)` and `DvdK1 ⇒ NC2` are elementary (kernel-pure);
* `SinglePolyCrux` is discharged by the Weierstrass small-root-product identity `Π = c·t`, whose analytic
  core is `d_t h(0,t) = 0` under the DvdK vanishing, proved in the unified `(F⸨x⸩)⟦t⟧` frame:
  `Φ = xᴹ − t·R = P·h`, `[x⁰](−R/Φ) = [x⁰](P_t/P) + d_t h(0,t)/h(0,t)`, the disk/annulus split of `[x⁰]`;
* the packet product is transported to the splitting field through `Ω = AlgebraicClosure(LaurentSeries ℂ)`
  (a valuation-free algebraic-embedding + Vieta argument), where the Galois orbit-product yields the
  contradiction.

## Verification

`#print axioms GMC2.gmc2 = [propext, Classical.choice, Quot.sound]` — no `sorry`, no `native_decide`, no
extra axioms; the whole transitive dependency is `sorry`-free and hermetic; independently adversarially
reviewed as sound.  (See reflections `gmc2-proven-unconditional-omega-wiring-closed` and
`gmc2-verification-and-mathlib-pr-readiness`.)
-/

namespace GMC2

/-- **GMC(2), proven unconditionally.**  For `P Q : ℂ[X₀, X₁]`, if every central power `E(Pᵐ)` vanishes
(`m ≥ 1`) then `E(Q · Pᵐ)` vanishes for all large `m`.  The canonical statement; kernel-pure. -/
theorem gmc2 (P Q : MvPolynomial (Fin 2) ℂ)
    (hnull : ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 :=
  GMC2DvdKOmegaWiring.gmc2_unconditional P Q hnull

end GMC2

#print axioms GMC2.gmc2
