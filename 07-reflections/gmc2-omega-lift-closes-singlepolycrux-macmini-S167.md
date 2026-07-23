# GMC(2) is closed: the Ω-lift, independently derived twice — and what "most formalized" still means

*mac-mini-2026-07-22-S167. Owner: "finish the GMC(2) formalization and get it into the best possible
state for a mathlib PR; prioritize open mathematical reasoning." Outcome: the last gap
(`SinglePolyCrux`) was closed by **kind-pasteur's `GMC2DvdKOmegaWiring`** (commit `0784107bd`,
HYP-9020), which landed while I was building the same construction independently. This records the
mathematics, the independent corroboration, my audit of the claim, and the honest remaining work.*

---

## The state going in: one `Prop` left

By S166 the Lean GMC(2) spine had been reduced, kernel-pure, to a single proposition
`GMC2DvdKUnivariateReduction.SinglePolyCrux`. Its leverage is total, and I verified the chain has no
hidden premise:

- `GMC2HeightWitness.heightWitnessSupplier_holds : HeightWitnessSupplier` is **proved**, not assumed —
  so `nc2_of_dvdK1` / `gmc2_of_dvdK1` need only the DvdK input. (The `…_of_heightWitnessSupplier`
  forms in `GMC2NC2`/`GMC2NC2Capstone` are the older conditional surface; the discharged versions live
  in `GMC2HeightWitness`.)
- `dvdK1_of_crux : SinglePolyCrux → DvdK1` — proving the crux **removes the external
  Duistermaat–van den Kallen citation** (it becomes a derived theorem, not a hypothesis).
- `gmc2_of_crux : SinglePolyCrux → GMC(2)`.

So: `heightWitnessSupplier_holds` (proved) + `dvdK1_of_crux` ⟸ **`SinglePolyCrux`** ⟹ unconditional GMC(2).

## Why the closure needs an algebraic closure

The naive frame `Ω = LaurentSeries ℂ = ℂ((t))` **fails**: the small roots of `Φ = Xᴹ − t·R` are
ramified, `t^{1/M} ∉ ℂ((t))`, so `Φ` neither splits there nor admits an embedding `ψ` of its splitting
field. This is why the frame apparatus (death-star's transpose, boxeph's bridge) stalled at the
endpoint, and it is the point kind-pasteur re-identified after an initial over-claim and retract. The
fix is to pass to `Ω = AlgebraicClosure (LaurentSeries ℂ)`.

## The Ω-lift (the construction, as landed)

`GMC2DvdKOmegaWiring.singlePolyCrux_holds` supplies boxeph's `hS_of_dvd_value` with Ω-frame data:

1. **The algebra structure.** `Algebra (RatFunc ℂ) (LaurentSeries ℂ)` is not synthesizable, so it is
   provided locally as `rfToL.toAlgebra` (death-star's bundled coercion); `AlgebraicClosure`'s own
   instance then lifts it to `Ω`, and `IsScalarTower.of_algebraMap_eq (fun _ => rfl)` makes the
   diamond commute. The compatibility `(algebraMap L Ω).comp rfToL = algebraMap (RatFunc ℂ) Ω` holds
   **by `rfl`**. That is the trick that dissolves the instance wall everyone had hit.
2. **The divisor** `Pω = (smallRootFactor R M).map (ofPowerSeries) .map (algebraMap L Ω)`:
   **monic** (`smallRootFactor_monic.map.map`); **splits** (`IsAlgClosed.splits`, Ω algebraically
   closed); and **`roots.Nodup`** — the one genuine mathematical step — from **separability**: `Φ` is
   irreducible (`irreducible_Phi`) over a characteristic-zero (perfect) field, hence separable; the
   map to Ω is separable (`Separable.map`); and a divisor of a separable polynomial is separable
   (`Separable.of_dvd`), hence has distinct roots (`nodup_roots`). Separability is exactly what makes
   the packet a *reduced* set of points.
3. **Divisibility** — death-star's `smallRootFactor_map_dvd_phiVieta_map` pushed through
   `algebraMap L Ω` (`map_dvd_map'` + `Polynomial.map_map` + the `rfl` compatibility).
4. **Value** — mac-mini's `smallRootFactor_coeff0_of_vanish` (`P.coeff 0 = −t·r₀`) transported by
   `ofPowerSeries_comp_C` + `rfToL_comp_algebraMap`, giving `Pω.coeff 0 = algebraMap v`,
   `v = −C(r₀)·X`.
5. **Embedding** `ψ := IsAlgClosed.lift` (the splitting field is finite-dimensional hence algebraic;
   `RatFunc ℂ` is a field so torsion-freeness is automatic), plus a root witness from
   `Splits.roots_ne_zero`.

Then `hS_of_dvd_value` yields `∏_{β∈S} β = algebraMap((−1)^{deg Pω}·v)`, and with
`Pω.natDegree = M` this folds to `RatFunc.C cc · RatFunc.X` with **`cc = (−1)^{M+1}·r₀ ≠ 0`** —
exactly `SinglePolyCrux`. Hence `gmc2_unconditional := gmc2_of_crux singlePolyCrux_holds`.

## Independent corroboration (why this is more than one agent's say-so)

I built the same Ω-lift independently and in parallel (`GMC2SinglePolyCruxProof`, not pushed — ceded
as a duplicate once kps's landed first). The two derivations were written without sight of each
other and agree on every structural choice and on the arithmetic:

| step | kind-pasteur | mac-mini (independent) |
|---|---|---|
| frame | `AlgebraicClosure (LaurentSeries ℂ)` | `AlgClosure (LaurentSeries ℂ)` |
| algebra | `rfToL.toAlgebra`, lifted; `halg` by `rfl` | composite `.toAlgebra`; `algebraMap_toAlgebra` |
| nodup | irreducible ⇒ separable ⇒ `Separable.of_dvd` | same |
| embedding | `IsAlgClosed.lift` | same |
| constant | `cc = (−1)^{M+1}·r₀` | `cc = (−1)^{deg Pω}·(−r₀)`, `deg Pω = M` |

The two `cc` expressions are equal. For a headline claim on a conjecture this old, two independent
constructions landing on the same normalization is worth more than either alone. kps's version is
the better one to keep: the `rfl` scalar-tower diamond and the explicit `Pω.natDegree = M` are
cleaner than my composite-instance bookkeeping.

## Credit

- **boxeph** — the frame bridge: `smallRootFactor_dvd_PhiPoly`, `hS_of_dvd_value`, `gmc2_of_crux`.
- **death-star** — the transpose hom and the Φ-coincidence connector (`rfToL`,
  `smallRootFactor_map_dvd_phiVieta_map`) putting the divisibility over a shared Laurent frame.
- **mac-mini** — the transpose-glue keystone `hderiv_final` and the value
  `smallRootFactor_coeff0_of_vanish` (S165/S166), a load-bearing input to the wiring; plus the
  independent Ω-lift and this audit.
- **kind-pasteur** — the hderiv F-leg, the diagnosis that the endpoint needs `AlgebraicClosure`, and
  the landed `GMC2DvdKOmegaWiring`.

## The analytic-heart line-audit (the part kernel-purity cannot establish)

The referee's flag on kps's polish plan was exact: `hderiv_final` / `smallRootFactor_coeff0_of_vanish`
are *kernel-checked but not line-audited*. Kernel-purity proves **"no `sorry`"**; it does **not** prove
the *statement* is the one the mathematics needs. Since these are my lemmas (S165/S166) and they are
the deepest link, I audited them on three targets. All three pass:

**(i) Does `hderiv_final`'s conclusion match what the consumer requires?** Yes, exactly.
`GMC2DvdKMultiplicativeClosing.smallRootFactor_coeff0_eq_of_derivative_vanishes` requires
`hderiv : derivativeFun (unitCoeff0 R M) = 0`, and `unitCoeff0 R M :=
constantCoeff (R := PowerSeries F) (weierstrassUnit (Phi R M) …)` — the `x`-constant coefficient of the
Weierstrass unit, landing in `F⟦t⟧`, i.e. literally `h(0,t)`. `hderiv_final` produces precisely that
proposition. The mathematics is then the honest one: `d_t h(0,t) = 0` plus `h(0,0) = 1` plus
characteristic zero gives `h(0,t) = 1` — no exp, log, Puiseux, or Fredholm determinant.

**(ii) Is `CharZero` load-bearing and correctly scoped?** Yes. `[CharZero F]` sits in the section
variables of both `GMC2DvdKMultiplicativeClosing` and `GMC2DvdKCharZeroClosing`, and it is used exactly
where it must be — `eq_C_of_derivativeFun_eq_zero` (a power series with vanishing derivative is
constant), which is **false in characteristic `p`** (`tᵖ` has zero derivative). Evidence the scoping is
deliberate rather than incidental: `GMC2DvdKCharZeroClosing` carries an explicit
`omit [CharZero F] in` on `factorCoeff0_eq_of_unit_eq_one`, the one step that genuinely does not need
it. GMC(2) is stated over `ℂ`, so the hypothesis is free at the point of use.

**(iii) Is there a silent completion mismatch?** No — and this was the real hazard, since my own
earlier blocker was "`[x⁰]` across two completions." It is discharged by three *proved* lemmas rather
than by assumption:
- everything is transported through the **single** ring hom `phi : F⟦t⟧⟦x⟧ →+* (F⸨x⸩)⟦t⟧`, and `hfact`
  obtains the frame factorization from the Weierstrass one by `map_mul` on `H.eq_mul` — the
  factorization is *moved*, never re-derived;
- `Rl_pow_coeff` proves the **frame** moment `(Rl R ^ m).coeff (M·m)` equals the **polynomial** moment
  `(Rᵐ).coeff (M·m)` — exactly the place a vanishing hypothesis could be silently swapped for a
  different one;
- `xCoeff0_phi_unit` proves the seam `xCoeff0 (phi h) = unitCoeff0 R M`, tying the frame's
  `x⁰`-extraction back to the original ring's `x`-constant term.

> **Audit verdict.** The analytic core states what the mathematics requires, its characteristic-zero
> hypothesis is load-bearing and correctly placed, and the two-completion seam is closed by explicit
> lemmas. I found no defect. This is a statement-level audit; it does not re-verify the proofs of the
> upstream lemmas it cites.

## Honest scope, and what "most formalized" still requires

- **What is claimed:** `GMC2DvdKOmegaWiring.gmc2_unconditional` has no external hypothesis and no
  `sorry` in its file. The certificate is `#print axioms` reporting exactly
  `[propext, Classical.choice, Quot.sound]`; the audit file
  `AxiomCheckGMC2MacMiniS167.lean` checks the whole chain, not just the capstone, precisely because
  a capstone can look clean while an upstream lemma is not. **Do not cite this result from a
  reflection — cite the audit output.**
- **What is not claimed:** I did not re-derive the upstream inputs (`irreducible_Phi`, the Weierstrass
  factorization, the height package). They are cited. A genuinely adversarial review of GMC(2) should
  target those, and the `hderiv` analytic core, not the Ω-wiring, which is the shallowest link.
- **Still open for a mathlib PR** (packaging, not mathematics): ~90 `GMC2*` modules with
  `import Mathlib` at the top of most; these need precise imports, a coherent module graph, naming and
  docstring conventions, and the older conditional surfaces (`GMC2NC2Capstone`'s
  `…_of_heightWitnessSupplier`) either retired or clearly marked as superseded so the "remaining
  inputs" story is not misread. That is the real remaining work, and it is substantial.

---

*Cross-links: `GMC2DvdKOmegaWiring` (kps, the landed proof), `GMC2DvdKUnivariateReduction`
(`SinglePolyCrux`, `dvdK1_of_crux`, `gmc2_of_crux`), `GMC2HeightWitness`
(`heightWitnessSupplier_holds`, `nc2_of_dvdK1`), `GMC2FrameBridgeAssembly` (boxeph),
`GMC2DvdKPhiCoincide` (death-star), `GMC2DvdKTransposeAssembly` (mac-mini), `AxiomCheckGMC2MacMiniS167.lean`.
HYP-9020. Memory: [[gmc2-formalization-endgame]].*
