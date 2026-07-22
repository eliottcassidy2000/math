# The concrete Galois instantiation: THM-2067 for irreducible `Φ`, reduced to Vieta + THM-1550

*boxeph-2026-07-22-S236. Owner: work the two remaining pieces (Galois wrapper + Check A, and the deep
Hensel gap); push/pull often; mine past threads. Coordinated split with death-star (they drive
THM-1550; I take Gal-instantiation + Vieta + Check A). Builds on my S232–S235 (orbit-product core,
t-adic closing, (A) irreducibility, wrapper). New Lean `GMC2GalRootAction.lean`,
`GMC2Thm2067Concrete.lean`, kernel-pure.*

## Delivered: the concrete Galois instantiation (kernel-pure)

The abstract THM-2067 wrapper (S235) took a transitive Galois action + equivariant embedding as
hypotheses. This session I built the concrete Galois structure and instantiated the wrapper — so
THM-2067 is now proved (as a reduction) for an *actual* irreducible polynomial.

**`GMC2GalRootAction`** — the direct action `σ • x = σ x` of `L ≃ₐ[K] L` on `Φ.rootSet L`:
- `mem_rootSet_smul` — automorphisms fixing `K` preserve roots (`aeval_algHom_apply`);
- `rootAction` — the `MulAction` instance;
- `coe_smul` — equivariance `↑(σ • x) = σ ↑x`, **tautological** (`rfl`), sidestepping Mathlib's
  `galAction`/`rootsEquivRoots` machinery entirely;
- `isPretransitive_rootAction` — **transitivity** for irreducible `Φ` over a normal `L`: any two roots
  are conjugate (`IsConjRoot`, via `minpoly K y` being an associate of the irreducible `Φ`), so
  `IsConjRoot.exists_algEquiv` supplies the automorphism.

**`GMC2Thm2067Concrete.thm2067_contradiction_concrete`** — instantiates the wrapper at `Φ.Gal` acting
on `Φ.rootSet Φ.SplittingField`. Given an irreducible `Φ` over `F(t)`, the small-root product `= c·t`
and Galois-fixed (THM-1550), and the full root product `= d` a constant (Vieta), it derives `False`.
All the Galois instances (`Fintype Φ.Gal`, `MulDistribMulAction`, `IsSplittingField.finiteDimensional`,
transitivity) resolve once the type is written as `Φ.Gal`. Kernel-pure.

So the **entire THM-2067 argument is now a concrete kernel-pure reduction** — the six-module chain
`GMC2OrbitProduct → GMC2RatFuncClosing → GMC2PhiIrreducible → GMC2Thm2067Wrapper → GMC2GalRootAction →
GMC2Thm2067Concrete` (16 theorems) — whose only open inputs are **Vieta** and **THM-1550**.

## Coordination with death-star (the split)

death-star is driving THM-1550 and has **landed obstacle (i)** — `HenselianLocalRing (PowerSeries F)`
kernel-pure (`GMC2Henselian.lean`) — plus the monic M-th-root Hensel step, and refined past the
missing degree-dropping factorization theorem via the `a_j·Y_j` reparametrization (build the `M` small
roots individually; the `M`-th-root step is monic Hensel, the remainder is a PowerSeries fixed-point
contraction). Split agreed: **death-star** takes the fixed-point convergence + Vieta-for-`Π` +
Wiener–Hopf (`D_m=0 ⟺ Π=c·t`); **I** take the Gal-instantiation (done), **Vieta** for the full product
`hΩ`, and **Check A** (`CT(Λ^m)=[u^{Mm}]R^m`). My wrapper takes THM-1550 as exactly the two
hypotheses `hS`/`hfix`, so death-star's output slots straight in.

## Remaining of my pieces

- **Vieta** (`hΩ`: `∏ roots = (−1)^d r₀/r_d`, a constant): non-monic Vieta (`Φ = C(leadingCoeff)·∏(X−root)`
  ⟹ `∏ roots = (−1)^d coeff₀/leadingCoeff`; for `Φ = X^M − tR` the `t` cancels) + separability
  (`rootSet` product = roots-multiset product). Bounded Mathlib.
- **Check A** (`CT(Λ^m) = [u^{Mm}]R^m`): connects my `constantTermRelation` to `Polynomial.coeff (R^m)`;
  no `Polynomial.coeff_pow` in Mathlib, so this needs a custom coefficient-of-power lemma over the
  `piAntidiag`/`totalCharge` indexing. Bounded but nontrivial.

## Scope

Honest: the concrete Galois instantiation of THM-2067 is **complete and kernel-pure** — the transitive
action, tautological equivariance, and the full reduction for an irreducible `Φ` over `F(t)`. This
closes the "Galois wrapper" piece the owner named. Not full GMC(2): Vieta and Check A (mine, bounded)
and THM-1550 (death-star, actively progressing — obstacle (i) done) remain. Two checkpoints + this
close-out pushed; coordinated split with death-star.

Links: HYP-8956, HYP-8951 (wrapper, S235), HYP-8946 ((A), S234), HYP-8941 (core, S232), THM-2067,
THM-1550, death-star `GMC2Henselian.lean`,
[[the-thm2067-wrapper-assembled-and-the-two-remaining-hard-pieces-boxeph-S235]].
