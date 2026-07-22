# THM-2067: the abstract contradiction assembled, the t-adic closing over ℂ(t), and the precise remaining interface

*boxeph-2026-07-22-S233. Owner: long session, keep going until the GMC(2) formalization is complete;
pull from other agents as they supply the needed ideas. Builds on my S232 (orbit-product core),
codex THM-2067, death-star-S106 (the THM-2067 map + unramified-Hensel gap) and S107 (the
manufactured-valuation reading). New Lean in `GMC2OrbitProduct.lean` and `GMC2RatFuncClosing.lean`,
kernel-pure.*

## State: GMC(2) ⟸ DvdK1 ⟸ THM-2067, and where THM-2067 stands in Lean

`GMC2NC2.gmc2_of_dvdK1` reduces GMC(2) to the single input `DvdK1` (HeightWitnessSupplier discharged,
death-star-S105). Uniform `DvdK1` = codex THM-2067. My S232 formalized its abstract orbit-product core.
This session I assembled the full abstract contradiction and proved the concrete closing, and pinned
the exact remaining interface — two hard inputs I asked death-star for, plus a Galois-on-roots wrapper
I will build.

## New this session (kernel-pure)

- **`GMC2OrbitProduct.orbit_product_contradiction`** — the abstract THM-2067 contradiction, complete:
  a transitive `G`-action, equivariant `f`, a `G`-fixed subset product `p`, and a valuation `v` with
  `v(∏_Ω f) = 0` but `v(p) ≠ 0`, is `False`. This packages `prod_pow_card_group_eq` +
  `valuation_zero_of_prod_fixed` into the endpoint THM-2067 invokes.
- **`GMC2RatFuncClosing.monomial_pow_ne_const`** — the concrete closing over `ℂ(t)`: `a·tᴺ ≠ (const)`
  for `a ≠ 0`, `N ≥ 1`, by pulling back along the injection `ℂ[t] ↪ ℂ(t)` and comparing `natDegree`.

## The refinement that matters: the valuation is on ℂ(t), not the splitting field

In S232 I flagged "a valuation on the splitting field of ℂ(t)" as part of the gap. Re-reading THM-2067
(and death-star-S107's manufactured-valuation reading), that is **not** needed. The orbit-product
equation `Π^{|G|} = C_Φ^{|S|·|Stab|}` has **both sides in ℂ(t)**: `C_Φ = (−1)^d r_0/r_d` (Vieta, a
constant) and `Π = c·t` (THM-1550). So the contradiction is realized entirely inside ℂ(t) by the
elementary degree argument `monomial_pow_ne_const` — no valuation on the ramified splitting field, no
valuation API at all. The equation itself is derived in the splitting field (where `G` acts on the
roots), then pulled back to ℂ(t) via injectivity. This shrinks the gap to exactly THM-1550.

## The remaining interface (precise), and who supplies what

The wrapper instantiates `prod_pow_card_group_eq` at `G = Φ.Gal`, `Ω = Φ.rootSet E`, `f =` the root
inclusion. Mathlib supplies every structural hook (confirmed): `Polynomial.Gal.galAction_isPretransitive`
(transitivity from `Irreducible Φ`), `Polynomial.Gal.smul_def` (equivariance of the inclusion),
`MulDistribMulAction Φ.Gal Φ.SplittingField`, and the `Fintype Φ.Gal` / `Fintype (Φ.rootSet E)`
instances. So the abstract core genuinely plugs in. What is left:

| piece | content | owner |
|---|---|---|
| **(A) irreducibility** of `Φ = X^M − t·R` over ℂ(t) | Gauss + `gcd(X^M,R)=1` (`R(0)≠0`) ⟹ transitivity | asked **death-star** (hard Mathlib: bivariate/RatFunc Gauss) |
| **(B) THM-1550** `Π = ∏_{small} = c·t` | the unramified-Hensel small-root product (S106) | **death-star** — *the one true gap*; supplies `hfix` (Π ∈ ℂ(t)) **and** `v(Π)=1` |
| (C) Vieta `C_Φ = (−1)^d r_0/r_d` | product of roots = coeff ratio; `t` cancels ⟹ constant | me (Mathlib `Polynomial` roots) |
| (D) equivariance + wrapper | wire `galAction`/`smul_def`, instantiate the core | me |
| (E) Check A `CT(Λ^m) = [u^{Mm}]R^m` | connect `constantTermRelation` to the `R^m` coefficient / generating function | me (combinatorial) |

(A) and (B) are the hard/gap inputs; I sent death-star a precise request (they own THM-2067/THM-1550).
(C)–(E) are mine — gap-free but substantial Mathlib. Once (A)+(B) land as lemmas, the assembly is
short: instantiate (D), apply (C) for `v(C_Φ)=0`, `monomial_pow_ne_const` for the closing, and (E) to
land on the `DvdK1` interface's `constantTermRelation`.

## Scope

Honest: the GMC(2) formalization is **not** complete. This session assembled the full *abstract*
THM-2067 contradiction and proved the concrete `ℂ(t)` closing, both kernel-pure; refined the gap to
exactly THM-1550 (no splitting-field valuation); confirmed Mathlib has every hook for the Galois
wrapper; and handed death-star the two hard inputs (irreducibility, THM-1550) with a precise
interface. Three pushes this session-arc. The remaining work is the wrapper (mine, gap-free) plus the
death-star THM-1550/Hensel gap — the completion is gated on that one analytic piece.

Links: HYP-8942, HYP-8941 (S232 core), HYP-8935 (death-star S106 map), THM-2067, THM-1550,
[[the-thm2067-orbit-product-core-in-lean-the-gap-independent-heart-of-general-dvdk-boxeph-S232]].
