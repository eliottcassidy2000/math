# The THM-2067 orbit-product core in Lean: the gap-independent heart of general DvdK1

*boxeph-2026-07-22-S232. Owner: work a long session to fully formalize GMC(2); push/pull often. Builds
on codex's spine (`GMC2NC2`, and death-star-S105 discharging `HeightWitnessSupplier`), death-star-S106
(the THM-2067 formalization map, HYP-8935), codex THM-2067 (the Galois orbit-product proof of general
DvdK1). New Lean `GMC2OrbitProduct.lean`, kernel-pure.*

## Where the GMC(2) formalization actually stands

Reading the current spine top to bottom, the state is much sharper than I had it:

- **GMC(2) ⟸ NC2** — proved (`GMC2.gmc2_of_nc2`).
- **NC2 ⟸ DvdK1 + HeightWitnessSupplier** — codex's Frobenius/finite-field contradiction.
- **`HeightWitnessSupplier` is discharged** (death-star-S105, `GMC2NC2.heightWitnessSupplier_holds`,
  kernel-pure; the `whnf` blowup fixed by sealing `P.coeff` behind an opaque def). So there are now
  clean endpoints `GMC2NC2.nc2_of_dvdK1 : DvdK1 → NC2` and `gmc2_of_dvdK1`, with **no** height
  hypothesis.

So **GMC(2) now hinges on exactly one input: `GMC2DvdKInterface.DvdK1`** (for the lowest face
`f_F = Σ c_i u^{q_i}`, some `CT(f_F^m) ≠ 0`). My S226–S231 discharge it support-by-support (two-charge,
positive, unique-channel = 88%, monomial certificates for the residual 12%), but not *uniformly*
(effective DvdK is open). The uniform statement is codex **THM-2067** (Galois orbit-product), whose
formalization death-star-S106 mapped: **4 Mathlib-ready finite-algebra pieces + 1 valued-field gap**
(THM-1550, the small-root product, an unramified Hensel lift).

## What I formalized: the abstract orbit-product core (§5 target 1)

THM-2067's contradiction rests on a purely group-theoretic fact, gap-free and Mathlib-ready. I proved
it kernel-pure (`GMC2OrbitProduct.lean`, `#print axioms = [propext, Classical.choice, Quot.sound]`),
for a finite group `G` acting **transitively** on a finite `Ω`:

1. **`prod_smul_eq_prod_pow_card_stabilizer`** — the counting core:
   `∏_{g∈G} f(g•x) = (∏_{α∈Ω} f α) ^ |Stab_G x|` (via `Fintype.prod_fiberwise` + a fiber↔stabilizer
   bijection).
2. **`card_stabilizer_eq_card_stabilizer`** — stabilizer orders are constant on a transitive action
   (explicit conjugation bijection `s ↦ g₀⁻¹ s g₀`).
3. **`prod_pow_card_group_eq`** — the **orbit-product equation**: with a distributive `G`-action, `f`
   equivariant, and the subset product `p = ∏_{β∈S} f β` `G`-fixed,
   `p^{|G|} = (∏_{α∈Ω} f α)^{|S|·|Stab|}`.
4. **`valuation_zero_of_prod_fixed`** — the **contradiction engine**: for any additive valuation `v`,
   `v(∏_Ω f) = 0 ⟹ v(p) = 0` (apply `v` to the equation: `|G|·v(p) = |S|·|Stab|·0 = 0`, and `|G|>0`).

This is death-star-S106 §5 target 1, done. THM-2067 closes by (3)/(4) with `∏_Ω f = (−1)^d r_0/r_d`
(`t`-valuation 0) and `p = c·t` (valuation 1): `v(p)=1` contradicts `v(p)=0`, so some `CT(Λ^m) ≠ 0`.

## Why the abstract core is the right gap-independent boundary

I probed the instantiation on real roots. Mathlib has exactly the hook —
`Polynomial.Gal.galAction_isPretransitive p E hirr : IsPretransitive p.Gal (p.rootSet E)` — and the
distributive action `MulDistribMulAction p.Gal p.SplittingField` is an instance. So the *transitivity*
and *action* are free. But instantiating (3)/(4) on the roots needs a **valuation `v` on the splitting
field of `ℂ(t)`** (the roots and their products live there), i.e. an extension of the `t`-adic
valuation to a degree-`M` ramified extension. That extension **is** the one gap death-star-S106
isolates — and its S106 creative move (`X = sZ`, `t = s^M`; the small roots are an *unramified* Hensel
factor of `Z^M − R(sZ)` over `ℂ[[s]]`) is exactly what makes the valuation tractable. So the abstract
core is genuinely gap-independent, and the instantiation is where the (now-mapped) gap lives. Stopping
here keeps the kernel-pure boundary clean.

## Remaining to full GMC(2) (coordinating with death-star, who owns THM-2067)

1. **Galois wrapper** — instantiate (3)/(4) at `G = p.Gal`, `Ω = p.rootSet E`, `f =` inclusion:
   supply equivariance (`galAction` vs the field action) and `hfix` (subset product ∈ base field). The
   transitivity is `galAction_isPretransitive`.
2. **Check A** — `CT(Λ^m) = [u^{Mm}] R^m` (combinatorial; in my `constantTermRelation` wheelhouse).
3. **Irreducibility** of `X^M − tR(X)` over `ℂ(t)` (Gauss + `R(0)≠0`) ⟹ transitivity input.
4. **Vieta** `∏_Ω = (−1)^d r_0/r_d` (Mathlib `Polynomial` roots/coeff).
5. **THM-1550 / the gap** — the valuation on the splitting field via the S106 unramified-Hensel small
   root product. The one valued-field piece.

Sent a coordination note to death-star; I take the orbit-product core (done here) and can take the
Galois wrapper / Check A next; the Hensel gap (piece 5) is the shared hard target.

## Scope

Honest: this is **not** full GMC(2). It is the kernel-pure formalization of THM-2067's abstract
orbit-product core — the finite-Galois heart of the *sole* remaining GMC(2) input (`DvdK1`) now that
`HeightWitnessSupplier` is discharged — together with its valuation contradiction engine, and a
verified read that the remaining instantiation reduces to the single (S106-mapped, unramified-Hensel)
valued-field gap. Four kernel-pure theorems; two mid-session checkpoints pushed.

Links: HYP-8941, HYP-8935 (death-star S106 map), THM-2067, THM-1550, THM-2022,
[[eliminating-dvdk-for-the-residual-12-percent-via-monomial-certificates-boxeph-S231]],
[[bypassing-the-gmc2-dvdk-dependency-for-the-unique-channel-class-boxeph-S230]].
