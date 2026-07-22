# (A) done: `X^M − t·R` is irreducible over `ℂ(t)`, kernel-pure — and (B) the Hensel lift, scoped

*boxeph-2026-07-22-S234. Owner: work the hard bivariate/RatFunc Gauss (A) and the unramified Hensel
lift (B); complete them; pull often. Builds on my S232/S233 (orbit-product core + t-adic closing),
death-star-S106 (the THM-2067 map). New Lean `GMC2PhiIrreducible.lean`, kernel-pure.*

## (A) COMPLETE: irreducibility of `Φ = X^M − t·R(X)` over `F(t)`

This was the harder-than-expected piece, and it is now fully proved kernel-pure
(`#print axioms = [propext, Classical.choice, Quot.sound]`). `GMC2PhiIrreducible.phi_irreducible_ratfunc`:
for `R ∈ F[X]` with `R(0) ≠ 0` and `M ≥ 1`, `Φ = X^M − t·R(X)` is irreducible over `F(t)`. The route:

1. **`phi_t_irreducible`** — `Φ` viewed as a *linear polynomial in `t`* over `F[X]`,
   `Φ = C(−R)·t + C(Xᴹ)`, has coprime coefficients `Xᴹ` and `−R` (`IsCoprime` from `X ∤ R`, i.e.
   `R(0) ≠ 0`, via `prime_X.coprime_iff_not_dvd` + `.pow_left`), hence irreducible in `F[X][t]`
   (`Polynomial.irreducible_of_degree_eq_one_of_isRelPrime_coeff`).
2. **`swap_phi_eq`** — transport across the bivariate swap `Polynomial.Bivariate.swap : F[X][Y] ≃ₐ
   F[X][Y]` (`swap_C`, `swap_Y`): `swap Φ = map C (Xᴹ) − map C R · C X`, now with the outer variable
   `Y` playing the DvdK variable and the inner `F[X]` playing `F[t]`.
3. **primitivity** — `swap Φ` is primitive over `F[t]`: its `Yᴹ`-coefficient is `1 − C(r_M)·t` and its
   `Y⁰`-coefficient is `−C(r₀)·t` (computed by `coeff_map` + `coeff_mul_C`); any `r` dividing both
   divides `t` (since `r₀ ≠ 0` makes `C(r₀)` a unit) and then divides `1` (Bézout `C(r_M)·t +
   (1 − C(r_M)·t) = 1`), so is a unit.
4. **Gauss** — `IsPrimitive.irreducible_iff_irreducible_map_fraction_map` sends irreducibility in
   `F[t][X]` to irreducibility over `F(t) = RatFunc F`.

`Irreducible Φ` over `F(t)` is exactly the hypothesis of `Polynomial.Gal.galAction_isPretransitive`,
so this delivers the **transitivity** input consumed by my S232 orbit-product core. (A) is closed.

## (B) the unramified Hensel lift (THM-1550): scoped, foundations assessed

THM-1550 (death-star-S106): assuming all constant terms `D_m = 0`, the small-root product `Π = c·t`
(`t`-valuation 1). The unramified route: substitute `X = sZ`, `t = sᴹ`, so `Φ(sZ) = sᴹ·ψ(Z)` with
`ψ(Z) = Zᴹ − R(sZ) ∈ ℂ[[s]][Z]`; `ψ ≡ Zᴹ − r₀ mod s` is **separable** (Mathlib
`Polynomial.separable_X_pow_sub_C`, char 0, `r₀ ≠ 0`). Hensel then lifts the small factor `A(Z)`
(monic degree `M`, `A ≡ Zᴹ − r₀ mod s`), and `Π = t·(−1)ᴹ A(0)`, with `A(0) ≡ −r₀` — so `Π = c·t`
iff `A(0)` is constant iff all `D_m = 0`.

**Honest feasibility (probed against Mathlib):** the building blocks exist — `HenselianLocalRing`,
`IsAdicComplete.henselianRing`, `hensels_lemma`, `PowerSeries.exists_isWeierstrassFactorization`,
`separable_X_pow_sub_C` — but three real obstacles make (B) a *large, multi-session* formalization,
not a single-session completion:
- `HenselianLocalRing (PowerSeries F)` is **not** a free instance (`infer_instance` fails); it needs
  an `IsAdicComplete`/local-ring derivation for the maximal ideal `(s)`.
- The **degree-dropping** Hensel factorization of `ψ` (degree `d` in `Z`, reducing to degree `M`
  mod `s`) into `A·B` with `A` monic degree `M` — this is a Weierstrass-type factorization, not the
  plain coprime-lift `hensels_lemma`.
- The **Wiener–Hopf / generating-function** bridge `D_m = 0 ∀m ⟺ Π = c·t` (the `log`-of-factorization
  argument, S106 §2) — a separate analytic-combinatorial formalization.

This is genuinely THM-1550's owner death-star's hard piece; I've mapped the exact Mathlib entry points
and the sub-lemma structure so it can be attacked incrementally, and flagged it for coordination.

## State after this session

With (A) done, the path to full GMC(2) is: `GMC(2) ⟸ DvdK1 ⟸ THM-2067`, and THM-2067 now needs only
the **Galois wrapper** — Vieta (`v(C_Φ)=0`), equivariance (`Gal.smul_def`), Check A
(`CT(Λ^m)=[u^{Mm}]R^m`), instantiating the orbit-product core with transitivity from
`phi_irreducible_ratfunc` — plus **(B) THM-1550**. The wrapper pieces are mine and tractable; (B) is
the one deep analytic gap.

## Scope

Honest: (A) is **complete** and kernel-pure — a genuine milestone (the bivariate/RatFunc Gauss
irreducibility, which turned out to be the standard-but-intricate transitivity input). (B) is **not**
complete; it is scoped, its Mathlib building blocks identified, and its three real obstacles named — a
multi-session formalization I've set up for incremental attack and death-star coordination. Four
kernel-pure theorems; two checkpoints + this close-out pushed.

Links: HYP-8946, HYP-8942 (S233), HYP-8941 (S232 core), HYP-8935 (death-star S106 map), THM-2067,
THM-1550,
[[thm2067-abstract-contradiction-assembled-and-the-t-adic-closing-over-Ct-boxeph-S233]].
