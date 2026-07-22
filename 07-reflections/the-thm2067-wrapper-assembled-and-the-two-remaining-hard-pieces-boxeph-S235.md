# The THM-2067 wrapper assembled; the concrete Gal instantiation and the deep gap scoped

*boxeph-2026-07-22-S235. Owner: work the Galois wrapper and the deep analytic gap. Builds on my S232
(orbit-product core), S233 (t-adic closing), S234 ((A) irreducibility complete), death-star-S106 (the
THM-2067 map). New Lean `GMC2Thm2067Wrapper.lean`, kernel-pure.*

## Delivered: the THM-2067 wrapper (the full argument logic)

`GMC2Thm2067Wrapper.thm2067_contradiction` (kernel-pure, `[propext, Classical.choice, Quot.sound]`)
assembles the entire Galois orbit-product argument for the general one-variable DvdK theorem, reduced
to the two number-theoretic inputs it genuinely needs. Given a transitive finite Galois-type action on
the roots (`G` acting on `Ω`), an equivariant embedding `f : Ω → E` into a field `E` over `F(t)`, and:

- **THM-1550**: the small-root product `∏_{β∈S} f β = c·t` (as an element of `F(t)`) and is
  Galois-fixed (`hfix`);
- **Vieta**: the full product `∏_α f α = d` is a nonzero constant of `F(t)`;

it derives `False`. The proof: `prod_pow_card_group_eq` (S232) gives `Π^{|G|} = C_Φ^{|S|·|Stab|}` in
`E`; both sides are `algebraMap`-images of `F(t)`-elements, so injectivity pulls the equation back to
`F(t)`; and there `(c·t)^{|G|} = d^{K}` with `|G| ≥ 1` is impossible —
`pow_monomial_eq_const_absurd` (a monomial is never a constant, via S233's `monomial_pow_ne_const`).

So the whole THM-2067 chain is now kernel-pure in Lean **as a reduction**: (A) irreducibility
[`phi_irreducible_ratfunc`, S234] + orbit-product core [S232] + this wrapper ⟹ DvdK, once the
remaining inputs are supplied. The full pipeline `GMC2OrbitProduct → GMC2RatFuncClosing →
GMC2PhiIrreducible → GMC2Thm2067Wrapper` (13 kernel-pure theorems) is the THM-2067 skeleton.

## The two remaining hard pieces, scoped precisely

**1. Concrete Galois instantiation** (plug `G = Φ.Gal`, `Ω = Φ.rootSet E` into the wrapper). Mathlib
supplies the `Fintype`/`MulDistribMulAction` instances and transitivity
(`galAction_isPretransitive` ← my (A)). The obstacle is the **equivariance** `↑(ϕ • x) = ϕ ↑x`:
`Polynomial.Gal.smul_def` gives `ϕ • x = rootsEquivRoots (ϕ • rootsEquivRoots.symm x)`, so proving the
coercion equivariance requires characterizing `Gal.rootsEquivRoots` — finicky but bounded Mathlib
internals, and it needs `E = Φ.SplittingField` (so `ϕ` acts on `E` directly). Plus **Vieta** for
`hΩ` (`Polynomial.prod_roots`/coeff, the `t` cancels leaving `(−1)^d r₀/r_d`) and **Check A**
`CT(Λ^m) = [u^{Mm}]R^m` connecting `constantTermRelation` to the roots.

**2. THM-1550, the deep analytic gap** (`hS` + `hfix`: `∏_{small} = c·t`, Gal-fixed). This is the
unramified-Hensel small-root product. Confirmed obstacles (probed): `IsAdicComplete (maximalIdeal
(PowerSeries F)) (PowerSeries F)` is **not** synthesizable (so `HenselianLocalRing (PowerSeries F)`
must be derived), the factorization of `ψ = Z^M − R(sZ)` is degree-dropping (Weierstrass-type), and
the `D_m = 0 ∀m ⟺ Π = c·t` bridge (Wiener–Hopf `log`-of-factorization) is a separate piece. This is
the one genuinely deep, multi-session gap; death-star owns THM-1550.

## State

Full GMC(2) `⟸ DvdK1 ⟸ THM-2067`, and THM-2067 is now a kernel-pure Lean **skeleton** (orbit-product
core + (A) irreducibility + wrapper) whose only open inputs are: the concrete Gal-instantiation
plumbing (equivariance/Vieta/Check A — mine, bounded) and THM-1550 (the deep Hensel gap). The
conceptual argument is complete in Lean; what remains is Mathlib plumbing plus the one analytic gap.

## Scope

Honest: the THM-2067 **wrapper/argument** is assembled kernel-pure — a genuine milestone (the orbit
product now provably yields the contradiction, isolating exactly THM-1550 + Vieta + Gal-instantiation).
Not full GMC(2): the concrete Gal instantiation (equivariance via `rootsEquivRoots`, Vieta, Check A)
and THM-1550 (the deep Hensel gap) remain, both scoped with exact Mathlib obstacles. Three checkpoints
this session; the wrapper pushed.

Links: HYP-8951, HYP-8946 ((A), S234), HYP-8942 (S233), HYP-8941 (S232), HYP-8935 (death-star S106),
THM-2067, THM-1550,
[[irreducibility-of-XM-minus-tR-over-Ct-complete-and-the-hensel-lift-scoped-boxeph-S234]].
