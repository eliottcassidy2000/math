# Closing the THM-1550 Hensel gap: obstacle (i) done kernel-pure, obstacle (ii) refined past the factorization wall

**death-star-2026-07-22-S111** (HYP-8960). Owner: work to close the Hensel gap (THM-1550), pushing/pulling often
and treating agent pulls as connection signal. **The full THM-1550 is not closed** (it is multi-session); this
session delivers two kernel-pure Lean foundations and *refines* the hardest obstacle so it no longer needs a
Mathlib theorem that does not exist. Honest throughout, in the spirit of MISTAKE-241 (which correctly flagged my
S106 numerics as unproven — the Lean here is actually checked).

## Context (coordinated with boxeph, live)

The GMC(2) DvdK residual = THM-2067; boxeph's `GMC2OrbitProduct` kernel-checks its abstract orbit-product
contradiction, and this session's pulls showed **boxeph closed piece (A)** — `Φ = X^M − t·R(X)` irreducible over
`F(t)`, kernel-pure (`GMC2PhiIrreducible`, HYP-8946) — the transitivity input. That leaves **(B) = THM-1550**
(the small-root product `Π = c·t`), which I own. boxeph and I independently mapped the same three sub-obstacles:
(i) `HenselianLocalRing (PowerSeries F)` is not a free instance; (ii) the small-factor extraction is
degree-dropping; (iii) the Wiener–Hopf product identity `Π = c·t` under CT-vanishing.

## Delivered (kernel-pure, `GMC2Henselian.lean`, both `[propext, Classical.choice, Quot.sound]`)

1. **Obstacle (i) — `powerSeries_henselianLocalRing : HenselianLocalRing (PowerSeries F)`.**
   `HenselianRing (F⟦X⟧) (span {X})` is free (`IsAdicComplete.henselianRing` + the `(X)`-adic-completeness
   instance); `maximalIdeal_eq_span_X` bridges to the maximal ideal; the derivative-unit hypothesis transfers
   along `Ideal.Quotient.mk`. This is the Henselian foundation the whole S106 lift stands on.

2. **The refined-route key — `exists_pow_eq_of_constantCoeff_pow`:** if `a₀^M = constantCoeff u` then `u` has an
   `M`-th root `Y` in `F⟦X⟧` with `constantCoeff Y = a₀`. Proof: `X^M − C u` is **monic**, `a₀` is a simple root
   mod `(X)` (derivative `M·a₀^{M−1}` a unit), so the Henselian simple-root property lifts it.

## The refinement of obstacle (ii) (the real content)

I first hoped the M small roots — being *simple* roots of the separable reduction `Z^M − r₀` — could be lifted
directly. **That is wrong, and I correct it:** the Henselian simple-root property requires a **monic** `f`, and
`ψ(Z) = Z^M − R(sZ)` is non-monic (genuine `Z`-degree `d`, leading `−r_d s^d`), with `ψ mod s = Z^M − r₀` dropping
the degree `d → M`. I verified `HenselianLocalRing.TFAE` has **no factorization item** — only the three
simple-root variants — so Mathlib genuinely lacks a Henselian *factorization* theorem, and obstacle (ii) as boxeph
framed it (extract the degree-`M` factor) needs that missing theorem (a real, ~multi-session Mathlib contribution).

**The fix removes that need.** Do not factor `ψ`; build the `M` small roots one at a time. Each is
`Z_j = a_j · Y_j`, where `a_j` ranges over the `M` distinct `M`-th roots of `r₀` (in `ℂ`, algebraically closed),
and `Y_j` is a *principal* unit (`constantCoeff 1`) solving `Y^M = R(s·a_j·Y)/r₀`. The `M`-th-root step is
**monic** Hensel (`Z^M − C u` is monic) — that is lemma 2, now done. So the refined THM-1550 route is:

```
  [ monic M-th roots of units ]     ← DONE (exists_pow_eq_of_constantCoeff_pow)
+ [ fixed point Y = (R(s·a·Y)/r₀)^{1/M} ]   ← PowerSeries contraction, converges by adic completeness
+ [ Vieta: Π = t·(∏ a_j)·(∏ Y_j),  ∏ a_j = (−1)^{M+1} r₀ ]
+ [ (iii) Wiener–Hopf: ∏ Y_j = ±1 (constant) ⟺ all D_m = 0 ⟺ Π = c·t ]
```

**No general Henselian factorization theorem is needed anymore** — that is the payoff of the `a_j·Y_j`
reparametrization. The remaining pieces are the fixed-point convergence (a manual `(s)`-adic contraction, tractable
next) and the genuinely analytic Wiener–Hopf identity (iii), which is THM-1550's core and stays the deep gap.

## Honest status

- **Closed this session:** obstacle (i) (`HenselianLocalRing (PowerSeries F)`) and the monic `M`-th-root lemma —
  both kernel-pure, built into the project.
- **Refined, not closed:** obstacle (ii) no longer needs a factorization theorem; it needs the fixed-point
  convergence (next) — a much smaller, well-scoped target.
- **Still open:** (iii) the Wiener–Hopf `Π = c·t` identity (my S106 §2 formal-log argument), the analytic core of
  THM-1550. And assembling all of it into boxeph's Galois wrapper. THM-1550 / general `DvdK1` remains open.
- **Correction:** my earlier "simple roots circumvent the degree drop" was wrong (needs monic); the working
  circumvention is `a_j·Y_j` + monic `M`-th roots.

Cross-links: S106 (the route this makes Lean-real — HYP-8935/MISTAKE-241, now with actually-checked pieces),
boxeph S232/S233 (`GMC2OrbitProduct` abstract core) and S234 (piece (A) irreducibility, `GMC2PhiIrreducible`),
THM-2067 (codex, the residual), THM-1550 (the target). Files:
`04-computation/lean/TournamentH7/TournamentH7/GMC2Henselian.lean` (2 kernel-pure theorems). HYP-8960.
