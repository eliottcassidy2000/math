---
id: THM-494
name: The transverse-resonance angles are exactly the rational-cosine rotations;
      the perfect 30° bisector ℤ[ζ₁₂] is OFF the Moser ladder and transverse-free
status: PROVED (from THM-434 + Kronecker + Niven) + VERIFIED exact-integer
        (densest-patch search, n=21..30)
date: 2026-06-13
session: monad-explorer-2026-06-13
depends_on:
  - THM-434   # #units(L_t)=12+r_E(t); transverse vectors α(1−ω_t), |1−ω_t|²=1/t
  - THM-493   # resonant-product decomposition; Δ_t bonus = the crossing
  - HYP-2461  # tie-vs-crossing dichotomy; the t=4 control; the ℤ[ζ₁₂] handoff question
resolves:
  - "HYP-2461 next-explorer question: does the exact-30° lattice ℤ[ζ₁₂] cross 3N at n=28?"
external:
  - "Niven, *Irrational Numbers* (Carus Monograph 11, 1956), Cor. 3.12: if θ/π∈ℚ and
     cosθ∈ℚ then cosθ∈{0,±½,±1}."
  - "Kronecker 1857: an algebraic integer all of whose conjugates have absolute value
     1 is a root of unity."
---

# THM-494: transverse resonance ⟺ rational cosine; the bisector is off-ladder

## Context

THM-434 counts the unit vectors of the Moser-ladder bridge lattice
`L_t = ℤ[ζ₆] ⊕ ω_t·ℤ[ζ₆]` (two triangular lattices glued at the Moser angle
`ω_t`, `cos = (2t−1)/2t`): `12 + r_E(t)`, the extra `r_E(t)` being the **transverse**
unit vectors `α(1−ω_t)` with Eisenstein norm `N(α)=t`. THM-493 shows those transverse
vectors are exactly the **resonance bonus** `Δ_t` that carries `u(28) ≥ 85 > 84`
across `3N`. HYP-2461's free-patch search + the `t=4` control localized the crossing
to `t=3` (`√−11`) and asked the natural sharpening:

> Does the **exact-30° bisector** lattice `ℤ[ζ₁₂]` — the geometrically *perfect*
> interleaving angle, which the `L_t` family brackets (`t=3`→33.6°, `t=4`→29.0°)
> but can never hit — cross `3N` at `n=28`? If even it fails, `√−11` is
> *arithmetically* singular, not merely off the geometric optimum.

This theorem answers it: **No.** And it explains why with a clean characterization
of *which* gluing angles can carry transverse unit vectors at all.

## Setup

For a unit `ω = e^{iθ} ∉ ℚ(√−3)` glue two Eisenstein copies into the rank-4 lattice
`L_ω = ℤ[ζ₆] ⊕ ω·ℤ[ζ₆]`. A **transverse** unit vector is `z = α + ωβ`,
`α,β ∈ ℤ[ζ₆]`, *both nonzero*, with `|z| = 1` (a unit vector using both copies — as
opposed to the two trivial rosettes `ζ₆^j` and `ζ₆^j ω`, always present).
`ℤ[ζ₁₂] = ℤ[ζ₆] ⊕ ζ₁₂·ℤ[ζ₆]` is the case `ω = ζ₁₂ = e^{iπ/6}` (basis
`1, ζ₁₂, ζ₁₂², ζ₁₂³`; `ζ₁₂² = ζ₆`), the perfect 30° bisector.

## Statement

**(A) The diagonal transverse family pins the angle.** The diagonal vector
`α(1−ω)` is a unit vector iff
```
        |1 − ω|²  =  2 − 2cosθ  =  1/N(α),     i.e.   cosθ = (2N(α)−1)/(2N(α)).
```
So the *only* relative rotation angles admitting a diagonal transverse unit vector
of "radius" `α` are the **Moser angles** `cosθ = (2t−1)/2t`, `t = N(α)` a Loeschian
integer. The Moser ladder is precisely the set of **rational-cosine** rotations
`{θ : cosθ = (2t−1)/2t}` — equivalently `{θ : |1−ω|² ∈ ℚ}` realised at `1/t`.

**(B) Niven obstruction.** For `t ≥ 2`, `cosθ_t = (2t−1)/2t ∉ {0, ±½, ±1}`, so by
**Niven's theorem** `θ_t / π` is **irrational**. Contrapositive: every rotation by a
*rational* multiple of `π` (other than the degenerate `θ = π/3`, the triangular
self-gluing `t=1`) has **irrational** cosine and is therefore **not a Moser angle**.
In particular the cyclotomic "clean-angle" gluings
```
   ℤ[ζ₁₂]  (θ = π/6 = 30°,  cosθ = √3/2),       ℤ[ζ₈]  (θ = π/4 = 45°,  cosθ = √2/2)
```
are **off the ladder**: `|1−ζ₁₂|² = 2−√3` and `|1−ζ₈|² = 2−√2` are irrational, so no
scalar multiple of `α(1−ω)` is ever a unit vector.

**(C) ℤ[ζ₁₂] is transverse-free: exactly 12 unit vectors.** Every `z ∈ ℤ[ζ₁₂]` with
`|z| = 1` is a 12th root of unity; the 12 split as `μ₆` (even powers, one hexagon)
and `ζ₁₂·μ₆` (odd powers, the 30°-rotated hexagon), with **no** vector mixing the
two. Hence `#units(ℤ[ζ₁₂]) = 12 = 6 + 6 + 0` — the THM-434 shape `12 + r_E` with
`r_E = 0`, structurally identical to the **non-Loeschian** rungs `t = 2, 5`
(`√7, √19`), even though those are *on* the ladder and `ℤ[ζ₁₂]` is *off* it.

**(D) [VERIFIED, exact integer] The bisector does not tie and does not cross.** The
identical exact-integer densest-patch search that reproduces Engel on `L_3` gives,
for `ℤ[ζ₁₂]`, the deficit profile (`n = 21…30`):
```
   n     21 22 23 24 25 26 27 28 29 30
   u≥    57 59 62 66 70 74 78 83 85 88
   3n    63 66 69 72 75 78 81 84 87 90
   u−3n  −6 −7 −7 −6 −5 −4 −3 −1 −2 −2
```
At `n = 27` it caps at **78** (deficit `−3`) — it **cannot even build the 81 tie**;
at `n = 28` it caps at **83 = P(28)** (the generic product cap, deficit `−1`) and
**does not cross**. This is *bit-for-bit identical* to the transverse-free `t = 2,5`
rungs (78@27, 83@28 in the HYP-2461 table).

## Proof

**(A)** `|α(1−ω)|² = N(α)·|1−ω|² = N(α)(2−2cosθ)`; set `= 1`. ∎
(Matches THM-434's `|1−ω_t|² = 1/t` at `cosθ = (2t−1)/2t`.)

**(B)** Niven's theorem directly; `(2t−1)/2t ∈ {0,±½,±1}` only at `t=1` (`=½`). ∎

**(C)** Two arguments, agreeing.
*Kronecker.* `ℤ[ζ₁₂]` is the ring of integers of the CM field `ℚ(ζ₁₂)`, real subfield
`ℚ(√3)`; complex conjugation `z ↦ z̄` is its nontrivial automorphism. If `|z|=1` then
`z z̄ ∈ ℤ[√3]` and maps to `1` under the chosen real embedding; since `√3` is
irrational, the only element of `ℤ[√3]` with that image is `1`, so `z z̄ = 1`
*exactly*. Then **both** archimedean places give `|z| = 1`, so every conjugate of `z`
has absolute value 1, and by Kronecker `z` is a root of unity, i.e. `z ∈ μ₁₂`. ∎
*Direct (the THM-434 Step-1 mechanism at `ω=ζ₁₂`).* Write `z = α + ζ₁₂ β`,
`γ := ᾱβ = p + qζ₆ ∈ ℤ[ζ₆]`. With `cos30° = √3/2`, `sin30° = 1/2`,
`Re(γ ζ₁₂) = (p+q/2)(√3/2) − (q√3/2)(1/2) = (√3/2)p`, so
`|z|² = N(α) + N(β) + p√3`. For `|z|² = 1 ∈ ℚ` and `√3` irrational, `p = 0`, leaving
`N(α) + N(β) = 1`. If both `α,β ≠ 0` then `N(α)+N(β) ≥ 2 > 1` — impossible. So one of
`α,β` vanishes: **no transverse unit vector**, and the units are the two rosettes,
`12` total. ∎ (Brute-force confirmed: the only solutions of the exact integer test in
a box are the 12 listed offsets.)

**(D)** `04-computation/unit_distance_zeta12_bisector_monad.py`
(`05-knowledge/results/…out`). Engine copied verbatim from
`unit_distance_bridge_lattice_family_monad.py`; adjacency is the exact integer test
`P₁Q₁+P₂Q₂ = 0 ∧ P₁²+3Q₁²+P₂²+3Q₂² = 4` (`P₁=2a+c, Q₁=b, P₂=b+2d, Q₂=c`); every
reported count independently exact-recounted; the triangular sub-hexagon calibrates
to the exact triangular maximum (63@27, 65@28). ∎

## Consequences (for OPEN-Q-057 / N*)

1. **The handoff is resolved, with a reframe.** `ℤ[ζ₁₂]` does **not** cross at `n=28`,
   so the `u(28)=85` crossing is **arithmetic, not geometric** — confirming
   HYP-2461. But the reason is *not* that 30° is the wrong angle in a band: it is that
   **the perfect bisector is off the resonance ladder entirely** (irrational cosine),
   so it carries **zero** transverse vectors and caps exactly like the transverse-free
   `t=2,5`. The "perfect bisector" is among the **worst** carriers, not the best.

2. **Geometric optimality and arithmetic admissibility are disjoint (Niven).** The
   crossing lives on **rational-cosine, irrational-angle** rotations (the Moser
   ladder). The geometrically clean **rational-angle** rotations (`π/6, π/4`, the
   cyclotomic lattices) are exactly the ones Niven forbids from the ladder. You cannot
   have both; the unit-distance crossing sits firmly on the rational-cosine side.

3. **Third independent confirmation of THM-493.** The bisector achieves precisely the
   generic product cap `P(28)=83` at `n=28` and falls short by **exactly** the
   resonance bonus `Δ₃ = 2`. So `ℤ[ζ₁₂]` is "the product at a non-resonant rational
   angle" — generic count, zero bonus — independently corroborating the
   THM-493 decomposition `U = (product) + Δ_t` from a carrier the `L_t` family cannot
   reach.

## Scope / honesty

- (A),(B),(C) are PROVED. (D) is an exact-integer *lower-bound* (annealing) search:
  a value below the transverse-free profile would be "not found," not a proved ceiling
  — but the agreement with `t=2,5` to the unit across all `n=21..30`, the
  calibration to the exact triangular maximum, and the structural transverse-freeness
  (C) make the cap robust. This does **not** prove `u(27)=81` (the AMP upper bound is
  still 90); it removes `ℤ[ζ₁₂]` as a crossing candidate and pins *why*.
- Open: is `√−11` (`t=3`) truly unique among **all** Loeschian rungs for the `n=28`
  crossing (HYP-2461 found only `t=3` crosses among `t=3,4,13,21,31,49`)? THM-493
  attributes it to `28 = 4·7` admitting a `√3`-bearing edge-dense factor pair; a
  clean arithmetic characterization of the crossing `t` (vs the merely transverse `t`)
  is the next target.
