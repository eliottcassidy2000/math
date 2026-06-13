---
source: monad-explorer-2026-06-13 (deep-research; dispatched seed = "prove u(21)",
  already settled by AMP — ranged to the live N* frontier and the THM-433/THM-434 seam)
status: REFLECTION on THM-493 (proved + exact-verified). One operation — the
  Minkowski product of two triangular patches — has a GENERIC face (THM-433: avgdeg
  additive) and a RESONANT face (THM-434: transverse unit vectors). The "non-product"
  3N-crossover is the resonant face, and the resonance bonus is literally the crossing.
tags: [unit-distance, moser-lattice, resonant-product, THM-493, THM-433, THM-434,
  N-star, 3N-crossover, eisenstein, transverse-vectors, everything-is-the-triangle,
  two-CM-tower, lrc-clock-shell, bonus-hostile]
---

# The resonance bonus IS the crossing — and 27 is bonus-hostile

## The seam between two theorems

Two of this project's unit-distance theorems sat next to each other without touching:

- **THM-433**: the *generic-angle* Minkowski product `G □ H` has `avgdeg` additive;
  the product family ties `3N` at `27, 30` and first beats it at `32`; hence the true
  crossover `N* ∈ [25,28]` is **"non-product"**.
- **THM-434**: the Moser-ladder lattice `L_t = ℤ[ζ₆] ⊕ ω_t·ℤ[ζ₆]` has exactly
  `12 + r_E(t)` unit vectors — `12` rosette + `r_E(t)` **transverse**.

THM-493 closes the seam. `L_t` is *literally* the Minkowski product of the triangular
lattice with a copy of itself rotated by the **Moser angle** `ω_t`. So the two
theorems describe **one operation at two angles**:

```
   product at a GENERIC angle   →  Cartesian product, Δ = 0          (THM-433)
   product at the RESONANT ω_t  →  Cartesian product + transverse Δ  (THM-434)
```

and the exact count is

```
   U(G ⊞_t H) = e(G)|H| + |G|e(H) + Δ_t(G,H),
   Δ_t = ½ Σ_{N(α)=t} m_α(G) m_α(H)   (m_α = # of √t-separations α in a factor).
```

The bonus `Δ_t` is a **correlation of the two factors' norm-`t` displacement
spectra**. It is nonzero only when *both* factors carry a matching `√t`-pair — so the
abstract "transverse unit vectors of `L_t`" become concrete extra edges you can see
and count.

## The crossing, made of two edges

The sharpest consequence is a single construction. Glue the Eisenstein rosette
`W₇` (7 pts, 12 edges) to the unit rhombus `R` (4 pts, 5 edges) at the Moser angle
`ω₃ = (5+i√11)/6`:

```
   U(W₇ ⊞₃ R) = 12·4 + 7·5 + Δ₃ = 48 + 35 + 2 = 85   on 28 points   (> 84 = 3·28).
```

The **same 28-point product graph** has only `83` unit distances at a generic angle —
exactly the THM-433 product cap `P(28)=83 < 84`. The resonant angle adds precisely
the `Δ₃ = 2` transverse edges, and **`83 < 84 < 85`**: those two edges are the *entire*
crossing of `3N`. This is `u(28) ≥ 85` (Engel's bound, THM-431's cited ceiling for
`N*`) reproduced as an explicit resonant product, exact-integer over `ℚ(√3,√11)`.

So THM-433's "the crossover is non-product/irreducible" sharpens to:

> The first construction to beat `3N` is the *same product* that realises the product
> cap, evaluated at the **resonant angle**; the resonance bonus is what carries it
> across the line. "Non-product" meant "product, plus the part that isn't additive."

## Why 27 holds the line: the edge-density / bonus-richness tension

THM-433 explained the `27 = 3³` tie as `avgdeg(K₃^□3) = 2+2+2 = 6` exactly. THM-493
explains why the tie is *robust against the resonance bonus that breaks 28* — a
strictly finer fact. Two resources compete inside a factor:

- **edge-density** — unit (norm-1) separations, the only thing the generic product
  rewards. The unit triangle `K₃` is the extremal atom: `avgdeg = 2 = κ/3`, the
  homomorphism extreme of THM-433.
- **bonus-richness** — norm-`t` (`√t`) separations, the only thing the resonant
  product rewards. `K₃` has **none** (all its separations are unit).

These pull against each other: a `√t`-rich patch spends structure on long
separations and loses unit edges. Now look at the factorizations:

```
   27 = 3·9 = 3·3·3   →   EVERY factorization forces a size-3 factor.
```

The densest 3-point unit-distance graph is `K₃`, which is **√t-free for every `t`**.
So at 27 the bonus-optimal route (`K₃`-free factors) doesn't exist and the
edge-optimal route (`K₃^□3`) earns **zero bonus** — it ties at 81 and cannot be
lifted. **27 is bonus-hostile by arithmetic**: `3³` admits no factor that is both
edge-dense and `√t`-bearing. The exact search confirms it — no two-factor resonant
product reaches 81 at 27 (best is 75), and `K₃^□3 = 81` carries bonus 0:

```
   n      24  25  26  27  | 28
   3n     72  75  78  81  | 84
   best   68  72  61  75  | 85   ← first composite whose factorization (4·7)
   bonus   2   2   0   0  |  2     admits a √3-bearing, edge-dense factor pair
```

`28 = 4·7` is the first composite past 27 with a factorization (`rhombus × W₇`) where
one factor (the rhombus) carries a `√3`-pair *while staying edge-dense* and the other
(`W₇`) is `√3`-rich. That is the structural reason the crossover lands at 28 and not
before — and it is strong (constructive) evidence for the standing conjecture
`u(27) = 81`, `N* = 28`.

## Where this points (honestly labelled)

- **PROVED:** the decomposition `U = e(G)|H|+|G|e(H)+Δ_t` (from THM-434), the
  injectivity, `u(28) ≥ 85` as `W₇ ⊞₃ R`, and the 27 search (a lower bound over a
  curated patch family — evidence, not a proof of `u(27)=81`).
- **OPEN:** an *upper* bound `u(27) ≤ 81` would settle `N* = 28`. THM-493 reframes the
  target: show no 27-point set beats the `K₃^□3` tie — and the bonus-hostility of `3³`
  says the obstruction is arithmetic (no edge-dense `√t`-factor at size 3), not merely
  geometric.
- **CONJECTURE (cross-domain transfer):** the S4 Moser reflection flagged that *two*
  lanes found "a product of two CM/cyclotomic pieces glued together" — triangular
  `⊕ √−11` here, clock `× shell` (`ℤ/n × ℤ/(2n−1)`, THM-427) in the LRC worry-set.
  THM-493 says the geometry-side gluing is a **resonant product with a correlation
  bonus**. The natural LRC question: does the worry-set count also split as
  `(clock count)(shell count) + (a correlation bonus that fires when clock and shell
  displacements are matched)`? If so, "irregularity beats the symmetric floor" (the
  LRC analogue of beating `3N`) would again be a *resonance bonus*, not a new
  mechanism. Worth a dedicated probe — the convergence is too clean to be nothing.

The mathematics handed us `u(28)=85` through a citation; the resonant-product lens
shows it is two edges of bonus on a product graph, and the same lens explains why 27
refuses the same trick. Follow the bonus.
