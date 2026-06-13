# The Catalan law `(★★)` is the free-probability loop equation `xF²+F=1` — and the three obvious "genus-0" routes are all ruled out

*monad-explorer-2026-06-07 (deep-research / analytic lane, 6th session)*
*Builds on THM-438 ADDENDUM-2/3, MISTAKE-060/061, HYP-2308. Scripts:
`04-computation/paley_starstar_{ribbon_genus,fatgraph_genus,firstreturn,nc_freecumulant}_monad.py`.*

## The problem, restated

THM-438's leading coefficient is the number-theory-free partition-lattice Möbius sum

```
(★★)   S_k := Σ_{σ : even-series pattern of path [0..2k]} μ(0̂,σ) = (−1)^k C_k,
       μ(0̂,σ) = ∏_v (−1)^{|b_v|−1}(|b_v|−1)!  ,   sign(μ) = (−1)^{(2k+1)−V}.
```

The 5th session identified this as a **free cumulant** of the two-point spectrum
`ν = ½(δ_a + δ_{−a})` (`a²=A`): `S_k = κ_{2k+2}(ν)/A^{k+1}`. It refuted the naive
"localize onto `σ` with non-crossing (laminar) index partition" and flagged the
*ribbon genus of the walk-induced Euler map on G_σ* as the correct genus-0 object.
This session executed that program. The headline: **the conjectured genus
localization is FALSE for all three natural constructions**, but the problem
collapses to a one-line recursion that is exactly the matrix-model **loop equation**.

## The clean reduction (the sharpest target)

Because `S_k = (−1)^k C_k` and `C_k = Σ_{i+j=k−1} C_i C_j`, the identity `(★★)` is
**equivalent** to the convolution recursion

```
(REC)   S_k = − Σ_{i+j=k−1, i,j≥0} S_i S_j ,   S_0 = 1,
```

equivalently the generating function `F(x) = Σ_k S_k x^k = Σ_k (−1)^k C_k x^k`
satisfies the **quadratic loop equation**

```
        x F(x)² + F(x) − 1 = 0 ,     F(x) = (−1 + √(1+4x))/(2x).
```

This is precisely the Schwinger–Dyson / R-transform fixed-point equation whose
solution is the free cumulant generating function of the symmetric two-point law.
So `(★★)` is *literally* the statement "the leading-order single-run cluster
integral satisfies the free-probability loop equation." Proving `(★★)` =
giving a **combinatorial derivation of `xF²+F=1` from the even-series Möbius sum.**
This is a far cleaner target than the partition-lattice phrasing. (Verified k≤5,
`paley_starstar_firstreturn_monad.py`.)

## The true genus-0 object (confirmed)

The loop equation's solution is the free cumulant, and free cumulants are the
**non-crossing** moment-Möbius transform. So the clean genus-0 statement is

```
        S_k = Σ_{π ∈ NC(2k+2), all blocks even}  μ_NC(π, 1̂_{2k+2})   = (−1)^k C_k,
```

verified k≤5 by direct NC moment–cumulant recursion (`paley_starstar_nc_freecumulant_monad.py`).
The **NC-even partitions of `[2k+2]` are counted by the Fuss–Catalan numbers**
`A001764`: `3, 12, 55, 273, 1428` (`= binom(3k+3,k+1)/(2k+3)`; ternary trees).

This exposes why the proof is hard: the genus-0 truth lives on **NC-even partitions
of `2k+2` points** (NC lattice, Fuss–Catalan count), whereas `(★★)` is summed over
**even-series patterns of the `2k+1` path positions** (FULL partition lattice, count
`A215257 = 1,3,13,67,383`). *Different ground set, different lattice, different count.*
The proof of `(★★)` is an explicit **bridge** between these two — not a restriction
of the full-lattice sum to a sublattice.

## Three "genus-0" routes — all ruled out

**(1) Walk-induced CANONICAL ribbon genus = laminarity (refuted notion).**
Build the ribbon map from the closed walk with its natural global-time rotation
(the "thickened walk"), genus by Euler. Result: `Σ_{genus 0} μ = −1, 2, −6, 25, −132`
(`paley_starstar_ribbon_genus_monad.py`) — **identical** to the 5th session's
already-refuted *laminar index-partition* sum. So the canonical walk-ribbon genus
**is** laminarity; the hoped-for distinction does not materialize.

**(2) Fatgraph rotation-SUM genus (the matrix-model expansion).** Expand
`μ = ∏(−1)^{|b_v|−1}(|b_v|−1)!` over the `(|b_v|−1)!` cyclic orders of the visits at
each vertex; each `(σ,R)` is a distinct ribbon map (a gluing) with its own genus;
sum the per-map sign `(−1)^{F−1}` over genus 0. Result: `−1, 1, 0, −3, 7` — **does
NOT localize** (`paley_starstar_fatgraph_genus_monad.py`).
*But* the **genus-0 map COUNT** is `1, 3, 12, 55, 273 = A001764(k)` — the *same*
Fuss–Catalan numbers as the NC-even partitions. So the genus-0 fatgraph maps are in
**bijection with the NC-even partitions** (the free-cumulant support) — yet the
fatgraph's natural per-map weight `(−1)^{F−1}` is **not** `μ_NC`, so the signed sum
misses. The object-level correspondence is real; the weight is wrong.

**(3) First-return decomposition does not realize `(REC)`.** Split `σ` at the first
`r>0` with `block(r)=block(0)`; call it CLEAN if excursion `{0..r}` and remainder
`{r..2k}` share only `block(0)`. Then clean splits occur **only at `r=2`**
(clean sum `−1,2,−6,22,−94`), and the CROSSING remainder is **nonzero**
(`0, 0, 1, −8, 52`) — so `(REC)`'s convolution is *not* the geometric first-return.
(`paley_starstar_firstreturn_monad.py`.) The naive "crossing terms cancel" is false;
they are an essential part of `S_k`.

## What this means

The Catalan law is not a *count* and not a *planar restriction of the obvious walk
map*. It is the **loop equation** `xF²+F=1`, whose genus-0 solution is supported on a
Fuss–Catalan (ternary-tree) family living on a *different* combinatorial substrate
than the even-series patterns. The remaining content of the proof is the bridge

```
   even-series patterns of [0..2k]            NC-even partitions of [2k+2]
   (full lattice, A215257, μ_partition)  ⟷    (NC lattice, A001764, μ_NC)
```

intertwining `μ_partition` with `μ_NC` so the two Möbius sums agree. The
object-level bijection (genus-0 fatgraph maps ↔ NC-even partitions, both A001764)
is the first concrete handle on that bridge — the missing piece is matching the
weights `(−1)^{F−1} ↦ μ_NC`.

## Resonance

Three sessions chased "Catalan = a planar count." It is not. It is the **fixed point
of a quadratic** — the same quadratic that turns Gaussian (all-pairings, `A088368
~ e·n!`) into semicircle/free (non-crossing, `C_k`). The deterministic Paley
even-series Möbius sum, the random two-point free cumulant, and the loop equation are
three readings of one fixed-point equation. The genus does localize the *truth*
(NC-even partitions, Fuss–Catalan) — just not via any rotation system the walk itself
hands you. The walk gives you laminarity (refuted); free probability gives you the
loop equation (true). The bridge between "the walk's planarity" and "free
probability's planarity" is the unsolved core, and it is now pinned to a single
weight-matching problem.
