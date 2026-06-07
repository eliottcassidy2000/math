# Genus-blindness: the Catalan cancellation in `(★★)` lives ACROSS cycle rank, not inside genus

*monad-explorer-2026-06-07 (deep-research / analytic lane, 7th session)*
*Builds on THM-438 ADDENDUM-3/4, the cycle-rank triangle, MISTAKE-060/061.
Scripts: `04-computation/paley_starstar_{rootmap,cyclerank_triangle,topological_factor}_monad.py`.*

## The one-line theorem that closes the whole genus program

Recall `(★★)`: `S_k := Σ_{σ : even-series pattern of [0..2k]} μ(0̂,σ) = (−1)^k C_k`,
with `μ(0̂,σ) = ∏_v (−1)^{|B_v|−1}(|B_v|−1)!`.

**Sign identity (PROVED).** The walk has `2k+1` positions `{0,…,2k}`, so
`Σ_v(|B_v|−1) = (2k+1) − V`, and the cycle rank of `G_σ` is `m = E−V+1 = 2k−V+1`.
These are EQUAL: `(2k+1)−V = 2k−V+1 = m`. Hence

```
        μ(0̂,σ) = (−1)^{m(σ)} · ∏_v (|B_v|−1)! ,        S_k = Σ_σ (−1)^{m(σ)} ∏_v (|B_v|−1)!.
```

The `∏_v(|B_v|−1)!` is positive: it counts the cyclic orders of the `|B_v|` visits to
vertex `v`. Choosing one cyclic order at every vertex is a **ribbon (rotation) system**
`R` on `G_σ`, i.e. a combinatorial map. So `S_k = Σ_{(σ,R)} (−1)^{m(σ)}`.

**Genus-blindness (PROVED).** For any connected map with `V` vertices, `E` edges,
`F` faces, genus `g`: `V−E+F = 2−2g`, so `F = (E−V+1)+1−2g = m+1−2g`, giving

```
        (−1)^{F(R)−1} = (−1)^{m−2g} = (−1)^{m(σ)}      for EVERY rotation system R.
```

The per-ribbon-map sign **does not depend on the genus**. Verified exhaustively over
all `403` `(σ,R)` pairs with `k≤3` (`paley_starstar_rootmap_monad.py`).

## Why this kills all three of ADDENDUM-4's routes at once

ADDENDUM-3 and -4 chased a "genus-0 localization": find a planar/ribbon structure on
each `G_σ` so that `Σ_{genus 0} (weight) = (−1)^k C_k`. Three constructions were tried
and all failed (laminar = canonical walk-ribbon; fatgraph rotation-sum; first-return).
Genus-blindness explains **why no such construction can exist**:

- `S_k = Σ_{(σ,R)} (−1)^{F(R)−1}` is an **all-genus** map sum, and the summand sign
  `(−1)^{F−1}=(−1)^{m(σ)}` is constant across the genus of a fixed `G_σ`. Restricting to
  genus 0 keeps `(−1)^{m(σ)}·#{genus-0 R}` — a *positive multiple* of the sign — which is
  exactly the (non-Catalan) fatgraph genus-0 number. The cancellation that produces
  `C_k` is **between different graphs `G_σ` of different cycle rank `m`**, NOT between
  rotation systems of one graph. A planarity filter can never see it.

This is the corrected mental model. The right home for the cancellation is the
**cycle-rank triangle** (ADDENDUM-3):
```
   t(k,m) = Σ_{even-series σ, rank m} ∏_v(|B_v|−1)!  (all positive)
   k=1: 1 | k=2: 1 3 | k=3: 1 9 13 | k=4: 1 18 72 69 | k=5: 1 30 230 580 421
   Σ_m (−1)^m t(k,m) = (−1)^k C_k.                          ← THIS is (★★)
```
`m=1` is the single `2k`-cycle (`t(k,1)=1`); `m=k` is the bigon-tree all-pairings
overcount (`t(k,k)=A088368 ~ e·k!`). The Catalan number is the signed collapse along
the `m`-axis.

## A refuted sub-conjecture (recorded so it is not retried)

Tempting idea: group even-series `σ` by `e =` number of series-classes (= edges `E_H`
of the topological reduction) and hope `G(k,e) := Σ_{e-classes}(−1)^m∏(b−1)!` factors as
`g_e·C(k−1,e−1)` with `g_e` `k`-independent. **FALSE** (`paley_starstar_topological_factor_monad.py`):
`G(k,2) = 3, 8, 15, 24 = k²−1` (not `3(k−1)`). Different topological types sharing the
same edge-count `e` have different `m` and weight, so `e` alone does not factor.

The binomial-inverse constants `g_e = −1, 3, −10, 36, −137, 543 = (−1)^e A002212(e)`
ARE real but **tautological**: `F(x)=Σ(−1)^kC_k x^k = G(x/(1−x))` with
`G(y)=F(y/(1+y))=Σ(−1)^e A002212(e) y^e` (sympy-verified). This is just the classical
fact that the loop equation `xF²+F=1` transported by `x↔y/(1+y)` is A002212's quadratic
`yG²+(1+y)(G−1)=0`. No new reduction.

## The sharpened handoff

The cancellation is graded by **cycle rank `m`**, not genus. So the loop equation must
come from a recursion/involution on the `m`-grading. Concretely:

```
   U(x,y) := Σ_{k,m} t(k,m) x^k y^m   (catalytic variable y marks cycle rank).
   (★★)  ⟺  U(x,−1) = F(x)  satisfies  x·U(x,−1)² + U(x,−1) − 1 = 0.
```

Find the algebraic/Tutte equation for `U(x,y)` (the marked Eulerian trail is the root),
then specialize `y=−1`. The catalytic variable is the right `y`; the three planarity
routes are closed off; the target is now a single bivariate functional equation whose
`y=−1` slice is the loop equation. A sign-reversing involution on even-series patterns
that changes `m` by `±1` (pairing rank-`m` with rank-`(m±1)` configurations) would do it
directly.

## Resonance

The first six sessions read the Catalan law as a *planar count* — a genus-0 object. It
is not. `(−1)^{F−1}=(−1)^{m}` is the precise reason: the surface a ribbon graph lives on
is invisible to the sign. What carries the cancellation is the **first Betti number of
the underlying graph** — a homotopy invariant, not a genus. Free probability's
"non-crossing dominance" here is realized not by suppressing crossings on one surface but
by an alternating sum over the homological complexity `m` of the whole family of walk
graphs. The planar answer (`C_k`) is genuine; its mechanism is homological, not
topological-surface. That is the lesson: when a signed sum refuses to localize on a
surface, look at the Betti grading underneath it.
