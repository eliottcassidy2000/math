---
source: kind-pasteur-2026-07-07-S61
status: PROVED (the arc bound + the per-V witness criterion; witnesses machine-verified) /
  TABULATED (explicit V0; the 2-torus GAP ledger) / delineation of what remains.
tags:
  - lonely-runner
  - LRC14
  - route-1
  - part-a
  - spread-residual
  - quantitative-factoring
---

# Part A factored, and the spread residual stratified

**kind-pasteur-2026-07-07-S61 (HYP-4857).** Owner: work the spread residual and the
quantitative factoring. Both turned out to hinge on the same move from S59/S60: **substitute
the ledger subset for the true good set.**

## 1. The o(Vmax) obstruction dissolves (Part A factored)

monad-HYP-4787's audit left Part A one missing written piece: the finite-`Vmax` correction
is `O(#arcs/Vmax)`, and the true good set's arc count grows (`~S^0.45`) for spread-`~Vmax`
shapes. The fix is that the witness argument never needed the true good set — **any
positive-measure few-arc subset works**, and the S60 intersection-ledger subset

> `R(P, D+1; Vmax) = G_P^{+3/(7Vmax)} ∩ {roof_{D+1} ≥ 1/7 + 12D/(7Vmax)}`

(widened cuts, drift-shifted threshold — both slacks conservative for the teeth/slow drift
across one ruler period) has a **proven absolute arc bound**:

> `#arcs(R) ≤ 13 + Σ_{p∈P}(p+1) ≤ 83`, independent of `E`, its spread, and `Vmax` —
> the roof superlevel has ≤ 13 components because adjacent Farey-6 rationals are separated
> by their mediant (denominator ≥ 7, where the roof dips to ≤ 1/7 < θ), and widened-`G_P`
> has ≤ `Σ(p+1)+1` components; an intersection of interval-unions adds components minus one.

With ruler centers spaced exactly `1/Vmax`: `#good centers ≥ Vmax·meas(R) − #arcs(R)`, so

> **`Vmax · meas(R(P,D+1;Vmax)) > A_P  ⟹  M(S) ≥ 1/14`, witnessed** (the fast phase
> sweeps the full window each period; put it at the widest center-config gap's center; the
> `0∈E` tooth keeps it inside `(1/14, 13/14)` automatically).

Since `θ` and the widening shrink as `Vmax` grows, `Vmax·meas(R)` is **monotone** — so with
the proven `A_P` the condition persists for all larger `Vmax`, giving the clean two-regime
factoring monad's target-2 asked for:

> **For every covering shape `(P, E)` with cluster diameter `D` inside its leg's S60 bite:
> `Vmax ≥ V₀abs(P,D)  ⟹  M(S) ≥ 1/14`.** Exact tables: sharp per-V `V₀ = 140..1064`;
> rigorous ∀-V `V₀abs ≤ 1106` at every bite edge (k=13: 1064; k=12: 784; k=11: 1047;
> k=10: 990; k=9: 695; k=8: 1106). **The finite check lives below height ~1106** — three
> orders of magnitude, not an asymptotic gap. Direct verification: 8/8 concrete covering
> sets at and above `V₀` produced exact witnesses `τ*` with min clearance `0.082–0.486 ≥ 1/14`.

Honest residue of the factoring: the finite check (`Vmax < V₀abs`, diam within bites) is
*specified, not executed* — its raw family count is large; per-`(P,D)` sharp `V₀` (as low
as 140) and the covering-sieve reductions shrink it; it is Tao-2018-flavored territory with
a dramatically smaller constant.

## 2. The spread residual stratifies (the GAP-grid superset ledger)

For the residual beyond the bites, the same subset lemma applies to **rank-2 supersets**
(mac-mini-S41's cover-or-decorrelate, handoff (b)): if `E ⊆ a + {i·d₁ + j·d₂}` (an
`(n₁,n₂)`-grid), then `μ_{1/7}(E) ≥ μ_{1/7}(grid)`, and the grid's tail is governed by the
`(d₁,d₂)`-geodesic average of a **2-torus function** — so it converges (as `d₁d₂` grows,
coprime) to a **step-independent constant** `μ₂(n₁,n₂)`:

- **Step-independence measured:** spreads `0.015–0.035` across steps from `(11,2)` to
  `(211,30)`; error `≈ C/(d₁d₂)`, `C ≤ ~5` empirically — the same Fourier mechanism as the
  S59 pair bound (modes `m₁d₁ + m₂d₂ = 0`), provable-shape.
- **The `μ₂` table clears the bar everywhere sampled:** `μ₂(n₁,n₂) ≥ m_P` for every grid
  with `N = n₁n₂ ≤ ~78` points (`n₂ ≤ 6`), with per-`n₂` crossings all at `N ≈ 78–80` —
  matching the rank-1 crossing (`AP_76`) at equal point-count: the `~4.3/N` Farey-window
  law is rank-independent.
- **Coverability stratifies the residual:** the record family is rank-1-covered (diam 20);
  the parity diam-80 family is `(11,6)`-grid-covered (steps `(8,3)`, `N=66`); the two-block
  `{0..5}∪{70..76}` is `(7,2)`-grid-covered (steps `(1,70)`, `μ₂ = 0.62` — my search's
  `n₂ ≤ 6`-cap missed the swapped orientation, found by hand). The families with **no**
  small grid cover — the deep interlacing `10·{0..10}∪{49,51}` and mac-mini's 3-adic
  cascade — are exactly the **decorrelation-lane** families, and their directly-measured
  `μ` is high (0.5–0.6): sparse in every GAP ⟹ weak additive structure ⟹ (S59 deficit
  frame) near-iid tails.

So the spread residual now reads:

```
k=13 residual (diam >= 76):
  rank-2 GAP-coverable (N <= ~78, steps equidistributed)  ->  mu2 ledger  [tabulated, >= m_P]
  not coverable                                            ->  sparse/decorrelation lane
                                                               (PZ/CE floors 2-4x bar; the
                                                               honest open analytic core)
```

The open mathematical items, sharply: (i) exact-rational `μ₂` (2-dim three-distance /
Farey-cell decomposition — the formalizable version of the table); (ii) the geodesic error
bound for the sharp indicator (smoothing or bounded-variation Koksma); (iii) the
non-coverable lane's floor (mac-mini PZ-on-U / monad CE — unchanged owners).

## The through-line

Three sessions, one lemma: *the point set of a subset family contains the point set of its
superset... backwards* — nested point sets, pointwise-dominated max-gaps, and every
functional, every threshold, every slow-part intersection, every rank of superset inherits
exact constants from it. S59 bounded the tail by diameter; S60 pushed it through `G_P`;
S61 pushed it through the ruler-period drift (Part A) and up a rank (GAP covers). Each time
the "hard analytic" obstruction (extremality, union-bound loss, arc growth, spread) was
bypassed rather than solved, and what remains is genuinely smaller: a bounded finite check
and one sparse-family analytic lane.

## Ledger

- Files: `lrc_quantitative_partA_kps_S61.py`, `lrc_gap_grid_ledger_kps_S61.py` (+outs).
- Builds on: kps-S59/S60 (subset lemma, ledgers), opus-S134+S135 (roof, now PROVED
  citation-free — THM-637), THM-527-A (ruler mechanics), boxeph-S1 (1/7 sharpness),
  monad-HYP-4787 (target 2), mac-mini-S41 (cover-or-decorrelate, handoff b).
- Does NOT prove LRC(14). Factors Part A over the ledger-covered domain with explicit
  constants; stratifies the spread residual; the sparse lane and the finite check remain.
