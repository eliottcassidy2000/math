---
id: THM-738
title: The cluster-clock lemma (additive pack-clock) — every far-cluster family B ∪ (W+J) has L ≥ Area(B,J) − C(B,J)/W with Area and C exactly computable; uniform tail closure of the additively-coherent slice of the j ≥ 7 seam (beyond kps THM-735's Bonferroni boundary j ≤ 6), with |J| = 1 degenerating to the far-element peel
status: PROVED (proof below) + VERIFIED-EXACT (Area exact by adaptive piecewise-linear Fraction integration with midpoint-linearity self-check; bound valid at all tested W; klein-S289's {1,90..101} sits at 99.9% of its Area limit) — below-threshold sweep in the companion .out
source: opus-2026-07-13-S273 (owner prompt: work the gcd-incoherent sector with the perspective frame)
depends_on: []   # self-contained counting/measure argument; no LRC citation (the base good set G_B is an exact finite computation)
related:
  - THM-735 (kps-S128c3, Bonferroni multi-peel)  # closes j ≤ 6 far slots via union bound (base 1 − j/7); THIS lemma takes the ADDITIVELY-COHERENT part of j ≥ 7, where the cluster's arcs overlap instead of union-bounding
  - THM-737 (opus-S272, pack-clock sampling)      # the MULTIPLICATIVE coherence twin (s = ct); this is the ADDITIVE twin (s = {Wt}); one move: ride the fastest coherent clock, freeze the slow base, count/integrate the detuned remainder
  - THM-731/732 (klein S287/288)                  # |J| = 1 degenerates to the far-element peel main term Area = (6/7)|G_B|
  - THM-668 (monad-S3)                            # the multiplicative point-witness dispatch; same family tree
  - HYP-6248 / opus-S270 frozen-fan lemma          # the base-freezing mechanism, here quantified as the b/W fattening
  - klein-S289 (HYP-6505 second use)               # {1,90..101}, {1,2,3,50..59}: the non-isolated counterexample shapes this lemma closes uniformly
---

# THM-738 — the cluster-clock lemma (additive pack-clock)

## Statement

Let `B` be a finite set of nonzero integers (the **base**), `J ⊆ {0, 1, …, J_m}` a finite offset
pattern with `max(J) = J_m` (the **cluster**), and for `W ≥ max(B) + 2` consider the 13-speed
family `V_W = B ∪ {W + j : j ∈ J}` (`|B| + |J| = 13`). Define

- `G_B = {u ∈ [0,1) : ‖bu‖ ≥ 1/14 ∀b ∈ B}` (base good set — exact intervals),
- `A_J(u) = 1 − |⋃_{j∈J} (−ju − 1/14, −ju + 1/14) mod 1|` (cluster admissibility — piecewise
  linear in `u`),
- `Area(B,J) = ∫_{G_B} A_J(u) du` (exactly computable rational),
- `C(B,J) = V(A_J) + 2·#comp(G_B) + 2·Σ_{b∈B} b + |J|·J_m` (`V` = total variation on `[0,1]`).

Then the `1/14`-safe measure satisfies

> **`L(V_W) ≥ Area(B,J) − C(B,J)/W`  for every `W`.**

In particular, if `Area(B,J) > 0` then every `W > W₀ = ⌈C/Area⌉` is closed (LRC(14) holds for
`V_W` with a positive-measure witness set), and `W ≤ W₀` is a finite list of exactly-decidable
bodies (rational interval arithmetic).

## Proof

Ride the cluster's clock: `s = {Wt}`, branches `t = (m+s)/W`, `m = 0..W−1`, `s ∈ [0,1)`,
`u_m = m/W`.

1. **Frozen fan (base).** For `b ∈ B`: `b·t ∈ [b u_m, b u_m + b/W]`, so if
   `‖b u_m‖ ≥ 1/14 + b/W` then branch `m` is `b`-safe for *every* `s`. The set of good branches
   `{m : u_m ∈ G_B^{fat}}` with `G_B^{fat} = {u : ‖bu‖ ≥ 1/14 + b/W ∀b}` satisfies
   `|G_B ∖ G_B^{fat}| ≤ Σ_b 2b/W` (each `b`-danger arc fattens by `b/W` per side; measure of
   `‖bu‖ < x` is `2x`).
2. **Cluster = fattened AP.** For `j ∈ J`: `(W+j)t ≡ s + jt (mod 1)` **exactly**, and
   `jt ∈ [j u_m, j u_m + j/W]`. So `s` avoiding
   `(−j u_m − 1/14 − J_m/W, −j u_m + 1/14)` for all `j` suffices; the branch-`m` safe-`s`
   measure is `≥ A_J(u_m) − |J|·J_m/W`.
3. **Riemann sum.** `L ≥ (1/W)Σ_{u_m ∈ G_B^{fat}} [A_J(u_m) − |J|J_m/W]`. Since `A_J` is
   piecewise linear (BV) and `G_B^{fat}` a finite interval union, the grid sum deviates from
   `∫_{G_B^{fat}} A_J` by at most `(V(A_J) + 2·#comp)/W`; combining with step 1's measure loss
   gives the stated constant. ∎

## The two computed shapes (klein-S289's counterexample family, exact)

| shape | `B` | `J` | `Area` (exact) | `V(A_J)` | `C` | `W₀` |
|---|---|---|---|---|---|---|
| far cluster `{1}∪{W..W+11}` | {1} | {0..11} | **`1039/13860` = 0.074964** | 10 | 146 | **1948** |
| two-scale `{1,2,3}∪{W..W+9}` | {1,2,3} | {0..9} | **`2767/61740` = 0.044817** | 346/35 | 4196/35 | **2676** |

Verification: the bound holds at every tested `W`; the true `L` converges far faster than
`C/W` — klein's actual body `{1,90..101}` (`W = 90`) has `L = 0.074877` = **99.9% of the Area
limit** (observed `|L − Area|·W ≈ 0.05`, i.e. the rigorous `C = 146` is ~3000× conservative;
the dominant crude term is the within-branch drift fattening `|J|·J_m`, which a second-order
equidistribution treatment would shrink — not needed for finiteness of `W₀`). At `W = 2` the
first shape degenerates to the AP `{1..13}` (`L = 0`, the extremal — the family boundary).

## Positivity mechanism — the collapse does positive work

`Area > 0` is carried by the **resonant windows**: near `u = p/q` (`q ≤ 6`, `p/q ∈ G_B`) the
cluster AP `{ju mod 1 : j ∈ J}` collapses onto `q` points, so the fattened union covers
`≤ q/7 < 1` and `A_J ≥ 1 − q/7 > 0` on a window. This is the S270 frame-collapse phenomenon
(the double-counted pair sector of the perspective decomposition degenerating) — here doing
*positive* work: the cluster's internal differences `δ = j′−j` (small, additively coherent)
are exactly what stack the arcs. The breakpoints of `A_J` are `(k ± 1/7)/δ` and `k/δ` — the
additive tile structure of the cluster's pair sector.

## Scope and position

1. **The seam at j = 7 (kps-S128 cont.3):** THM-735's Bonferroni closes `j ≤ 6` far slots
   (base `1 − j/7`); for `j ≥ 7` the union bound is void. This lemma replaces the union bound
   by the *overlap structure* whenever the far elements are **additively coherent** (a cluster
   with bounded pattern `J`), closing `|J|` up to 12. Remaining in the `j ≥ 7` seam: families
   whose ≥ 7 far elements are neither multiplicatively (THM-668/737) nor additively (this)
   coherent — genuinely spread two-scale families, the density route's home turf, plus
   multi-cluster shapes (extension: hierarchical clocks — ride the fastest cluster's clock,
   the lower cluster joins the frozen base with its own `A`-factor; not written here).
2. **Unification (the perspective arc S270→S273):** multiplicative pack ⟹ THM-737 (`s = ct`);
   additive cluster ⟹ this (`s = {Wt}`); single far element ⟹ `|J| = 1` degeneration
   (`A ≡ 6/7`, `Area = (6/7)|G_B|` — the far-element peel main term, THM-731's regime). One
   move: **ride the fastest coherent clock; the slow base freezes (frozen fan); the detuned
   remainder is counted (THM-737) or integrated (this) exactly.**
3. **Lean shape:** exact interval arithmetic + a piecewise-linear integral (finite rational
   breakpoint partition with midpoint-linearity certificates) + the two fattening bookkeeping
   steps. No analysis at runtime; same infrastructure class as THM-737/LRCDetunedDispatch.

## Verification & files

`04-computation/lrc14_cluster_clock_lemma_thm738_opus_S273.py` (+ `.out`): exact Area/V/C/W₀
for both shapes (adaptive exact integration, linearity self-checked at every piece), bound
validity at `W ∈ {50,…,1600}`. Companion sweep `lrc14_cluster_clock_sweep_opus_S273.out`:
float `L > 0` for all `W ∈ [max(B)+2, W₀]` with exact re-verification at the family minimum.
