---
source: opus-2026-07-09-S169
status: made the arc-count good-period route partly a-priori/Lean. (1) LEAN: formalized the PIGEONHOLE
  HEART -- good_period_of_arccount (#arcs < rho*.V => a grid point j/V is Good = a good period),
  exists_gridpoint_Ico (a long interval catches a grid point), exists_long_Ico (N intervals summing
  to >N/V => one exceeds 1/V), + the contrapositive. Kernel-pure [propext,Classical.choice,Quot.sound],
  built (8476 jobs). bounded-arc-count enters as the clean hypothesis. (2) SHARPENED TARGET: the correct
  closing inequality is the rho*-pigeonhole #arcs<rho*(spread+1) (UNIVERSAL, huge margin), NOT mac-mini's
  c<D3 (which FAILS at small spread: k=13 spread=40 c=0.70>D3=0.59). The min-arc-length route is DEAD
  (arcs shrink ~0.008/spread). The trivial a-priori #arcs bound O(k^2 spread) is 200-1300x too loose --
  the small-constant #arcs<=c.spread is the ONE open analytic item. (3) LEMNISCATE (owner): the node =
  the collision (e_i-e_j)x in Z = the maxgap breakpoint = the measure-zero exact resonance (S168 clock)
  = the arc boundary -- one object, four views; and the elliptic arc-length reparam DESINGULARIZES the
  near-tangential 1/7-crossings (the short arcs that make the bound hard), a concrete proof direction.
tags:
  - lrc14
  - good-period
  - arc-count
  - lean
  - lemniscate
  - desingularization
---

# The lemniscate node is the collision; the arc-count pigeonhole goes to Lean

**opus-2026-07-09-S169.** Owner: make bounded-arc-count + `c<D3` fully a-priori/Lean; consider the
lemniscate `(x²+y²)² = x²-y²` as a strange source of inspiration.  Three results.

## (1) The pigeonhole heart is now Lean (kernel-pure)

The dissociated good-period branch closes by EXISTENCE of a good period `j` (any `j`, not small-`j`)
from the pigeonhole `#arcs(Good ∩ [0,1)) < rho*·V`.  I formalized its mathematical heart in
`TournamentH7.LRCArcCount` (built, `[propext, Classical.choice, Quot.sound]`):

- `exists_gridpoint_Ico`: if `Ico a b` is longer than `1/V`, then `⌈a·V⌉/V ∈ [a,b)` -- a long
  interval catches a grid point.  (`Int.le_ceil` + `Int.ceil_lt_add_one`, `1 < V(b-a)`.)
- `exists_long_Ico`: of `N` intervals with `∑ lengths > N/V`, one exceeds `1/V` (pigeonhole).
- `good_period_of_arccount`: if `⋃ᵢ Ico(aᵢ)(bᵢ) ⊆ Good θ E` and `∑(bᵢ-aᵢ) > N/V`, then
  `∃ j : ℤ, (j/V) ∈ Good θ E` -- **a good period exists.**
- `arccount_le_of_no_good_period` (contrapositive): no good period ⟹ `∑ lengths ≤ N/V` -- the exact
  obstruction, realized only by the tight complete-residue AP at `V` prime (opus-S164), the extremal.

The two INPUTS stay hypotheses (as they should): `rho* ≥ D3 ≥ bar > 0` is the density floor (THM-661 /
klein's finite `D3_c` table, done), and `#arcs ≤ c·spread` with `c < rho*` is bounded-arc-count.  The
REDUCTION -- "those two ⟹ good period" -- is now unconditional Lean.  What was an informal
measure/pigeonhole hand-wave (mac-mini-S58) is a checked theorem.

## (2) The a-priori target, sharpened and de-risked

Testing the candidate closing inequalities on genuinely dissociated clusters (longest-AP `≤ k-6`, the
only branch arc-count serves), k=11,13, spread `40 → 1280`:

| inequality | holds? | note |
|---|---|---|
| **(i) `#arcs < rho*·(spread+1)`** (the pigeonhole) | **YES, universally** | margin huge: `c ~ .04-.7` vs `rho* ~ .94-.998` |
| (ii) `c = #arcs/spread < D3(E)` (mac-mini-S62) | **NO at small spread** | k=13 spread 40: `c=.70 > D3=.59` |
| (iii) `min-arc-length ≥ const/spread` | **NO** | min-arc·spread `.29 → .008`, arcs shrink |

So the CORRECT closing inequality is the `rho*`-pigeonhole (i), not `c < D3`.  `D3 ~ .59-.85` is the
a-priori *floor* on `rho*`, but the *actual* `rho* ~ .94-.998` is far larger, and (ii) trades that
margin away -- failing below `spread ≈ 120` at k=13 (exactly mac-mini's `spread < 200` finite-check
carve-out).  Using `rho*` directly (or the honest floor `rho* ≥ D3`) needs no carve-out.  The min-arc
route (iii) is DEAD: arcs get arbitrarily short (`~0.008/spread`), so `#arcs ≤ rho*/min-len` is useless
-- the short arcs are the near-tangential `1/7`-crossings (see §3).

**The one open item, made crisp.** All that remains a-priori is `#arcs ≤ c·spread` with `c < rho*`.
The trivial bound `#arcs ≤ 2∑|eᵢ-eⱼ| = O(k² spread)` is **200-1300× too loose** (measured), giving
`c = O(k²) ≈ 121 ≫ rho*`.  The truth is `#arcs ~ spread^{0.92}` (opus-S168), so `c → 0`.  Closing the
`~10³×` gap -- controlling the near-tangential crossings -- is the sole remaining analytic obstruction
of the covering capstone.

## (3) The lemniscate: the node is the collision (one object, four views)

The owner's lemniscate `(x²+y²)²=x²-y²` (`r²=cos 2θ`) is a figure-eight whose SELF-INTERSECTION at the
origin makes two circle-far runners momentarily coincide -- the "crossing shortcut."  The owner is
right that Euclidean lemniscate distance BREAKS LRC (non-isometric) and one must pull back the circle
metric.  But the node itself is the gift.  **The lemniscate node is, simultaneously:**

- the **runner collision** `(eᵢ - eⱼ)x ∈ ℤ` (two orbit points coincide);
- the **breakpoint** of the piecewise-linear `maxgap(x)` (where the cyclic order, hence the gap slopes,
  changes) -- the only places `#arcs` boundaries can sit;
- the **measure-zero exact resonance** -- the broken clock being *exactly* right (opus-S168), negligible
  vs the positive-measure near-resonance lingering;
- the **arc boundary** candidate of `Good ∩ [0,1)`.

One object, four views -- collision = node = exact resonance = arc edge.  This unifies the S168 clock
picture with the arc-count combinatorics: `#arcs` is controlled by the collisions, and the collisions
are the lemniscate's nodes.

**The genuine proof lead: desingularization.** The short arcs that kill the min-arc route (§2, iii) are
near-*tangential* `1/7`-crossings of `maxgap` -- passages that graze the level `1/7` and retreat, the
analog of the lemniscate's near-origin passage where arc-length (the elliptic integral
`∫dθ/√(cos 2θ)`) has its singular, curvature-blowup structure.  The owner's "use the arc-length metric,
runners become variable-speed (elliptic integrals)" is exactly a **reparametrization that resolves the
tangency into a transversal crossing.**  `#arcs` is reparametrization-INVARIANT (a topological count),
so this changes no answer -- but a coordinate in which every `1/7`-crossing is transversal (slope
bounded below) makes `#arcs` *countable by* the crossing number, which is what an a-priori bound needs.
Desingularizing the near-tangencies via a lemniscate-type arc-length coordinate is the concrete (hard)
route to the open `#arcs ≤ c·spread`.  The variable-speed non-uniformity the owner flags IS the
variable arc-length of the good set -- the short-arc problem is intrinsic, and the elliptic coordinate
is where it might become uniform.

## Ledger

- LEAN: `TournamentH7.LRCArcCount` -- `good_period_of_arccount`, `exists_gridpoint_Ico`,
  `exists_long_Ico`, `arccount_le_of_no_good_period`; kernel-pure, built (8476 jobs).  The arc-count
  REDUCTION is now unconditional; the two inputs (`rho*≥D3`, bounded-arc-count) are the hypotheses.
- TARGET: the closing inequality is the `rho*`-pigeonhole `#arcs<rho*(spread+1)` (universal), NOT
  `c<D3` (fails small-spread k=13).  min-arc route DEAD.  Open item: `#arcs≤c·spread`, `c<rho*`;
  trivial bound `O(k² spread)` is `10³×` loose; truth `~spread^{0.92}`.
- LEMNISCATE: node = collision = maxgap breakpoint = exact resonance (S168) = arc edge; elliptic
  arc-length DESINGULARIZES the near-tangential `1/7`-crossings (the short arcs) -- the concrete lead
  for the open bound.  Do NOT use Euclidean lemniscate distance (breaks LRC).
- CONVERGENCE (kps-2026-07-09-S94 / LEM-013, independent, same day): kps hit the SAME `c≥D3` sliver
  (k=13 small-spread: `c=0.675 ≥ D3=0.659`) and resolved it exactly as §2 predicts -- **skip the `c<D3`
  proxy, test good-period EXISTENCE directly**: exhaustive `spread ≤ 22` (621,455 primitive dissociated
  clusters, ZERO failures, min margin `7·maxgap/Vmax = 1.1053 = 21/19`), and an adversarial hill-climb
  shows the margin INCREASES with spread (`1.36` at `s=27`, `1.72-2.31` at `s∈[50,200]`).  So the
  `c≥D3` sliver is a certificate artifact, not a covering gap.  This is the DIRECT-EXISTENCE twin of my
  `rho*`-pigeonhole (i): my `good_period_of_arccount` is the Lean-checked REDUCTION behind LEM-013's
  existence, and my `#arcs ~ spread^{0.92}` (S168) is the mechanism behind kps's margin-grows-with-spread
  (fewer arcs per unit `V` ⟹ bigger gaps ⟹ bigger margin).  Together: [Lean reduction: arcs+measure ⟹
  existence] + [kps: inputs verified exhaustive `s≤22` + adversarial-robust] = the dissociated branch,
  modulo a clean a-priori `μ>1` on the band `s∈[22,~100]` (kps: "not a genuine risk").
- Files: `lrc14_arccount_apriori_target_opus_S169` (+out), `LRCArcCount.lean`.  -> mac-mini-S58/S61/S62
  (arc-count, `c<D3`), kps-S94 (LEM-013 direct existence), opus-S167/S168 (interlock, clock),
  klein-S187 (`D3_c` table), THM-661.
