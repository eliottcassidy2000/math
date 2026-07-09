# The continuum bridge: where the grid AND the drift desingularize together

*kind-pasteur-2026-07-09-S112. Owner: formalize the continuum reformulation bridge on the smooth
surrogate. This note records why the continuum is the natural home of the bridge — it is the single
limit in which BOTH finite-scale obstructions that stalled the fleet vanish at once — and what the four
new Lean lemmas (`LRCSmoothBridge`) actually pin down.*

---

## Two obstructions, one limit

The reformulation `ρ* > 0 ⟹ ∃ lonely τ` had been blocked by two *different-looking* obstructions, each
found by mac-mini on the finite-`Vmax` ruler:

1. **Grid-invisible pinches** (MISTAKE-130, kps-S107). The sharp good set `G* = {x : maxgap{frac(e_i
   x)} > 1/7}` is bounded by rational pinches `x = m/d` (`d` a cluster difference), a measure-zero
   nodal set (the lemniscate node, `LRCPinch`). *Sampling* `G*` on any ruler grid `x_j = (j+½)/Vmax`
   over- or under-counts arcs — the grid is blind to a measure-zero boundary.
2. **The drift** (mac-mini's local hembed refutation). Even given a good period, the single-`φ`
   embedding `nearInt((Vmax−e_i)τ) = nearInt(φ − e_i·τ)` carries a drift `e_i·φ/Vmax`
   (`LRCSlowFast.drift_eq`) that is `O(1)` in the good-period window, so the good period does **not**
   certify loneliness at its own `(j,φ)` — the observer wobbles off the gap.

These look like separate walls. They are the **same wall seen from two sides**: both are artifacts of
working at a *finite scale* `Vmax`. And they have a **common desingularization** — the `Vmax → ∞`
continuum:

- The grid obstruction dissolves because you stop sampling. A positive **measure** `∫₀¹ W > 0` hands
  you a positive **point** `W(x) > 0` *directly* (`exists_pos_of_integral_pos`) — the pinches are
  measure-zero and simply do not move the integral. No grid, nothing to be blind to.
- The drift obstruction dissolves because `e_i·φ/Vmax → 0`. In the limit the embedding is
  `nearInt(φ − e_i·τ)` with **no** drift term, so the gap-midpoint observer is exact
  (`observer_of_confined`).

The continuum is not one convenience among several; it is the *unique* limit in which both finite-scale
artifacts vanish simultaneously. That is why the bridge belongs there.

## What the four lemmas pin (`LRCSmoothBridge`, sorry-free)

The bridge factors into a clean chain, each link a named lemma:

| step | lemma | content |
|---|---|---|
| density floor → good point | `exists_pos_of_integral_pos` | `W ≥ 0`, `∫₀¹W > 0` ⟹ `∃ x, W(x) > 0`. The **desingularization**: measure → point, no grid. |
| good point → wide gap | `exists_gap_gt_of_smoothW_pos` | `W = Σ(gap_i − 1/7)_+ > 0` ⟹ `∃ i, gap_i > 1/7`. A genuine clearing gap at the continuum. |
| wide gap → observer | `observer_of_confined` / `observer_at_threshold` | cluster in an arc `≤ 6/7` (gap `≥ 1/7`) ⟹ gap-midpoint `φ` has `nearInt(φ − q_i) ≥ 1/14` ∀ i. **Drift-free**, exact LRC margin. |
| observer → reach | `mreach_ge_of_smooth_surrogate` / `_density_floor` | composes with `Mreach_ge_of_lonely_instant` (kps-S99b) to `Mreach ≥ 1/14`. |

The **middle two** are the new mathematical content — they are the smooth surrogate *earning its keep*.
`exists_pos_of_integral_pos` is the S107 prescription made rigorous: the reason the smooth `W` (Fourier
`1/m²`, `C⁰` through pinches) works where the sharp indicator (`1/m`, jump at pinches) fails is that its
positive integral is *robust* to the measure-zero nodal set. `observer_of_confined` is the reconstruction
`good period ⟹ lonely phase` **exactly**, in the drift-free limit — the clean core that mac-mini's finite
refutation could not touch, because the thing mac-mini refuted (finite drift) is gone.

## The residue is Diophantine, not analytic

After these four lemmas the only thing still *hypothesised* (the `hrefl` argument of the bridge) is the
**Kronecker realization**: that the cluster `{e_i·τ mod 1}` can be simultaneously confined to a short arc
for a genuine `τ`, and that the observer phase `φ` lifts back to a real time. This is a statement of
**simultaneous Diophantine approximation** (Weyl/Kronecker on the co-offsets `e_i = Vmax − v_i`), not an
analytic or grid statement. The bridge has thereby sorted its own difficulty: everything *analytic*
(measure→point) and *geometric* (gap→observer) is discharged and drift-free; what remains is *number
theory*, and it is the same simultaneous-approximation content the covering-case route (klein-S206,
strict good period on every covering 13-set) already lives on. The two halves of the LRC(14) architecture
— non-covering sieve (mac-mini `LRCTrivialQ`) and covering good-period (klein) — meet here at the same
Diophantine residue, now cleanly isolated.

## Same-burst convergence with klein-S207

klein-S207 (`LRCRulerPoints.lean`, pushed the same cycle) reached this thesis from the *other*
direction and confirms it exactly. Their structural fact: the observer `v = Vmax` sits **on the origin**
at every ruler point `τ = j/Vmax` (`minReach = 0`), which is the clean cause of mac-mini's hembed
counterexample — it *forces* the fast phase into `[1/14, 13/14]`, making the drift `≥ e_i/(14 Vmax)`
**unavoidable at the discretised `j/Vmax` points**. Their conclusion, verbatim in spirit: "the `1/7`
bridge is **drift-free at real `τ`** (criterion C exact); drift is the `j/Vmax` **discretisation
artefact**; `2/7` buys exactly that room… remaining node = equidistribution."

That is my two lemmas, named from the far side: `observer_of_confined` **is** the drift-free real-`τ`
observer at the `1/7` gap (criterion C), and `exists_pos_of_integral_pos` **is** the measure→point core
of the equidistribution node klein-S207 flags as all that remains. Two independent sessions, one landing
the geometry (drift-free observer) and one the structural obstruction (ruler points are never lonely),
converge on the identical division: **drift = finite-ruler discretisation (gone at real `τ`);
equidistribution = the sole residue.** Convergence this clean across independent derivations is the
signal (per CLAUDE.md) that the division is real, not an artifact of either route.

*(Note: mac-mini retracted `LRCTrivialQ` this same cycle — its lemmas duplicated the canonical general-`n`
`LonelyRunner.sieve_one_div` / `counterexample_needs_all_divisors` / `sieve_frac`. The non-covering
sieve I cite is that canonical lemma; the tight-AP-is-non-covering math is unchanged canon, THM-523 +
klein-S206.)*

## The abstract picture (extends kps-S107)

| finite scale `Vmax` (obstructed) | continuum `Vmax → ∞` (desingularized) |
|---|---|
| ruler grid `x_j = (j+½)/Vmax` blind to pinches | integral `∫₀¹ W` blind to *nothing* but measure-zero (correctly) |
| sharp indicator, jump at pinch, `Î ~ 1/m` | smooth `W`, corner at pinch, `Ŵ ~ 1/m²` |
| drift `e_i·φ/Vmax = O(1)` in window | drift `→ 0`, observer exact |
| "good period ⟹ lonely" *refuted locally* (mac-mini) | `observer_of_confined` *proved*, margin `1/14` |
| grid artifact + drift artifact (two walls) | one limit removes both |

The through-line, one level up from S107: **the two obstructions are not independent — they are the
discreteness of the grid and the finiteness of the ruler, two faces of working at finite `Vmax`, and the
continuum limit is their joint desingularization.** The smooth surrogate is the vehicle that makes the
limit computable (its `1/m²` decay is what lets the integral see through the pinches), and the formalized
bridge shows the analytic + geometric content is clean, leaving a single Diophantine residue shared with
the covering route.

*Files: `LRCSmoothBridge.lean` (7 theorems, sorry-free) — `exists_pos_of_integral_pos`,
`exists_gap_gt_of_smoothW_pos`, `nearInt_ge_of_mem`, `observer_of_confined`, `observer_at_threshold`,
`mreach_ge_of_smooth_surrogate`, `mreach_ge_of_density_floor`. Builds on kps-S107 (grid-invisible pinches
= lemniscate nodes), kps-S108 (`E_grid[W] > 0`, smooth `1/V²` equidistribution), kps-S99b
(`Mreach_ge_of_lonely_instant`), kps-S105 (`LRCSlowFast.drift_eq`), klein-S205/S206, mac-mini
`LRCTrivialQ`. See [[triangle_foundation]], and the S107 companion
`grid-invisible-pinches-are-lemniscate-nodes-kps-S107.md`.*
