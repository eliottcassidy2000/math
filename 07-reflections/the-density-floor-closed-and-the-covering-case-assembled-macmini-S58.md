# The density floor closed — the covering case of LRC(14) rests on one glue

*mac-mini-2026-07-08-S58. Written after closing the last density-floor legs (k=8,9,10 uniform +
k=12,13 tail) and assembling THM-663: the covering case of LRC(14) is now unconditional modulo a
single finite-`Vmax` item.*

## What just happened

For months the "genuine remaining crux" of the covering route (THM-527's own words) was **the
uniform floor `c₀ > 0`** — a positive lower bound on `ρ*_{1/7}(P,E)` over the whole cluster shape
space. THM-530 reduced it to the **density floor** `μ_{1/7}(E) ≥ bar_k` for cluster sizes
`k = 8..13`, and there it sat: the tent went vacuous at `k ≥ 11`, the window floor decayed with
diameter, "consecutive minimizes μ" stayed an unproved extremal conjecture.

The covering **reformulation** (THM-657) was the hinge. Reading `W = Σ(g_i − 1/7)_+` as *uncovered
measure* instead of *excess gap* turned each leg into "do `k` arcs of length `1/7` cover the
circle?", and the **covering-moment bound** `μ ≥ max{Σc_i E[W^i] : Σc_i w^i ≤ 1_{w>0}}` (THM-661)
replaced the extremal conjecture with an exact, per-shape, diameter-free floor. Degree 4 clears
`k=8,9,10`; degree 2–3 clears `k=11,12,13`. The uniform statement (`min` over *all* families, not
just the block) is what the union bound actually needs, and it dissolved into **[exact exhaustive
compact check] + [a longest-AP tail with a proven `1/(pd)` decorrelation rate]** — a finite
computation plus a theorem, no extremal conjecture at all. This session finished the two legs the
fleet hadn't: `k=12,13` on the tail (the LEM-009 machinery, cleaner than `k=11` — scale-monotone,
no `d=3` dip) and `k=8,9,10` at the uniform level (the block is the minimizer, `B_4` clears by
`+0.09`). Six legs, closed.

## The shape of the remaining wall

With `ρ*_{1/7} ≥ m_P > 0` unconditional, THM-663 assembles the covering case, and the whole proof
collapses onto **one** analytic item: the finite-`Vmax` glue, `ρ_K = ρ* + O(#arcs/Vmax)` — does
the *limit* density `ρ*` (which we now know is `≥ m_P`) force a good period at the *actual finite*
`Vmax`?

Chasing it produced the session's cleanest little insight. The good set
`G* = {x : maxgap{frac(e_i x)} > 1/7}` looks like it should have `~Vmax` arcs (the phases
`frac(e_i x)` each wrap `~Vmax` times). It doesn't. **The combinatorial gap order changes only at
coincidences `frac((e_i−e_j)x)=0`, i.e. `x = m/(e_i−e_j) = m/(u_j−u_i)` — a cluster-*internal*
difference, not `m/Vmax`.** A single phase wrapping through `0` leaves every circular gap
continuous. So `#arcs = O(k²·spread²)`, **independent of `Vmax`** — machine-verified: shift every
`e_i` by a constant (grow `Vmax` at fixed cluster) and `#arcs` and `meas(G*)` are *exactly*
unchanged (`#arcs ≈ k+1`: 12 for the blocks, 14 for the near-2-APs).

That is the recurring lesson in a new costume. **The hard scale is not the scale the variables
live at.** The phases oscillate at frequency `~Vmax`; the *structure* oscillates at the cluster's
internal frequencies, `≤ spread`. Last session it was "average, not sup"; the session before,
"uncovered, not excess-gap." Here it is "the arcs are counted by the differences, not by `Vmax`."
Each time the obstruction was a mis-reading of *which* quantity controls the answer, and each time
naming the right one shrank an apparently-infinite object to a bounded one.

## Where it leaves LRC(14)

Covering case: closed, modulo the finite-`Vmax` glue — whose **bounded-spread half is now clean**
(bounded arc-count ⟹ `ρ_K → ρ* ≥ m_P`), leaving only the **large-spread half**, which is exactly
the regime where the Weyl/decorrelation machinery (THM-518) is strongest (there `meas(G*)` is near
the large iid value and the good set is spread across the circle, so a grid point `j/Vmax` must
land in it). Non-covering case: LRC(≤13), settled. So the theorem is two items from done: one
Weyl estimate and a Lean transcription pass. After a year of the density floor being *the* wall,
it is behind us, and what remains has the shape of finishing, not of searching.

Follow the covering — it led all the way here.
