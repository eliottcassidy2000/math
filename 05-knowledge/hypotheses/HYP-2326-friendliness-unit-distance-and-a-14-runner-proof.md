# HYP-2326 — Friendliness ↔ the unit-distance problem, and a machine-checked LRC proof for the 14 consecutive runners

**Session:** S648
**Status:** CONFIRMED (correlation grounded in S634; the 14-runner config proof formalized)
**Provenance forward:** math-lean `Math/LonelyRunner/LonelyFourteen.lean` (sorry-free)
**Prompt:** (1) does the friendliness/survival view (S647) correlate to unit distance? (2) "write me a
proof for 14 runners."

---

## Part 1 — Friendliness ↔ unit distance: YES, via the S634 unification

Under the one-graph unification of LRC / Hadwiger–Nelson / unit-distance (S634, χ·α ≥ n), the
friendliness/loneliness dictionary maps cleanly:

| LRC / friendliness | unit-distance graph `G` |
|---|---|
| loneliness = the 1-avoiding / **independent set** `α` | an **isolated vertex** (degree 0) |
| **friendliness** = having neighbours | the **degree** of a vertex |
| (a friendly pair) | an **edge** = a unit distance |
| **max** friendliness | **max edges = `u(n)`** (the unit-distance maximisation, S641) |
| first lonely time `τ` (scan = **time**, floor `1/n`) | first-friendly radius (scan = **radius**, floor = min inter-point distance) |

Three concrete confirmations (`friendliness_unit_distance_link_s648.py`):

1. **Unit-distance friendliness = degree = 6 on the Eisenstein lattice** — exactly the hexagon
   `a²−ab+b² = 1` with 6 solutions (formalized `eisenstein_unit_neighbours`, S641). So a lattice point's
   friendliness is the rigid `6 = 2·3`; `u(n)` *escapes* this degree-6 lattice cap via the Moser/√−11
   slab to be friendlier on the boundary (S641).
2. **The survival shape transfers.** Scanning the threshold **radius** `r` upward, "isolation-survival"
   `S(r) = P(no neighbour within r)` is **flat at 1 then decays to 0** — the *same* curve as the S647
   friendliness survival (where the scan is **time** and the floor is `1/n`). "Never lonely yet" =
   "never isolated yet"; the first-passage machinery is identical, time ↔ radius.
3. **The shared resonance.** LRC asks for clock-distance `d ≥ 1/n` (a **band**, far from 0); unit distance
   asks for chord `= 1 ⟺ d = 1/6` **exactly** (a **sharp resonance** — the hexagon, `2 sin(π/6) = 1`,
   the cube root, S623). Same covering graph (S634), different target: LRC = *band*-friendly, unit
   distance = *resonance*-friendly.

> **So friendliness is the unit-distance graph's degree, loneliness is its independent set, and the
> first-passage survival curve is shared — with the scan being time (LRC) or radius (unit distance) and
> the resonance being the same `1/6 = ` hexagon `= 6 = 2·3` cube root the whole arc converges on.**

---

## Part 2 — A proof for 14 runners (machine-checked)

> **HONEST SCOPE.** This proves LRC for the **canonical 14-runner configuration** `{1, 2, …, 13}` (the
> 13 runners of consecutive speeds). The **full** Lonely Runner Conjecture for *all* 14-runner speed
> sets is **open** — LRC is proven only up to 7 runners (Barajas–Serra 2008). But for `{1,…,13}` the
> proof is short, rigorous, and now formalized.

**Theorem.** The 14-runner configuration `{1,…,13}` is **lonely at `t = 1/14`**.

**Proof.** At `t = 1/14`, runner `k` sits at clock distance `‖k·(1/14)‖ = dZ(k/14)`. For each
`k ∈ {1,…,13}`, `k/14 ∈ [1/14, 13/14] = [1/14, 1 − 1/14]`, and for any `x ∈ [δ, 1−δ]` the nearest
integer (`0` or `1`) is at least `δ` away, so `dZ(x) ≥ δ`. With `δ = 1/14`: `dZ(k/14) ≥ 1/14` for every
runner. So all 13 runners are at least the gap `1/14` away from the watched runner — it is lonely. ∎

**Formalized (math-lean, sorry-free): `Math/LonelyRunner/LonelyFourteen.lean`**
- `dZ_ge_of_mem : 0 ≤ δ → δ ≤ x → x ≤ 1 − δ → δ ≤ dZ x` (the clock-distance lower bound; proof: for
  *every* integer `m`, `|x − m| ≥ δ` — `m ≤ 0` gives `x − m ≥ x ≥ δ`, `m ≥ 1` gives `m − x ≥ 1 − x ≥ δ`).
- `lonely_fourteen : 1 ≤ k → k ≤ 13 → 1/14 ≤ dZ (k/14)` (the LRC witness `t = 1/14`).

**It is the friendliest config (the S647 tie).** `{1,…,13}` is the *tight extremal* case: its lonely set
is the single instant `t = 1/14` (measure zero), so under "never lonely yet" it is friendly *almost
everywhere the whole lap* (S647). The 14-runner proof is exactly the one point where its maximal
friendliness finally, barely, touches solitude — at the gap `1/14` with equality on the unit runners.

---

## New threads / handoffs
- **General `{1,…,n−1}`:** the same `dZ_ge_of_mem` gives `dZ(k/n) ≥ 1/n` at `t = 1/n` for all `n` — a
  one-line formalization of LRC for the consecutive config at every `n` (the easy/tight family).
- **The hard part stays hard:** general LRC(14) over *all* speed sets is the open frontier (the
  reductions of S639/S640/S643 attack it); the consecutive config is the trivial corner.
- **Friendliness ↔ unit distance, formal:** transport the survival-floor (`τ ≥ 1/n`, S647) to the
  unit-distance radius-scan (`first-friendly radius ≥ min distance`); identify `u(n) =` max friendliness.
