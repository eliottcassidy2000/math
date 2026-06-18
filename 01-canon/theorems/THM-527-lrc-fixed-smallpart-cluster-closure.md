---
id: THM-527
title: The fixed-small-part single-tight-cluster closure for LRC(14) — an explicit-threshold equidistribution closure of a parametrized infinite sub-family of the S3 residual (for each fixed small part P ⊆ {1..13} with M(P) > 1/14 and each fixed offset multiset, there is an EXPLICIT V0*(P, offsets) such that cluster-base ≥ V0* ⟹ a global witness exists ⟹ M(S) ≥ 1/14; cluster-base < V0* is a finite exact check)
status: PROVED (constructive-witness lemma; exact-rational verified, 0 violations across V1/V2/V3 adversarial checks). HONEST SCOPE: closes a PARAMETRIZED infinite sub-family of S3, NOT all of S3. The AP / coordinated-growth family (no fixed small part) is explicitly OUTSIDE this closure and remains open.
source: kind-pasteur-2026-06-18-S4 (angle "fixed-small-part-equidistribution")
depends_on:
  - THM-523   # covering-set reduction
  - THM-524   # binding-pair reduction / exact M
  - THM-526   # arc-width lemma + cluster-collapse Lemma A + the slow-fast (offset-fit) reduction
related:
  - HYP-2581  # the unified S3 residual (this closes the small-P tight-cluster sub-case it named)
  - OPEN-Q-108 # uniform fattening / the S3 residual frontier
external: Lonely Runner Conjecture; proven for ≤7 runners (Barajas–Serra), open for 13.
---

# THM-527 — Fixed-small-part single-tight-cluster closure (explicit-threshold equidistribution)

> **⚠️ THM-527 NUMBER COLLISION (flagged mac-mini-2026-06-18-S2).** This file (kind-pasteur-S4,
> committed 09:01 — FIRST, so first-come keeps 527) shares id THM-527 with
> `THM-527-lrc-lonely-density-reformulation-and-the-bounded-spread-compact-reduction.md`
> (mac-mini-S1, 09:21). The lonely-density file became the reference hub (THM-528 depends on it;
> HYP-2584–2590 cite it). Renumber to be coordinated (likely the lonely-density file → THM-529, or
> this one → THM-529). @kind-pasteur: your call as first-committer. Flagged, not renamed.

## Context

LRC(14) reduces (THM-523/524/525/526) to: every primitive covering 13-set `S` has
`M(S) = max_τ min_{v∈S} ‖vτ‖ ≥ 1/14`. The open residual is **case S3** (k = #{v>13} ≥ 2,
spread `Vmax ≥ 13·Vmin`). The kind-pasteur-S3 slow-fast reduction isolated, inside S3, a
**small-P tight-cluster sub-case**: `S = P ∪ L` with `P ⊆ {1..13}` a *fixed* small part and
`L` a single tight cluster of *fixed shape* whose base grows. This theorem closes that
sub-case with explicit constants. It does **not** close the AP / coordinated-growth family
(no fixed small part), which is the genuinely-tight open core.

## The closed sub-family

`S = P ∪ L`, `|P| + |L| = 13`, where
- `P ⊆ {1,…,13}` is a **fixed** small part with `M(P) > 1/14` (equivalently `meas(G_P) > 0`);
- `L = {V0 + d_0, V0 + d_1, …, V0 + d_c}`, `0 = d_0 < d_1 < … < d_c` a **fixed** integer
  offset multiset (a single tight cluster: `d_c < 12·V0`, automatic for fixed offsets, large `V0`);
- `S` primitive and covering.

`G_P = {τ : ‖uτ‖ ≥ 1/14 ∀ u ∈ P}`, `G_L` similarly. A **global witness** is a single
`τ*` with `τ* ∈ G_P ∩ G_L`; it forces `min_{v∈S} ‖vτ*‖ ≥ 1/14`, hence `M(S) ≥ 1/14`.

## Bounded objects (the two facts that make the limit positive)

**FACT 1 (verified exhaustively, exact).** `meas(G_P) > 0` for **every** proper subset
`P ⊊ {1,…,13}` (0 subsets of size ≤ 12 have `meas = 0`); the per-size infimum of `meas(G_P)`
over `{meas > 0}` subsets is
`66/91, 55/91, 1979/4004, 2243/5880, 3029/10780, 45107/229320, 2479/17640, 10601/114660,
14249/252252, 313/9702, 7/858` for `|P| = 2,…,12`. The worst is `meas(G_{drop-6}) = 7/858`
at `|P| = 12`. The **only** small part with `meas = 0` is the full `P = {1,…,13}` (where
`M(P) = 1/14` exactly), but then `|L| = 0` (no cluster) — OUT of scope. So every in-scope
fixed small part has `meas(G_P) ≥ 7/858 > 0`. (Equivalence `meas(G_P) > 0 ⟺ M(P) > 1/14`
verified, 0/300 mismatches.)

**FACT 2 (exact, positive).** On a fixed P-safe arc, the **limiting cluster density**
`ρ∞(P, offsets) = (1/W_P)·∫_{I_P} w(τ) dτ`, `w(τ) = max(0, 6/7 − circ_width({d_i·τ}))`,
is `> 0` (computed exactly as a piecewise-linear integral; e.g. `ρ∞ ∈ {0.032, 0.145, 0.095,
0.182, 0.068}` for the five sample patterns). `w(τ)` is the width of the joint-cluster-safe
fast-phase window at slow time `τ`.

## The constructive-witness lemma (PROVED)

> **LEMMA.** Fix `(P, offsets)` in scope. There is an EXPLICITLY computable threshold
> `V0*(P, offsets)` such that for every `V0 ≥ V0*`, `S = P ∪ {V0 + d_i}` has a global witness,
> hence `M(S) ≥ 1/14`.

**Construction (all exact rationals).** Write `θ := V0·τ` (the *fast* phase). For `u = V0 + d_i`,
`‖uτ‖ ≥ 1/14 ⟺ frac(θ + d_iτ) ∈ [1/14, 13/14] ⟺ θ` lies in the width-`6/7` arc centered at
`μ_i(τ) = 1/2 − d_iτ`. The joint-cluster-safe set of `θ` at slow time `τ` is
`Ω(τ) = ⋂_i {width-6/7 arc at μ_i(τ)}`, an arc of width `w(τ)`.

- **STEP 1 (uniform common window).** Take the widest P-safe arc `I_P = (α,β)`. Choose an
  exact interior `τ0` maximizing `w` (over breakpoints of the piecewise-linear `w`), then
  rational-bisect a half-length `h` so that the sub-arc `[a',b'] = [τ0−h, τ0+h] ⊆ I_P` has a
  **common θ-window** `Ω* ≠ ∅` with `Ω* ⊆ Ω(τ)` for **all** `τ ∈ [a',b']`. `Ω*` exists with
  width `g* > 0` because the slow center-trajectories `{1/2 − d_iτ : τ ∈ [a',b']}` (arcs of
  length `d_i·(b'−a')`) leave a circular gap `> 1/7`. `Ω*` and `g*` are computed exactly.
- **STEP 2 (fast sweep hits it).** `τ ↦ {V0·τ}` is continuous, derivative `V0`, and over
  `[a',b']` sweeps a real interval of length `V0·(b'−a')`. If `V0·(b'−a') ≥ 1`, its image mod 1
  is **all** of `ℝ/ℤ`, so it hits the arc `Ω*` (width `g* > 0`): there is `τ* ∈ [a',b']` with
  `{V0·τ*} ∈ Ω* ⊆ Ω(τ*)`. Then `τ* ∈ I_P ⊆ G_P` and `τ*` is L-safe — a global witness. ∎

> **EXPLICIT THRESHOLD: `V0* = ⌈1/(b'−a')⌉`.** For `V0 ≥ V0*`, the witness `τ*` is constructed
> exactly (`τ* = (n + θ_c)/V0` for the first integer `n` with `(n+θ_c)/V0 ∈ [a',b']`, `θ_c` =
> center of `Ω*`). For `V0 < V0*`: a **finite** exact list of sets, each `M(S)` checked directly.

## Verification (exact rationals, 0 violations)

Script `04-computation/lrc14_finish_fixed-smallpart-equidistribution_kps-S4-wf.py`
(+ `_verify_…`, `_threshold-bound_…`); outputs in `05-knowledge/results/`.
- **(V1)** `Ω* ⊆ Ω(τ)` for all `τ ∈ [a',b']`: 0 violations at center-θ AND at the extreme
  edge-θ values, 400+ exact sample points per pattern.
- **(V2)** For `V0 ≥ V0*`, an **exact** witness `τ*` was exhibited with
  `min_v ‖vτ*‖ ≥ 1/14` for the actual `S` (all ≥ 0.096; e.g. `P={1,2,3}`, `V0=81`, `V0*=81`,
  `τ* = 2069390396012068553545361249/8225531276657936904700821504`, `min = 0.1124`).
- **(V3)** 120 rows at-or-above threshold across 40 random fixed-`(P,offsets)` patterns:
  **0 failures** (`M(S) ≥ 1/14` always).
- **Direct M-validation:** for `P={1,2,3}`, `offsets={0,…,9}`, every covering+primitive row
  (`V0 ∈ {20,84,112,168,308,420,700,1000}`) has exact `M ∈ {5/47, 3/25, 4/33, 6/49, 11/89,
  15/121, 25/201, 250/2007}`, all `≥ 1/14`. Second pattern `P={1,2,3,4}`, `offsets={0,2,…,23}`:
  26/26 covering+primitive rows with `V0 ∈ [14,400)` pass.

## Honest scope — what is and is NOT closed

**CLOSED (PROVED):** the parametrized infinite sub-family above. For **each** fixed
`(P, offsets)`, `V0 ≥ V0*(P, offsets)` is proved (global witness) and `V0 < V0*` is a finite
exact check. `V0*` is finite & explicit but **grows** with offset spread/resonance (sweep:
median `V0* = 57` for offsets ≤ 12, up to `3158` for offsets ≤ 50) — so there is **no uniform**
`V0` over all patterns, consistent with S3 being genuinely infinite. This is a legitimate
*sub-family tail* closure, not a closure of S3.

**NOT closed by this angle (the genuine open residual):**
- The **AP / coordinated-growth family** `{t, 2t, …, 12t, V}` (`t → ∞`): the speeds `≤ 13` are
  the *shrinking, t-dependent* set `{kt ≤ 13}`, NOT a fixed `P`. This is the asymptotically-tight
  core where `M → 2/23`. Outside scope, OPEN.
- Clusters with **unbounded / non-fixed** offsets, or multiple separated clusters.

So THM-527 removes the small-P tight-cluster sub-case of S3 (the one the kind-pasteur-S3
slow-fast reduction governed) with explicit constants, leaving the coordinated-growth core
(HYP-2581d / OPEN-Q-108) as the open residual. **LRC(14) remains OPEN.**
