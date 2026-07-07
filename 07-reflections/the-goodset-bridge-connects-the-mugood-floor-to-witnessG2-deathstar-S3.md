---
source: death-star-2026-07-07-S3
status: FORMALIZATION WIN. A complete kernel-pure chain from the bounded-diameter muGood floor to
  the exact witnessG2 quantity: `(slowμ(goodSet E)).toReal ≥ m_P` for any nonempty diameter-≤30 E —
  the `hlarge` density-floor obligation in the pure-cluster (P-empty) case.
tags:
  - lonely-runner
  - LRC14
  - lean
  - formalization
  - density-floor
  - witnessG2
  - bridge
---

# The goodSet bridge connects the muGood floor to `witnessG2`

Owner: *work the natural next step (the `witnessG2 = muGood` wiring), continuing next steps as
they come in.* I closed the bridge between the two "good set" formalisms and carried the floor all
the way to the exact real quantity `witnessG2` is.

## The two good-sets, and why they were disconnected

- `LRCTailDiameter.Good (1/7) E` (real-anchored: `∃ real a, ∀ e, frac(e·x−a) > 1/7`) — what the
  muGood diameter floors (death-star-S2 `muGood_diam_floor` `D≤30`; boxeph `AP30`/`AP44`/`AP76`) bound.
- `GoodSet.goodSet E` (phase-anchored: `∃ c∈E, ∀ b, frac((b−c)x) ∉ (0,1/7]`) — the CONCRETE
  measurable event whose `slowμ`-measure IS the skeleton's `witnessG2` (via `LRCEventMeasureBridge`).

Both mean "max-gap > 1/7", but the muGood machinery lived on `Good` and `witnessG2` on `goodSet`.

## What I proved (`TournamentH7/LRCGoodSetBridge.lean`, kernel-pure, in root import)

1. **`Good_subset_goodSet`** — `Good (1/7) E.toFinset ⊆ goodSet E`. Proof: take
   `c = argmax_{e∈E} frac(e·x−a)` (`Finset.exists_max_image`); then for every `b`,
   `u_b − u_c ≤ 0` with `u_b,u_c ∈ (1/7,1)`, so `frac((b−c)x) = frac(u_b−u_c)` is `0` (if
   `u_b=u_c`) or `u_b−u_c+1 ∈ (1/7,1)` — never in `(0,1/7]`. (Uses `GoodSet`'s own
   `fract_sub_eq_fract_fract_sub_fract`.)

2. **`slowμ = volume.restrict (Ico 0 1)`** (found in `LRCDenseCovers`) — so
   `slowμ(goodSet E) = volume(goodSet E ∩ Ico 0 1)`, and (modulo the null endpoint `{1}`,
   handled via `Ico_ae_eq_Icc`) the muGood floor transfers:
   **`slowμ_goodSet_ge_mP`** — `ofReal(m_P) ≤ slowμ(goodSet E)`, and
   **`slowμ_goodSet_toReal_ge_mP`** — `m_P ≤ (slowμ(goodSet E)).toReal`, for any nonempty `E`
   whose translate lies in `{0,…,D}`, `2 ≤ D ≤ 30`.

**Axiom audit:** all four theorems depend on exactly `[propext, Classical.choice, Quot.sound]`.
Full project builds green (8479 jobs).

## Where this plugs in

The skeleton splits the density floor by `clusterSize (shapeOf v)`: `k ≤ 7` is pigeonhole-elementary
(`maxgap ≥ 1/k ≥ 1/7`); the open leg is **`hlarge` (`k = 8..13`)** — `witnessMP ≤ witnessG2`.
A census family like `{8,…,20}` has co-offset cluster `{0,…,12}` (`k=13`, diameter `≤ 19`), so the
`hlarge` leg is exactly the bounded-big-cluster regime my floor covers. And
`witnessG2 s = (slowμ(goodSet(cluster s) ∩ GP s)).toReal` (`LRCEventMeasureBridge.hwitness`), so in
the **pure-cluster case** (`P = ∅`, `GP = univ`), `witnessG2 = (slowμ(goodSet cluster)).toReal`, and
`slowμ_goodSet_toReal_ge_mP` **IS** the `hlarge` floor `witnessMP ≤ witnessG2` for a nonempty cluster
of diameter `≤ 30`. (NB: no file currently consumes `hlarge`/`AP76Certificate` — this leg is open;
my chain supplies the measure fact it needs.)

**Remaining wiring (mac-mini/kps shape layer):** (a) the `shapeOf → cluster` co-offset extraction
and the `hwitness` identification for concrete shapes; (b) the `GP` intersection for non-pure-cluster
shapes (Bonferroni: `witnessG2 ≥ nuShape + measGP − 1`, needs `nuShape = slowμ(goodSet)` large — for
small clusters it's the pigeonhole `1`; my `m_P` floor is the pure-cluster / `GP=univ` slice).

## Ledger

- **Files:** `TournamentH7/LRCGoodSetBridge.lean` (kernel-pure). Builds on death-star-S2
  `LRCAP20Certificate` (`muGood_diam_floor` `D≤30`), mac-mini `LRCTailDiameter`/`LRCGoodSet`,
  kps `LRCDenseCovers` (`slowμ`), the `LRCEventMeasureBridge` `hwitness` interface.
- **Complements:** boxeph `LRCFareyRoofBridge`/`AP44` (wide regime, `D≤43`, toward `AP76`).
- **Does NOT prove LRC(14).** Connects the (unconditional, bounded-diameter) density floor to the
  exact `witnessG2` measure; the final `census/hlarge` discharge needs the shape/GP layer.
