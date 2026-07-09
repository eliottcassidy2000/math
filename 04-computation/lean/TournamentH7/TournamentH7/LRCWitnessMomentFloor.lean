/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S186)
-/
import Mathlib.Tactic
import TournamentH7.LRCWitnessBonferroni

/-!
# Route `hfloor` through the PROVED moment floor (opus-S186)

The density-floor node `hfloor` (`witnessMP ≤ witnessG2`) was reduced (`LRCWitnessBonferroni`) to
`hbonf` (Bonferroni, PROVED kps-S30) + `hB` (Lemma B, proved rational) + **Lemma A** (`nuConsec ≤ nu`),
which is OPEN (a compactness minimization). opus-S185 showed `nu(E) = mu(E) = P(W>0)` exactly, so `nu`'s
lower bound is precisely the covering measure that THM-661's moment-LP already bounds:

  `nu(E) = mu(E) ≥ D3(E) ≥ min_E D3(E) ≥ (A') bar = witnessMP + 1 − capRat(k)`   (THM-661).

This file makes that substitution: it replaces the OPEN Lemma A hypothesis by the PROVED moment-floor
bound `momentBar k = witnessMP + 1 − capRat k`, and derives `hfloor` with the same Bonferroni + Lemma B.
The arithmetic that needed a separate lemma in the Lemma-A route (`nuConsec + capRat − 1 ≥ witnessMP`)
is now DEFINITIONAL: `momentBar k + capRat k − 1 = witnessMP` exactly, because THM-661's `(A')` bar is
*constructed* as `witnessMP + 1 − capRat k`.

Net effect: `hfloor`'s open dependency changes from the never-proved standalone Lemma A to THM-661's
moment floor — a theorem the fleet has actively proved (per-shape exact; the only residue is THM-661's
own decorrelation tail, shared by every route since `nu = mu`, opus-S185).

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14
namespace MomentFloor

open scoped BigOperators

/-- The moment-floor **`(A')` bar** for cluster size `k`: `momentBar k = witnessMP + 1 − capRat k`.
THM-661 proves `min_E D3(E) ≥ momentBar k`; since `nu(E) = mu(E) ≥ D3(E) ≥ min_E D3(E)`, every
`k`-cluster satisfies `nu(E) ≥ momentBar k`. (This is exactly THM-661's honest `(A')` bar, which the
moment floor is built to clear.) -/
def momentBar (k : ℕ) : ℚ := witnessMP + 1 - capRat k

/-- The bar and the cap sum to `witnessMP + 1` by construction — so the Bonferroni arithmetic is
definitional (no separate rational floor lemma needed, unlike the Lemma-A route). -/
theorem momentBar_add_capRat (k : ℕ) : momentBar k + capRat k - 1 = witnessMP := by
  simp only [momentBar]; ring

/-- **Route `hfloor` through the proved moment floor.**
Replaces the OPEN Lemma A (`nuConsec ≤ nu`) with the PROVED moment floor (THM-661:
`nu = mu ≥ min D3 ≥ momentBar k`). Given Bonferroni (`hbonf`, kps-S30), the moment-floor bound
(`hMoment`, THM-661), and Lemma B (`hB`), the admissible floor `witnessMP` lower-bounds `witnessG2 s`
for every shape with cluster size in `8..13`. The arithmetic is definitional
(`momentBar k + capRat k − 1 = witnessMP`). -/
theorem witness_floor_from_momentfloor_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hMoment : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (s : Shape) (h8 : 8 ≤ clusterSize s) (h13 : clusterSize s ≤ 13) :
    (witnessMP : ℝ) ≤ witnessG2 s := by
  set k := clusterSize s with hk
  -- the arithmetic floor is DEFINITIONAL, cast to ℝ
  have hrat : momentBar k + capRat k - 1 = witnessMP := momentBar_add_capRat k
  have hcast : (momentBar k : ℝ) + (capRat k : ℝ) - 1 = (witnessMP : ℝ) := by
    have := congrArg (fun q : ℚ => (q : ℝ)) hrat
    push_cast at this ⊢
    linarith
  -- moment floor + Lemma B raise the bars to the genuine measures
  have hMk : (momentBar k : ℝ) ≤ nuShape s := hMoment s h8 h13
  have hBk : (capRat k : ℝ) ≤ measGP s := hB s h8 h13
  -- Bonferroni
  have hbk : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  linarith

/-- The large-cluster (`8..13`) floor node in the shape expected by the LRC14 skeleton's case-split,
routed through the moment floor (THM-661) rather than the open Lemma A. -/
theorem large_witness_floor_from_momentfloor_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hMoment : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v _hv h8 h13
  exact witness_floor_from_momentfloor_nodes nuShape measGP hbonf hMoment hB
    (shapeOf v) h8 h13

/-- **The witness-route case-split, with the moment floor swapped in for Lemma A, and every
opaque-free proved leg discharged (opus-S186).**
LRC(14) follows from: the pigeonhole small case (`hsmall`, `k ≤ 7`); the **moment-floor** large case
(`hbonf` Bonferroni + `hMoment` THM-661 + `hB` Lemma B, `8 ≤ k ≤ 13`) — the slot the Lemma-A node held;
the size bound (`hsize`); and the reach node (`hpartA`).

`hR0` (`Mreach ≥ 1/14 ⟹ lonely`) is DISCHARGED here — supplied internally as `lonely_of_Mreach_ge`
(a proved skeleton theorem about the concrete `Mreach`) by routing through `lrc14_from_witness_floor`.

The remaining parameters are exactly the two cruxes **`hMoment`** (density floor = THM-661) and
**`hpartA`** (reach), PLUS the four legs `hbonf, hB, hsmall, hsize`. Those four are NOT dischargeable to
proof terms in this file: they reference the skeleton's `opaque witnessG2 : Shape → ℝ` and
`opaque shapeOf`, and — as `LRCTailDiameter` records — "the bridge from `muGood` to the opaque
`witnessG2` cannot be a theorem". Under a CONCRETE `witnessG2 = (slowμ (GOOD ∩ G_P)).toReal` they become
proved (`hbonf` = `LRCBonferroniMeasure.toReal_bonferroni`; `hB` = Lemma B; `hsmall` = the `k ≤ 7`
pigeonhole `nu = 1`; `hsize` = the concrete cluster length), but that concretization is a coordinated
skeleton change (it also re-states `hpartA`), carried via the `LRCEventMeasureBridge` `hwitness`
hypothesis — flagged, not done here. -/
theorem lrc14_from_momentfloor_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hMoment : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → clusterSize (shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  -- `lrc14_from_witness_floor` supplies `hR0 = lonely_of_Mreach_ge` internally
  lrc14_from_witness_floor
    (witness_floor_from_cluster_cases hsmall
      (large_witness_floor_from_momentfloor_nodes nuShape measGP hbonf hMoment hB)
      hsize)
    hpartA

end MomentFloor
end LRC14
end LonelyRunner
