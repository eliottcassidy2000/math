/-
  TournamentH7.LRCDeepCertificate — THE DECIDABLE CENSUS PIPELINE
  (death-star-2026-07-17-S43, HYP-7204; the end-to-end integration of the census
  funnel: THM-945's capped identities + THM-950's census criterion + the
  Mreach consumer, all behind DECIDABLE gates).

  `CoverageCapped`, `liveCount`, `bandCount`, and the deep count are all computable
  over the finite multiplier range — so the entire census funnel is DECIDABLE for
  concrete `(v, q)`:

  * `instDecidableCoverageCapped` — the gate: `CoverageCapped v q C` decides.
  * `lonely_of_census` — THE PIPELINE: for nonzero speeds,
    `CoverageCapped v q 6` + `deepSix < liveCount` ⟹ `∃ t, Lonely 14 v t`
    (THM-945's exact identity `B5 = live − deepSix` + `mreach_ge_of_B5_pos` +
    `lonely_of_Mreach_ge`).  With the cap and the census both decidable, concrete
    families acquire END-TO-END machine-checkable loneliness certificates through
    the B5 funnel — no witness search, no analysis.
  * `lonely_of_census_capped5` — the cap-5 instance: `CoverageCapped v q 5` +
    `0 < liveCount` suffices (`B5 = liveCount` there).
  * `census_demo` — a concrete family certified end-to-end by `decide` through the
    pipeline (the first loneliness proof in the corpus produced by the census
    funnel rather than a direct witness or a pack row).

  Kernel-pure except the demo's `decide` (finite, kernel-executed).
-/
import Mathlib
import TournamentH7.LRCDeepCount

namespace LonelyRunner
namespace LRC14Concrete

open Finset

instance instDecidableCoverageCapped (v : Fin 13 → ℤ) (q C : ℕ) :
    Decidable (CoverageCapped v q C) := by
  unfold CoverageCapped
  infer_instance

/-- **THE CENSUS PIPELINE**: cap at 6 plus a live majority over depth-six forces a
lonely instant — end-to-end, decidably. -/
theorem lonely_of_census (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hv : ∀ i, v i ≠ 0)
    (hcap : CoverageCapped v q 6)
    (hrace : (((Finset.Ioo 0 q).filter fun p => bandCount v q p = 6).card)
      < liveCount v q) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hid := B5_eq_live_sub_deepSix v q hcap
  have hpos : 0 < B5 v q := by
    rw [hid]
    have h : (((Finset.Ioo 0 q).filter fun p => bandCount v q p = 6).card : ℤ)
        < (liveCount v q : ℤ) := by exact_mod_cast hrace
    linarith
  exact lonely_of_Mreach_ge v hv (mreach_ge_of_B5_pos v q hq hpos)

/-- The cap-5 instance: `B5 = liveCount` there, so one live multiplier suffices. -/
theorem lonely_of_census_capped5 (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hv : ∀ i, v i ≠ 0)
    (hcap : CoverageCapped v q 5)
    (hlive : 0 < liveCount v q) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hid := B5_eq_live_capped5 v q hcap
  have hpos : 0 < B5 v q := by
    rw [hid]
    exact_mod_cast hlive
  exact lonely_of_Mreach_ge v hv (mreach_ge_of_B5_pos v q hq hpos)

/-- **The unconditional weighted census pipeline.**  No coverage cap is
needed: THM-950 charges every multiplier of depth at least six at the sharp
universal pointwise cost `792`.  Thus one finite, decidable live/deep race
already produces an actual lonely time. -/
theorem lonely_of_weighted_census (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hv : ∀ i, v i ≠ 0)
    (hrace : 792 * ((Finset.Ioo 0 q).filter
      (fun p => 6 ≤ bandCount v q p)).card < liveCount v q) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hraceZ :
      792 * (((Finset.Ioo 0 q).filter
        (fun p => 6 ≤ bandCount v q p)).card : ℤ) <
        (liveCount v q : ℤ) := by
    exact_mod_cast hrace
  exact lonely_of_Mreach_ge v hv
    (mreach_ge_of_B5_pos v q hq (B5_pos_of_live_beats_deep v q hraceZ))

/-- **The demo**: the arithmetic family `v_i = i + 2` is certified lonely
END-TO-END through the census funnel at `q = 31` — cap and census both by
`decide`; the first loneliness proof in the corpus produced by the B5 pipeline
rather than a direct witness or a pack row. -/
theorem census_demo :
    ∃ t : ℝ, Lonely 14 (fun i : Fin 13 => (i : ℤ) + 2) t := by
  apply lonely_of_census (fun i : Fin 13 => (i : ℤ) + 2) 31 (by norm_num)
    (fun i => by
      show (i : ℤ) + 2 ≠ 0
      have h := Int.natCast_nonneg (i : ℕ)
      omega)
  · decide
  · decide

/-! ## Axiom audit -/
#print axioms lonely_of_census
#print axioms lonely_of_census_capped5
#print axioms lonely_of_weighted_census
#print axioms census_demo

end LRC14Concrete
end LonelyRunner
