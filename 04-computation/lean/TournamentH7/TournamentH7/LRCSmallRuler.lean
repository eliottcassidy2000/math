/-
# LRCSmallRuler — the covering hypothesis calibrates the ruler window exactly at 15

Attacking `StrictlyLiveSupply` (the wall).  Two facts, one proved here and one measured.

**Proved.**  A zero residue kills a modulus outright: if `q ∣ vᵢ` then `(vᵢ·p) % q = 0` for every `p`, and
`0` never lies strictly inside the band `(q/14, 13q/14)`.  Since `CoveringFamily v` supplies, for *every*
`q ∈ [2,14]`, some `i` with `q ∣ vᵢ`, **no covering family has a strictly-live ruler below 15**:

  `strictlyLive_modulus_ge_15 : CoveringFamily v → StrictlyLive v q p → 15 ≤ q`.

So the covering branch of the assembly does not merely dispatch small-`q` families — it *calibrates the
search window*, closing exactly `[2,14]` and no more.

**Measured** (`lrc14_small_ruler_law_kps_S127`).  Over covering, ratio-`> 13` families the minimal
strictly-live modulus is `15 ≤ min q ≤ 26`, **independent of scale** (`Vmax` from 45 to 1200, zero
failures, mean ≈ 17).  Adversarial hill-climbing reaches only 30.  The one adversary that pushes it
higher — a speed divisible by many moduli, killing each by a zero residue — makes the family *detuned*,
and the residual's `hdiv` hypothesis excludes exactly that.  On the residual class (covering, ratio > 13,
compressed, not detuned, difference-primitive) the observed maximum of `min q` is **25**.

This suggests the sharpened conjecture recorded below as `BoundedStrictlyLiveSupply`: the wall's witness
lives in a *bounded* window `[15, B]`, so the remaining obligation is a small-modulus, residue-level
statement rather than an asymptotic one.  It is a conjecture, not a theorem.

kind-pasteur-2026-07-10-S127.
-/
import Mathlib
import TournamentH7.LRCStrictRuler

namespace LonelyRunner
namespace LRC14Grand

/-- **A zero residue kills the modulus.**  If `q ∣ vᵢ` then `(vᵢ·p) % q = 0` for every `p`, and `0` is
never strictly above `q/14`.  So no multiplier is live at `q`. -/
theorem not_strictlyLive_of_dvd {v : Fin 13 → ℤ} {q p : ℤ} {i : Fin 13}
    (hq : 0 < q) (hdvd : q ∣ v i) : ¬ StrictlyLive v q p := by
  rintro ⟨hp0, hpq, hband⟩
  have hz : (v i * p) % q = 0 := Int.emod_eq_zero_of_dvd (hdvd.mul_right p)
  have hb := (hband i).1
  rw [hz] at hb
  omega

/-- **Covering calibrates the ruler window: every strictly-live modulus of a covering family is `≥ 15`.**
Covering supplies a zero residue at every `q ∈ [2,14]`, and `0` is never inside the band. -/
theorem strictlyLive_modulus_ge_15 {v : Fin 13 → ℤ} {q p : ℤ}
    (hcov : LRC14.CoveringFamily v) (h : StrictlyLive v q p) : 15 ≤ q := by
  have hp0 : 0 < p := h.1
  have hpq : p < q := h.2.1
  have hq : 0 < q := lt_trans hp0 hpq
  have hq2 : 2 ≤ q := by omega
  by_contra hlt
  push_neg at hlt
  obtain ⟨i, hi⟩ := hcov q.toNat (by omega) (by omega)
  have hdvd : q ∣ v i := by rwa [Int.toNat_of_nonneg (le_of_lt hq)] at hi
  exact not_strictlyLive_of_dvd hq hdvd h

/-- Every residual family is covering, so its strictly-live rulers all lie at `q ≥ 15`. -/
theorem residual_strictlyLive_ge_15 {v : Fin 13 → ℤ} {q p : ℤ}
    (hres : IsResidual v) (h : StrictlyLive v q p) : 15 ≤ q :=
  strictlyLive_modulus_ge_15 hres.2.1 h

/-! ### The sharpened conjecture: the witness lives in a bounded window. -/

/-- **`BoundedStrictlyLiveSupply B`** — every residual family has a strictly-live ruler at some modulus
`q ≤ B`.  Combined with `residual_strictlyLive_ge_15`, the witness is confined to `[15, B]`.

Measured: `B = 25` suffices on every residual family sampled (scale-independent, `Vmax ≤ 1200`).  This is
a CONJECTURE; `StrictlyLiveSupply` (its `B = ∞` form) is the open case of LRC(14). -/
def BoundedStrictlyLiveSupply (B : ℤ) : Prop :=
  ∀ v : Fin 13 → ℤ, IsResidual v → ∃ q p : ℤ, q ≤ B ∧ StrictlyLive v q p

/-- A bounded supply is a supply. -/
theorem strictlyLiveSupply_of_bounded {B : ℤ} (h : BoundedStrictlyLiveSupply B) :
    StrictlyLiveSupply := by
  intro v hres
  obtain ⟨q, p, _, hl⟩ := h v hres
  exact ⟨q, p, hl⟩

/-- **LRC(14) from a bounded strictly-live ruler window.**  Kernel-pure.  With
`residual_strictlyLive_ge_15` this says: to finish LRC(14) it suffices to exhibit, for each residual
family, a modulus `q ∈ [15, B]` and a multiplier `p` with all thirteen residues strictly inside the band. -/
theorem lrc14_of_boundedStrictlyLiveSupply (cite : LRCUpTo13) {B : ℤ}
    (h : BoundedStrictlyLiveSupply B) : LRC14.LRC14Statement :=
  lrc14_of_strictlyLiveSupply cite (strictlyLiveSupply_of_bounded h)

end LRC14Grand
end LonelyRunner
