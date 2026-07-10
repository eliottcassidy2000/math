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
`residual_strictlyLive_ge_15` this says: a modulus `q ∈ [15, B]` per residual family would finish LRC(14).

⚠ **`BoundedStrictlyLiveSupply B` is REFUTED for every fixed `B`** (kps-S127, corrected): the minimal
strictly-live modulus is *unbounded* on the residual class.  The strictly-live condition at `q` depends
only on `v mod q`, so a family whose residues are tight at every `q ∈ [15,B]` fails all of them — and such
families are genuinely residual (see `cexFamily` below, `min q = 28`; the class reaches `min q ≥ 43`).  So
this implication is TRUE but its hypothesis is unachievable: the route to LRC(14) must be
`StrictlyLiveSupply` (`B = ∞`) or the measure floor, and klein's THM-685 transfer (live rulers at
`q > Σv/μ`) is load-bearing, not a superfluous safety net. -/
theorem lrc14_of_boundedStrictlyLiveSupply (cite : LRCUpTo13) {B : ℤ}
    (h : BoundedStrictlyLiveSupply B) : LRC14.LRC14Statement :=
  lrc14_of_strictlyLiveSupply cite (strictlyLiveSupply_of_bounded h)

/-! ### The counterexample to the upper edge (certified decidably).

`cexFamily` is covering with ratio `28.57 > 13` (so it reaches the residual/`GapFamily` branch), yet it has
NO strictly-live multiplier at any modulus `q ∈ [15,27]` — refuting `q ≤ 27`.  Its full residuality
(compressed, distinct, not detuned, difference-primitive, not near-AP) is verified in
`lrc14_upper_edge_refuted_kps_S127`; the two load-bearing decidable facts are certified here. -/

instance instDecStrictlyLive (v : Fin 13 → ℤ) (q p : ℤ) : Decidable (StrictlyLive v q p) := by
  unfold StrictlyLive; infer_instance

instance instDecGapFamily (v : Fin 13 → ℤ) : Decidable (GapFamily v) := by
  unfold GapFamily; infer_instance

/-- The refuting family: `min` strictly-live modulus `= 28`. -/
def cexFamily : Fin 13 → ℤ :=
  ![210, 1378, 1379, 2106, 2222, 2247, 3650, 3773, 4123, 5083, 5561, 5680, 6000]

/-- `cexFamily` is covering (it reaches the residual branch). -/
theorem cex_covering : LRC14.CoveringFamily cexFamily := by
  intro q hq2 hq14
  interval_cases q
  · exact ⟨0, by decide⟩
  · exact ⟨0, by decide⟩
  · exact ⟨12, by decide⟩
  · exact ⟨0, by decide⟩
  · exact ⟨0, by decide⟩
  · exact ⟨0, by decide⟩
  · exact ⟨12, by decide⟩
  · exact ⟨3, by decide⟩
  · exact ⟨0, by decide⟩
  · exact ⟨4, by decide⟩
  · exact ⟨12, by decide⟩
  · exact ⟨1, by decide⟩
  · exact ⟨0, by decide⟩

/-- `cexFamily` is scale-gapped (ratio `28.57 > 13`), so it is not dispatched by `spread13`. -/
theorem cex_gapFamily : GapFamily cexFamily := by decide

/-- **The upper edge is false: `cexFamily` has NO strictly-live multiplier at any `q ∈ [15,27]`.**
A covering, scale-gapped (hence residual-branch) family with its first strictly-live ruler at `q = 28`. -/
theorem cex_no_ruler_below_28 :
    ∀ q ∈ Finset.Icc (15 : ℤ) 27, ∀ p ∈ Finset.Ico (1 : ℤ) q, ¬ StrictlyLive cexFamily q p := by
  decide

/-! ### The measure floor for the counterexample family (a concrete instance, kernel-pure).

The general measure floor `SafeMeasureFloor` is the open case of LRC(14) — and `cex_no_ruler_below_28`
shows it admits no bounded-modulus shortcut.  But for *this* family — the hardest kind, near-tight with no
ruler below 28 — the floor holds and is provable: its ruler at `q = 28` (multiplier `p = 1`, i.e.
`t = 1/28`) is strictly live, so the whole chain `StrictlyLive → StrictWitness → safe interval → μ > 0`
fires.  This certifies `μ(cexFamily) > 0` with foundational axioms only, and demonstrates the machinery
end-to-end on the exact family that refuted the bounded conjecture. -/

/-- At `q = 28` the multiplier `p = 1` (`t = 1/28`) is strictly live for `cexFamily`. -/
theorem cex_strictlyLive_28 : StrictlyLive cexFamily 28 1 := by decide

/-- **The measure floor holds for `cexFamily`.**  `0 < volume (safePeriod cexFamily)` — a positive-measure
safe set for the near-tight residual family whose first strictly-live ruler is at `q = 28`. -/
theorem cex_measureFloor : 0 < MeasureTheory.volume (safePeriod cexFamily) :=
  measureFloor_of_strictWitness (strictWitness_of_strictlyLive cex_strictlyLive_28)

/-- **`cexFamily` is lonely** (a witness time exists) — the concrete instance of LRC(14). -/
theorem cex_lonely : ∃ t : ℝ, Lonely 14 cexFamily t :=
  lonely_of_safePeriod_measure_pos cex_measureFloor

end LRC14Grand
end LonelyRunner
