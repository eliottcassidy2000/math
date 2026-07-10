/-
  TournamentH7.LRCLadderFattening  (boxeph-2026-07-09-S4)

  THM-667 Lemma A (klein-S208), the Lipschitz fattening of the LRC(≤13)
  citation witness — interval form.

  The adaptive-split ladder needs: any band `T` of ≤ 12 distinct nonzero
  speeds has a SAFE SET `G_T = {t : ∀ v ∈ T, ‖v·t‖ ≥ 1/14}` of positive
  measure, quantitatively `≥ 1/(91·maxT)`.  Proof shape (klein-S208, "3
  lines given the citation"): the LRC(≤13) citation gives a witness `t₀`
  clearing every band member at `1/(k+1) ≥ 1/13`; each residue `|w·t − m|`
  is `|w|`-Lipschitz in `t`, so every `t` within `1/(182·B)` of `t₀`
  (`B ≥ max |w|`) still clears at `1/13 − 1/182 = 1/14` exactly.  The safe
  set therefore CONTAINS the interval `[t₀ − 1/(182B), t₀ + 1/(182B)]` of
  length `1/(91B)` — the measure floor follows by monotonicity.

  Two theorems, both elementary:
    * `lonely_band_transport` — the Lipschitz transport of a `Lonely n`
      band to a nearby time at an explicitly degraded clearance.
    * `safe_interval_of_citation` — the citation ⟹ an explicit closed
      interval of `Lonely 14` times for any ≤ 12-speed band.

  Precedent mechanism: `LRCIntervalTransport.lonely_interval_of_margin`
  (same Lipschitz transport, typed at repeat 13-tuples for the coarse
  induction). This file states the BAND-LEVEL version the adaptive-split
  ladder (THM-667) consumes: an arbitrary stratum T with |T| = k <= 12,
  i.e. the G_T floor for every rung, not the repeat-tuple corner.

  Kernel-pure: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

namespace LonelyRunner

/-- **Lipschitz transport of a lonely band.**  If every runner of `w` clears
`1/n` at `t₀` and every speed is bounded by `B`, then at any nearby time `t`
every runner still clears `1/n − B·|t − t₀|`: each residue `|w·t − m|` moves
at rate `≤ |w| ≤ B`.  Pure triangle inequality. -/
theorem lonely_band_transport {ι : Type*} (n : ℕ) (w : ι → ℤ) (t0 t B c : ℝ)
    (hlon : Lonely n w t0) (hB : ∀ i, |(w i : ℝ)| ≤ B)
    (hc : c ≤ 1 / n - B * |t - t0|) :
    ∀ i, ∀ m : ℤ, c ≤ |(w i : ℝ) * t - m| := by
  intro i m
  have h0 : (1 : ℝ) / n ≤ |(w i : ℝ) * t0 - m| := hlon i m
  -- |w t − m| ≥ |w t₀ − m| − |w t₀ − w t|
  have htri : |(w i : ℝ) * t0 - m| - |(w i : ℝ) * t0 - (w i : ℝ) * t|
      ≤ |(w i : ℝ) * t - m| := by
    have h := abs_sub_abs_le_abs_sub ((w i : ℝ) * t0 - m) ((w i : ℝ) * t - m)
    have hrw : ((w i : ℝ) * t0 - m) - ((w i : ℝ) * t - m)
        = (w i : ℝ) * t0 - (w i : ℝ) * t := by ring
    rw [hrw] at h
    linarith
  -- the drift term: |w t₀ − w t| = |w|·|t₀ − t| ≤ B·|t − t₀|
  have hdrift : |(w i : ℝ) * t0 - (w i : ℝ) * t| ≤ B * |t - t0| := by
    have h1 : (w i : ℝ) * t0 - (w i : ℝ) * t = (w i : ℝ) * (t0 - t) := by ring
    rw [h1, abs_mul]
    have h2 : |t0 - t| = |t - t0| := abs_sub_comm t0 t
    rw [h2]
    exact mul_le_mul_of_nonneg_right (hB i) (abs_nonneg _)
  linarith

/-- **THM-667 Lemma A, interval form.**  Any band of `k ≤ 12` nonzero integer
speeds bounded by `B > 0` admits an explicit closed interval of length
`1/(91·B)` on which EVERY band member clears the LRC(14) band `1/14`:
the LRC(≤13) citation witness, fattened by `1/(182·B)` on each side
(`1/13 − 1/182 = 1/14` exactly).  The ladder's measure floor
`meas(G_T) ≥ 1/(91·maxT)` follows by measure monotonicity on this interval. -/
theorem safe_interval_of_citation (cite : LRCUpTo13) (k : ℕ) (hk : k ≤ 12)
    (w : Fin k → ℤ) (hw : ∀ i, w i ≠ 0) (B : ℝ)
    (hB : ∀ i, |(w i : ℝ)| ≤ B) (hB0 : 0 < B) :
    ∃ t0 : ℝ, ∀ t ∈ Set.Icc (t0 - 1 / (182 * B)) (t0 + 1 / (182 * B)),
      Lonely 14 w t := by
  obtain ⟨t0, hlon⟩ := cite k hk w hw
  refine ⟨t0, fun t ht => ?_⟩
  have habs : |t - t0| ≤ 1 / (182 * B) := by
    rcases ht with ⟨h1, h2⟩
    rw [abs_le]
    constructor <;> linarith
  intro i m
  refine lonely_band_transport (k + 1) w t0 t B (1 / 14) hlon hB ?_ i m
  -- 1/14 ≤ 1/(k+1) − B·|t − t₀|, since 1/(k+1) ≥ 1/13 and B·|t−t₀| ≤ 1/182
  have hk13 : ((k : ℝ) + 1) ≤ 13 := by
    have : (k : ℝ) ≤ 12 := by exact_mod_cast hk
    linarith
  have hkpos : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have h13 : (1 : ℝ) / 13 ≤ 1 / ((k : ℝ) + 1) :=
    one_div_le_one_div_of_le hkpos hk13
  have hcast : (1 : ℝ) / ((k + 1 : ℕ) : ℝ) = 1 / ((k : ℝ) + 1) := by
    push_cast
    ring_nf
  have hdrift : B * |t - t0| ≤ 1 / 182 := by
    have h1 : B * |t - t0| ≤ B * (1 / (182 * B)) :=
      mul_le_mul_of_nonneg_left habs (le_of_lt hB0)
    have h2 : B * (1 / (182 * B)) = 1 / 182 := by
      field_simp
    linarith
  rw [hcast]
  linarith

end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.lonely_band_transport
#print axioms LonelyRunner.safe_interval_of_citation
