/-
  TournamentH7.LRCSieveDispatch — the non-covering ⟹ sieve dispatch of LRC(14).
  boxeph-2026-07-18-S106.

  A 13-family `V` is COVERING (divisibility-sieve notion) if every modulus
  `2 ≤ n ≤ 14` divides some speed.  boxeph-S101 proved `M<1/13 ⟹` every `n ∈ {2..13}`
  divides some speed; the `≤ 14` form matches the `1/14` threshold.

  If `V` is NOT covering, some `n ∈ {2..14}` divides no speed, so `t = 1/n` is
  `n`-lonely (`lonely_of_no_multiple`), hence `1/14`-lonely (`n ≤ 14`).  This is the
  LRC(≤13) sieve witness for the non-covering class.

  `sieve_dispatch`   (PROVED, kernel-pure): `¬ Covering v ⟹ ∃ t, Lonely 14 v t`.
  `lonely14_dispatch`(PROVED): the full dichotomy — LRC(14) for `v` follows from the
  COVERING case alone; the non-covering case is discharged here by the sieve.

  With `LRCAPCoreBridge` (S105) the covering case reduces (via `INV` + LRC(≤13)) to
  `ap_core_bridge`, so:
      LRC(14)  ⟸  LRC(≤13) [cited]  +  INV [open]  +  {sieve dispatch, descent} [Lean].
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- **Covering (divisibility-sieve notion).** Every modulus `2 ≤ n ≤ 14` divides some
speed of the 13-family. -/
def Covering (v : Fin 13 → ℤ) : Prop :=
  ∀ n : ℕ, 2 ≤ n → n ≤ 14 → ∃ i, (n : ℤ) ∣ v i

/-- **Lonely-constant monotonicity.** A time lonely at gap `1/n` with `0 < n ≤ 14` is
lonely at `1/14` (the band only shrinks). -/
theorem lonely14_of_lonely_le {v : Fin 13 → ℤ} {t : ℝ} {n : ℕ}
    (hn0 : 0 < n) (hn : n ≤ 14) (h : Lonely n v t) : Lonely 14 v t := by
  intro i m
  have hpos : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn0
  have hle : (1 : ℝ) / 14 ≤ 1 / (n : ℝ) :=
    one_div_le_one_div_of_le hpos (by exact_mod_cast hn)
  exact le_trans hle (h i m)

/-- **The non-covering ⟹ sieve dispatch (PROVED).** If `v` is not covering, some
`n ∈ {2..14}` divides no speed; then `t = 1/n` is `1/14`-lonely by the empty-circle
sieve `lonely_of_no_multiple`. -/
theorem sieve_dispatch (v : Fin 13 → ℤ) (hnc : ¬ Covering v) :
    ∃ t : ℝ, Lonely 14 v t := by
  simp only [Covering] at hnc
  push_neg at hnc
  obtain ⟨n, hn2, hn14, hdiv⟩ := hnc
  exact ⟨1 / n,
    lonely14_of_lonely_le (by omega) hn14 (lonely_of_no_multiple n (by omega) v hdiv)⟩

/-- **The dispatch dichotomy (PROVED).** LRC(14) for `v` follows from the COVERING
case alone: if `v` is covering, apply the covering hypothesis; otherwise the sieve. -/
theorem lonely14_dispatch (v : Fin 13 → ℤ)
    (covering_case : Covering v → ∃ t : ℝ, Lonely 14 v t) :
    ∃ t : ℝ, Lonely 14 v t := by
  by_cases hc : Covering v
  · exact covering_case hc
  · exact sieve_dispatch v hc

/-- **The covering case of LRC(14) — the open crux.** Every covering 13-family is
`1/14`-lonely.  (For `M ≥ 1/14` families this is immediate; the hard sub-case `M < 1/14`
needs the inverse theorem — boxeph-S94: `M<1/13 ⟹ ρ ≥ 13`, then `ap_core_bridge` (S105).
Entered as a NAMED HYPOTHESIS, never a `sorry`.  NB: this is the honest covering crux,
*not* the too-strong "covering ⟹ ρ≥13" — e.g. `{2,…,14}` is covering yet has ρ < 13,
but is trivially `1/14`-lonely.) -/
def CoveringCase : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → Covering v → ∃ t : ℝ, Lonely 14 v t

/-- **LRC(14) reduces to the covering case (PROVED).**  Given the covering crux, every
13-family of positive speeds is `1/14`-lonely: the covering families by the crux, the
non-covering families by `sieve_dispatch`.  This records
`LRC(14) ⟸ CoveringCase + {sieve dispatch} [Lean]`, and (with S105)
`CoveringCase ⟸ LRC(≤13) + INV + descent`. -/
theorem lrc14_of_covering (hcov : CoveringCase)
    (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) : ∃ t : ℝ, Lonely 14 v t :=
  lonely14_dispatch v (hcov v hpos)

end LonelyRunner

#print axioms LonelyRunner.lonely14_of_lonely_le
#print axioms LonelyRunner.lrc14_of_covering
#print axioms LonelyRunner.sieve_dispatch
#print axioms LonelyRunner.lonely14_dispatch
