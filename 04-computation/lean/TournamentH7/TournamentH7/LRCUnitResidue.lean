/-
  TournamentH7.LRCUnitResidue  (mac-mini-2026-07-01-S95)

  The UNIT-RESIDUE LEMMA (THM-593 Part A + stability addendum), Lean-checked.

  Paper statement: if `M(S) = 1/q` exactly (a tight set) and `S` has no
  multiple of `q`, then for every `a` coprime to `q` the residue set
  `{a·v mod q : v ∈ S}` contains `+1` (and, applying it at `−a`, also `−1`).

  Lean form (the constructive strict-improvement half, which carries the whole
  mathematical content): if no speed is `≡ 0 (mod q)` and no speed has
  `a·v ≡ 1 (mod q)`, then the shifted witness `t = a/q − ε` with the explicit
  `ε = 1/(q(V+1))` (`V` = a bound on the speeds) is lonely with STRICT margin
  `1/q + ε`.  Hence a speed set missing a unit residue is not `1/q`-tight.
  Corollaries (paper, THM-593): tight sets represent every unit residue;
  dropped residues in the duplication+drop classification are non-units; at
  prime `q` tight residue systems are permutations of `{1,…,q−1}` — the
  structural explanation of the klein-S48 census zeros at n = 4, 6.

  Proof idea (the ε-shift refinement of `LonelyRunner.sieve_frac`): the
  numerator `D = v·a − m·q` has `D % q ∉ {0, 1}`, hence `D ∉ {0, 1}`; so
  either `D = −1` (distance `= 1/q + vε ≥ 1/q + ε`: the residue-(q−1) runner
  sits on the far side of the shift) or `|D| ≥ 2` (distance
  `≥ 2/q − Vε = 1/q + ε` by the choice of `ε`).

  Sorry-free (checked with `lake build TournamentH7.LRCUnitResidue`).
-/
import Mathlib.Tactic
import TournamentH7.LonelyRunner

namespace TournamentH7.LRCUnitResidue

open LonelyRunner

variable {ι : Type*}

/-- **The unit-residue strict-improvement lemma.**
If `q ≥ 3`, `a ⟂ q`, the speeds satisfy `1 ≤ v i ≤ V`, `q ∤ v i`, and the unit
residue `+1` is missed at `a` (`(v i * a) % q ≠ 1` for all `i`), then the time
`t = a/q − 1/(q(V+1))` keeps every runner at distance `≥ 1/q + 1/(q(V+1))`
from every integer — a lonely time with margin strictly better than `1/q`. -/
theorem unit_residue_improvement (q : ℕ) (a : ℤ) (V : ℤ) (v : ι → ℤ)
    (hq : 3 ≤ q) (hV : 1 ≤ V)
    (hpos : ∀ i, 1 ≤ v i) (hbd : ∀ i, v i ≤ V)
    (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i))
    (hcop : IsCoprime (q : ℤ) a)
    (hmiss : ∀ i, (v i * a) % q ≠ 1) :
    ∀ i, ∀ m : ℤ,
      1 / q + 1 / (q * (V + 1)) ≤
        |(v i : ℝ) * ((a : ℝ) / q - 1 / (q * (V + 1))) - m| := by
  intro i m
  have hq0 : (0 : ℝ) < (q : ℝ) := by
    have : 0 < q := by omega
    exact_mod_cast this
  have hqne : (q : ℝ) ≠ 0 := ne_of_gt hq0
  have hV1 : (0 : ℝ) < (V : ℝ) + 1 := by
    have : (1 : ℝ) ≤ (V : ℝ) := by exact_mod_cast hV
    linarith
  have hqVne : (q : ℝ) * ((V : ℝ) + 1) ≠ 0 := ne_of_gt (mul_pos hq0 hV1)
  have hv1 : (1 : ℝ) ≤ (v i : ℝ) := by exact_mod_cast hpos i
  have hv0 : (0 : ℝ) < (v i : ℝ) := lt_of_lt_of_le zero_lt_one hv1
  have hvV : (v i : ℝ) ≤ (V : ℝ) := by exact_mod_cast hbd i
  have hε0 : (0 : ℝ) < 1 / ((q : ℝ) * ((V : ℝ) + 1)) := by positivity
  have hvε0 : (0 : ℝ) < (v i : ℝ) * (1 / ((q : ℝ) * ((V : ℝ) + 1))) :=
    mul_pos hv0 hε0
  -- the numerator integer and the key algebraic identity
  set D : ℤ := v i * a - m * q with hD
  have key : (v i : ℝ) * ((a : ℝ) / q - 1 / (q * (V + 1))) - m
      = (D : ℝ) / q - (v i : ℝ) * (1 / (q * (V + 1))) := by
    rw [hD]
    push_cast
    field_simp
    ring
  -- residue bookkeeping:  D % q = (v i * a) % q ∉ {0, 1}, hence D ∉ {0, 1}
  have hmodD : D % q = (v i * a) % q := by
    rw [hD, Int.sub_emod, Int.mul_emod_left, sub_zero, Int.emod_emod_of_dvd _ dvd_rfl]
  have hnotdvd : ¬ ((q : ℤ) ∣ v i * a) := fun h => hdiv i (hcop.dvd_of_dvd_mul_right h)
  have hD0 : D ≠ 0 := by
    intro h
    apply hnotdvd
    apply Int.dvd_of_emod_eq_zero
    rw [← hmodD, h]
    simp
  have hD1 : D ≠ 1 := by
    intro h
    apply hmiss i
    rw [← hmodD, h]
    exact Int.emod_eq_of_lt (by omega) (by exact_mod_cast (by omega : 1 < q))
  rw [key]
  rcases eq_or_ne D (-1) with hDneg1 | hDneg1
  · -- CASE D = −1 : distance = 1/q + v·ε ≥ 1/q + ε
    rw [hDneg1]
    have hval : ((-1 : ℤ) : ℝ) / q - (v i : ℝ) * (1 / (q * (V + 1)))
        = -(1 / q + (v i : ℝ) * (1 / (q * (V + 1)))) := by
      push_cast
      ring
    rw [hval, abs_neg, abs_of_pos (by positivity)]
    have hstep : (1 : ℝ) * (1 / (q * (V + 1))) ≤ (v i : ℝ) * (1 / (q * (V + 1))) :=
      mul_le_mul_of_nonneg_right hv1 (le_of_lt hε0)
    linarith
  · -- CASE |D| ≥ 2 : distance ≥ 2/q − V·ε = 1/q + ε
    have habs2 : (2 : ℤ) ≤ |D| := by
      rcases lt_trichotomy D 0 with h | h | h
      · have h2 : D ≤ -2 := by omega
        rw [abs_of_neg h]; omega
      · exact absurd h hD0
      · have h2 : 2 ≤ D := by omega
        rw [abs_of_pos h]; omega
    have habs2R : (2 : ℝ) ≤ |(D : ℝ)| := by
      rw [← Int.cast_abs]
      exact_mod_cast habs2
    have hvle : (v i : ℝ) * (1 / (q * (V + 1))) ≤ (V : ℝ) * (1 / (q * (V + 1))) :=
      mul_le_mul_of_nonneg_right hvV (le_of_lt hε0)
    calc 1 / (q : ℝ) + 1 / (q * (V + 1))
        = 2 / q - (V : ℝ) * (1 / (q * (V + 1))) := by
          field_simp
          ring
      _ ≤ 2 / q - (v i : ℝ) * (1 / (q * (V + 1))) := by linarith
      _ ≤ |(D : ℝ)| / q - (v i : ℝ) * (1 / (q * (V + 1))) := by
          have h2q : (2 : ℝ) / q ≤ |(D : ℝ)| / q := by gcongr
          linarith
      _ = |(D : ℝ) / q| - |(v i : ℝ) * (1 / (q * (V + 1)))| := by
          rw [abs_div, abs_of_pos hq0, abs_of_pos hvε0]
      _ ≤ |(D : ℝ) / q - (v i : ℝ) * (1 / (q * (V + 1)))| :=
          abs_sub_abs_le_abs_sub _ _

/-- Packaging: under the same hypotheses the shifted witness is `Lonely` in the
`LonelyRunner` sense for the threshold `n = q`, with room to spare. -/
theorem lonely_of_missing_unit_residue (q : ℕ) (a : ℤ) (V : ℤ) (v : ι → ℤ)
    (hq : 3 ≤ q) (hV : 1 ≤ V)
    (hpos : ∀ i, 1 ≤ v i) (hbd : ∀ i, v i ≤ V)
    (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i))
    (hcop : IsCoprime (q : ℤ) a)
    (hmiss : ∀ i, (v i * a) % q ≠ 1) :
    Lonely q v ((a : ℝ) / q - 1 / (q * (V + 1))) := by
  intro i m
  have h := unit_residue_improvement q a V v hq hV hpos hbd hdiv hcop hmiss i m
  have hε0 : (0 : ℝ) < 1 / ((q : ℝ) * ((V : ℝ) + 1)) := by
    have hq0 : (0 : ℝ) < (q : ℝ) := by
      have : 0 < q := by omega
      exact_mod_cast this
    have : (1 : ℝ) ≤ (V : ℝ) := by exact_mod_cast hV
    positivity
  calc (1 : ℝ) / q ≤ 1 / q + 1 / (q * (V + 1)) := by linarith
    _ ≤ _ := h

end TournamentH7.LRCUnitResidue
