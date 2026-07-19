/-
  TournamentH7.LRCMod19Spread — the mod-19 antipodal-spread lemma (boxeph-2026-07-19-S126).

  Formalizes the necessary condition proved in S126 for any family in the gap regime `M(C) < 2/19`
  (which contains the depth-minimal gap value `3/38`, the edge `2/25`, and the floor `1/13`):

    if `M(C) < 2/19` then either `19 ∣` some speed, or the residues `{v mod 19}` cover ALL nine
    antipodal unit-pairs `{±1,…,±9}` of `ℤ/19`.

  The engine is the mod-19 middle-band sieve witness (`sieve19_middle_witness`): if at scale `b`
  every speed sits in the MIDDLE band `[2,17]` mod 19, then `t = b/19` puts every runner at
  circle-distance `≥ 2/19` from the integers, so `M(C) ≥ 2/19`.  Contrapositive: a family with
  `M(C) < 2/19` has, at every scale `b` with `19 ∤ b`, some UNIT speed hitting residue `±1`
  (`antipodal_spread`).  Since `b ↦ b⁻¹` is a bijection on the units of `ℤ/19`, this per-scale
  `±1`-hit is exactly the statement that the residues cover every antipodal unit-pair — the
  antipodal spread.

  Direct mod-19 analogue of `LRCMod13Blocking`; `2·19+... `— the modulus `19` is the apex of
  `38 = 2·19`, the natural alphabet of the `3/38` question (the `3/38` comb spectrum vanishes at
  `k ≡ 0 mod 19`, boxeph-S125).  Translation-sensitive, per opus's THM-1185/1220 triage.
-/
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- Integer core: a residue `r ∈ [2,17]` is at distance `≥ 2` from every multiple of 19.
`19k + r` never lands in `{-1,0,1}` when `2 ≤ r ≤ 17`. -/
theorem mod19_middle_far (r : ℤ) (h2 : 2 ≤ r) (h17 : r ≤ 17) (k : ℤ) :
    2 ≤ |19 * k + r| := by
  rcases le_or_gt 0 (19 * k + r) with h | h
  · rw [abs_of_nonneg h]; omega
  · rw [abs_of_neg h]; omega

/-- **Single-runner mod-19 middle-band witness.**  If `(v·b) mod 19 ∈ [2,17]`, then at `t = b/19`
the runner `v` is at distance `≥ 2/19` from every integer. -/
theorem sieve19_single (v b : ℤ)
    (h2 : 2 ≤ (v * b) % 19) (h17 : (v * b) % 19 ≤ 17) :
    ∀ m : ℤ, (2 : ℝ) / 19 ≤ |(v : ℝ) * ((b : ℝ) / 19) - m| := by
  intro m
  have hreal : (v : ℝ) * ((b : ℝ) / 19) - m
      = ((v * b - 19 * m : ℤ) : ℝ) / 19 := by push_cast; ring
  have hdecomp : v * b - 19 * m = 19 * (v * b / 19 - m) + (v * b) % 19 := by omega
  have hint : (2 : ℤ) ≤ |v * b - 19 * m| := by
    rw [hdecomp]; exact mod19_middle_far _ h2 h17 _
  have hintR : (2 : ℝ) ≤ |((v * b - 19 * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hint
  rw [hreal, abs_div, show |(19 : ℝ)| = 19 by norm_num]
  gcongr

/-- **The mod-19 middle-band sieve witness (family form).**  If at scale `b` every speed's residue
`(c_i·b) mod 19` lies in the middle band `[2,17]`, then at `t = b/19` every runner is at distance
`≥ 2/19` from every integer — a witness that `M(C) ≥ 2/19`. -/
theorem sieve19_middle_witness {ι : Type*} (c : ι → ℤ) (b : ℤ)
    (hmid : ∀ i, 2 ≤ (c i * b) % 19 ∧ (c i * b) % 19 ≤ 17) :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 19 ≤ |(c i : ℝ) * ((b : ℝ) / 19) - m| :=
  fun i => sieve19_single (c i) b (hmid i).1 (hmid i).2

/-- **Contrapositive — the middle band is emptied by a close runner.**  If at scale `b` some runner
is strictly within `2/19` of the integers, then not every residue can lie in the middle band. -/
theorem no_middle_band_of_close {ι : Type*} (c : ι → ℤ) (b : ℤ)
    (hclose : ∃ i, ∃ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 19) - m| < 2 / 19) :
    ¬ (∀ i, 2 ≤ (c i * b) % 19 ∧ (c i * b) % 19 ≤ 17) := by
  intro hmid
  obtain ⟨i, m, hi⟩ := hclose
  exact absurd (sieve19_middle_witness c b hmid i m) (not_le.mpr hi)

/-- **The mod-19 antipodal-spread lemma (per-scale form, PROVED).**  Suppose no speed is divisible
by 19, and at every scale `b` some runner is strictly within `2/19` of the integers (which holds
whenever `M(C) < 2/19`).  Then for every scale `b` with `19 ∤ b`, some runner has residue `±1`
mod 19:
    `∃ i, (c_i · b) mod 19 = 1 ∨ (c_i · b) mod 19 = 18`.
Via the bijection `b ↦ b⁻¹` on the units of `ℤ/19`, this says the residues `{c_i mod 19}` cover
every antipodal unit-pair `{±u}` — the antipodal spread. -/
theorem antipodal_spread {ι : Type*} (c : ι → ℤ)
    (hunit : ∀ i, ¬ ((19 : ℤ) ∣ c i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 19) - m| < 2 / 19)
    (b : ℤ) (hb : ¬ ((19 : ℤ) ∣ b)) :
    ∃ i, (c i * b) % 19 = 1 ∨ (c i * b) % 19 = 18 := by
  obtain ⟨i, mc, hi⟩ := hclose b
  refine ⟨i, ?_⟩
  by_contra hcon
  simp only [not_or] at hcon
  obtain ⟨hne1, hne18⟩ := hcon
  -- the residue is a genuine unit: not 0 (else 19 ∣ c i or 19 ∣ b, both excluded)
  have hne0 : (c i * b) % 19 ≠ 0 := by
    intro h0
    have hdvd : (19 : ℤ) ∣ c i * b := Int.dvd_of_emod_eq_zero h0
    have h19 : Prime (19 : ℤ) := by norm_num
    rcases h19.dvd_or_dvd hdvd with h | h
    · exact hunit i h
    · exact hb h
  -- so the residue lies in the middle band [2,17]
  have h0le : 0 ≤ (c i * b) % 19 := Int.emod_nonneg _ (by norm_num)
  have hlt19 : (c i * b) % 19 < 19 := Int.emod_lt_of_pos _ (by norm_num)
  have h2 : 2 ≤ (c i * b) % 19 := by omega
  have h17 : (c i * b) % 19 ≤ 17 := by omega
  -- but then the single-runner witness contradicts closeness at the witnessed integer `mc`
  exact absurd (sieve19_single (c i) b h2 h17 mc) (not_le.mpr hi)

/-- **The antipodal-covering corollary (in `ZMod 19`).**  Under the hypotheses of `antipodal_spread`,
the residues cover EVERY antipodal unit-pair of `ℤ/19`: for every nonzero `u : ZMod 19`, some runner
satisfies `c_i ≡ u` or `c_i ≡ -u (mod 19)`.  This is the `b ↦ b⁻¹` packaging of the per-scale
`±1`-hit — take `b = u⁻¹` and read the `±1` residue back as `±u`. -/
theorem antipodal_cover {ι : Type*} (c : ι → ℤ)
    (hunit : ∀ i, ¬ ((19 : ℤ) ∣ c i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 19) - m| < 2 / 19)
    (u : ZMod 19) (hu : u ≠ 0) :
    ∃ i, (c i : ZMod 19) = u ∨ (c i : ZMod 19) = -u := by
  haveI : Fact (Nat.Prime 19) := ⟨by norm_num⟩
  -- lift `u⁻¹` to an integer scale `bz`
  set bz : ℤ := ((u⁻¹).val : ℤ) with hbz
  have hvpos : 0 < (u⁻¹).val := ZMod.val_pos.mpr (inv_ne_zero hu)
  have hvlt : (u⁻¹).val < 19 := ZMod.val_lt _
  have hbne : ¬ ((19 : ℤ) ∣ bz) := by
    intro hdvd
    have hpos : (0 : ℤ) < bz := by rw [hbz]; exact_mod_cast hvpos
    have hle := Int.le_of_dvd hpos hdvd
    have hlt : bz < 19 := by rw [hbz]; exact_mod_cast hvlt
    omega
  have hbcast : (bz : ZMod 19) = u⁻¹ := by
    rw [hbz]; push_cast; exact ZMod.natCast_rightInverse u⁻¹
  obtain ⟨i, hi⟩ := antipodal_spread c hunit hclose bz hbne
  refine ⟨i, ?_⟩
  have hinv : u⁻¹ * u = 1 := inv_mul_cancel₀ hu
  rcases hi with h | h
  · left
    have hmod : (c i * bz) ≡ 1 [ZMOD (19 : ℤ)] := by
      show (c i * bz) % 19 = (1 : ℤ) % 19; omega
    have hz : ((c i * bz : ℤ) : ZMod 19) = ((1 : ℤ) : ZMod 19) :=
      (ZMod.intCast_eq_intCast_iff _ _ _).mpr hmod
    push_cast at hz
    rw [hbcast] at hz
    calc (c i : ZMod 19) = (c i : ZMod 19) * (u⁻¹ * u) := by rw [hinv, mul_one]
      _ = ((c i : ZMod 19) * u⁻¹) * u := by ring
      _ = 1 * u := by rw [hz]
      _ = u := one_mul u
  · right
    have hmod : (c i * bz) ≡ (-1) [ZMOD (19 : ℤ)] := by
      show (c i * bz) % 19 = (-1 : ℤ) % 19; omega
    have hz : ((c i * bz : ℤ) : ZMod 19) = ((-1 : ℤ) : ZMod 19) :=
      (ZMod.intCast_eq_intCast_iff _ _ _).mpr hmod
    push_cast at hz
    rw [hbcast] at hz
    calc (c i : ZMod 19) = (c i : ZMod 19) * (u⁻¹ * u) := by rw [hinv, mul_one]
      _ = ((c i : ZMod 19) * u⁻¹) * u := by ring
      _ = (-1) * u := by rw [hz]
      _ = -u := by ring

end LonelyRunner

#print axioms LonelyRunner.antipodal_cover
#print axioms LonelyRunner.mod19_middle_far
#print axioms LonelyRunner.sieve19_single
#print axioms LonelyRunner.sieve19_middle_witness
#print axioms LonelyRunner.no_middle_band_of_close
#print axioms LonelyRunner.antipodal_spread
