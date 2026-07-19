/-
  TournamentH7.LRCMod23Spread — the mod-23 antipodal-spread lemma (boxeph-2026-07-19-S132).

  Mod-23 analogue of LRCMod19Spread (same one-line engine, prime 23): for any family in
  `M(C) < 2/23` (contains the whole n=12 uniqueness gap `(1/13, 2/25)` and the live n=13
  window `(1/14, 3/41)`), either `23 ∣` some speed or the residues cover ALL ELEVEN
  antipodal unit-pairs of `ℤ/23`.  Load-bearing since boxeph-S131/S132 (HYP-7895/7900):
  the no-23-blocker coincidence costs 1-in-1191 on uniform slots, and the S132 exhaustive
  in-band census shows the 23-rung kills 80% of the S126 mod-19 evaders; with 12 slots the
  covering is a NEAR-BIJECTION (slack 1).  Original 19-header follows.

  Formalizes the necessary condition proved in S126 for any family in the gap regime `M(C) < 2/23`
  (which contains the depth-minimal gap value `3/38`, the edge `2/25`, and the floor `1/13`):

    if `M(C) < 2/23` then either `23 ∣` some speed, or the residues `{v mod 23}` cover ALL nine
    antipodal unit-pairs `{±1,…,±11}` of `ℤ/23`.

  The engine is the mod-23 middle-band sieve witness (`sieve23_middle_witness`): if at scale `b`
  every speed sits in the MIDDLE band `[2,21]` mod 23, then `t = b/23` puts every runner at
  circle-distance `≥ 2/23` from the integers, so `M(C) ≥ 2/23`.  Contrapositive: a family with
  `M(C) < 2/23` has, at every scale `b` with `23 ∤ b`, some UNIT speed hitting residue `±1`
  (`antipodal_spread23`).  Since `b ↦ b⁻¹` is a bijection on the units of `ℤ/23`, this per-scale
  `±1`-hit is exactly the statement that the residues cover every antipodal unit-pair — the
  antipodal spread.

  Direct mod-19 analogue of `LRCMod13Blocking`; `2·19+... `— the modulus `23` is the apex of
  `38 = 2·19`, the natural alphabet of the `3/38` question (the `3/38` comb spectrum vanishes at
  `k ≡ 0 mod 23`, boxeph-S125).  Translation-sensitive, per opus's THM-1185/1220 triage.
-/
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- Integer core: a residue `r ∈ [2,21]` is at distance `≥ 2` from every multiple of 23.
`23k + r` never lands in `{-1,0,1}` when `2 ≤ r ≤ 21`. -/
theorem mod23_middle_far (r : ℤ) (h2 : 2 ≤ r) (h21 : r ≤ 21) (k : ℤ) :
    2 ≤ |23 * k + r| := by
  rcases le_or_gt 0 (23 * k + r) with h | h
  · rw [abs_of_nonneg h]; omega
  · rw [abs_of_neg h]; omega

/-- **Single-runner mod-23 middle-band witness.**  If `(v·b) mod 23 ∈ [2,21]`, then at `t = b/23`
the runner `v` is at distance `≥ 2/23` from every integer. -/
theorem sieve23_single (v b : ℤ)
    (h2 : 2 ≤ (v * b) % 23) (h21 : (v * b) % 23 ≤ 21) :
    ∀ m : ℤ, (2 : ℝ) / 23 ≤ |(v : ℝ) * ((b : ℝ) / 23) - m| := by
  intro m
  have hreal : (v : ℝ) * ((b : ℝ) / 23) - m
      = ((v * b - 23 * m : ℤ) : ℝ) / 23 := by push_cast; ring
  have hdecomp : v * b - 23 * m = 23 * (v * b / 23 - m) + (v * b) % 23 := by omega
  have hint : (2 : ℤ) ≤ |v * b - 23 * m| := by
    rw [hdecomp]; exact mod23_middle_far _ h2 h21 _
  have hintR : (2 : ℝ) ≤ |((v * b - 23 * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hint
  rw [hreal, abs_div, show |(23 : ℝ)| = 23 by norm_num]
  gcongr

/-- **The mod-23 middle-band sieve witness (family form).**  If at scale `b` every speed's residue
`(c_i·b) mod 23` lies in the middle band `[2,21]`, then at `t = b/23` every runner is at distance
`≥ 2/23` from every integer — a witness that `M(C) ≥ 2/23`. -/
theorem sieve23_middle_witness {ι : Type*} (c : ι → ℤ) (b : ℤ)
    (hmid : ∀ i, 2 ≤ (c i * b) % 23 ∧ (c i * b) % 23 ≤ 21) :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 23 ≤ |(c i : ℝ) * ((b : ℝ) / 23) - m| :=
  fun i => sieve23_single (c i) b (hmid i).1 (hmid i).2

/-- **Contrapositive — the middle band is emptied by a close runner.**  If at scale `b` some runner
is strictly within `2/23` of the integers, then not every residue can lie in the middle band. -/
theorem no_middle_band_of_close23 {ι : Type*} (c : ι → ℤ) (b : ℤ)
    (hclose : ∃ i, ∃ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23) :
    ¬ (∀ i, 2 ≤ (c i * b) % 23 ∧ (c i * b) % 23 ≤ 21) := by
  intro hmid
  obtain ⟨i, m, hi⟩ := hclose
  exact absurd (sieve23_middle_witness c b hmid i m) (not_le.mpr hi)

/-- **The mod-23 antipodal-spread lemma (per-scale form, PROVED).**  Suppose no speed is divisible
by 23, and at every scale `b` some runner is strictly within `2/23` of the integers (which holds
whenever `M(C) < 2/23`).  Then for every scale `b` with `23 ∤ b`, some runner has residue `±1`
mod 23:
    `∃ i, (c_i · b) mod 23 = 1 ∨ (c_i · b) mod 23 = 22`.
Via the bijection `b ↦ b⁻¹` on the units of `ℤ/23`, this says the residues `{c_i mod 23}` cover
every antipodal unit-pair `{±u}` — the antipodal spread. -/
theorem antipodal_spread23 {ι : Type*} (c : ι → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ c i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23)
    (b : ℤ) (hb : ¬ ((23 : ℤ) ∣ b)) :
    ∃ i, (c i * b) % 23 = 1 ∨ (c i * b) % 23 = 22 := by
  obtain ⟨i, mc, hi⟩ := hclose b
  refine ⟨i, ?_⟩
  by_contra hcon
  simp only [not_or] at hcon
  obtain ⟨hne1, hne22⟩ := hcon
  -- the residue is a genuine unit: not 0 (else 23 ∣ c i or 23 ∣ b, both excluded)
  have hne0 : (c i * b) % 23 ≠ 0 := by
    intro h0
    have hdvd : (23 : ℤ) ∣ c i * b := Int.dvd_of_emod_eq_zero h0
    have h19 : Prime (23 : ℤ) := by norm_num
    rcases h19.dvd_or_dvd hdvd with h | h
    · exact hunit i h
    · exact hb h
  -- so the residue lies in the middle band [2,21]
  have h0le : 0 ≤ (c i * b) % 23 := Int.emod_nonneg _ (by norm_num)
  have hlt19 : (c i * b) % 23 < 23 := Int.emod_lt_of_pos _ (by norm_num)
  have h2 : 2 ≤ (c i * b) % 23 := by omega
  have h21 : (c i * b) % 23 ≤ 21 := by omega
  -- but then the single-runner witness contradicts closeness at the witnessed integer `mc`
  exact absurd (sieve23_single (c i) b h2 h21 mc) (not_le.mpr hi)

/-- **The antipodal-covering corollary (in `ZMod 23`).**  Under the hypotheses of `antipodal_spread23`,
the residues cover EVERY antipodal unit-pair of `ℤ/23`: for every nonzero `u : ZMod 23`, some runner
satisfies `c_i ≡ u` or `c_i ≡ -u (mod 23)`.  This is the `b ↦ b⁻¹` packaging of the per-scale
`±1`-hit — take `b = u⁻¹` and read the `±1` residue back as `±u`. -/
theorem antipodal_cover23 {ι : Type*} (c : ι → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ c i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23)
    (u : ZMod 23) (hu : u ≠ 0) :
    ∃ i, (c i : ZMod 23) = u ∨ (c i : ZMod 23) = -u := by
  haveI : Fact (Nat.Prime 23) := ⟨by norm_num⟩
  -- lift `u⁻¹` to an integer scale `bz`
  set bz : ℤ := ((u⁻¹).val : ℤ) with hbz
  have hvpos : 0 < (u⁻¹).val := ZMod.val_pos.mpr (inv_ne_zero hu)
  have hvlt : (u⁻¹).val < 23 := ZMod.val_lt _
  have hbne : ¬ ((23 : ℤ) ∣ bz) := by
    intro hdvd
    have hpos : (0 : ℤ) < bz := by rw [hbz]; exact_mod_cast hvpos
    have hle := Int.le_of_dvd hpos hdvd
    have hlt : bz < 23 := by rw [hbz]; exact_mod_cast hvlt
    omega
  have hbcast : (bz : ZMod 23) = u⁻¹ := by
    rw [hbz]; push_cast; exact ZMod.natCast_rightInverse u⁻¹
  obtain ⟨i, hi⟩ := antipodal_spread23 c hunit hclose bz hbne
  refine ⟨i, ?_⟩
  have hinv : u⁻¹ * u = 1 := inv_mul_cancel₀ hu
  rcases hi with h | h
  · left
    have hmod : (c i * bz) ≡ 1 [ZMOD (23 : ℤ)] := by
      show (c i * bz) % 23 = (1 : ℤ) % 23; omega
    have hz : ((c i * bz : ℤ) : ZMod 23) = ((1 : ℤ) : ZMod 23) :=
      (ZMod.intCast_eq_intCast_iff _ _ _).mpr hmod
    push_cast at hz
    rw [hbcast] at hz
    calc (c i : ZMod 23) = (c i : ZMod 23) * (u⁻¹ * u) := by rw [hinv, mul_one]
      _ = ((c i : ZMod 23) * u⁻¹) * u := by ring
      _ = 1 * u := by rw [hz]
      _ = u := one_mul u
  · right
    have hmod : (c i * bz) ≡ (-1) [ZMOD (23 : ℤ)] := by
      show (c i * bz) % 23 = (-1 : ℤ) % 23; omega
    have hz : ((c i * bz : ℤ) : ZMod 23) = ((-1 : ℤ) : ZMod 23) :=
      (ZMod.intCast_eq_intCast_iff _ _ _).mpr hmod
    push_cast at hz
    rw [hbcast] at hz
    calc (c i : ZMod 23) = (c i : ZMod 23) * (u⁻¹ * u) := by rw [hinv, mul_one]
      _ = ((c i : ZMod 23) * u⁻¹) * u := by ring
      _ = (-1) * u := by rw [hz]
      _ = -u := by ring

end LonelyRunner

#print axioms LonelyRunner.antipodal_cover23
#print axioms LonelyRunner.mod23_middle_far
#print axioms LonelyRunner.sieve23_single
#print axioms LonelyRunner.sieve23_middle_witness
#print axioms LonelyRunner.no_middle_band_of_close23
#print axioms LonelyRunner.antipodal_spread23
