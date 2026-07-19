import TournamentH7.LRCWitnessAttainment

/-!
# The mod-23 near-bijection rung for twelve speeds (THM-1263)

If twelve integer speeds avoid the blocking residue `0 mod 23` and the
loneliness margin is strictly below `2/23` at every time `b/23`, then their
nonzero residues cover all eleven antipodal pairs of `ZMod 23`.  Since twelve
runners cover eleven pairs, exactly one pair occurs twice and every other pair
occurs once.

The closeness hypothesis deliberately has shape `∃ i, ∃ m`, not the vacuous
`∃ i, ∀ m` shape corrected in MISTAKE-186.
-/

open TournamentH7.LRCWitness

namespace LRC14.Mod23NearBijection

theorem threshold_chain :
    (1 : ℝ) / 13 < 2 / 25 ∧ (2 : ℝ) / 25 < 2 / 23 := by
  norm_num

/-- Integer-periodicity in the form needed to reduce arbitrary `b/23` into
the compact interval `[0,1]`. -/
theorem margin_add_int (v : Fin 12 → ℤ) (t : ℝ) (q : ℤ) :
    margin v (t + (q : ℝ)) = margin v t := by
  unfold margin
  have hfun :
      (fun i => distZ ((v i : ℝ) * (t + (q : ℝ)))) =
        (fun i => distZ ((v i : ℝ) * t)) := by
    funext i
    rw [show (v i : ℝ) * (t + (q : ℝ)) =
        (v i : ℝ) * t + ((v i * q : ℤ) : ℝ) by push_cast; ring,
      distZ_add_int]
  rw [hfun]

/-- Integer core for the mod-23 middle band. -/
theorem mod23_middle_far (r : ℤ) (h2 : 2 ≤ r) (h21 : r ≤ 21) (k : ℤ) :
    2 ≤ |23 * k + r| := by
  rcases le_or_gt 0 (23 * k + r) with h | h
  · rw [abs_of_nonneg h]
    omega
  · rw [abs_of_neg h]
    omega

/-- A middle-band residue is at circle distance at least `2/23`. -/
theorem sieve23_single (v b : ℤ)
    (h2 : 2 ≤ (v * b) % 23) (h21 : (v * b) % 23 ≤ 21) :
    ∀ m : ℤ, (2 : ℝ) / 23 ≤ |(v : ℝ) * ((b : ℝ) / 23) - m| := by
  intro m
  have hreal : (v : ℝ) * ((b : ℝ) / 23) - m =
      ((v * b - 23 * m : ℤ) : ℝ) / 23 := by
    push_cast
    ring
  have hdecomp : v * b - 23 * m = 23 * (v * b / 23 - m) + (v * b) % 23 := by
    omega
  have hint : (2 : ℤ) ≤ |v * b - 23 * m| := by
    rw [hdecomp]
    exact mod23_middle_far _ h2 h21 _
  have hintR : (2 : ℝ) ≤ |((v * b - 23 * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]
    exact_mod_cast hint
  rw [hreal, abs_div, show |(23 : ℝ)| = 23 by norm_num]
  gcongr

/-- Correct contrapositive shape: one runner is close to one witnessed integer. -/
theorem no_middle_band_of_close {ι : Type*} (v : ι → ℤ) (b : ℤ)
    (hclose : ∃ i, ∃ m : ℤ,
      |(v i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23) :
    ¬ (∀ i, 2 ≤ (v i * b) % 23 ∧ (v i * b) % 23 ≤ 21) := by
  intro hmid
  obtain ⟨i, m, hi⟩ := hclose
  exact absurd (sieve23_single (v i) b (hmid i).1 (hmid i).2 m)
    (not_le.mpr hi)

/-- Per-scale form: every unit scale sees a residue `±1`. -/
theorem antipodal_spread {ι : Type*} (v : ι → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ v i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ,
      |(v i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23)
    (b : ℤ) (hb : ¬ ((23 : ℤ) ∣ b)) :
    ∃ i, (v i * b) % 23 = 1 ∨ (v i * b) % 23 = 22 := by
  obtain ⟨i, m, hi⟩ := hclose b
  refine ⟨i, ?_⟩
  by_contra hcon
  simp only [not_or] at hcon
  obtain ⟨hne1, hne22⟩ := hcon
  have hne0 : (v i * b) % 23 ≠ 0 := by
    intro h0
    have hdvd : (23 : ℤ) ∣ v i * b := Int.dvd_of_emod_eq_zero h0
    have h23 : Prime (23 : ℤ) := by norm_num
    rcases h23.dvd_or_dvd hdvd with h | h
    · exact hunit i h
    · exact hb h
  have h0le : 0 ≤ (v i * b) % 23 := Int.emod_nonneg _ (by norm_num)
  have hlt23 : (v i * b) % 23 < 23 := Int.emod_lt_of_pos _ (by norm_num)
  have h2 : 2 ≤ (v i * b) % 23 := by omega
  have h21 : (v i * b) % 23 ≤ 21 := by omega
  exact absurd (sieve23_single (v i) b h2 h21 m) (not_le.mpr hi)

/-- Every nonzero antipodal pair is covered. -/
theorem antipodal_cover {ι : Type*} (v : ι → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ v i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ,
      |(v i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23)
    (u : ZMod 23) (hu : u ≠ 0) :
    ∃ i, (v i : ZMod 23) = u ∨ (v i : ZMod 23) = -u := by
  haveI : Fact (Nat.Prime 23) := ⟨by norm_num⟩
  set bz : ℤ := ((u⁻¹).val : ℤ) with hbz
  have hvpos : 0 < (u⁻¹).val := ZMod.val_pos.mpr (inv_ne_zero hu)
  have hvlt : (u⁻¹).val < 23 := ZMod.val_lt _
  have hbne : ¬ ((23 : ℤ) ∣ bz) := by
    intro hdvd
    have hpos : (0 : ℤ) < bz := by
      rw [hbz]
      exact_mod_cast hvpos
    have hle := Int.le_of_dvd hpos hdvd
    have hlt : bz < 23 := by
      rw [hbz]
      exact_mod_cast hvlt
    omega
  have hbcast : (bz : ZMod 23) = u⁻¹ := by
    rw [hbz]
    push_cast
    exact ZMod.natCast_rightInverse u⁻¹
  obtain ⟨i, hi⟩ := antipodal_spread v hunit hclose bz hbne
  refine ⟨i, ?_⟩
  have hinv : u⁻¹ * u = 1 := inv_mul_cancel₀ hu
  rcases hi with hi | hi
  · left
    have hmod : (v i * bz) ≡ 1 [ZMOD (23 : ℤ)] := by
      show (v i * bz) % 23 = (1 : ℤ) % 23
      omega
    have hz : ((v i * bz : ℤ) : ZMod 23) = ((1 : ℤ) : ZMod 23) :=
      (ZMod.intCast_eq_intCast_iff _ _ _).mpr hmod
    push_cast at hz
    rw [hbcast] at hz
    calc
      (v i : ZMod 23) = (v i : ZMod 23) * (u⁻¹ * u) := by rw [hinv, mul_one]
      _ = ((v i : ZMod 23) * u⁻¹) * u := by ring
      _ = 1 * u := by rw [hz]
      _ = u := one_mul u
  · right
    have hmod : (v i * bz) ≡ (-1) [ZMOD (23 : ℤ)] := by
      show (v i * bz) % 23 = (-1 : ℤ) % 23
      omega
    have hz : ((v i * bz : ℤ) : ZMod 23) = ((-1 : ℤ) : ZMod 23) :=
      (ZMod.intCast_eq_intCast_iff _ _ _).mpr hmod
    push_cast at hz
    rw [hbcast] at hz
    calc
      (v i : ZMod 23) = (v i : ZMod 23) * (u⁻¹ * u) := by rw [hinv, mul_one]
      _ = ((v i : ZMod 23) * u⁻¹) * u := by ring
      _ = (-1) * u := by rw [hz]
      _ = -u := by ring

/-- The canonical representative `1,...,11` of an antipodal pair. -/
def antipodalRep (r : Fin 11) : ZMod 23 := ((r.val + 1 : ℕ) : ZMod 23)

theorem val_antipodalRep (r : Fin 11) :
    (antipodalRep r).val = r.val + 1 := by
  rw [antipodalRep, ZMod.val_natCast]
  exact Nat.mod_eq_of_lt (by omega)

theorem antipodalRep_ne_zero (r : Fin 11) : antipodalRep r ≠ 0 := by
  intro h
  have hval : (antipodalRep r).val = 0 := (ZMod.val_eq_zero _).mpr h
  rw [val_antipodalRep] at hval
  omega

theorem val_neg_antipodalRep (r : Fin 11) :
    (-antipodalRep r).val = 22 - r.val := by
  rw [ZMod.neg_val, if_neg (antipodalRep_ne_zero r), val_antipodalRep]
  omega

/-- Every nonzero residue belongs to a unique canonical antipodal pair. -/
theorem existsUnique_antipodalRep (x : ZMod 23) (hx : x ≠ 0) :
    ∃! r : Fin 11, x = antipodalRep r ∨ x = -antipodalRep r := by
  have hxpos : 0 < x.val := ZMod.val_pos.mpr hx
  have hxlt : x.val < 23 := ZMod.val_lt x
  by_cases hlo : x.val ≤ 11
  · let r : Fin 11 := ⟨x.val - 1, by omega⟩
    refine ⟨r, Or.inl ?_, ?_⟩
    · apply ZMod.val_injective 23
      rw [val_antipodalRep]
      simp only [r]
      omega
    · intro s hs
      apply Fin.ext
      rcases hs with hs | hs
      · have hval := congrArg ZMod.val hs
        rw [val_antipodalRep] at hval
        simp only [r]
        omega
      · have hval := congrArg ZMod.val hs
        rw [val_neg_antipodalRep] at hval
        simp only [r]
        omega
  · let r : Fin 11 := ⟨22 - x.val, by omega⟩
    refine ⟨r, Or.inr ?_, ?_⟩
    · apply ZMod.val_injective 23
      rw [val_neg_antipodalRep]
      simp only [r]
      omega
    · intro s hs
      apply Fin.ext
      rcases hs with hs | hs
      · have hval := congrArg ZMod.val hs
        rw [val_antipodalRep] at hval
        simp only [r]
        omega
      · have hval := congrArg ZMod.val hs
        rw [val_neg_antipodalRep] at hval
        simp only [r]
        omega

/-- Cardinality of a fiber of a pair-labeling map. -/
def fiberCard (pair : Fin 12 → Fin 11) (r : Fin 11) : ℕ :=
  (Finset.univ.filter fun i => pair i = r).card

theorem fiberCard_pos_of_surjective (pair : Fin 12 → Fin 11)
    (hpair : Function.Surjective pair) (r : Fin 11) :
    1 ≤ fiberCard pair r := by
  obtain ⟨i, hi⟩ := hpair r
  rw [fiberCard, Finset.one_le_card]
  exact ⟨i, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hi⟩⟩

theorem sum_fiberCard (pair : Fin 12 → Fin 11) :
    ∑ r : Fin 11, fiberCard pair r = 12 := by
  have h := Finset.card_eq_sum_card_fiberwise
    (f := pair) (s := (Finset.univ : Finset (Fin 12)))
    (t := (Finset.univ : Finset (Fin 11)))
    (fun i _ => Finset.mem_univ (pair i))
  simpa [fiberCard] using h.symm

/-- A surjection from twelve objects to eleven labels has exactly one double
fiber and ten singleton fibers. -/
theorem nearBijection_of_surjective (pair : Fin 12 → Fin 11)
    (hpair : Function.Surjective pair) :
    ∃! r : Fin 11,
      fiberCard pair r = 2 ∧
      ∀ s : Fin 11, s ≠ r → fiberCard pair s = 1 := by
  let card : Fin 11 → ℕ := fiberCard pair
  have hpos : ∀ r, 1 ≤ card r := fun r => fiberCard_pos_of_surjective pair hpair r
  have hsum : ∑ r, card r = 12 := by simpa [card] using sum_fiberCard pair
  have hex : ∃ r, 2 ≤ card r := by
    by_contra h
    push Not at h
    have hone : ∀ r, card r = 1 := by
      intro r
      have hposr := hpos r
      have hr := h r
      omega
    have : ∑ r, card r = 11 := by simp [hone]
    omega
  obtain ⟨r, hr2⟩ := hex
  have heraseCard : (Finset.univ.erase r : Finset (Fin 11)).card = 10 := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ r)]
    simp
  have hrestLower : 10 ≤ ∑ s ∈ Finset.univ.erase r, card s := by
    rw [← heraseCard, Finset.card_eq_sum_ones]
    exact Finset.sum_le_sum fun s hs => hpos s
  have hsplit := Finset.sum_erase_add (Finset.univ : Finset (Fin 11)) card
    (Finset.mem_univ r)
  have hr : card r = 2 := by
    have hsplit' : (∑ s ∈ Finset.univ.erase r, card s) + card r = 12 := by
      simpa using hsplit.trans hsum
    omega
  have hrestSum : ∑ s ∈ Finset.univ.erase r, card s = 10 := by
    have hsplit' : (∑ s ∈ Finset.univ.erase r, card s) + card r = 12 := by
      simpa using hsplit.trans hsum
    omega
  have hones : ∀ s ∈ Finset.univ.erase r, card s = 1 := by
    have heq :
        (∑ s ∈ Finset.univ.erase r, (1 : ℕ)) =
          ∑ s ∈ Finset.univ.erase r, card s := by
      simpa [heraseCard] using hrestSum.symm
    have hall := (Finset.sum_eq_sum_iff_of_le
      (s := Finset.univ.erase r) (f := fun _ => (1 : ℕ))
      (g := card) (fun s _ => hpos s)).mp heq
    exact fun s hs => (hall s hs).symm
  refine ⟨r, ⟨?_, ?_⟩, ?_⟩
  · simpa [card] using hr
  · intro s hsr
    have hs : s ∈ Finset.univ.erase r := Finset.mem_erase.mpr ⟨hsr, Finset.mem_univ s⟩
    simpa [card] using hones s hs
  · intro s hs
    by_contra hsr
    have hsone : fiberCard pair s = 1 := by
      have hs' : s ∈ Finset.univ.erase r :=
        Finset.mem_erase.mpr ⟨hsr, Finset.mem_univ s⟩
      simpa [card] using hones s hs'
    omega

/-- Full near-bijection package from the corrected existential closeness
hypothesis. -/
theorem mod23_nearBijection (v : Fin 12 → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ v i))
    (hclose : ∀ b : ℤ, ∃ i, ∃ m : ℤ,
      |(v i : ℝ) * ((b : ℝ) / 23) - m| < 2 / 23) :
    ∃ pair : Fin 12 → Fin 11,
      (∀ i, (v i : ZMod 23) = antipodalRep (pair i) ∨
        (v i : ZMod 23) = -antipodalRep (pair i)) ∧
      Function.Surjective pair ∧
      ∃! r : Fin 11,
        fiberCard pair r = 2 ∧
        ∀ s : Fin 11, s ≠ r → fiberCard pair s = 1 := by
  have hnonzero : ∀ i, (v i : ZMod 23) ≠ 0 := by
    intro i hzero
    apply hunit i
    exact (ZMod.intCast_zmod_eq_zero_iff_dvd (v i) 23).mp hzero
  have hlabel : ∀ i, ∃! r : Fin 11,
      (v i : ZMod 23) = antipodalRep r ∨
        (v i : ZMod 23) = -antipodalRep r :=
    fun i => existsUnique_antipodalRep (v i) (hnonzero i)
  choose pair hpair using fun i => (hlabel i).exists
  have hsurj : Function.Surjective pair := by
    intro r
    obtain ⟨i, hi⟩ := antipodal_cover v hunit hclose
      (antipodalRep r) (antipodalRep_ne_zero r)
    refine ⟨i, ?_⟩
    exact (hlabel i).unique (hpair i) hi
  exact ⟨pair, hpair, hsurj, nearBijection_of_surjective pair hsurj⟩

/-- Margin-form bridge.  This is the form consumed by the LRC ledger; a
global bound `M(v)<2/23` implies the stated pointwise hypothesis. -/
theorem mod23_nearBijection_of_margin (v : Fin 12 → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ v i))
    (hgap : ∀ b : ℤ, margin v ((b : ℝ) / 23) < 2 / 23) :
    ∃ pair : Fin 12 → Fin 11,
      (∀ i, (v i : ZMod 23) = antipodalRep (pair i) ∨
        (v i : ZMod 23) = -antipodalRep (pair i)) ∧
      Function.Surjective pair ∧
      ∃! r : Fin 11,
        fiberCard pair r = 2 ∧
        ∀ s : Fin 11, s ≠ r → fiberCard pair s = 1 := by
  apply mod23_nearBijection v hunit
  intro b
  have hnot : ¬ (∀ i, ∀ m : ℤ,
      (2 : ℝ) / 23 ≤ |(v i : ℝ) * ((b : ℝ) / 23) - m|) := by
    rw [← le_margin_iff v (2 / 23) ((b : ℝ) / 23)]
    exact not_le.mpr (hgap b)
  simp only [not_forall, not_le] at hnot
  exact hnot

/-- Global-reach form of THM-1263.  Here the global margin is the supremum
of `margin v` on one compact period. -/
theorem mod23_nearBijection_of_globalMargin (v : Fin 12 → ℤ)
    (hunit : ∀ i, ¬ ((23 : ℤ) ∣ v i))
    (hreach : sSup (margin v '' Set.Icc (0 : ℝ) 1) < 2 / 23) :
    ∃ pair : Fin 12 → Fin 11,
      (∀ i, (v i : ZMod 23) = antipodalRep (pair i) ∨
        (v i : ZMod 23) = -antipodalRep (pair i)) ∧
      Function.Surjective pair ∧
      ∃! r : Fin 11,
        fiberCard pair r = 2 ∧
        ∀ s : Fin 11, s ≠ r → fiberCard pair s = 1 := by
  apply mod23_nearBijection_of_margin v hunit
  intro b
  have hr0 : (0 : ℤ) ≤ b % 23 := Int.emod_nonneg _ (by norm_num)
  have hr23 : b % 23 < 23 := Int.emod_lt_of_pos _ (by norm_num)
  have hr0R : (0 : ℝ) ≤ ((b % 23 : ℤ) : ℝ) := by exact_mod_cast hr0
  have hr23R : ((b % 23 : ℤ) : ℝ) < 23 := by exact_mod_cast hr23
  have htmem : (((b % 23 : ℤ) : ℝ) / 23) ∈ Set.Icc (0 : ℝ) 1 := by
    constructor
    · positivity
    · exact (div_lt_one (by norm_num)).mpr hr23R |>.le
  obtain ⟨t0, ht0, hmax⟩ := exists_max_margin v
  have hbdd : BddAbove (margin v '' Set.Icc (0 : ℝ) 1) := by
    refine ⟨margin v t0, ?_⟩
    rintro y ⟨t, ht, rfl⟩
    exact hmax t ht
  have hle : margin v (((b % 23 : ℤ) : ℝ) / 23) ≤
      sSup (margin v '' Set.Icc (0 : ℝ) 1) :=
    le_csSup hbdd ⟨_, htmem, rfl⟩
  have hbint : b = 23 * (b / 23) + b % 23 := by omega
  have hbdec : (b : ℝ) / 23 =
      ((b % 23 : ℤ) : ℝ) / 23 + ((b / 23 : ℤ) : ℝ) := by
    conv_lhs => rw [hbint]
    push_cast
    ring
  rw [hbdec, margin_add_int]
  exact hle.trans_lt hreach

#print axioms threshold_chain
#print axioms margin_add_int
#print axioms mod23_middle_far
#print axioms sieve23_single
#print axioms no_middle_band_of_close
#print axioms antipodal_spread
#print axioms antipodal_cover
#print axioms existsUnique_antipodalRep
#print axioms nearBijection_of_surjective
#print axioms mod23_nearBijection
#print axioms mod23_nearBijection_of_margin
#print axioms mod23_nearBijection_of_globalMargin

end LRC14.Mod23NearBijection
