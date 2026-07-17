/-
  TournamentH7.LRCEndgameUniformThreePhase

  Exact cyclic normal form for the saturated uniform `(3,3,3)` branch left by
  `LRCEndgameParameterDischargeTwoThree`.

  If `g / gcd(delta,g) = 3`, write `delta = d*p` and `g = 3*d`.
  The reduced speed `p` is a unit modulo three, hence has residue `+1` or
  `-1`.  After the common factor `d` is removed, the ten harmonic speeds are
  multiples of three and the three exceptional speeds are units.  There are
  therefore exactly three relevant time branches, separated by `1/3`.

  A unit speed can kill at most one of those branches.  Consequently failure
  of all three branches is not merely a union-bound saturation: it is exactly
  a permutation matching the three exceptional coordinates to the three
  cyclic branch classes.  Summing the three matched phase inequalities gives
  the minimal scalar escape condition recorded below: the normalized phase
  frequency must avoid the `3/14` neighborhood of every integer.

  The sanctioned `LRC(<=13)` citation controls only ten harmonic frequencies
  here.  Integer transport of its witness rotates all three branch labels
  together and does not change the permutation obstruction.  Thus the scalar
  `3/14` phase-frequency inequality is stated honestly rather than inferred
  from the citation's `1/11` harmonic clearance.

  Zarankiewicz / relation audit: the new support-at-most-four tiny-scale
  relation forces an additive collision below scale `40`, but its quotient
  does not encode the real witness phase.  It therefore does not by itself
  rule out the cyclic matching proved here.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCDetunedOverlap

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- Exact arithmetic normal form of a reduced denominator-three coordinate. -/
theorem reducedDenominator_three_normalForm {δ g : ℤ} (hg : 0 < g)
    (hδ : ¬ g ∣ δ) (hq : g / (Int.gcd δ g : ℤ) = 3) :
    ∃ d p : ℤ, 0 < d ∧ δ = d * p ∧ g = d * 3 ∧ ¬ 3 ∣ p ∧
      (p % 3 = 1 ∨ p % 3 = 2) := by
  let d : ℤ := (Int.gcd δ g : ℤ)
  let p : ℤ := δ / d
  have hdpos : 0 < d := by
    dsimp [d]
    have : 0 < Int.gcd δ g := by
      rw [Int.gcd_pos_iff]
      right
      omega
    exact_mod_cast this
  have hddδ : d ∣ δ := Int.gcd_dvd_left δ g
  have hddg : d ∣ g := Int.gcd_dvd_right δ g
  have hδdp : δ = d * p := (Int.mul_ediv_cancel' hddδ).symm
  have hgd : g = d * 3 := by
    have hfactor := Int.mul_ediv_cancel' hddg
    change d * (g / d) = g at hfactor
    rw [hq] at hfactor
    exact hfactor.symm
  have hunit : ¬ 3 ∣ p := by
    intro hp
    obtain ⟨k, hk⟩ := hp
    apply hδ
    refine ⟨k, ?_⟩
    rw [hδdp, hk, hgd]
    ring
  have hremNonneg : 0 ≤ p % 3 := Int.emod_nonneg p (by norm_num)
  have hremLt : p % 3 < 3 := Int.emod_lt_of_pos p (by norm_num)
  have hremNe : p % 3 ≠ 0 := by
    intro hzero
    exact hunit (Int.dvd_iff_emod_eq_zero.mpr hzero)
  exact ⟨d, p, hdpos, hδdp, hgd, hunit, by omega⟩

/-- Every unit modulo three can be independently signed into residue `+1`.
This is the precise `±1 mod 3` normalization used by the phase sum. -/
theorem exists_sign_mul_residue_one (p : ℤ) (hunit : ¬ 3 ∣ p) :
    ∃ ε : ℤ, (ε = 1 ∨ ε = -1) ∧ (3 : ℤ) ∣ ε * p - 1 := by
  have hremNonneg : 0 ≤ p % 3 := Int.emod_nonneg p (by norm_num)
  have hremLt : p % 3 < 3 := Int.emod_lt_of_pos p (by norm_num)
  have hremNe : p % 3 ≠ 0 := by
    intro hzero
    exact hunit (Int.dvd_iff_emod_eq_zero.mpr hzero)
  have hdecomp := Int.emod_add_mul_ediv p 3
  rcases (show p % 3 = 1 ∨ p % 3 = 2 by omega) with hrem | hrem
  · refine ⟨1, Or.inl rfl, p / 3, ?_⟩
    rw [hrem] at hdecomp
    omega
  · refine ⟨-1, Or.inr rfl, -(p / 3) - 1, ?_⟩
    rw [hrem] at hdecomp
    omega

/-- Uniform denominator three forces the explicit common divisor `g/3`. -/
theorem uniformThree_commonDivisor (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hall : AllDetuningDenominatorsThree v g) :
    0 < g / 3 ∧ (g / 3) * 3 = g ∧ ∀ i, g / 3 ∣ v i := by
  have hg0 : (0 : ℤ) < g := by omega
  rw [nonMultCard] at hcard
  obtain ⟨i₁, i₂, i₃, _h12, _h13, _h23, hfilter⟩ :=
    Finset.card_eq_three.mp hcard
  have hδ1 : ¬ g ∣ v i₁ := by
    have : i₁ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
      rw [hfilter]
      simp
    simpa using this
  obtain ⟨d, p, hdpos, hδdp, hgd, _hunit, _hrem⟩ :=
    reducedDenominator_three_normalForm hg0 hδ1 (hall i₁ hδ1)
  have hdEq : g / 3 = d := by
    rw [hgd]
    omega
  refine ⟨by omega, by omega, ?_⟩
  intro i
  by_cases hi : g ∣ v i
  · exact dvd_trans ⟨3, by omega⟩ hi
  · obtain ⟨di, pi, _hdipos, hδi, hgdi, _hiunit, _hirem⟩ :=
      reducedDenominator_three_normalForm hg0 hi (hall i hi)
    have hdiEq : g / 3 = di := by
      rw [hgdi]
      omega
    rw [hdiEq]
    exact ⟨pi, hδi⟩

/-- On the primitive residual, the uniform denominator-three modulus is
literally `g = 3`; no larger three-adic scale survives. -/
theorem uniformThree_primitive_modulus_eq_three
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hall : AllDetuningDenominatorsThree v g)
    (hprimitive : LRC14.tupleGcd v = 1) :
    g = 3 := by
  obtain ⟨hdpos, hfactor, hdvd⟩ :=
    uniformThree_commonDivisor v g hg hcard hall
  have hdvdGcd : (g / 3).natAbs ∣ LRC14.tupleGcd v := by
    rw [LRC14.tupleGcd]
    apply Finset.dvd_gcd
    intro i _hi
    exact Int.natAbs_dvd_natAbs.mpr (hdvd i)
  rw [hprimitive] at hdvdGcd
  have habs : (g / 3).natAbs = 1 := Nat.eq_one_of_dvd_one hdvdGcd
  have hdEq : g / 3 = 1 := by
    rcases Int.natAbs_eq_iff.mp habs with h | h
    · exact h
    · omega
  omega

/-- The phase of reduced speed `p` on branch `c ∈ {0,1,2}`. -/
noncomputable def threeClassPhase (p : ℤ) (u : ℝ) (c : Fin 3) : ℝ :=
  (p : ℝ) * ((u + ((c : ℕ) : ℝ)) / 3)

/-- The reduced speed kills branch `c` at witness parameter `u`. -/
def ThreeClassBad (p : ℤ) (u : ℝ) (c : Fin 3) : Prop :=
  ∃ n : ℤ, |threeClassPhase p u c - n| < 1 / 14

/-- Negating a speed does not change which branch class it kills. -/
theorem threeClassBad_neg_iff (p : ℤ) (u : ℝ) (c : Fin 3) :
    ThreeClassBad (-p) u c ↔ ThreeClassBad p u c := by
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨-n, ?_⟩
    have heq : threeClassPhase (-p) u c - (n : ℝ) =
        -(threeClassPhase p u c - ((-n : ℤ) : ℝ)) := by
      dsimp [threeClassPhase]
      push_cast
      ring
    rw [heq, abs_neg] at hn
    simpa using hn
  · rintro ⟨n, hn⟩
    refine ⟨-n, ?_⟩
    have heq : threeClassPhase (-p) u c - ((-n : ℤ) : ℝ) =
        -(threeClassPhase p u c - (n : ℝ)) := by
      dsimp [threeClassPhase]
      push_cast
      ring
    rw [heq, abs_neg]
    exact hn

/-- Multiplication by a sign preserves the killed branch class. -/
theorem threeClassBad_sign_iff (p ε : ℤ) (u : ℝ) (c : Fin 3)
    (hε : ε = 1 ∨ ε = -1) :
    ThreeClassBad (ε * p) u c ↔ ThreeClassBad p u c := by
  rcases hε with rfl | rfl
  · simp
  · simpa using threeClassBad_neg_iff p u c

/-- One unit speed modulo three cannot kill two distinct branch classes. -/
theorem threeClassBad_branch_unique {p : ℤ} {u : ℝ} {c c' : Fin 3}
    (hunit : ¬ 3 ∣ p) (hc : ThreeClassBad p u c)
    (hc' : ThreeClassBad p u c') :
    c = c' := by
  obtain ⟨n, hn⟩ := hc
  obtain ⟨n', hn'⟩ := hc'
  rw [threeClassPhase] at hn hn'
  let k : ℤ := p * (((c : ℕ) : ℤ) - ((c' : ℕ) : ℤ)) - 3 * (n - n')
  have hklt : |(k : ℝ)| < 3 / 7 := by
    have heq : (k : ℝ) =
        3 * (((p : ℝ) * ((u + ((c : ℕ) : ℝ)) / 3) - n) -
          ((p : ℝ) * ((u + ((c' : ℕ) : ℝ)) / 3) - n')) := by
      dsimp [k]
      push_cast
      ring
    rw [heq, abs_lt] at ⊢
    rw [abs_lt] at hn hn'
    constructor <;> linarith
  have hkzero : k = 0 := by
    rw [abs_lt] at hklt
    have hloR : (-1 : ℝ) < (k : ℝ) := by linarith
    have hhiR : (k : ℝ) < 1 := by linarith
    have hlo : (-1 : ℤ) < k := by exact_mod_cast hloR
    have hhi : k < (1 : ℤ) := by exact_mod_cast hhiR
    omega
  have hdvdProd : (3 : ℤ) ∣
      p * (((c : ℕ) : ℤ) - ((c' : ℕ) : ℤ)) := by
    refine ⟨n - n', ?_⟩
    dsimp [k] at hkzero
    omega
  have hdvdDiff : (3 : ℤ) ∣
      (((c : ℕ) : ℤ) - ((c' : ℕ) : ℤ)) :=
    (Int.prime_three.dvd_mul.mp hdvdProd).resolve_left hunit
  obtain ⟨m, hm⟩ := hdvdDiff
  apply Fin.ext
  have hcLt : (c : ℕ) < 3 := c.isLt
  have hc'Lt : (c' : ℕ) < 3 := c'.isLt
  omega

/-- The exact saturated obstruction: a permutation matches the three speeds to
the three branch classes they kill. -/
def ThreeClassCyclicObstruction (p : Fin 3 → ℤ) (u : ℝ) : Prop :=
  ∃ σ : Equiv.Perm (Fin 3), ∀ c, ThreeClassBad (p (σ c)) u c

/-- Coordinatewise sign normalization preserves the exact cyclic
obstruction. -/
theorem cyclicObstruction_sign_iff (p : Fin 3 → ℤ) (ε : Fin 3 → ℤ)
    (u : ℝ) (hε : ∀ i, ε i = 1 ∨ ε i = -1) :
    ThreeClassCyclicObstruction (fun i => ε i * p i) u ↔
      ThreeClassCyclicObstruction p u := by
  constructor
  · rintro ⟨σ, hσ⟩
    refine ⟨σ, fun c => ?_⟩
    exact (threeClassBad_sign_iff (p (σ c)) (ε (σ c)) u c
      (hε (σ c))).mp (hσ c)
  · rintro ⟨σ, hσ⟩
    refine ⟨σ, fun c => ?_⟩
    exact (threeClassBad_sign_iff (p (σ c)) (ε (σ c)) u c
      (hε (σ c))).mpr (hσ c)

/-- Exact phase-intersection normal form.  For three units modulo three, some
branch clears all three speeds iff the cyclic permutation obstruction is
absent. -/
theorem exists_threeClassClear_iff_no_cyclicObstruction
    (p : Fin 3 → ℤ) (u : ℝ) (hunit : ∀ i, ¬ 3 ∣ p i) :
    (∃ c : Fin 3, ∀ i (n : ℤ),
      1 / 14 ≤ |threeClassPhase (p i) u c - n|) ↔
      ¬ ThreeClassCyclicObstruction p u := by
  constructor
  · rintro ⟨c, hclear⟩ ⟨σ, hσ⟩
    obtain ⟨n, hn⟩ := hσ c
    exact (not_lt_of_ge (hclear (σ c) n)) hn
  · intro hnobs
    by_contra hnclear
    have hkill : ∀ c : Fin 3, ∃ i, ThreeClassBad (p i) u c := by
      intro c
      by_contra hnone
      apply hnclear
      refine ⟨c, ?_⟩
      intro i n
      exact not_lt.mp (fun hlt => hnone ⟨i, n, hlt⟩)
    choose killer hkiller using hkill
    have hinj : Function.Injective killer := by
      intro c c' heq
      apply threeClassBad_branch_unique (hunit (killer c)) (hkiller c)
      simpa [heq] using hkiller c'
    have hbij : Function.Bijective killer :=
      (Finite.injective_iff_bijective).mp hinj
    let σ : Equiv.Perm (Fin 3) := Equiv.ofBijective killer hbij
    apply hnobs
    refine ⟨σ, ?_⟩
    intro c
    simpa [σ] using hkiller c

/-- The `Fin 3` cyclic phase model is exactly the concrete integer branch
model used by the three-detuned endgame at primitive modulus three. -/
theorem hasThreeDetunedGoodBranch_three_iff_no_cyclicObstruction
    (a : Fin 3 → ℤ) (u : ℝ) (hunit : ∀ i, ¬ 3 ∣ a i) :
    HasThreeDetunedGoodBranch (a 0) (a 1) (a 2) 3 u ↔
      ¬ ThreeClassCyclicObstruction a u := by
  constructor
  · rintro ⟨c, hcIco, hc0, hc1, hc2⟩
    rw [← exists_threeClassClear_iff_no_cyclicObstruction a u hunit]
    have hcBounds := Finset.mem_Ico.mp hcIco
    let branch : Fin 3 := ⟨c.toNat, by omega⟩
    have hbranchInt : ((branch : ℕ) : ℤ) = c := by
      dsimp [branch]
      rw [Int.toNat_of_nonneg hcBounds.1]
    have hbranchReal : ((branch : ℕ) : ℝ) = (c : ℝ) := by
      exact_mod_cast hbranchInt
    refine ⟨branch, ?_⟩
    intro i n
    have clear_of_not_bad (p : ℤ)
        (hnot : c ∉ detunedBadBranches p 3 u) :
        (1 : ℝ) / 14 ≤ |threeClassPhase p u branch - n| := by
      exact not_lt.mp fun hlt => hnot (by
        rw [detunedBadBranches, Finset.mem_filter]
        refine ⟨hcIco, n, ?_⟩
        rw [threeClassPhase, hbranchReal] at hlt
        norm_num at hlt ⊢
        exact hlt)
    fin_cases i
    · exact clear_of_not_bad (a 0) hc0
    · exact clear_of_not_bad (a 1) hc1
    · exact clear_of_not_bad (a 2) hc2
  · intro hnobs
    obtain ⟨branch, hclear⟩ :=
      (exists_threeClassClear_iff_no_cyclicObstruction a u hunit).mpr hnobs
    let c : ℤ := (branch : ℕ)
    have hcIco : c ∈ Finset.Ico (0 : ℤ) 3 := by
      rw [Finset.mem_Ico]
      constructor
      · exact Int.natCast_nonneg _
      · dsimp [c]
        exact_mod_cast branch.isLt
    refine ⟨c, hcIco, ?_, ?_, ?_⟩
    · intro hbad
      obtain ⟨n, hn⟩ := (Finset.mem_filter.mp hbad).2
      apply (not_lt_of_ge (hclear 0 n))
      simpa [threeClassPhase, c] using hn
    · intro hbad
      obtain ⟨n, hn⟩ := (Finset.mem_filter.mp hbad).2
      apply (not_lt_of_ge (hclear 1 n))
      simpa [threeClassPhase, c] using hn
    · intro hbad
      obtain ⟨n, hn⟩ := (Finset.mem_filter.mp hbad).2
      apply (not_lt_of_ge (hclear 2 n))
      simpa [threeClassPhase, c] using hn

/-- Failure of the concrete three-branch clearing is equivalent to the exact
cyclic matching obstruction. -/
theorem noThreeDetunedGoodBranch_three_iff_cyclicObstruction
    (a : Fin 3 → ℤ) (u : ℝ) (hunit : ∀ i, ¬ 3 ∣ a i) :
    ¬ HasThreeDetunedGoodBranch (a 0) (a 1) (a 2) 3 u ↔
      ThreeClassCyclicObstruction a u := by
  constructor
  · intro hno
    by_contra hnobs
    exact hno ((hasThreeDetunedGoodBranch_three_iff_no_cyclicObstruction
      a u hunit).mpr hnobs)
  · intro hobs hgood
    exact ((hasThreeDetunedGoodBranch_three_iff_no_cyclicObstruction
      a u hunit).mp hgood) hobs

/-- Exact synthesis of the Zarankiewicz and phase quotients at primitive
uniform denominator three: a cyclic permutation obstruction is equivalent to
the three concrete bad rows forming their full pairwise-disjoint partition of
the parallel-class circle. -/
theorem uniformThreeBadPartition_iff_cyclicObstruction
    (a : Fin 3 → ℤ) (u : ℝ) (hunit : ∀ i, ¬ 3 ∣ a i)
    (hq0 : (3 : ℤ) / (Int.gcd (a 0) 3 : ℤ) = 3)
    (hq1 : (3 : ℤ) / (Int.gcd (a 1) 3 : ℤ) = 3)
    (hq2 : (3 : ℤ) / (Int.gcd (a 2) 3 : ℤ) = 3) :
    UniformThreeBadPartition (a 0) (a 1) (a 2) 3 u ↔
      ThreeClassCyclicObstruction a u := by
  constructor
  · intro hpartition
    apply (noThreeDetunedGoodBranch_three_iff_cyclicObstruction
      a u hunit).mp
    rintro ⟨c, hcIco, hc0, hc1, hc2⟩
    have hcUnion : c ∈
        detunedBadBranches (a 0) 3 u ∪ detunedBadBranches (a 1) 3 u ∪
          detunedBadBranches (a 2) 3 u := by
      rw [hpartition.cover]
      exact hcIco
    simp [hc0, hc1, hc2] at hcUnion
  · intro hobs
    apply uniformThreeBadPartition_of_noGoodBranch
      (a 0) (a 1) (a 2) 3 u (by norm_num) hq0 hq1 hq2
    exact (noThreeDetunedGoodBranch_three_iff_cyclicObstruction
      a u hunit).mpr hobs

/-- Residue-one representatives make the scalar phase frequency integral. -/
def threePhaseFrequency (a : Fin 3 → ℤ) : ℤ :=
  (∑ i, a i) / 3

/-- A saturated cyclic matching forces the normalized sum frequency into the
open `3/14` neighborhood of an integer. -/
theorem frequency_bad_of_cyclicObstruction
    (a : Fin 3 → ℤ) (u : ℝ)
    (hresidue : ∀ i, (3 : ℤ) ∣ a i - 1)
    (hobs : ThreeClassCyclicObstruction a u) :
    ∃ n : ℤ, |(threePhaseFrequency a : ℝ) * u - n| < 3 / 14 := by
  obtain ⟨σ, hσ⟩ := hobs
  obtain ⟨n0, hn0⟩ := hσ 0
  obtain ⟨n1, hn1⟩ := hσ 1
  obtain ⟨n2, hn2⟩ := hσ 2
  let a0 : ℤ := a (σ 0)
  let a1 : ℤ := a (σ 1)
  let a2 : ℤ := a (σ 2)
  have ha0 : (3 : ℤ) ∣ a0 - 1 := hresidue (σ 0)
  have ha1 : (3 : ℤ) ∣ a1 - 1 := hresidue (σ 1)
  have ha2 : (3 : ℤ) ∣ a2 - 1 := hresidue (σ 2)
  have hsumDvd : (3 : ℤ) ∣ a0 + a1 + a2 := by
    obtain ⟨k0, hk0⟩ := ha0
    obtain ⟨k1, hk1⟩ := ha1
    obtain ⟨k2, hk2⟩ := ha2
    refine ⟨k0 + k1 + k2 + 1, ?_⟩
    omega
  have hweightDvd : (3 : ℤ) ∣ a1 + 2 * a2 := by
    obtain ⟨k1, hk1⟩ := ha1
    obtain ⟨k2, hk2⟩ := ha2
    refine ⟨k1 + 2 * k2 + 1, ?_⟩
    omega
  let B : ℤ := (a0 + a1 + a2) / 3
  let K : ℤ := (a1 + 2 * a2) / 3
  have hB : B * 3 = a0 + a1 + a2 := Int.ediv_mul_cancel hsumDvd
  have hK : K * 3 = a1 + 2 * a2 := Int.ediv_mul_cancel hweightDvd
  have hfreq : threePhaseFrequency a = B := by
    dsimp [threePhaseFrequency, B]
    rw [← Equiv.sum_comp σ a]
    simp [Fin.sum_univ_succ, a0, a1, a2, add_assoc]
  refine ⟨n0 + n1 + n2 - K, ?_⟩
  have heq : (B : ℝ) * u - ((n0 + n1 + n2 - K : ℤ) : ℝ) =
      (threeClassPhase a0 u 0 - n0) +
      (threeClassPhase a1 u 1 - n1) +
      (threeClassPhase a2 u 2 - n2) := by
    dsimp [threeClassPhase]
    push_cast
    have hBR : (B : ℝ) * 3 = (a0 : ℝ) + a1 + a2 := by exact_mod_cast hB
    have hKR : (K : ℝ) * 3 = (a1 : ℝ) + 2 * a2 := by exact_mod_cast hK
    have hBdiv : (B : ℝ) = ((a0 : ℝ) + a1 + a2) / 3 := by
      nlinarith
    have hKdiv : (K : ℝ) = ((a1 : ℝ) + 2 * a2) / 3 := by
      nlinarith
    rw [hBdiv, hKdiv]
    ring
  rw [hfreq, heq, abs_lt]
  change |threeClassPhase a0 u 0 - n0| < 1 / 14 at hn0
  change |threeClassPhase a1 u 1 - n1| < 1 / 14 at hn1
  change |threeClassPhase a2 u 2 - n2| < 1 / 14 at hn2
  rw [threeClassPhase, abs_lt] at hn0 hn1 hn2
  simp only [threeClassPhase] at ⊢
  constructor <;> linarith

/-- The minimal scalar escape inequality: if the normalized sum frequency is
`3/14`-clear at `u`, then one of the exact three branches clears every unit
speed. -/
theorem exists_threeClassClear_of_frequency_clear
    (a : Fin 3 → ℤ) (u : ℝ)
    (hresidue : ∀ i, (3 : ℤ) ∣ a i - 1)
    (hfrequency : ∀ n : ℤ,
      3 / 14 ≤ |(threePhaseFrequency a : ℝ) * u - n|) :
    ∃ c : Fin 3, ∀ i (n : ℤ),
      1 / 14 ≤ |threeClassPhase (a i) u c - n| := by
  have hunit : ∀ i, ¬ (3 : ℤ) ∣ a i := by
    intro i hai
    obtain ⟨k, hk⟩ := hai
    obtain ⟨m, hm⟩ := hresidue i
    omega
  rw [exists_threeClassClear_iff_no_cyclicObstruction a u hunit]
  intro hobs
  obtain ⟨n, hn⟩ := frequency_bad_of_cyclicObstruction a u hresidue hobs
  exact (not_lt_of_ge (hfrequency n)) hn

/-! ## A sharp witness-selection counterexample -/

/-- Three primitive denominator-three speeds for which the harmonic witness
`u = 1/7` realizes the cyclic obstruction. -/
def phaseCounterexampleSpeeds : Fin 3 → ℤ := ![1, 29, 28]

/-- The unique coordinatewise signs normalizing the example to residue one
modulo three. -/
def phaseCounterexampleSigns : Fin 3 → ℤ := ![1, -1, 1]

def phaseCounterexampleNormalized : Fin 3 → ℤ :=
  fun i => phaseCounterexampleSigns i * phaseCounterexampleSpeeds i

/-- Ten harmonic quotient speeds simultaneously clear `1/11` at `u = 1/7`.
Together with the cyclic obstruction below, this shows that an arbitrary
`LRC(10)` witness need not be the phase selected by the q-three finish. -/
def phaseCounterexampleHarmonics : Fin 10 → ℤ :=
  ![1, 2, 3, 4, 5, 6, 8, 10, 11, 13]

theorem phaseCounterexample_harmonic_clear :
    Lonely 11 phaseCounterexampleHarmonics ((1 : ℝ) / 7) := by
  have hdiv : ∀ i, ¬ ((7 : ℤ) ∣ phaseCounterexampleHarmonics i) := by
    intro i
    fin_cases i <;> norm_num [phaseCounterexampleHarmonics]
  simpa using sieve_frac 11 7 1 phaseCounterexampleHarmonics
    (by norm_num) (by norm_num) (by norm_num) hdiv

theorem phaseCounterexample_units :
    ∀ i, ¬ (3 : ℤ) ∣ phaseCounterexampleSpeeds i := by
  intro i
  fin_cases i <;> norm_num [phaseCounterexampleSpeeds]

/-- The identity matching kills branch classes zero, one, and two with
speeds `1`, `29`, and `28`, respectively. -/
theorem phaseCounterexample_cyclicObstruction :
    ThreeClassCyclicObstruction phaseCounterexampleSpeeds ((1 : ℝ) / 7) := by
  refine ⟨Equiv.refl (Fin 3), ?_⟩
  intro c
  fin_cases c
  · refine ⟨0, ?_⟩
    norm_num [threeClassPhase, phaseCounterexampleSpeeds]
  · refine ⟨11, ?_⟩
    norm_num [threeClassPhase, phaseCounterexampleSpeeds]
  · refine ⟨20, ?_⟩
    norm_num [threeClassPhase, phaseCounterexampleSpeeds]

theorem phaseCounterexample_residue_one :
    ∀ i, (3 : ℤ) ∣ phaseCounterexampleNormalized i - 1 := by
  intro i
  fin_cases i <;>
    norm_num [phaseCounterexampleNormalized, phaseCounterexampleSigns,
      phaseCounterexampleSpeeds]

theorem phaseCounterexample_normalized_cyclicObstruction :
    ThreeClassCyclicObstruction phaseCounterexampleNormalized ((1 : ℝ) / 7) := by
  apply (cyclicObstruction_sign_iff phaseCounterexampleSpeeds
    phaseCounterexampleSigns ((1 : ℝ) / 7) (by
      intro i
      fin_cases i <;> simp [phaseCounterexampleSigns])).mpr
  exact phaseCounterexample_cyclicObstruction

/-- The normalized signed sum frequency vanishes exactly because
`1 + 28 = 29`. -/
theorem phaseCounterexample_frequency_zero :
    threePhaseFrequency phaseCounterexampleNormalized = 0 := by
  norm_num [threePhaseFrequency, phaseCounterexampleNormalized,
    phaseCounterexampleSigns, phaseCounterexampleSpeeds, Fin.sum_univ_succ]

/-- Hence the sufficient `3/14` scalar gate is not forced by harmonic
clearance plus the structural q-three hypotheses.  The remaining task is
witness selection, not a pointwise contradiction of the cyclic normal form. -/
theorem phaseCounterexample_not_frequency_clear :
    ¬ ∀ n : ℤ,
      3 / 14 ≤
        |(threePhaseFrequency phaseCounterexampleNormalized : ℝ) *
          ((1 : ℝ) / 7) - n| := by
  intro hclear
  have hzero := hclear 0
  rw [phaseCounterexample_frequency_zero] at hzero
  norm_num at hzero

/-! ## Axiom audit -/

#print axioms reducedDenominator_three_normalForm
#print axioms exists_sign_mul_residue_one
#print axioms uniformThree_commonDivisor
#print axioms uniformThree_primitive_modulus_eq_three
#print axioms threeClassBad_sign_iff
#print axioms threeClassBad_branch_unique
#print axioms cyclicObstruction_sign_iff
#print axioms exists_threeClassClear_iff_no_cyclicObstruction
#print axioms hasThreeDetunedGoodBranch_three_iff_no_cyclicObstruction
#print axioms noThreeDetunedGoodBranch_three_iff_cyclicObstruction
#print axioms uniformThreeBadPartition_iff_cyclicObstruction
#print axioms frequency_bad_of_cyclicObstruction
#print axioms exists_threeClassClear_of_frequency_clear
#print axioms phaseCounterexample_harmonic_clear
#print axioms phaseCounterexample_cyclicObstruction
#print axioms phaseCounterexample_normalized_cyclicObstruction
#print axioms phaseCounterexample_frequency_zero
#print axioms phaseCounterexample_not_frequency_clear

end LRC14Grand
end LonelyRunner
