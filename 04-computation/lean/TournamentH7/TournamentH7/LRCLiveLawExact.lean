import TournamentH7.LRCLiveLaw

namespace LonelyRunner
namespace LRC14Concrete

open Finset

theorem scale_dvd_of_resonant_live (m p : ℕ) (hm : 0 < m)
    (hpq : p < 14 * m) (hlive : bandCount canonical (14 * m) p = 0) :
    m ∣ p := by
  classical
  have hq : 0 < 14 * m := by omega
  have hceil : (14 * m + 13) / 14 = m := by omega
  let f : ℕ → ℕ := fun index => (index * p) % (14 * m)
  have hf_lt (index : ℕ) : f index < 14 * m := by
    exact Nat.mod_lt _ hq
  have hgap_sep : ∀ first second : ℕ, first < second → second ≤ 13 →
      f first / m ≠ f second / m := by
    intro first second hfirstSecond hsecond13 hEq
    obtain ⟨hgap, -, -, -⟩ :=
      live_gap (14 * m) p hq hlive first second hfirstSecond hsecond13
    rw [hceil] at hgap
    change (m ≤ f second - f first ∨ m ≤ f first - f second) at hgap
    have hfirstDiv := Nat.div_add_mod (f first) m
    have hsecondDiv := Nat.div_add_mod (f second) m
    have hfirstMod : f first % m < m := Nat.mod_lt _ hm
    have hsecondMod : f second % m < m := Nat.mod_lt _ hm
    rw [← hEq] at hsecondDiv
    have hcancel : f first + f second % m = f second + f first % m := by
      linarith [hfirstDiv, hsecondDiv]
    rcases hgap with hgap | hgap <;> omega
  have hquotientInjective :
      Function.Injective (fun index : Fin 14 => f index / m) := by
    intro first second hEq
    apply Fin.ext
    rcases lt_trichotomy (first : ℕ) (second : ℕ) with hlt | heq | hgt
    · exact absurd hEq (hgap_sep first second hlt (by omega))
    · exact heq
    · exact absurd hEq.symm (hgap_sep second first hgt (by omega))
  let block : Fin 14 → Fin 14 := fun index =>
    ⟨f index / m, (Nat.div_lt_iff_lt_mul hm).2 (hf_lt index)⟩
  have hblockInjective : Function.Injective block := by
    intro first second hEq
    apply hquotientInjective
    simpa [block] using congrArg Fin.val hEq
  have hblockSurjective : Function.Surjective block :=
    Finite.surjective_of_injective hblockInjective
  let owner : Fin 14 → Fin 14 := Function.invFun block
  have hownerBlock : Function.RightInverse owner block := by
    simpa [owner] using Function.rightInverse_invFun hblockSurjective
  have hownerInjective : Function.Injective owner := hownerBlock.injective
  let residue : Fin 14 → ℕ := fun blockIndex => f (owner blockIndex)
  let offset : Fin 14 → ℕ := fun blockIndex => residue blockIndex % m
  have hresidueBlock (blockIndex : Fin 14) :
      residue blockIndex / m = (blockIndex : ℕ) := by
    have h := congrArg Fin.val (hownerBlock blockIndex)
    simpa [block, residue] using h
  have hpairGap (first second : Fin 14) (hne : first ≠ second) :
      m ≤ f first - f second ∨ m ≤ f second - f first := by
    rcases lt_trichotomy (first : ℕ) (second : ℕ) with hlt | heq | hgt
    · have hgap :=
        (live_gap (14 * m) p hq hlive first second hlt (by omega)).1
      rw [hceil] at hgap
      change (m ≤ f second - f first ∨ m ≤ f first - f second) at hgap
      exact hgap.elim Or.inr Or.inl
    · exact (hne (Fin.ext heq)).elim
    · have hgap :=
        (live_gap (14 * m) p hq hlive second first hgt (by omega)).1
      rw [hceil] at hgap
      change (m ≤ f first - f second ∨ m ≤ f second - f first) at hgap
      exact hgap
  have hpairCogap (first second : Fin 14) (hne : first ≠ second) :
      m ≤ 14 * m - (f first - f second) ∧
        m ≤ 14 * m - (f second - f first) := by
    rcases lt_trichotomy (first : ℕ) (second : ℕ) with hlt | heq | hgt
    · obtain ⟨-, -, hforward, hbackward⟩ :=
        live_gap (14 * m) p hq hlive first second hlt (by omega)
      rw [hceil] at hforward hbackward
      change m ≤ 14 * m - (f second - f first) at hforward
      change m ≤ 14 * m - (f first - f second) at hbackward
      exact ⟨hbackward, hforward⟩
    · exact (hne (Fin.ext heq)).elim
    · obtain ⟨-, -, hforward, hbackward⟩ :=
        live_gap (14 * m) p hq hlive second first hgt (by omega)
      rw [hceil] at hforward hbackward
      change m ≤ 14 * m - (f first - f second) at hforward
      change m ≤ 14 * m - (f second - f first) at hbackward
      exact ⟨hforward, hbackward⟩
  have hresidue_lt (first second : Fin 14) (hlt : first < second) :
      residue first < residue second := by
    by_contra hnot
    have hle : residue second ≤ residue first := Nat.le_of_not_gt hnot
    have hdiv := Nat.div_le_div_right (c := m) hle
    rw [hresidueBlock second, hresidueBlock first] at hdiv
    omega
  have hoffset_mono (index : ℕ) (hindex : index < 13) :
      offset ⟨index, by omega⟩ ≤ offset ⟨index + 1, by omega⟩ := by
    let first : Fin 14 := ⟨index, by omega⟩
    let second : Fin 14 := ⟨index + 1, by omega⟩
    have hfirstSecond : first < second := by simp [first, second]
    have hownerNe : owner first ≠ owner second := by
      intro hEq
      exact (ne_of_lt hfirstSecond) (hownerInjective hEq)
    have hresidueOrder : residue first < residue second :=
      hresidue_lt first second hfirstSecond
    have hgap : m ≤ residue second - residue first := by
      rcases hpairGap (owner first) (owner second) hownerNe with hgap | hgap
      · change m ≤ residue first - residue second at hgap
        omega
      · change m ≤ residue second - residue first at hgap
        exact hgap
    have hfirstDecomp :
        m * index + residue first % m = residue first := by
      have h := Nat.div_add_mod (residue first) m
      rw [hresidueBlock first] at h
      simpa [first] using h
    have hsecondDecomp :
        m * (index + 1) + residue second % m = residue second := by
      have h := Nat.div_add_mod (residue second) m
      rw [hresidueBlock second] at h
      simpa [second] using h
    have hsuccessor : m * (index + 1) = m * index + m := by ring
    change residue first % m ≤ residue second % m
    omega
  have hblockZero : block (0 : Fin 14) = 0 := by
    apply Fin.ext
    simp [block, f]
  have hownerZero : owner (0 : Fin 14) = 0 := by
    apply hblockInjective
    exact (hownerBlock 0).trans hblockZero.symm
  have hresidueZero : residue (0 : Fin 14) = 0 := by
    simp [residue, hownerZero, f]
  have hownerTopNe : owner (0 : Fin 14) ≠ owner (13 : Fin 14) := by
    intro hEq
    exact (by decide : (0 : Fin 14) ≠ 13) (hownerInjective hEq)
  have htopWrap :=
    (hpairCogap (owner 0) (owner 13) hownerTopNe).2
  change m ≤ 14 * m - (residue 13 - residue 0) at htopWrap
  rw [hresidueZero, Nat.sub_zero] at htopWrap
  have htopDecomp := Nat.div_add_mod (residue (13 : Fin 14)) m
  rw [hresidueBlock 13] at htopDecomp
  have htopModLt : residue (13 : Fin 14) % m < m := Nat.mod_lt _ hm
  have hfourteen : 14 * m = 13 * m + m := by ring
  have htopOffset : offset (13 : Fin 14) = 0 := by
    change residue (13 : Fin 14) % m = 0
    omega
  have h01 := hoffset_mono 0 (by omega)
  have h12 := hoffset_mono 1 (by omega)
  have h23 := hoffset_mono 2 (by omega)
  have h34 := hoffset_mono 3 (by omega)
  have h45 := hoffset_mono 4 (by omega)
  have h56 := hoffset_mono 5 (by omega)
  have h67 := hoffset_mono 6 (by omega)
  have h78 := hoffset_mono 7 (by omega)
  have h89 := hoffset_mono 8 (by omega)
  have h910 := hoffset_mono 9 (by omega)
  have h1011 := hoffset_mono 10 (by omega)
  have h1112 := hoffset_mono 11 (by omega)
  have h1213 := hoffset_mono 12 (by omega)
  have h01' : offset (0 : Fin 14) ≤ offset 1 := by simpa using h01
  have h12' : offset (1 : Fin 14) ≤ offset 2 := by simpa using h12
  have h23' : offset (2 : Fin 14) ≤ offset 3 := by simpa using h23
  have h34' : offset (3 : Fin 14) ≤ offset 4 := by simpa using h34
  have h45' : offset (4 : Fin 14) ≤ offset 5 := by simpa using h45
  have h56' : offset (5 : Fin 14) ≤ offset 6 := by simpa using h56
  have h67' : offset (6 : Fin 14) ≤ offset 7 := by simpa using h67
  have h78' : offset (7 : Fin 14) ≤ offset 8 := by simpa using h78
  have h89' : offset (8 : Fin 14) ≤ offset 9 := by simpa using h89
  have h910' : offset (9 : Fin 14) ≤ offset 10 := by simpa using h910
  have h1011' : offset (10 : Fin 14) ≤ offset 11 := by simpa using h1011
  have h1112' : offset (11 : Fin 14) ≤ offset 12 := by simpa using h1112
  have h1213' : offset (12 : Fin 14) ≤ offset 13 := by simpa using h1213
  have h0zero : offset (0 : Fin 14) = 0 := by omega
  have h1zero : offset (1 : Fin 14) = 0 := by omega
  have h2zero : offset (2 : Fin 14) = 0 := by omega
  have h3zero : offset (3 : Fin 14) = 0 := by omega
  have h4zero : offset (4 : Fin 14) = 0 := by omega
  have h5zero : offset (5 : Fin 14) = 0 := by omega
  have h6zero : offset (6 : Fin 14) = 0 := by omega
  have h7zero : offset (7 : Fin 14) = 0 := by omega
  have h8zero : offset (8 : Fin 14) = 0 := by omega
  have h9zero : offset (9 : Fin 14) = 0 := by omega
  have h10zero : offset (10 : Fin 14) = 0 := by omega
  have h11zero : offset (11 : Fin 14) = 0 := by omega
  have h12zero : offset (12 : Fin 14) = 0 := by omega
  have hoffsetZero (blockIndex : Fin 14) : offset blockIndex = 0 := by
    fin_cases blockIndex
    · exact h0zero
    · exact h1zero
    · exact h2zero
    · exact h3zero
    · exact h4zero
    · exact h5zero
    · exact h6zero
    · exact h7zero
    · exact h8zero
    · exact h9zero
    · exact h10zero
    · exact h11zero
    · exact h12zero
    · exact htopOffset
  have hownerOne : owner (block (1 : Fin 14)) = 1 := by
    apply hblockInjective
    exact hownerBlock (block 1)
  have hpMod : p % m = 0 := by
    have h := hoffsetZero (block (1 : Fin 14))
    change residue (block (1 : Fin 14)) % m = 0 at h
    simp only [residue, hownerOne] at h
    have hfOne : f 1 = p := by
      simp [f, Nat.mod_eq_of_lt hpq]
    change f 1 % m = 0 at h
    rwa [hfOne] at h
  exact Nat.dvd_of_mod_eq_zero hpMod

theorem live_multiplier_mem_unit_scales (m p : ℕ) (hm : 0 < m)
    (hp : 0 < p) (hpq : p < 14 * m)
    (hlive : bandCount canonical (14 * m) p = 0) :
    p ∈ (({1, 3, 5, 9, 11, 13} : Finset ℕ).image fun unit => unit * m) := by
  obtain ⟨unit, hunit⟩ := scale_dvd_of_resonant_live m p hm hpq hlive
  have hunitLt : unit < 14 := by
    rw [← Nat.mul_lt_mul_left hm]
    simpa [hunit, Nat.mul_comm] using hpq
  have hunitPos : 0 < unit := by
    by_contra hnot
    have hzero : unit = 0 := by omega
    rw [hzero] at hunit
    simp at hunit
    omega
  have hnonzero : ∀ speed : ℕ, 1 ≤ speed → speed ≤ 13 →
      (speed * unit) % 14 ≠ 0 := by
    intro speed hspeed1 hspeed13 hzero
    obtain ⟨hlower, -⟩ := live_safe (14 * m) p hlive speed hspeed1 hspeed13
    have hscale : (speed * p) % (14 * m) = ((speed * unit) % 14) * m := by
      rw [hunit]
      have hreassoc : speed * (m * unit) = (speed * unit) * m := by ring
      rw [hreassoc, Nat.mul_mod_mul_right]
    rw [hscale, hzero] at hlower
    simp at hlower
    omega
  have hcoprime : Nat.gcd unit 14 = 1 := by
    by_contra hnot
    let divisor := Nat.gcd unit 14
    have hdivisorPos : 0 < divisor := by
      dsimp [divisor]
      exact Nat.gcd_pos_of_pos_right unit (by norm_num)
    have hdivisorTwo : 2 ≤ divisor := by
      have hneOne : divisor ≠ 1 := by simpa [divisor] using hnot
      omega
    have hdivisorFourteen : divisor ∣ 14 := by
      dsimp [divisor]
      exact Nat.gcd_dvd_right unit 14
    obtain ⟨speed, hfourteenFactor⟩ := hdivisorFourteen
    have hspeedPos : 0 < speed := by
      by_contra hspeed
      have : speed = 0 := by omega
      simp [this] at hfourteenFactor
    have hspeed13 : speed ≤ 13 := by
      nlinarith
    have hdivisorUnit : divisor ∣ unit := by
      dsimp [divisor]
      exact Nat.gcd_dvd_left unit 14
    obtain ⟨factor, hunitFactor⟩ := hdivisorUnit
    apply hnonzero speed (by omega) hspeed13
    apply Nat.mod_eq_zero_of_dvd
    refine ⟨factor, ?_⟩
    rw [hfourteenFactor, hunitFactor]
    ring
  have hunitMem : unit ∈ ({1, 3, 5, 9, 11, 13} : Finset ℕ) := by
    interval_cases unit
    · norm_num
    · norm_num at hcoprime
    · norm_num
    · norm_num at hcoprime
    · norm_num
    · norm_num at hcoprime
    · norm_num at hcoprime
    · norm_num at hcoprime
    · norm_num
    · norm_num at hcoprime
    · norm_num
    · norm_num at hcoprime
    · norm_num
  exact Finset.mem_image.mpr ⟨unit, hunitMem, by simpa [Nat.mul_comm] using hunit.symm⟩

theorem liveCount_le_six_of_dvd (m : ℕ) (hm : 0 < m) :
    liveCount canonical (14 * m) ≤ 6 := by
  unfold liveCount
  have hsubset :
      ((Finset.Ioo 0 (14 * m)).filter fun p =>
          bandCount canonical (14 * m) p = 0) ⊆
        (({1, 3, 5, 9, 11, 13} : Finset ℕ).image fun unit => unit * m) := by
    intro p hp
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp
    exact live_multiplier_mem_unit_scales m p hm hp.1.1 hp.1.2 hp.2
  calc
    ((Finset.Ioo 0 (14 * m)).filter fun p =>
        bandCount canonical (14 * m) p = 0).card
        ≤ (({1, 3, 5, 9, 11, 13} : Finset ℕ).image fun unit => unit * m).card :=
          Finset.card_le_card hsubset
    _ ≤ ({1, 3, 5, 9, 11, 13} : Finset ℕ).card := Finset.card_image_le
    _ = 6 := by decide

theorem liveCount_eq_six_of_dvd (m : ℕ) (hm : 0 < m) :
    liveCount canonical (14 * m) = 6 :=
  le_antisymm (liveCount_le_six_of_dvd m hm) (liveCount_ge_six_of_dvd m hm)

theorem canonical_liveCount_exact (q : ℕ) (hq : 0 < q) :
    liveCount canonical q = if 14 ∣ q then 6 else 0 := by
  split_ifs with hdivisible
  · obtain ⟨m, rfl⟩ := hdivisible
    exact liveCount_eq_six_of_dvd m (by omega)
  · exact liveCount_eq_zero_of_not_dvd q hq hdivisible

#print axioms scale_dvd_of_resonant_live
#print axioms live_multiplier_mem_unit_scales
#print axioms liveCount_le_six_of_dvd
#print axioms liveCount_eq_six_of_dvd
#print axioms canonical_liveCount_exact

end LRC14Concrete
end LonelyRunner
