/-
  TournamentH7.LRCEndgameTwoEight

  Parity-sliced closure of the q-eight edge of the detuned triple census.
  A q-eight bad row may occupy two mod-eight classes globally, but inside a
  fixed parity class it occupies at most one.  After deleting a nonempty
  q-two row, the remaining branch circle is one parity class, so q-eight costs
  only g/8 there rather than its global g/4 bound.

  This exact local-density/Zarankiewicz value closes the exceptional triples
  (2,4,8) and (2,8,8).  Together with LRCEndgameTwoThreeSix, the residue-method
  frontier is reduced to (2,2,q) phase opposition and (2,4,4).

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCEndgameTwoThreeSix

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- Modulo two there is a unique residue class different from a fixed class. -/
theorem modEq_two_of_not_modEq {a b c : ℤ}
    (ha : ¬ a ≡ c [ZMOD 2]) (hb : ¬ b ≡ c [ZMOD 2]) :
    a ≡ b [ZMOD 2] := by
  change a % 2 = b % 2
  change ¬ a % 2 = c % 2 at ha
  change ¬ b % 2 = c % 2 at hb
  rcases Int.emod_two_eq_zero_or_one a with ha0 | ha1 <;>
    rcases Int.emod_two_eq_zero_or_one b with hb0 | hb1 <;>
      rcases Int.emod_two_eq_zero_or_one c with hc0 | hc1 <;> omega

theorem four_mul_badCount_eq_of_reducedDenominator_four
    (δ g : ℤ) (hg : 0 < g)
    (hq4 : g / (Int.gcd δ g : ℤ) = 4) :
    4 * DetunedD3.badCount δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 4 := by
    rw [hq4]
    rfl
  have hbad : DetunedD3.badCount δ g = Int.gcd δ g := by
    rw [DetunedD3.badCount, htoNat]
    norm_num
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq4] at hfactor
  rw [hbad]
  omega

theorem four_mul_badCount_eq_of_reducedDenominator_eight
    (δ g : ℤ) (hg : 0 < g)
    (hq8 : g / (Int.gcd δ g : ℤ) = 8) :
    4 * DetunedD3.badCount δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 8 := by
    rw [hq8]
    rfl
  have hbad : DetunedD3.badCount δ g = 2 * Int.gcd δ g := by
    rw [DetunedD3.badCount, htoNat]
    norm_num [Nat.mul_comm]
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq8] at hfactor
  rw [hbad]
  omega

/-- q-eight rigidity inside one parity sheet. -/
theorem detunedBadBranches_pair_modEq_of_q_eight_of_modEq_two
    (δ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq8 : g / (Int.gcd δ g : ℤ) = 8)
    {c c' : ℤ}
    (hc : c ∈ detunedBadBranches δ g u)
    (hc' : c' ∈ detunedBadBranches δ g u)
    (hpar : c ≡ c' [ZMOD 2]) :
    c ≡ c' [ZMOD 8] := by
  have hg0 : (0 : ℤ) < g := hg
  set d : ℤ := (Int.gcd δ g : ℤ) with hd
  have hdpos : 0 < d := by
    rw [hd]
    have : Int.gcd δ g ≠ 0 := by
      intro hzero
      rw [Int.gcd_eq_zero_iff] at hzero
      omega
    exact_mod_cast Nat.pos_of_ne_zero this
  have hd0 : d ≠ 0 := ne_of_gt hdpos
  have hddδ : d ∣ δ := by
    rw [hd]
    exact Int.gcd_dvd_left δ g
  have hddg : d ∣ g := by
    rw [hd]
    exact Int.gcd_dvd_right δ g
  set p : ℤ := δ / d with hp
  set q : ℤ := g / d with hq
  have hδdp : δ = d * p := (Int.mul_ediv_cancel' hddδ).symm
  have hgdq : g = d * q := (Int.mul_ediv_cancel' hddg).symm
  have hqpos : 0 < q := by
    by_contra hnot
    push Not at hnot
    have : d * q ≤ 0 := mul_nonpos_of_nonneg_of_nonpos (le_of_lt hdpos) hnot
    rw [← hgdq] at this
    omega
  have hq8' : q = 8 := by simpa [q, d, hd] using hq8
  have hcop : IsCoprime p q := by
    set A := Int.gcdA δ g
    set B := Int.gcdB δ g
    have hbez : d = δ * A + g * B := by
      rw [hd]
      exact Int.gcd_eq_gcd_ab δ g
    refine ⟨A, B, ?_⟩
    have hfac : d * (A * p + B * q) = d * 1 := by
      have h1 : d * (A * p + B * q) = (d * p) * A + (d * q) * B := by ring
      rw [h1, ← hδdp, ← hgdq, ← hbez, mul_one]
    exact mul_left_cancel₀ hd0 hfac
  obtain ⟨n, hn⟩ := (Finset.mem_filter.mp hc).2
  obtain ⟨n', hn'⟩ := (Finset.mem_filter.mp hc').2
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hqpos
  have hphase : (δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n =
      (p : ℝ) * (u + (c : ℝ)) / (q : ℝ) - n := by
    rw [hδdp, hgdq]
    push_cast
    field_simp
  have hphase' : (δ : ℝ) * ((u + (c' : ℝ)) / (g : ℝ)) - n' =
      (p : ℝ) * (u + (c' : ℝ)) / (q : ℝ) - n' := by
    rw [hδdp, hgdq]
    push_cast
    field_simp
  rw [hphase] at hn
  rw [hphase'] at hn'
  set z : ℤ := p * (c - c') - q * (n - n') with hz
  have hdiff : |(z : ℝ) / (q : ℝ)| < (1 : ℝ) / 7 := by
    have hid : (z : ℝ) / (q : ℝ) =
        ((p : ℝ) * (u + (c : ℝ)) / (q : ℝ) - n) -
          ((p : ℝ) * (u + (c' : ℝ)) / (q : ℝ) - n') := by
      rw [hz]
      push_cast
      field_simp
      ring
    rw [hid]
    calc
      |((p : ℝ) * (u + (c : ℝ)) / (q : ℝ) - n) -
          ((p : ℝ) * (u + (c' : ℝ)) / (q : ℝ) - n')|
          ≤ |(p : ℝ) * (u + (c : ℝ)) / (q : ℝ) - n| +
              |(p : ℝ) * (u + (c' : ℝ)) / (q : ℝ) - n'| := abs_sub _ _
      _ < (1 : ℝ) / 14 + 1 / 14 := add_lt_add hn hn'
      _ = (1 : ℝ) / 7 := by norm_num
  have hzlt : |(z : ℤ)| < 2 := by
    have hzR : |(z : ℝ)| < 2 := by
      have hmul := mul_lt_mul_of_pos_left hdiff hqR
      have habsEq : (q : ℝ) * |(z : ℝ) / (q : ℝ)| = |(z : ℝ)| := by
        calc
          (q : ℝ) * |(z : ℝ) / (q : ℝ)| =
              |(q : ℝ)| * |(z : ℝ) / (q : ℝ)| := by rw [abs_of_pos hqR]
          _ = |(q : ℝ) * ((z : ℝ) / (q : ℝ))| := (abs_mul _ _).symm
          _ = |(z : ℝ)| := by
            congr 1
            field_simp
      rw [habsEq] at hmul
      have hq8R : (q : ℝ) = 8 := by exact_mod_cast hq8'
      rw [hq8R] at hmul
      nlinarith
    exact_mod_cast hzR
  obtain ⟨k, hk⟩ := hpar.dvd
  have hcc : c - c' = -2 * k := by omega
  have hzdvd : (2 : ℤ) ∣ z := by
    refine ⟨-(p * k + 4 * (n - n')), ?_⟩
    rw [hz, hcc, hq8']
    ring
  obtain ⟨w, hw⟩ := hzdvd
  have hzBounds := abs_lt.mp hzlt
  have hz0 : z = 0 := by omega
  have hqdvdpc : q ∣ p * (c - c') := by
    refine ⟨n - n', ?_⟩
    have : p * (c - c') = q * (n - n') := by linarith [hz, hz0]
    exact this
  have hqdvdc : q ∣ c - c' := hcop.symm.dvd_of_dvd_mul_left hqdvdpc
  rw [Int.modEq_iff_dvd]
  rw [← hq8']
  have := dvd_neg.mpr hqdvdc
  simpa only [neg_sub] using this

/-- A nonempty q-two row leaves exactly half of the branch circle. -/
theorem two_mul_qTwoComplement_card_eq
    (δ₂ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hrow2 : (detunedBadBranches δ₂ g u).Nonempty) :
    2 * (Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u).card = g.toNat := by
  have hrowCard := two_mul_detunedBadBranches_card_eq_of_nonempty
    δ₂ g u hg hq2 hrow2
  have hbranches : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
    rw [Int.card_Ico]
    congr 1
    omega
  have hsum := Finset.card_sdiff_add_card_eq_card
    (detunedBadBranches_subset_Ico δ₂ g u)
  rw [hbranches] at hsum
  omega

/-- Exact parity-sliced Zarankiewicz value: a q-eight row occupies at most
`g/8` classes in the parity sheet complementary to a nonempty q-two row. -/
theorem eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
    (δ₂ δ₈ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq8 : g / (Int.gcd δ₈ g : ℤ) = 8)
    (hrow2 : (detunedBadBranches δ₂ g u).Nonempty) :
    8 * ((Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u) ∩
      detunedBadBranches δ₈ g u).card ≤ g.toNat := by
  let available := Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u
  let slice := available ∩ detunedBadBranches δ₈ g u
  by_cases hslice : slice.Nonempty
  · obtain ⟨c, hc⟩ := hslice
    have hcParts := Finset.mem_inter.mp hc
    have hcAvail := Finset.mem_sdiff.mp hcParts.1
    obtain ⟨a, ha2⟩ := hrow2
    have hsubset : slice ⊆
        {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c [ZMOD 8]} := by
      intro x hx
      have hxParts := Finset.mem_inter.mp hx
      have hxAvail := Finset.mem_sdiff.mp hxParts.1
      have hxaNot : ¬ x ≡ a [ZMOD 2] := by
        intro hmod
        apply hxAvail.2
        apply detunedBadBranches_mem_of_modEq_reducedDenominator
          δ₂ g u hg ha2 hxAvail.1
        simpa [hq2] using hmod.symm
      have hcaNot : ¬ c ≡ a [ZMOD 2] := by
        intro hmod
        apply hcAvail.2
        apply detunedBadBranches_mem_of_modEq_reducedDenominator
          δ₂ g u hg ha2 hcAvail.1
        simpa [hq2] using hmod.symm
      have hpar : x ≡ c [ZMOD 2] := modEq_two_of_not_modEq hxaNot hcaNot
      rw [Finset.mem_filter]
      exact ⟨hxAvail.1,
        detunedBadBranches_pair_modEq_of_q_eight_of_modEq_two
          δ₈ g u hg hq8 hxParts.2 hcParts.2 hpar⟩
    have h8dvd : (8 : ℤ) ∣ g := by
      have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ₈ g)
      rw [hq8] at hfactor
      exact ⟨(Int.gcd δ₈ g : ℤ), by simpa [mul_comm] using hfactor.symm⟩
    have hfilterCard := Ico_zero_filter_modEq_card_of_dvd
      g 8 c (by omega) (by norm_num) h8dvd
    have hcardLe := Finset.card_le_card hsubset
    rw [hfilterCard] at hcardLe
    have hfactorZ := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ₈ g)
    rw [hq8] at hfactorZ
    have hdiv : g / 8 = (Int.gcd δ₈ g : ℤ) := by
      calc
        g / 8 = ((Int.gcd δ₈ g : ℤ) * 8) / 8 :=
          congrArg (fun x : ℤ => x / 8) hfactorZ.symm
        _ = (Int.gcd δ₈ g : ℤ) := by omega
    rw [hdiv] at hcardLe
    have hfactorNat : 8 * Int.gcd δ₈ g = g.toNat := by omega
    change 8 * slice.card ≤ g.toNat
    simpa using (show 8 * slice.card ≤ g.toNat by omega)
  · have hzero : slice.card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hslice
    change 8 * slice.card ≤ g.toNat
    rw [hzero]
    omega

/-- The parity-sliced q-eight value closes every `(2,4,8)` phase. -/
theorem hasThreeDetunedGoodBranch_two_four_eight
    (δ₂ δ₄ δ₈ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq4 : g / (Int.gcd δ₄ g : ℤ) = 4)
    (hq8 : g / (Int.gcd δ₈ g : ℤ) = 8) :
    HasThreeDetunedGoodBranch δ₂ δ₄ δ₈ g u := by
  have hg0 : (0 : ℤ) < g := hg
  have hbad4 := four_mul_badCount_eq_of_reducedDenominator_four δ₄ g hg0 hq4
  have hbad8 := four_mul_badCount_eq_of_reducedDenominator_eight δ₈ g hg0 hq8
  have hcard2 := detunedBadBranches_card_le δ₂ g hg u
  have hcard4 := detunedBadBranches_card_le δ₄ g hg u
  have hcard8 := detunedBadBranches_card_le δ₈ g hg u
  by_cases hrow2 : (detunedBadBranches δ₂ g u).Nonempty
  · let available := Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u
    let slice4 := available ∩ detunedBadBranches δ₄ g u
    let slice8 := available ∩ detunedBadBranches δ₈ g u
    have havailable := two_mul_qTwoComplement_card_eq δ₂ g u hg hq2 hrow2
    change 2 * available.card = g.toNat at havailable
    have hslice4Card : slice4.card ≤ (detunedBadBranches δ₄ g u).card :=
      Finset.card_le_card Finset.inter_subset_right
    have hslice4Scaled : 4 * slice4.card ≤ g.toNat := by omega
    have hslice8Scaled :=
      eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
        δ₂ δ₈ g u hg hq2 hq8 hrow2
    change 8 * slice8.card ≤ g.toNat at hslice8Scaled
    have hunionCard := Finset.card_union_le slice4 slice8
    have hunionLt : (slice4 ∪ slice8).card < available.card := by omega
    have hnotsub : ¬ available ⊆ slice4 ∪ slice8 := by
      intro hsub
      exact (Nat.not_le_of_lt hunionLt) (Finset.card_le_card hsub)
    rw [Finset.not_subset] at hnotsub
    obtain ⟨c, hcAvailable, hcOutside⟩ := hnotsub
    simp only [Finset.mem_union, not_or] at hcOutside
    have hcNot4 : c ∉ detunedBadBranches δ₄ g u := by
      intro hc4
      exact hcOutside.1 (Finset.mem_inter.mpr ⟨hcAvailable, hc4⟩)
    have hcNot8 : c ∉ detunedBadBranches δ₈ g u := by
      intro hc8
      exact hcOutside.2 (Finset.mem_inter.mpr ⟨hcAvailable, hc8⟩)
    have hcAvailable' := hcAvailable
    dsimp [available] at hcAvailable'
    have hcParts := Finset.mem_sdiff.mp hcAvailable'
    exact ⟨c, hcParts.1, hcParts.2, hcNot4, hcNot8⟩
  · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₄ δ₈ g u
    have hzero : (detunedBadBranches δ₂ g u).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hrow2
    omega

/-- Two q-eight rows together occupy at most one quarter of the parity sheet
complementary to a q-two row, so every `(2,8,8)` phase clears. -/
theorem hasThreeDetunedGoodBranch_two_eight_eight
    (δ₂ δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8) :
    HasThreeDetunedGoodBranch δ₂ δ₈a δ₈b g u := by
  have hg0 : (0 : ℤ) < g := hg
  have hbad8a := four_mul_badCount_eq_of_reducedDenominator_eight δ₈a g hg0 hq8a
  have hbad8b := four_mul_badCount_eq_of_reducedDenominator_eight δ₈b g hg0 hq8b
  have hcard2 := detunedBadBranches_card_le δ₂ g hg u
  have hcard8a := detunedBadBranches_card_le δ₈a g hg u
  have hcard8b := detunedBadBranches_card_le δ₈b g hg u
  by_cases hrow2 : (detunedBadBranches δ₂ g u).Nonempty
  · let available := Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u
    let slice8a := available ∩ detunedBadBranches δ₈a g u
    let slice8b := available ∩ detunedBadBranches δ₈b g u
    have havailable := two_mul_qTwoComplement_card_eq δ₂ g u hg hq2 hrow2
    change 2 * available.card = g.toNat at havailable
    have hslice8aScaled :=
      eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
        δ₂ δ₈a g u hg hq2 hq8a hrow2
    change 8 * slice8a.card ≤ g.toNat at hslice8aScaled
    have hslice8bScaled :=
      eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
        δ₂ δ₈b g u hg hq2 hq8b hrow2
    change 8 * slice8b.card ≤ g.toNat at hslice8bScaled
    have hunionCard := Finset.card_union_le slice8a slice8b
    have hunionLt : (slice8a ∪ slice8b).card < available.card := by omega
    have hnotsub : ¬ available ⊆ slice8a ∪ slice8b := by
      intro hsub
      exact (Nat.not_le_of_lt hunionLt) (Finset.card_le_card hsub)
    rw [Finset.not_subset] at hnotsub
    obtain ⟨c, hcAvailable, hcOutside⟩ := hnotsub
    simp only [Finset.mem_union, not_or] at hcOutside
    have hcNot8a : c ∉ detunedBadBranches δ₈a g u := by
      intro hc8a
      exact hcOutside.1 (Finset.mem_inter.mpr ⟨hcAvailable, hc8a⟩)
    have hcNot8b : c ∉ detunedBadBranches δ₈b g u := by
      intro hc8b
      exact hcOutside.2 (Finset.mem_inter.mpr ⟨hcAvailable, hc8b⟩)
    have hcAvailable' := hcAvailable
    dsimp [available] at hcAvailable'
    have hcParts := Finset.mem_sdiff.mp hcAvailable'
    exact ⟨c, hcParts.1, hcParts.2, hcNot8a, hcNot8b⟩
  · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₈a δ₈b g u
    have hzero : (detunedBadBranches δ₂ g u).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hrow2
    omega

theorem threeDetunedInstanceClearing_two_four_eight
    (δ₂ δ₄ δ₈ g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq4 : g / (Int.gcd δ₄ g : ℤ) = 4)
    (hq8 : g / (Int.gcd δ₈ g : ℤ) = 8) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₄ δ₈ g := by
  intro u
  exact (hasThreeDetunedGoodBranch_two_four_eight
    δ₂ δ₄ δ₈ g u hg hq2 hq4 hq8).clearances

theorem threeDetunedInstanceClearing_two_eight_eight
    (δ₂ δ₈a δ₈b g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₈a δ₈b g := by
  intro u
  exact (hasThreeDetunedGoodBranch_two_eight_eight
    δ₂ δ₈a δ₈b g u hg hq2 hq8a hq8b).clearances

theorem lonely14_of_three_detuned_two_four_eight
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₄ i₈ : Fin 13) (h24 : i₂ ≠ i₄) (h28 : i₂ ≠ i₈) (h48 : i₄ ≠ i₈)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₄ → j ≠ i₈ → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq4 : g / (Int.gcd (v i₄) g : ℤ) = 4)
    (hq8 : g / (Int.gcd (v i₈) g : ℤ) = 8) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ i₄ i₈ h24 h28 h48 hdvd
    (threeDetunedInstanceClearing_two_four_eight
      (v i₂) (v i₄) (v i₈) g (by omega) hq2 hq4 hq8)

theorem lonely14_of_three_detuned_two_eight_eight
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₈a i₈b : Fin 13)
    (h2a : i₂ ≠ i₈a) (h2b : i₂ ≠ i₈b) (hab : i₈a ≠ i₈b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₈a → j ≠ i₈b → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq8a : g / (Int.gcd (v i₈a) g : ℤ) = 8)
    (hq8b : g / (Int.gcd (v i₈b) g : ℤ) = 8) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ i₈a i₈b h2a h2b hab hdvd
    (threeDetunedInstanceClearing_two_eight_eight
      (v i₂) (v i₈a) (v i₈b) g (by omega) hq2 hq8a hq8b)

/-- The remaining `(2,4,4)` boundary is genuine for the phase-uniform branch
interface: at `g=4`, rows `2`, `7`, and `9` cover all four branches. -/
theorem no_three_detuned_goodBranch_two_seven_nine :
    ¬ HasThreeDetunedGoodBranch 2 7 9 4 ((11 : ℝ) / 100) := by
  rintro ⟨c, hcIco, hc1, hc2, hc3⟩
  have hcBounds := Finset.mem_Ico.mp hcIco
  have hc0 : 0 ≤ c := hcBounds.1
  have hc4 : c < 4 := hcBounds.2
  interval_cases c
  · apply hc1
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 0, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc2
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 2, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc1
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 1, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc3
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 7, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]

/-- Selected nonmultiples exhaust an exactly-three detuned set. -/
theorem dvd_of_nonMultCard_three_of_selected
    (v : Fin 13 → ℤ) (g : ℤ)
    (hcard : nonMultCard v g = 3)
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) (hδ3 : ¬ g ∣ v i₃) :
    ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j := by
  have hselectedCard : ({i₁, i₂, i₃} : Finset (Fin 13)).card = 3 :=
    Finset.card_eq_three.mpr ⟨i₁, i₂, i₃, h12, h13, h23, rfl⟩
  have hsubset : ({i₁, i₂, i₃} : Finset (Fin 13)) ⊆
      Finset.univ.filter (fun i => ¬ g ∣ v i) := by
    intro i hi
    simp only [Finset.mem_insert, Finset.mem_singleton] at hi
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    rcases hi with rfl | rfl | rfl
    · exact hδ1
    · exact hδ2
    · exact hδ3
  have hfilterCard :
      (Finset.univ.filter (fun i => ¬ g ∣ v i)).card = 3 := by
    simpa [nonMultCard] using hcard
  have hselected : ({i₁, i₂, i₃} : Finset (Fin 13)) =
      Finset.univ.filter (fun i => ¬ g ∣ v i) := by
    apply Finset.eq_of_subset_of_card_le hsubset
    rw [hfilterCard, hselectedCard]
  intro j hj1 hj2 hj3
  by_contra hj
  have hjFilter : j ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by simp [hj]
  have hjSelected : j ∈ ({i₁, i₂, i₃} : Finset (Fin 13)) := by
    rw [hselected]
    exact hjFilter
  simp only [Finset.mem_insert, Finset.mem_singleton] at hjSelected
  rcases hjSelected with h | h | h
  · exact hj1 h
  · exact hj2 h
  · exact hj3 h

def HasTwoFourEightPattern (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i₂ i₄ i₈ : Fin 13,
    i₂ ≠ i₄ ∧ i₂ ≠ i₈ ∧ i₄ ≠ i₈ ∧
    ¬ g ∣ v i₂ ∧ ¬ g ∣ v i₄ ∧ ¬ g ∣ v i₈ ∧
    g / (Int.gcd (v i₂) g : ℤ) = 2 ∧
    g / (Int.gcd (v i₄) g : ℤ) = 4 ∧
    g / (Int.gcd (v i₈) g : ℤ) = 8

def HasTwoEightEightPattern (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i₂ i₈a i₈b : Fin 13,
    i₂ ≠ i₈a ∧ i₂ ≠ i₈b ∧ i₈a ≠ i₈b ∧
    ¬ g ∣ v i₂ ∧ ¬ g ∣ v i₈a ∧ ¬ g ∣ v i₈b ∧
    g / (Int.gcd (v i₂) g : ℤ) = 2 ∧
    g / (Int.gcd (v i₈a) g : ℤ) = 8 ∧
    g / (Int.gcd (v i₈b) g : ℤ) = 8

def HasTwoEightClosedPattern (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  HasTwoFourEightPattern v g ∨ HasTwoEightEightPattern v g

/-- The exceptional dispatch after removing both q-two/q-three CRT patterns
and the parity-sliced q-eight patterns. -/
def DeepExceptionalDetunedDispatchAfterTwoEight : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    ((nonMultCard v g = 2 ∧ ¬ TwoAdicLiftTerminates v g) ∨
      (nonMultCard v g = 3 ∧
        (HasTwoWithSubNineCompanion v g ∨
          AllDetuningDenominatorsThree v g) ∧
        ¬ HasTwoThreeCRTClosedPattern v g ∧
        ¬ HasTwoEightClosedPattern v g)) →
    ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

theorem deepExceptionalDetunedDispatchAfterTwoThreeCRT_of_afterTwoEight
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchAfterTwoEight) :
    DeepExceptionalDetunedDispatchAfterTwoThreeCRT := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with hpair | ⟨hcard, hpattern, hnoCRT⟩
  · exact hdeep v hv g hg (Or.inl hpair) hnongeneric
  · by_cases hclosed8 : HasTwoEightClosedPattern v g
    · rcases hclosed8 with h248 | h288
      · obtain ⟨i₂, i₄, i₈, h24, h28, h48, hδ2, hδ4, hδ8,
          hq2, hq4, hq8⟩ := h248
        exact lonely14_of_three_detuned_two_four_eight cite v hv g hg
          i₂ i₄ i₈ h24 h28 h48
          (dvd_of_nonMultCard_three_of_selected
            v g hcard i₂ i₄ i₈ h24 h28 h48 hδ2 hδ4 hδ8)
          hq2 hq4 hq8
      · obtain ⟨i₂, i₈a, i₈b, h2a, h2b, hab, hδ2, hδ8a, hδ8b,
          hq2, hq8a, hq8b⟩ := h288
        exact lonely14_of_three_detuned_two_eight_eight cite v hv g hg
          i₂ i₈a i₈b h2a h2b hab
          (dvd_of_nonMultCard_three_of_selected
            v g hcard i₂ i₈a i₈b h2a h2b hab hδ2 hδ8a hδ8b)
          hq2 hq8a hq8b
    · exact hdeep v hv g hg
        (Or.inr ⟨hcard, hpattern, hnoCRT, hclosed8⟩) hnongeneric

theorem lrc14_from_afterTwoEight_and_relationBudget
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchAfterTwoEight)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_afterTwoThreeCRT_and_relationBudget cite
    (deepExceptionalDetunedDispatchAfterTwoThreeCRT_of_afterTwoEight cite hdeep)
    hsupply

/-! ## Axiom audit -/

#print axioms detunedBadBranches_pair_modEq_of_q_eight_of_modEq_two
#print axioms eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
#print axioms threeDetunedInstanceClearing_two_four_eight
#print axioms threeDetunedInstanceClearing_two_eight_eight
#print axioms lonely14_of_three_detuned_two_four_eight
#print axioms lonely14_of_three_detuned_two_eight_eight
#print axioms no_three_detuned_goodBranch_two_seven_nine
#print axioms deepExceptionalDetunedDispatchAfterTwoThreeCRT_of_afterTwoEight
#print axioms lrc14_from_afterTwoEight_and_relationBudget

end
end LRC14Grand
end LonelyRunner
