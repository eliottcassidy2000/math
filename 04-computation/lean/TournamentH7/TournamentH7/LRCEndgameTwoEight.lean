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

/-- For reduced denominator at most seven, the observed bad row occupies at
most one residue class, improving the universal q-seven count from `2g/7` to
the sharp observed value `g/7`. -/
theorem reducedDenominator_mul_detunedBadBranches_card_le_of_le_seven
    (δ g q : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hqpos : 0 < q) (hq7 : q ≤ 7)
    (hq : g / (Int.gcd δ g : ℤ) = q) :
    q.toNat * (detunedBadBranches δ g u).card ≤ g.toNat := by
  by_cases hrow : (detunedBadBranches δ g u).Nonempty
  · obtain ⟨c, hc⟩ := hrow
    have hrowEq : detunedBadBranches δ g u =
        {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c [ZMOD q]} := by
      ext x
      constructor
      · intro hx
        rw [Finset.mem_filter]
        refine ⟨detunedBadBranches_subset_Ico δ g u hx, ?_⟩
        simpa [hq] using detunedBadBranches_pair_modEq_of_q_le_seven
          δ g u hg (by simpa [hq] using hq7) hx hc
      · intro hx
        have hx' := Finset.mem_filter.mp hx
        apply detunedBadBranches_mem_of_modEq_reducedDenominator
          δ g u hg hc hx'.1
        simpa [hq] using hx'.2.symm
    have hfactorZ := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ g)
    rw [hq] at hfactorZ
    have hqdvd : q ∣ g :=
      ⟨(Int.gcd δ g : ℤ), by simpa [mul_comm] using hfactorZ.symm⟩
    have hcard := Ico_zero_filter_modEq_card_of_dvd
      g q c (by omega) hqpos hqdvd
    have hdiv : g / q = (Int.gcd δ g : ℤ) := by
      calc
        g / q = ((Int.gcd δ g : ℤ) * q) / q :=
          congrArg (fun x : ℤ => x / q) hfactorZ.symm
        _ = (Int.gcd δ g : ℤ) := by
          exact Int.mul_ediv_cancel (Int.gcd δ g : ℤ) (ne_of_gt hqpos)
    have hfactorNat : q.toNat * Int.gcd δ g = g.toNat := by
      have h := congrArg Int.toNat hfactorZ
      rw [Int.toNat_mul (by positivity) (le_of_lt hqpos)] at h
      simpa [Nat.mul_comm] using h
    rw [hrowEq, hcard, hdiv]
    simpa using le_of_eq hfactorNat
  · have hzero : (detunedBadBranches δ g u).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hrow
    rw [hzero]
    omega

/-- A nonempty observed row of reduced denominator at most seven is exactly
one reduced residue class. -/
theorem detunedBadBranches_eq_filter_modEq_of_nonempty_q_le_seven
    (δ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq7 : g / (Int.gcd δ g : ℤ) ≤ 7)
    (hbad : (detunedBadBranches δ g u).Nonempty) :
    ∃ c ∈ detunedBadBranches δ g u,
      detunedBadBranches δ g u =
        {x ∈ Finset.Ico (0 : ℤ) g |
          x ≡ c [ZMOD g / (Int.gcd δ g : ℤ)]} := by
  obtain ⟨c, hc⟩ := hbad
  refine ⟨c, hc, ?_⟩
  ext x
  constructor
  · intro hx
    exact Finset.mem_filter.mpr ⟨(Finset.mem_filter.mp hx).1,
      detunedBadBranches_pair_modEq_of_q_le_seven δ g u hg hq7 hx hc⟩
  · intro hx
    have hx' := Finset.mem_filter.mp hx
    exact detunedBadBranches_mem_of_modEq_reducedDenominator
      δ g u hg hc hx'.1 hx'.2.symm

/-- Every `q≥4` row costs at most one quarter of the q-two complement sheet,
and the inequality is strict unless `q=4`.  The cases are: one-class rigidity
for q=4..7, the parity-sliced q-eight theorem, and the universal q≥9 bound. -/
theorem four_mul_companion_slice_card_le_and_lt_unless_four
    (δ₂ δₓ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hqx4 : 4 ≤ g / (Int.gcd δₓ g : ℤ))
    (hrow2 : (detunedBadBranches δ₂ g u).Nonempty) :
    4 * ((Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u) ∩
      detunedBadBranches δₓ g u).card ≤ g.toNat ∧
    (g / (Int.gcd δₓ g : ℤ) ≠ 4 →
      4 * ((Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u) ∩
        detunedBadBranches δₓ g u).card < g.toNat) := by
  let slice := (Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u) ∩
    detunedBadBranches δₓ g u
  have hsliceCard : slice.card ≤ (detunedBadBranches δₓ g u).card :=
    Finset.card_le_card Finset.inter_subset_right
  change 4 * slice.card ≤ g.toNat ∧
    (g / (Int.gcd δₓ g : ℤ) ≠ 4 → 4 * slice.card < g.toNat)
  by_cases hq7 : g / (Int.gcd δₓ g : ℤ) ≤ 7
  · have hcases : g / (Int.gcd δₓ g : ℤ) = 4 ∨
        g / (Int.gcd δₓ g : ℤ) = 5 ∨
        g / (Int.gcd δₓ g : ℤ) = 6 ∨
        g / (Int.gcd δₓ g : ℤ) = 7 := by omega
    rcases hcases with hq4 | hq5 | hq6 | hq7eq
    · have hrowBound :=
        reducedDenominator_mul_detunedBadBranches_card_le_of_le_seven
          δₓ g 4 u hg (by norm_num) (by norm_num) hq4
      change 4 * (detunedBadBranches δₓ g u).card ≤ g.toNat at hrowBound
      constructor <;> omega
    · have hrowBound :=
        reducedDenominator_mul_detunedBadBranches_card_le_of_le_seven
          δₓ g 5 u hg (by norm_num) (by norm_num) hq5
      change 5 * (detunedBadBranches δₓ g u).card ≤ g.toNat at hrowBound
      constructor <;> omega
    · have hrowBound :=
        reducedDenominator_mul_detunedBadBranches_card_le_of_le_seven
          δₓ g 6 u hg (by norm_num) (by norm_num) hq6
      change 6 * (detunedBadBranches δₓ g u).card ≤ g.toNat at hrowBound
      constructor <;> omega
    · have hrowBound :=
        reducedDenominator_mul_detunedBadBranches_card_le_of_le_seven
          δₓ g 7 u hg (by norm_num) (by norm_num) hq7eq
      change 7 * (detunedBadBranches δₓ g u).card ≤ g.toNat at hrowBound
      constructor <;> omega
  · have hq8le : 8 ≤ g / (Int.gcd δₓ g : ℤ) := by omega
    by_cases hq8 : g / (Int.gcd δₓ g : ℤ) = 8
    · have hsliceBound :=
        eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
          δ₂ δₓ g u hg hq2 hq8 hrow2
      change 8 * slice.card ≤ g.toNat at hsliceBound
      constructor <;> omega
    · have hq9 : 9 ≤ g / (Int.gcd δₓ g : ℤ) := by omega
      have hbadBound := nine_mul_badCount_le_two_mul δₓ g (by omega) hq9
      have hrowBound := detunedBadBranches_card_le δₓ g hg u
      constructor <;> omega

/-- Corrected observed-degree classification: with one q-two row and two
companions of denominator at least four, every phase clears unless both
companions have denominator exactly four. -/
theorem hasThreeDetunedGoodBranch_two_companions_four_le_unless_four_four
    (δ₂ δa δb g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hqa4 : 4 ≤ g / (Int.gcd δa g : ℤ))
    (hqb4 : 4 ≤ g / (Int.gcd δb g : ℤ))
    (hnot44 : g / (Int.gcd δa g : ℤ) ≠ 4 ∨
      g / (Int.gcd δb g : ℤ) ≠ 4) :
    HasThreeDetunedGoodBranch δ₂ δa δb g u := by
  by_cases hrow2 : (detunedBadBranches δ₂ g u).Nonempty
  · let available := Finset.Ico (0 : ℤ) g \ detunedBadBranches δ₂ g u
    let sliceA := available ∩ detunedBadBranches δa g u
    let sliceB := available ∩ detunedBadBranches δb g u
    have havailable := two_mul_qTwoComplement_card_eq δ₂ g u hg hq2 hrow2
    change 2 * available.card = g.toNat at havailable
    have hboundA := four_mul_companion_slice_card_le_and_lt_unless_four
      δ₂ δa g u hg hq2 hqa4 hrow2
    change 4 * sliceA.card ≤ g.toNat ∧
      (g / (Int.gcd δa g : ℤ) ≠ 4 → 4 * sliceA.card < g.toNat) at hboundA
    have hboundB := four_mul_companion_slice_card_le_and_lt_unless_four
      δ₂ δb g u hg hq2 hqb4 hrow2
    change 4 * sliceB.card ≤ g.toNat ∧
      (g / (Int.gcd δb g : ℤ) ≠ 4 → 4 * sliceB.card < g.toNat) at hboundB
    have hstrict : 4 * sliceA.card < g.toNat ∨ 4 * sliceB.card < g.toNat := by
      rcases hnot44 with ha | hb
      · exact Or.inl (hboundA.2 ha)
      · exact Or.inr (hboundB.2 hb)
    have hunionCard := Finset.card_union_le sliceA sliceB
    have hunionLt : (sliceA ∪ sliceB).card < available.card := by omega
    have hnotsub : ¬ available ⊆ sliceA ∪ sliceB := by
      intro hsub
      exact (Nat.not_le_of_lt hunionLt) (Finset.card_le_card hsub)
    rw [Finset.not_subset] at hnotsub
    obtain ⟨c, hcAvailable, hcOutside⟩ := hnotsub
    simp only [Finset.mem_union, not_or] at hcOutside
    have hcNotA : c ∉ detunedBadBranches δa g u := by
      intro hca
      exact hcOutside.1 (Finset.mem_inter.mpr ⟨hcAvailable, hca⟩)
    have hcNotB : c ∉ detunedBadBranches δb g u := by
      intro hcb
      exact hcOutside.2 (Finset.mem_inter.mpr ⟨hcAvailable, hcb⟩)
    have hcAvailable' := hcAvailable
    dsimp [available] at hcAvailable'
    have hcParts := Finset.mem_sdiff.mp hcAvailable'
    exact ⟨c, hcParts.1, hcParts.2, hcNotA, hcNotB⟩
  · have hzero : (detunedBadBranches δ₂ g u).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hrow2
    have hbadA := seven_mul_badCount_le_two_mul δa g (by omega) hqa4
    have hbadB := seven_mul_badCount_le_two_mul δb g (by omega) hqb4
    have hcardA := detunedBadBranches_card_le δa g hg u
    have hcardB := detunedBadBranches_card_le δb g hg u
    apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δa δb g u
    omega

theorem threeDetunedInstanceClearing_two_companions_four_le_unless_four_four
    (δ₂ δa δb g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hqa4 : 4 ≤ g / (Int.gcd δa g : ℤ))
    (hqb4 : 4 ≤ g / (Int.gcd δb g : ℤ))
    (hnot44 : g / (Int.gcd δa g : ℤ) ≠ 4 ∨
      g / (Int.gcd δb g : ℤ) ≠ 4) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δa δb g := by
  intro u
  exact (hasThreeDetunedGoodBranch_two_companions_four_le_unless_four_four
    δ₂ δa δb g u hg hq2 hqa4 hqb4 hnot44).clearances

theorem lonely14_of_three_detuned_two_companions_four_le_unless_four_four
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ ia ib : Fin 13) (h2a : i₂ ≠ ia) (h2b : i₂ ≠ ib) (hab : ia ≠ ib)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ ia → j ≠ ib → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hqa4 : 4 ≤ g / (Int.gcd (v ia) g : ℤ))
    (hqb4 : 4 ≤ g / (Int.gcd (v ib) g : ℤ))
    (hnot44 : g / (Int.gcd (v ia) g : ℤ) ≠ 4 ∨
      g / (Int.gcd (v ib) g : ℤ) ≠ 4) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ ia ib h2a h2b hab hdvd
    (threeDetunedInstanceClearing_two_companions_four_le_unless_four_four
      (v i₂) (v ia) (v ib) g (by omega) hq2 hqa4 hqb4 hnot44)

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

/-- Exact parallel-class normal form for a failing `(2,4,4)` phase.  The rows
are one mod-two class plus two distinct mod-four classes, both of parity
opposite to the q-two row.  Hence the degrees are `g/2,g/4,g/4`, all pair
intersections vanish, and the elementary incidence bound is attained. -/
theorem qTwo_four_four_failure_normal_form
    (δ₂ δ₄a δ₄b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hnogood : ¬ HasThreeDetunedGoodBranch δ₂ δ₄a δ₄b g u) :
    ∃ c₂ c₄a c₄b,
      c₂ ∈ detunedBadBranches δ₂ g u ∧
      c₄a ∈ detunedBadBranches δ₄a g u ∧
      c₄b ∈ detunedBadBranches δ₄b g u ∧
      detunedBadBranches δ₂ g u =
        {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c₂ [ZMOD 2]} ∧
      detunedBadBranches δ₄a g u =
        {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c₄a [ZMOD 4]} ∧
      detunedBadBranches δ₄b g u =
        {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c₄b [ZMOD 4]} ∧
      ¬ c₂ ≡ c₄a [ZMOD 2] ∧
      ¬ c₂ ≡ c₄b [ZMOD 2] ∧
      ¬ c₄a ≡ c₄b [ZMOD 4] := by
  have hg0 : (0 : ℤ) < g := by omega
  have hg1 : (1 : ℤ) ≤ g := by omega
  have hbad2 := two_mul_badCount_eq δ₂ g hg0 hq2
  have hbad4a :=
    four_mul_badCount_eq_of_reducedDenominator_four δ₄a g hg0 hq4a
  have hbad4b :=
    four_mul_badCount_eq_of_reducedDenominator_four δ₄b g hg0 hq4b
  have hbudget : DetunedD3.badCount δ₂ g + DetunedD3.badCount δ₄a g +
      DetunedD3.badCount δ₄b g ≤ g.toNat := by omega
  have hno2a : ¬ (detunedBadBranches δ₂ g u ∩
      detunedBadBranches δ₄a g u).Nonempty := by
    intro hinter
    exact hnogood (hasThreeDetunedGoodBranch_of_pairOverlap
      δ₂ δ₄a δ₄b g u hg1 hbudget (Or.inl hinter))
  have hno2b : ¬ (detunedBadBranches δ₂ g u ∩
      detunedBadBranches δ₄b g u).Nonempty := by
    intro hinter
    exact hnogood (hasThreeDetunedGoodBranch_of_pairOverlap
      δ₂ δ₄a δ₄b g u hg1 hbudget (Or.inr (Or.inl hinter)))
  have hno4ab : ¬ (detunedBadBranches δ₄a g u ∩
      detunedBadBranches δ₄b g u).Nonempty := by
    intro hinter
    exact hnogood (hasThreeDetunedGoodBranch_of_pairOverlap
      δ₂ δ₄a δ₄b g u hg1 hbudget (Or.inr (Or.inr hinter)))
  have hnotCardLt : ¬
      (detunedBadBranches δ₂ g u).card +
        (detunedBadBranches δ₄a g u).card +
        (detunedBadBranches δ₄b g u).card < g.toNat := by
    intro hlt
    exact hnogood (hasThreeDetunedGoodBranch_of_card_sum_lt
      δ₂ δ₄a δ₄b g u hlt)
  have hcard2 := detunedBadBranches_card_le δ₂ g hg1 u
  have hcard4a := detunedBadBranches_card_le δ₄a g hg1 u
  have hcard4b := detunedBadBranches_card_le δ₄b g hg1 u
  have hsatur2 : 2 * (detunedBadBranches δ₂ g u).card = g.toNat := by omega
  have hsatur4a : 4 * (detunedBadBranches δ₄a g u).card = g.toNat := by omega
  have hsatur4b : 4 * (detunedBadBranches δ₄b g u).card = g.toNat := by omega
  have hnonempty2 : (detunedBadBranches δ₂ g u).Nonempty := by
    rw [← Finset.card_pos]
    omega
  have hnonempty4a : (detunedBadBranches δ₄a g u).Nonempty := by
    rw [← Finset.card_pos]
    omega
  have hnonempty4b : (detunedBadBranches δ₄b g u).Nonempty := by
    rw [← Finset.card_pos]
    omega
  obtain ⟨c₂, hc₂, hrow2⟩ :=
    detunedBadBranches_eq_filter_modEq_of_nonempty_q_le_seven
      δ₂ g u hg1 (by omega) hnonempty2
  obtain ⟨c₄a, hc₄a, hrow4a⟩ :=
    detunedBadBranches_eq_filter_modEq_of_nonempty_q_le_seven
      δ₄a g u hg1 (by omega) hnonempty4a
  obtain ⟨c₄b, hc₄b, hrow4b⟩ :=
    detunedBadBranches_eq_filter_modEq_of_nonempty_q_le_seven
      δ₄b g u hg1 (by omega) hnonempty4b
  have hpar2a : ¬ c₂ ≡ c₄a [ZMOD 2] := by
    intro hmod
    apply hno2a
    refine ⟨c₄a, Finset.mem_inter.mpr ⟨?_, hc₄a⟩⟩
    apply detunedBadBranches_mem_of_modEq_reducedDenominator
      δ₂ g u hg1 hc₂ (detunedBadBranches_subset_Ico δ₄a g u hc₄a)
    simpa [hq2] using hmod
  have hpar2b : ¬ c₂ ≡ c₄b [ZMOD 2] := by
    intro hmod
    apply hno2b
    refine ⟨c₄b, Finset.mem_inter.mpr ⟨?_, hc₄b⟩⟩
    apply detunedBadBranches_mem_of_modEq_reducedDenominator
      δ₂ g u hg1 hc₂ (detunedBadBranches_subset_Ico δ₄b g u hc₄b)
    simpa [hq2] using hmod
  have hmod4ab : ¬ c₄a ≡ c₄b [ZMOD 4] := by
    intro hmod
    apply hno4ab
    refine ⟨c₄b, Finset.mem_inter.mpr ⟨?_, hc₄b⟩⟩
    apply detunedBadBranches_mem_of_modEq_reducedDenominator
      δ₄a g u hg1 hc₄a (detunedBadBranches_subset_Ico δ₄b g u hc₄b)
    simpa [hq4a] using hmod
  refine ⟨c₂, c₄a, c₄b, hc₂, hc₄a, hc₄b, ?_, ?_, ?_,
    hpar2a, hpar2b, hmod4ab⟩
  · simpa [hq2] using hrow2
  · simpa [hq4a] using hrow4a
  · simpa [hq4b] using hrow4b

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

/-- In an exactly-three detuned set, two selected distinct nonmultiples have a
unique third companion, and every other coordinate is divisible by `g`. -/
theorem exists_third_nonmultiple_of_nonMultCard_three
    (v : Fin 13 → ℤ) (g : ℤ)
    (hcard : nonMultCard v g = 3)
    (i j : Fin 13) (hij : i ≠ j)
    (hδi : ¬ g ∣ v i) (hδj : ¬ g ∣ v j) :
    ∃ k : Fin 13, i ≠ k ∧ j ≠ k ∧ ¬ g ∣ v k ∧
      ∀ x, x ≠ i → x ≠ j → x ≠ k → g ∣ v x := by
  let detuned := Finset.univ.filter (fun x => ¬ g ∣ v x)
  have hi : i ∈ detuned := by simp [detuned, hδi]
  have hj : j ∈ detuned := by simp [detuned, hδj]
  have hpairSub : ({i, j} : Finset (Fin 13)) ⊆ detuned := by
    intro x hx
    simp only [Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with rfl | rfl
    · exact hi
    · exact hj
  have hdetunedCard : detuned.card = 3 := by
    simpa [detuned, nonMultCard] using hcard
  have hpairCard : ({i, j} : Finset (Fin 13)).card = 2 := Finset.card_pair hij
  have hsum := Finset.card_sdiff_add_card_eq_card hpairSub
  have hdiffCard : (detuned \ ({i, j} : Finset (Fin 13))).card = 1 := by
    rw [hdetunedCard, hpairCard] at hsum
    omega
  obtain ⟨k, hkEq⟩ := Finset.card_eq_one.mp hdiffCard
  have hkDiff : k ∈ detuned \ ({i, j} : Finset (Fin 13)) := by
    rw [hkEq]
    simp
  have hkParts := Finset.mem_sdiff.mp hkDiff
  have hkNe : k ≠ i ∧ k ≠ j := by
    simpa only [Finset.mem_insert, Finset.mem_singleton, not_or] using hkParts.2
  have hδk : ¬ g ∣ v k := by simpa [detuned] using hkParts.1
  refine ⟨k, hkNe.1.symm, hkNe.2.symm, hδk, ?_⟩
  exact dvd_of_nonMultCard_three_of_selected
    v g hcard i j k hij hkNe.1.symm hkNe.2.symm hδi hδj hδk

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

/-- One selected q-two row has a distinct q-two companion; the third detuned
row is retained explicitly because its phase still matters. -/
def HasTwoTwoPattern (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i₂a i₂b iₓ : Fin 13,
    i₂a ≠ i₂b ∧ i₂a ≠ iₓ ∧ i₂b ≠ iₓ ∧
    ¬ g ∣ v i₂a ∧ ¬ g ∣ v i₂b ∧ ¬ g ∣ v iₓ ∧
    g / (Int.gcd (v i₂a) g : ℤ) = 2 ∧
    g / (Int.gcd (v i₂b) g : ℤ) = 2

def HasTwoFourFourPattern (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i₂ i₄a i₄b : Fin 13,
    i₂ ≠ i₄a ∧ i₂ ≠ i₄b ∧ i₄a ≠ i₄b ∧
    ¬ g ∣ v i₂ ∧ ¬ g ∣ v i₄a ∧ ¬ g ∣ v i₄b ∧
    g / (Int.gcd (v i₂) g : ℤ) = 2 ∧
    g / (Int.gcd (v i₄a) g : ℤ) = 4 ∧
    g / (Int.gcd (v i₄b) g : ℤ) = 4

/-- Exact phase-uniform residue frontier.  The q-two branch has been reduced
to `(2,2,q)` and `(2,4,4)`; the no-q-two branch is uniform `(3,3,3)`. -/
def DeepExceptionalDetunedDispatchFinalResidues : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    ((nonMultCard v g = 2 ∧ ¬ TwoAdicLiftTerminates v g) ∨
      (nonMultCard v g = 3 ∧
        (HasTwoTwoPattern v g ∨ HasTwoFourFourPattern v g ∨
          AllDetuningDenominatorsThree v g))) →
    ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

/-- The observed-degree theorem turns the final residue interface into the
earlier q-two/q-three dispatch.  This is the exact machine-checked denominator
classification; no raw `badCount` classification is claimed. -/
theorem deepExceptionalDetunedDispatchTwoThree_of_finalResidues
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchFinalResidues) :
    DeepExceptionalDetunedDispatchTwoThree := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with hpair | ⟨hcard, hpattern⟩
  · exact hdeep v hv g hg (Or.inl hpair) hnongeneric
  · rcases hpattern with hcompanion | hallThree
    · obtain ⟨i₂, j, h2j, hδ2, hδj, hq2, -⟩ := hcompanion
      obtain ⟨k, h2k, hjk, hδk, hdvd⟩ :=
        exists_third_nonmultiple_of_nonMultCard_three
          v g hcard i₂ j h2j hδ2 hδj
      have hqj2 : 2 ≤ g / (Int.gcd (v j) g : ℤ) :=
        reducedDetuningDenominator_ge_two hg hδj
      have hqk2 : 2 ≤ g / (Int.gcd (v k) g : ℤ) :=
        reducedDetuningDenominator_ge_two hg hδk
      by_cases hqjEq2 : g / (Int.gcd (v j) g : ℤ) = 2
      · exact hdeep v hv g hg
          (Or.inr ⟨hcard, Or.inl
            ⟨i₂, j, k, h2j, h2k, hjk, hδ2, hδj, hδk, hq2, hqjEq2⟩⟩)
          hnongeneric
      · by_cases hqkEq2 : g / (Int.gcd (v k) g : ℤ) = 2
        · exact hdeep v hv g hg
            (Or.inr ⟨hcard, Or.inl
              ⟨i₂, k, j, h2k, h2j, hjk.symm, hδ2, hδk, hδj, hq2, hqkEq2⟩⟩)
            hnongeneric
        · have hqj3 : 3 ≤ g / (Int.gcd (v j) g : ℤ) := by omega
          have hqk3 : 3 ≤ g / (Int.gcd (v k) g : ℤ) := by omega
          by_cases hqjEq3 : g / (Int.gcd (v j) g : ℤ) = 3
          · by_cases hqkEq3 : g / (Int.gcd (v k) g : ℤ) = 3
            · exact lonely14_of_three_detuned_two_three_three cite v hv g hg
                i₂ j k h2j h2k hjk hdvd hq2 hqjEq3 hqkEq3
            · exact lonely14_of_three_detuned_two_three_four_le cite v hv g hg
                i₂ j k h2j h2k hjk hdvd hq2 hqjEq3 (by omega)
          · by_cases hqkEq3 : g / (Int.gcd (v k) g : ℤ) = 3
            · exact lonely14_of_three_detuned_two_three_four_le cite v hv g hg
                i₂ k j h2k h2j hjk.symm
                (fun x hx2 hxk hxj => hdvd x hx2 hxj hxk)
                hq2 hqkEq3 (by omega)
            · have hqj4 : 4 ≤ g / (Int.gcd (v j) g : ℤ) := by omega
              have hqk4 : 4 ≤ g / (Int.gcd (v k) g : ℤ) := by omega
              by_cases hqjEq4 : g / (Int.gcd (v j) g : ℤ) = 4
              · by_cases hqkEq4 : g / (Int.gcd (v k) g : ℤ) = 4
                · exact hdeep v hv g hg
                    (Or.inr ⟨hcard, Or.inr (Or.inl
                      ⟨i₂, j, k, h2j, h2k, hjk, hδ2, hδj, hδk,
                        hq2, hqjEq4, hqkEq4⟩)⟩) hnongeneric
                · exact lonely14_of_three_detuned_two_companions_four_le_unless_four_four
                    cite v hv g hg i₂ j k h2j h2k hjk hdvd hq2 hqj4 hqk4
                    (Or.inr hqkEq4)
              · exact lonely14_of_three_detuned_two_companions_four_le_unless_four_four
                  cite v hv g hg i₂ j k h2j h2k hjk hdvd hq2 hqj4 hqk4
                  (Or.inl hqjEq4)
    · exact hdeep v hv g hg
        (Or.inr ⟨hcard, Or.inr (Or.inr hallThree)⟩) hnongeneric

theorem lrc14_from_finalResidues_and_relationBudget
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchFinalResidues)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_relationBudget cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite hdeep) hsupply

/-- The genuinely open nonterminating two-adic pair supplier. -/
def NonterminatingPairTowerSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    nonMultCard v g = 2 → ¬ TwoAdicLiftTerminates v g →
    ¬ genericCount v g → ∃ t : ℝ, Lonely 14 v t

/-- Joint harmonic/detuned phase selection for the `(2,2,q)` residue. -/
def TwoTwoSelectedWitnessSupply : Prop :=
  ∀ (v : Fin 13 → ℤ) (g : ℤ), 2 ≤ g →
    ∀ i₂a i₂b iₓ : Fin 13,
    i₂a ≠ i₂b → i₂a ≠ iₓ → i₂b ≠ iₓ →
    (∀ j, j ≠ i₂a → j ≠ i₂b → j ≠ iₓ → g ∣ v j) →
    g / (Int.gcd (v i₂a) g : ℤ) = 2 →
    g / (Int.gcd (v i₂b) g : ℤ) = 2 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂a i₂b iₓ u ∧
      HasThreeDetunedGoodBranch (v i₂a) (v i₂b) (v iₓ) g u

/-- Joint harmonic/detuned phase selection for the extremal `(2,4,4)` residue. -/
def TwoFourFourSelectedWitnessSupply : Prop :=
  ∀ (v : Fin 13 → ℤ) (g : ℤ), 2 ≤ g →
    ∀ i₂ i₄a i₄b : Fin 13,
    i₂ ≠ i₄a → i₂ ≠ i₄b → i₄a ≠ i₄b →
    (∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j) →
    g / (Int.gcd (v i₂) g : ℤ) = 2 →
    g / (Int.gcd (v i₄a) g : ℤ) = 4 →
    g / (Int.gcd (v i₄b) g : ℤ) = 4 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂ i₄a i₄b u ∧
      HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u

/-- Joint harmonic/detuned phase selection for the uniform q-three residue. -/
def UniformThreeSelectedWitnessSupply : Prop :=
  ∀ (v : Fin 13 → ℤ) (g : ℤ), 2 ≤ g →
    ∀ i₁ i₂ i₃ : Fin 13,
    i₁ ≠ i₂ → i₁ ≠ i₃ → i₂ ≠ i₃ →
    (∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j) →
    g / (Int.gcd (v i₁) g : ℤ) = 3 →
    g / (Int.gcd (v i₂) g : ℤ) = 3 →
    g / (Int.gcd (v i₃) g : ℤ) = 3 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u

/-- The four mathematically genuine endgame suppliers discharge the abstract
final-residue dispatch.  The selected-witness consumer avoids any universal
phase-clearing hypothesis. -/
theorem deepExceptionalDetunedDispatchFinalResidues_of_selectedWitnessSupplies
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply) :
    DeepExceptionalDetunedDispatchFinalResidues := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with ⟨hcard, hnonterm⟩ | ⟨hcard, hpattern⟩
  · exact hpairs v hv g hg hcard hnonterm hnongeneric
  · rcases hpattern with h22pattern | h244pattern | hallThree
    · obtain ⟨i₂a, i₂b, iₓ, h2ab, h2ax, h2bx, hδ2a, hδ2b, hδx,
        hq2a, hq2b⟩ := h22pattern
      have hdvd := dvd_of_nonMultCard_three_of_selected
        v g hcard i₂a i₂b iₓ h2ab h2ax h2bx hδ2a hδ2b hδx
      apply lonely14_of_three_detuned_selectedWitness v g hg i₂a i₂b iₓ hdvd
      exact h22 v g hg i₂a i₂b iₓ h2ab h2ax h2bx hdvd hq2a hq2b
    · obtain ⟨i₂, i₄a, i₄b, h24a, h24b, h4ab, hδ2, hδ4a, hδ4b,
        hq2, hq4a, hq4b⟩ := h244pattern
      have hdvd := dvd_of_nonMultCard_three_of_selected
        v g hcard i₂ i₄a i₄b h24a h24b h4ab hδ2 hδ4a hδ4b
      apply lonely14_of_three_detuned_selectedWitness v g hg i₂ i₄a i₄b hdvd
      exact h244 v g hg i₂ i₄a i₄b h24a h24b h4ab hdvd hq2 hq4a hq4b
    · have hcard' := hcard
      rw [nonMultCard] at hcard'
      obtain ⟨i₁, i₂, i₃, h12, h13, h23, hfilter⟩ :=
        Finset.card_eq_three.mp hcard'
      have hδ1 : ¬ g ∣ v i₁ := by
        have : i₁ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          rw [hfilter]
          simp
        simpa using this
      have hδ2 : ¬ g ∣ v i₂ := by
        have : i₂ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          rw [hfilter]
          simp
        simpa using this
      have hδ3 : ¬ g ∣ v i₃ := by
        have : i₃ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          rw [hfilter]
          simp
        simpa using this
      have hdvd := dvd_of_nonMultCard_three_of_selected
        v g hcard i₁ i₂ i₃ h12 h13 h23 hδ1 hδ2 hδ3
      apply lonely14_of_three_detuned_selectedWitness v g hg i₁ i₂ i₃ hdvd
      exact h333 v g hg i₁ i₂ i₃ h12 h13 h23 hdvd
        (hallThree i₁ hδ1) (hallThree i₂ hδ2) (hallThree i₃ hδ3)

/-- Sharp supplier-level capstone: after `LRCUpTo13`, only the pair tower,
three selected-witness theorems, and the trapped B5 relation budget remain. -/
theorem lrc14_from_selectedWitnessSupplies_and_relationBudget
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_finalResidues_and_relationBudget cite
    (deepExceptionalDetunedDispatchFinalResidues_of_selectedWitnessSupplies
      hpairs h22 h244 h333) hsupply

/-! ## Axiom audit -/

#print axioms detunedBadBranches_pair_modEq_of_q_eight_of_modEq_two
#print axioms eight_mul_detunedBadBranches_inter_qTwoComplement_card_le
#print axioms threeDetunedInstanceClearing_two_four_eight
#print axioms threeDetunedInstanceClearing_two_eight_eight
#print axioms lonely14_of_three_detuned_two_four_eight
#print axioms lonely14_of_three_detuned_two_eight_eight
#print axioms no_three_detuned_goodBranch_two_seven_nine
#print axioms qTwo_four_four_failure_normal_form
#print axioms deepExceptionalDetunedDispatchAfterTwoThreeCRT_of_afterTwoEight
#print axioms lrc14_from_afterTwoEight_and_relationBudget
#print axioms reducedDenominator_mul_detunedBadBranches_card_le_of_le_seven
#print axioms threeDetunedInstanceClearing_two_companions_four_le_unless_four_four
#print axioms lonely14_of_three_detuned_two_companions_four_le_unless_four_four
#print axioms deepExceptionalDetunedDispatchTwoThree_of_finalResidues
#print axioms lrc14_from_finalResidues_and_relationBudget
#print axioms deepExceptionalDetunedDispatchFinalResidues_of_selectedWitnessSupplies
#print axioms lrc14_from_selectedWitnessSupplies_and_relationBudget

end
end LRC14Grand
end LonelyRunner
