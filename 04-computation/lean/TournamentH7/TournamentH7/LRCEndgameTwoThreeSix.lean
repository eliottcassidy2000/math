/-
  TournamentH7.LRCEndgameTwoThreeSix

  First phase-sensitive closure inside the q-two exceptional triple alphabet.
  For reduced denominator at most seven, two bad branches lie in the same
  reduced residue class; membership is invariant along that class.  Therefore
  nonempty q-two and q-three bad rows meet by the concrete CRT modulo six.

  In the saturated degree pattern `(2,3,6)`, the three universal row bounds
  sum exactly to `g`.  If either of the first two rows is empty, the observed
  union is strictly subcritical; if both are nonempty, their CRT intersection
  pays one unit of overlap debt.  Either way a common good branch exists for
  every harmonic phase, and the standard three-detuned construction produces
  an actual LRC(14) witness.

  Parallel-class audit: vertices are branch residue classes, the pairwise
  observable is exact row intersection, and CRT is the gauge converting the
  q-two/q-three pair into a binary collision.  The quotient preserves the
  common-good-branch predicate here and discards the internal ordering of each
  residue class.  No tournament orientation or K_{2,t}-free assumption is used.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCDetunedOverlap
import TournamentH7.LRCB5RelationEndgame
import Mathlib.Data.Int.CardIntervalMod

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

theorem detunedBadBranches_pair_modEq_of_q_le_seven
    (δ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq7 : g / (Int.gcd δ g : ℤ) ≤ 7)
    {c c' : ℤ}
    (hc : c ∈ detunedBadBranches δ g u)
    (hc' : c' ∈ detunedBadBranches δ g u) :
    c ≡ c' [ZMOD g / (Int.gcd δ g : ℤ)] := by
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
  have hq7' : q ≤ 7 := by simpa [q, d, hd] using hq7
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
  have hgR : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  have hdR : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hdpos
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
  have hzlt : |(z : ℤ)| < 1 := by
    have hzR : |(z : ℝ)| < 1 := by
      have hq7R : (q : ℝ) ≤ 7 := by exact_mod_cast hq7'
      have hmul := mul_lt_mul_of_pos_left hdiff hqR
      have habsEq : (q : ℝ) * |(z : ℝ) / (q : ℝ)| = |(z : ℝ)| := by
        calc
          (q : ℝ) * |(z : ℝ) / (q : ℝ)| = |(q : ℝ)| * |(z : ℝ) / (q : ℝ)| := by
            rw [abs_of_pos hqR]
          _ = |(q : ℝ) * ((z : ℝ) / (q : ℝ))| := (abs_mul _ _).symm
          _ = |(z : ℝ)| := by
            congr 1
            field_simp
      rw [habsEq] at hmul
      nlinarith
    exact_mod_cast hzR
  have hz0 : z = 0 := Int.abs_lt_one_iff.mp hzlt
  have hqdvdpc : q ∣ p * (c - c') := by
    refine ⟨n - n', ?_⟩
    have : p * (c - c') = q * (n - n') := by linarith [hz, hz0]
    exact this
  have hqdvdc : q ∣ c - c' := hcop.symm.dvd_of_dvd_mul_left hqdvdpc
  rw [Int.modEq_iff_dvd]
  have := dvd_neg.mpr hqdvdc
  simpa only [neg_sub] using this

theorem detunedBadBranches_mem_of_modEq_reducedDenominator
    (δ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    {c c' : ℤ}
    (hc : c ∈ detunedBadBranches δ g u)
    (hc'Ico : c' ∈ Finset.Ico (0 : ℤ) g)
    (hmod : c ≡ c' [ZMOD g / (Int.gcd δ g : ℤ)]) :
    c' ∈ detunedBadBranches δ g u := by
  have hg0 : (0 : ℤ) < g := hg
  set d : ℤ := (Int.gcd δ g : ℤ) with hd
  have hdpos : 0 < d := by
    rw [hd]
    have : Int.gcd δ g ≠ 0 := by
      intro hzero
      rw [Int.gcd_eq_zero_iff] at hzero
      omega
    exact_mod_cast Nat.pos_of_ne_zero this
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
  obtain ⟨k, hk⟩ := hmod.dvd
  obtain ⟨n, hn⟩ := (Finset.mem_filter.mp hc).2
  rw [detunedBadBranches, Finset.mem_filter]
  refine ⟨hc'Ico, n + p * k, ?_⟩
  have heq : (δ : ℝ) * ((u + (c' : ℝ)) / (g : ℝ)) - ((n + p * k : ℤ) : ℝ) =
      (δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n := by
    have hkR : (c' : ℝ) - (c : ℝ) = (q : ℝ) * (k : ℝ) := by exact_mod_cast hk
    rw [hδdp, hgdq]
    push_cast
    have hdR : (d : ℝ) ≠ 0 := by exact_mod_cast (ne_of_gt hdpos)
    have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast (ne_of_gt hqpos)
    field_simp
    linear_combination (p : ℝ) * hkR
  rwa [heq]

/-- Two low-denominator bad rows with the same reduced denominator are equal
as soon as they meet.  Thus their only possible phase relation is equality or
disjointness; for q-two this is the exact same-parity/opposite-parity split. -/
theorem detunedBadBranches_eq_of_overlap_same_reducedDenominator
    (δ₁ δ₂ g q : ℤ) (u : ℝ) (hg : 1 ≤ g) (hq7 : q ≤ 7)
    (hq1 : g / (Int.gcd δ₁ g : ℤ) = q)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = q)
    (hoverlap :
      (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₂ g u).Nonempty) :
    detunedBadBranches δ₁ g u = detunedBadBranches δ₂ g u := by
  obtain ⟨c, hc⟩ := hoverlap
  have ⟨hc1, hc2⟩ := Finset.mem_inter.mp hc
  ext x
  constructor
  · intro hx
    have hxIco := detunedBadBranches_subset_Ico δ₁ g u hx
    have hxmod : x ≡ c [ZMOD q] := by
      simpa [hq1] using detunedBadBranches_pair_modEq_of_q_le_seven
        δ₁ g u hg (by simpa [hq1] using hq7) hx hc1
    apply detunedBadBranches_mem_of_modEq_reducedDenominator
      δ₂ g u hg hc2 hxIco
    simpa [hq2] using hxmod.symm
  · intro hx
    have hxIco := detunedBadBranches_subset_Ico δ₂ g u hx
    have hxmod : x ≡ c [ZMOD q] := by
      simpa [hq2] using detunedBadBranches_pair_modEq_of_q_le_seven
        δ₂ g u hg (by simpa [hq2] using hq7) hx hc2
    apply detunedBadBranches_mem_of_modEq_reducedDenominator
      δ₁ g u hg hc1 hxIco
    simpa [hq1] using hxmod.symm

theorem exists_Ico_modEq_two_three (g a b : ℤ) (hg6 : 6 ≤ g) :
    ∃ c ∈ Finset.Ico (0 : ℤ) g, c ≡ a [ZMOD 2] ∧ c ≡ b [ZMOD 3] := by
  let c : ℤ := (3 * a + 4 * b) % 6
  have hc0 : 0 ≤ c := Int.emod_nonneg _ (by norm_num)
  have hc6 : c < 6 := Int.emod_lt_of_pos _ (by norm_num)
  have hcmod6 : c ≡ 3 * a + 4 * b [ZMOD 6] := Int.mod_modEq _ _
  have hcmod2 : c ≡ a [ZMOD 2] := by
    apply (hcmod6.of_dvd (by norm_num : (2 : ℤ) ∣ 6)).trans
    rw [Int.modEq_iff_dvd]
    refine ⟨-(a + 2 * b), ?_⟩
    ring
  have hcmod3 : c ≡ b [ZMOD 3] := by
    apply (hcmod6.of_dvd (by norm_num : (3 : ℤ) ∣ 6)).trans
    rw [Int.modEq_iff_dvd]
    refine ⟨-(a + b), ?_⟩
    ring
  exact ⟨c, Finset.mem_Ico.mpr ⟨hc0, by omega⟩, hcmod2, hcmod3⟩

theorem detunedBadBranches_two_three_overlap
    (δ₂ δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g) (hg6 : 6 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hbad2 : (detunedBadBranches δ₂ g u).Nonempty)
    (hbad3 : (detunedBadBranches δ₃ g u).Nonempty) :
    (detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u).Nonempty := by
  obtain ⟨c₂, hc₂⟩ := hbad2
  obtain ⟨c₃, hc₃⟩ := hbad3
  obtain ⟨c, hcIco, hc2mod, hc3mod⟩ := exists_Ico_modEq_two_three g c₂ c₃ hg6
  have hc₂' : c ∈ detunedBadBranches δ₂ g u := by
    apply detunedBadBranches_mem_of_modEq_reducedDenominator δ₂ g u hg hc₂ hcIco
    simpa [hq2] using hc2mod.symm
  have hc₃' : c ∈ detunedBadBranches δ₃ g u := by
    apply detunedBadBranches_mem_of_modEq_reducedDenominator δ₃ g u hg hc₃ hcIco
    simpa [hq3] using hc3mod.symm
  exact ⟨c, Finset.mem_inter.mpr ⟨hc₂', hc₃'⟩⟩

theorem Ico_zero_filter_modEq_card_of_dvd
    (g r v : ℤ) (hg : 0 ≤ g) (hr : 0 < r) (hrg : r ∣ g) :
    ({x ∈ Finset.Ico (0 : ℤ) g | x ≡ v [ZMOD r]}).card = (g / r).toNat := by
  rw [← Int.ofNat_inj]
  rw [Int.Ico_filter_modEq_card (r := r) 0 g hr v]
  obtain ⟨k, rfl⟩ := hrg
  have hk : 0 ≤ k := by nlinarith
  push_cast
  have hrat : (((r : ℚ) * (k : ℚ) - (v : ℚ)) / (r : ℚ)) =
      (-(v : ℚ) / (r : ℚ)) + (k : ℚ) := by
    have hrQ : (r : ℚ) ≠ 0 := by exact_mod_cast hr.ne'
    field_simp
    ring
  rw [hrat, Int.ceil_add_intCast]
  have hzero : (((0 : ℚ) - (v : ℚ)) / (r : ℚ)) =
      (-(v : ℚ) / (r : ℚ)) := by ring
  rw [hzero]
  have hdiv : r * k / r = k := Int.mul_ediv_cancel_left k hr.ne'
  rw [hdiv]
  omega

/-- A nonempty q-two bad row attains its universal half-circle bound exactly. -/
theorem two_mul_detunedBadBranches_card_eq_of_nonempty
    (δ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ g : ℤ) = 2)
    (hbad : (detunedBadBranches δ g u).Nonempty) :
    2 * (detunedBadBranches δ g u).card = g.toNat := by
  obtain ⟨c, hc⟩ := hbad
  have hrow : detunedBadBranches δ g u =
      {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c [ZMOD 2]} := by
    ext x
    constructor
    · intro hx
      rw [Finset.mem_filter]
      refine ⟨detunedBadBranches_subset_Ico δ g u hx, ?_⟩
      simpa [hq2] using detunedBadBranches_pair_modEq_of_q_le_seven
        δ g u hg (by omega) hx hc
    · intro hx
      have hx' := Finset.mem_filter.mp hx
      apply detunedBadBranches_mem_of_modEq_reducedDenominator δ g u hg
        hc hx'.1
      simpa [hq2] using hx'.2.symm
  have h2dvd : (2 : ℤ) ∣ g := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ g)
    rw [hq2] at hfactor
    exact ⟨(Int.gcd δ g : ℤ), by simpa [mul_comm] using hfactor.symm⟩
  have hcard := Ico_zero_filter_modEq_card_of_dvd
    g 2 c (by omega) (by norm_num) h2dvd
  have hfactor := (badCount_of_q_two (δ := δ) (g := g) (by omega) hq2).2
  have hfactorZ := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ g)
  rw [hq2] at hfactorZ
  have hdiv : g / 2 = (Int.gcd δ g : ℤ) := by
    calc
      g / 2 = ((Int.gcd δ g : ℤ) * 2) / 2 :=
        congrArg (fun x : ℤ => x / 2) hfactorZ.symm
      _ = (Int.gcd δ g : ℤ) := by omega
  rw [hrow, hcard, hdiv]
  simpa using hfactor

theorem detunedBadBranches_two_three_inter_eq_modEq
    (δ₂ δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g) (hg6 : 6 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hbad2 : (detunedBadBranches δ₂ g u).Nonempty)
    (hbad3 : (detunedBadBranches δ₃ g u).Nonempty) :
    ∃ c ∈ detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u,
      detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u =
        {x ∈ Finset.Ico (0 : ℤ) g | x ≡ c [ZMOD 6]} := by
  obtain ⟨c, hc⟩ := detunedBadBranches_two_three_overlap
    δ₂ δ₃ g u hg hg6 hq2 hq3 hbad2 hbad3
  refine ⟨c, hc, ?_⟩
  ext x
  constructor
  · intro hx
    have hxrows := Finset.mem_inter.mp hx
    have hcrows := Finset.mem_inter.mp hc
    have hx2 : x ≡ c [ZMOD 2] := by
      simpa [hq2] using detunedBadBranches_pair_modEq_of_q_le_seven
        δ₂ g u hg (by omega) hxrows.1 hcrows.1
    have hx3 : x ≡ c [ZMOD 3] := by
      simpa [hq3] using detunedBadBranches_pair_modEq_of_q_le_seven
        δ₃ g u hg (by omega) hxrows.2 hcrows.2
    rw [Finset.mem_filter]
    refine ⟨(Finset.mem_filter.mp hxrows.1).1, ?_⟩
    simpa using (Int.modEq_and_modEq_iff_modEq_mul
      (by norm_num : (2 : ℤ).natAbs.Coprime (3 : ℤ).natAbs)).mp ⟨hx2, hx3⟩
  · intro hx
    have hx' := Finset.mem_filter.mp hx
    have hx2 : x ≡ c [ZMOD 2] := hx'.2.of_dvd (by norm_num)
    have hx3 : x ≡ c [ZMOD 3] := hx'.2.of_dvd (by norm_num)
    have hcrows := Finset.mem_inter.mp hc
    rw [Finset.mem_inter]
    constructor
    · apply detunedBadBranches_mem_of_modEq_reducedDenominator δ₂ g u hg
        hcrows.1 hx'.1
      simpa [hq2] using hx2.symm
    · apply detunedBadBranches_mem_of_modEq_reducedDenominator δ₃ g u hg
        hcrows.2 hx'.1
      simpa [hq3] using hx3.symm

theorem detunedBadBranches_two_three_inter_card_eq
    (δ₂ δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g) (hg6 : 6 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hbad2 : (detunedBadBranches δ₂ g u).Nonempty)
    (hbad3 : (detunedBadBranches δ₃ g u).Nonempty) :
    (detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u).card =
      (g / 6).toNat := by
  obtain ⟨c, -, hinter⟩ := detunedBadBranches_two_three_inter_eq_modEq
    δ₂ δ₃ g u hg hg6 hq2 hq3 hbad2 hbad3
  rw [hinter]
  apply Ico_zero_filter_modEq_card_of_dvd g 6 c (by omega) (by norm_num)
  have h2dvd : (2 : ℤ) ∣ g := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ₂ g)
    rw [hq2] at hfactor
    exact ⟨(Int.gcd δ₂ g : ℤ), by simpa [mul_comm] using hfactor.symm⟩
  have h3dvd : (3 : ℤ) ∣ g := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ₃ g)
    rw [hq3] at hfactor
    exact ⟨(Int.gcd δ₃ g : ℤ), by simpa [mul_comm] using hfactor.symm⟩
  simpa using (show IsCoprime (2 : ℤ) 3 by norm_num).mul_dvd h2dvd h3dvd


theorem hasThreeDetunedGoodBranch_of_card_sum_lt
    (δ₁ δ₂ δ₃ g : ℤ) (u : ℝ)
    (hcards : (detunedBadBranches δ₁ g u).card +
      (detunedBadBranches δ₂ g u).card +
      (detunedBadBranches δ₃ g u).card < g.toNat) :
    HasThreeDetunedGoodBranch δ₁ δ₂ δ₃ g u := by
  let first := detunedBadBranches δ₁ g u
  let second := detunedBadBranches δ₂ g u
  let third := detunedBadBranches δ₃ g u
  let branches := Finset.Ico (0 : ℤ) g
  have hbranches : branches.card = g.toNat := by
    dsimp [branches]
    rw [Int.card_Ico]
    congr 1
    omega
  have hunion_lt : (first ∪ second ∪ third).card < branches.card := by
    calc
      (first ∪ second ∪ third).card ≤ (first ∪ second).card + third.card :=
        Finset.card_union_le _ _
      _ ≤ first.card + second.card + third.card := by
        exact Nat.add_le_add_right (Finset.card_union_le first second) third.card
      _ < g.toNat := hcards
      _ = branches.card := hbranches.symm
  have hnotsub : ¬ branches ⊆ first ∪ second ∪ third := by
    intro hsub
    exact (Nat.not_le_of_lt hunion_lt) (Finset.card_le_card hsub)
  rw [Finset.not_subset] at hnotsub
  obtain ⟨c, hcBranches, hcUnion⟩ := hnotsub
  simp only [Finset.mem_union, not_or] at hcUnion
  exact ⟨c, hcBranches, hcUnion.1.1, hcUnion.1.2, hcUnion.2⟩

theorem six_mul_badCount_eq (δ g : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd δ g : ℤ) = 6) :
    6 * DetunedD3.badCount δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 6 := by
    rw [hq]
    rfl
  have hbad : DetunedD3.badCount δ g = Int.gcd δ g := by
    rw [DetunedD3.badCount, htoNat]
    norm_num
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq] at hfactor
  rw [hbad]
  omega

theorem six_le_of_reducedDenominator_six (δ g : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd δ g : ℤ) = 6) : 6 ≤ g := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq] at hfactor
  have hdpos : (0 : ℤ) < (Int.gcd δ g : ℤ) := by
    have : 0 < Int.gcd δ g := by
      rw [Int.gcd_pos_iff]
      right
      omega
    exact_mod_cast this
  omega

/-- Simultaneous reduced denominators two and three force at least six branch
classes, independently of the third row. -/
theorem six_le_of_reducedDenominators_two_three
    (δ₂ δ₃ g : ℤ) (hg : 0 < g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3) :
    6 ≤ g := by
  have hdvd2 : ((Int.gcd δ₂ g : ℤ)) ∣ g := Int.gcd_dvd_right δ₂ g
  have hdvd3 : ((Int.gcd δ₃ g : ℤ)) ∣ g := Int.gcd_dvd_right δ₃ g
  have hfactor2 := Int.mul_ediv_cancel' hdvd2
  have hfactor3 := Int.mul_ediv_cancel' hdvd3
  rw [hq2] at hfactor2
  rw [hq3] at hfactor3
  omega

/-- The exact Zarankiewicz value of a nonempty q-two/q-three row
intersection, scaled back to the full branch circle. -/
theorem six_mul_detunedBadBranches_two_three_inter_card_eq
    (δ₂ δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hbad2 : (detunedBadBranches δ₂ g u).Nonempty)
    (hbad3 : (detunedBadBranches δ₃ g u).Nonempty) :
    6 * (detunedBadBranches δ₂ g u ∩
      detunedBadBranches δ₃ g u).card = g.toNat := by
  have hg6 := six_le_of_reducedDenominators_two_three δ₂ δ₃ g hg hq2 hq3
  rw [detunedBadBranches_two_three_inter_card_eq
    δ₂ δ₃ g u hg hg6 hq2 hq3 hbad2 hbad3]
  have hfactor2 := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ₂ g)
  have hfactor3 := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ₃ g)
  rw [hq2] at hfactor2
  rw [hq3] at hfactor3
  omega

/-- The only local obstruction left by a `(2,2,3)` denominator pattern: both
q-two rows are nonempty and occupy opposite parity classes. -/
def TwoTwoThreePhaseOpposition (δ₂a δ₂b g : ℤ) (u : ℝ) : Prop :=
  (detunedBadBranches δ₂a g u).Nonempty ∧
  (detunedBadBranches δ₂b g u).Nonempty ∧
  Disjoint (detunedBadBranches δ₂a g u) (detunedBadBranches δ₂b g u)

/-- Away from opposite q-two parity rows, every `(2,2,3)` phase has a common
good branch.  If the q-two rows meet, they coincide and contribute an exact
half-circle overlap credit; if either is empty, the observed degree sum is
already subcritical. -/
theorem hasThreeDetunedGoodBranch_two_two_three_of_not_opposition
    (δ₂a δ₂b δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hnot : ¬ TwoTwoThreePhaseOpposition δ₂a δ₂b g u) :
    HasThreeDetunedGoodBranch δ₂a δ₂b δ₃ g u := by
  have hg0 : (0 : ℤ) < g := hg
  have hbad2a := two_mul_badCount_eq δ₂a g hg0 hq2a
  have hbad2b := two_mul_badCount_eq δ₂b g hg0 hq2b
  have hbad3 := three_mul_badCount_eq δ₃ g hg0 hq3
  have hcard2a := detunedBadBranches_card_le δ₂a g hg u
  have hcard2b := detunedBadBranches_card_le δ₂b g hg u
  have hcard3 := detunedBadBranches_card_le δ₃ g hg u
  by_cases hrow2a : (detunedBadBranches δ₂a g u).Nonempty
  · by_cases hrow2b : (detunedBadBranches δ₂b g u).Nonempty
    · by_cases hoverlap :
        (detunedBadBranches δ₂a g u ∩
          detunedBadBranches δ₂b g u).Nonempty
      · have hroweq := detunedBadBranches_eq_of_overlap_same_reducedDenominator
          δ₂a δ₂b g 2 u hg (by norm_num) hq2a hq2b hoverlap
        have hinter : detunedBadBranches δ₂a g u ∩
            detunedBadBranches δ₂b g u = detunedBadBranches δ₂a g u := by
          rw [hroweq]
          simp
        have hintercard := two_mul_detunedBadBranches_card_eq_of_nonempty
          δ₂a g u hg hq2a hrow2a
        apply hasThreeDetunedGoodBranch_of_overlapDebt δ₂a δ₂b δ₃ g u hg
        unfold ThreeDetunedOverlapDebtPaid
        left
        rw [hinter]
        omega
      · exfalso
        apply hnot
        refine ⟨hrow2a, hrow2b, Finset.disjoint_left.mpr ?_⟩
        intro c hc2a hc2b
        exact hoverlap ⟨c, Finset.mem_inter.mpr ⟨hc2a, hc2b⟩⟩
    · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂a δ₂b δ₃ g u
      have hzero : (detunedBadBranches δ₂b g u).card = 0 := by
        rw [Finset.card_eq_zero]
        exact Finset.not_nonempty_iff_eq_empty.mp hrow2b
      omega
  · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂a δ₂b δ₃ g u
    have hzero : (detunedBadBranches δ₂a g u).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hrow2a
    omega

/-- Opposite nonempty q-two rows cover the whole branch circle, independently
of the q-three row. -/
theorem not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
    (δ₂a δ₂b δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hopposite : TwoTwoThreePhaseOpposition δ₂a δ₂b g u) :
    ¬ HasThreeDetunedGoodBranch δ₂a δ₂b δ₃ g u := by
  intro hgood
  obtain ⟨c, hcIco, hc2a, hc2b, -⟩ := hgood
  have hcard2a := two_mul_detunedBadBranches_card_eq_of_nonempty
    δ₂a g u hg hq2a hopposite.1
  have hcard2b := two_mul_detunedBadBranches_card_eq_of_nonempty
    δ₂b g u hg hq2b hopposite.2.1
  have hbranches : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
    rw [Int.card_Ico]
    congr 1
    omega
  have hsubset : detunedBadBranches δ₂a g u ∪
      detunedBadBranches δ₂b g u ⊆ Finset.Ico (0 : ℤ) g :=
    Finset.union_subset
      (detunedBadBranches_subset_Ico δ₂a g u)
      (detunedBadBranches_subset_Ico δ₂b g u)
  have hcover : detunedBadBranches δ₂a g u ∪
      detunedBadBranches δ₂b g u = Finset.Ico (0 : ℤ) g := by
    apply Finset.eq_of_subset_of_card_le hsubset
    rw [hbranches, Finset.card_union_of_disjoint hopposite.2.2]
    omega
  have hcUnion : c ∈ detunedBadBranches δ₂a g u ∪
      detunedBadBranches δ₂b g u := by
    rw [hcover]
    exact hcIco
  rcases Finset.mem_union.mp hcUnion with hc | hc
  · exact hc2a hc
  · exact hc2b hc

/-- Exact local normal form for `(2,2,3)`: failure is equivalent to opposite
q-two parity rows. -/
theorem not_hasThreeDetunedGoodBranch_two_two_three_iff_opposition
    (δ₂a δ₂b δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3) :
    ¬ HasThreeDetunedGoodBranch δ₂a δ₂b δ₃ g u ↔
      TwoTwoThreePhaseOpposition δ₂a δ₂b g u := by
  constructor
  · intro hfail
    by_contra hnot
    exact hfail (hasThreeDetunedGoodBranch_two_two_three_of_not_opposition
      δ₂a δ₂b δ₃ g u hg hq2a hq2b hq3 hnot)
  · exact not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
      δ₂a δ₂b δ₃ g u hg hq2a hq2b

/-- Kernel-checked warning that the `(2,2,3)` phase obstruction is genuine:
at `g=6`, the q-two rows of `3` and `27` cover alternating branch classes at
`u=11/100`, irrespective of the third row.  Taking `δ₃=2` gives the literal
reduced-denominator triple `(2,2,3)`. -/
theorem no_three_detuned_goodBranch_three_twentySeven
    (δ₃ : ℤ) :
    ¬ HasThreeDetunedGoodBranch 3 27 δ₃ 6 ((11 : ℝ) / 100) := by
  rintro ⟨c, hcIco, hc1, hc2, -⟩
  have hcBounds := Finset.mem_Ico.mp hcIco
  have hc0 : 0 ≤ c := hcBounds.1
  have hc6 : c < 6 := hcBounds.2
  interval_cases c
  · apply hc1
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 0, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc2
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 5, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc1
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 1, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc2
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 14, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc1
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 2, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]
  · apply hc2
    rw [detunedBadBranches, Finset.mem_filter]
    refine ⟨by norm_num, 23, ?_⟩
    norm_num [abs_of_nonneg, abs_of_nonpos]

theorem twoTwoThree_phaseOpposition_counterexample :
    TwoTwoThreePhaseOpposition 3 27 6 ((11 : ℝ) / 100) :=
  (not_hasThreeDetunedGoodBranch_two_two_three_iff_opposition
    3 27 2 6 ((11 : ℝ) / 100) (by norm_num) (by norm_num)
      (by norm_num) (by norm_num)).mp
    (no_three_detuned_goodBranch_three_twentySeven 2)

/-- Global `(2,2,3)` clearing is now reduced to the exact phase-selection
obligation that the two q-two rows are not opposite at every harmonic phase. -/
theorem threeDetunedInstanceClearing_two_two_three_of_no_opposition
    (δ₂a δ₂b δ₃ g : ℤ) (hg : 1 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hphase : ∀ u : ℝ, ¬ TwoTwoThreePhaseOpposition δ₂a δ₂b g u) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂a δ₂b δ₃ g := by
  intro u
  exact (hasThreeDetunedGoodBranch_two_two_three_of_not_opposition
    δ₂a δ₂b δ₃ g u hg hq2a hq2b hq3 (hphase u)).clearances

theorem lonely14_of_three_detuned_two_two_three_of_no_opposition
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂a i₂b i₃ : Fin 13)
    (h22 : i₂a ≠ i₂b) (h2a3 : i₂a ≠ i₃) (h2b3 : i₂b ≠ i₃)
    (hdvd : ∀ j, j ≠ i₂a → j ≠ i₂b → j ≠ i₃ → g ∣ v j)
    (hq2a : g / (Int.gcd (v i₂a) g : ℤ) = 2)
    (hq2b : g / (Int.gcd (v i₂b) g : ℤ) = 2)
    (hq3 : g / (Int.gcd (v i₃) g : ℤ) = 3)
    (hphase : ∀ u : ℝ,
      ¬ TwoTwoThreePhaseOpposition (v i₂a) (v i₂b) g u) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂a i₂b i₃ h22 h2a3 h2b3 hdvd
    (threeDetunedInstanceClearing_two_two_three_of_no_opposition
      (v i₂a) (v i₂b) (v i₃) g (by omega) hq2a hq2b hq3 hphase)

/-- Two q-three rows cannot saturate together with one q-two row: the q-two
row earns an exact `g/6` overlap credit from each of them. -/
theorem threeDetunedInstanceClearing_two_three_three
    (δ₂ δ₃a δ₃b g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3a : g / (Int.gcd δ₃a g : ℤ) = 3)
    (hq3b : g / (Int.gcd δ₃b g : ℤ) = 3) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₃a δ₃b g := by
  have hg0 : (0 : ℤ) < g := hg
  have hbad2 := two_mul_badCount_eq δ₂ g hg0 hq2
  have hbad3a := three_mul_badCount_eq δ₃a g hg0 hq3a
  have hbad3b := three_mul_badCount_eq δ₃b g hg0 hq3b
  intro u
  have hcard2 := detunedBadBranches_card_le δ₂ g hg u
  have hcard3a := detunedBadBranches_card_le δ₃a g hg u
  have hcard3b := detunedBadBranches_card_le δ₃b g hg u
  have hgood : HasThreeDetunedGoodBranch δ₂ δ₃a δ₃b g u := by
    by_cases hrow2 : (detunedBadBranches δ₂ g u).Nonempty
    · by_cases hrow3a : (detunedBadBranches δ₃a g u).Nonempty
      · by_cases hrow3b : (detunedBadBranches δ₃b g u).Nonempty
        · have hinter3a := six_mul_detunedBadBranches_two_three_inter_card_eq
            δ₂ δ₃a g u hg hq2 hq3a hrow2 hrow3a
          have hinter3b := six_mul_detunedBadBranches_two_three_inter_card_eq
            δ₂ δ₃b g u hg hq2 hq3b hrow2 hrow3b
          have hbranches : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
            rw [Int.card_Ico]
            congr 1
            omega
          apply exists_outside_three_of_two_first_overlapDebt
            (Finset.Ico (0 : ℤ) g)
            (detunedBadBranches δ₂ g u)
            (detunedBadBranches δ₃a g u)
            (detunedBadBranches δ₃b g u)
          rw [hbranches]
          omega
        · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₃a δ₃b g u
          have hzero : (detunedBadBranches δ₃b g u).card = 0 := by
            rw [Finset.card_eq_zero]
            exact Finset.not_nonempty_iff_eq_empty.mp hrow3b
          omega
      · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₃a δ₃b g u
        have hzero : (detunedBadBranches δ₃a g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrow3a
        omega
    · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₃a δ₃b g u
      have hzero : (detunedBadBranches δ₂ g u).card = 0 := by
        rw [Finset.card_eq_zero]
        exact Finset.not_nonempty_iff_eq_empty.mp hrow2
      omega
  exact hgood.clearances

theorem lonely14_of_three_detuned_two_three_three
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₃a i₃b : Fin 13)
    (h23a : i₂ ≠ i₃a) (h23b : i₂ ≠ i₃b) (h33 : i₃a ≠ i₃b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₃a → j ≠ i₃b → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq3a : g / (Int.gcd (v i₃a) g : ℤ) = 3)
    (hq3b : g / (Int.gcd (v i₃b) g : ℤ) = 3) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ i₃a i₃b h23a h23b h33 hdvd
    (threeDetunedInstanceClearing_two_three_three
      (v i₂) (v i₃a) (v i₃b) g (by omega) hq2 hq3a hq3b)

/-- Full-credit CRT closure: once q-two and q-three rows are present, every
third row of bad degree strictly below `g/3` is discharged. -/
theorem threeDetunedInstanceClearing_two_three_of_three_mul_badCount_lt
    (δ₂ δ₃ δₓ g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hsmall : 3 * DetunedD3.badCount δₓ g < g.toNat) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₃ δₓ g := by
  have hg0 : (0 : ℤ) < g := hg
  have hbad2 := two_mul_badCount_eq δ₂ g hg0 hq2
  have hbad3 := three_mul_badCount_eq δ₃ g hg0 hq3
  intro u
  have hcard2 := detunedBadBranches_card_le δ₂ g hg u
  have hcard3 := detunedBadBranches_card_le δ₃ g hg u
  have hcardx := detunedBadBranches_card_le δₓ g hg u
  have hgood : HasThreeDetunedGoodBranch δ₂ δ₃ δₓ g u := by
    by_cases hrow2 : (detunedBadBranches δ₂ g u).Nonempty
    · by_cases hrow3 : (detunedBadBranches δ₃ g u).Nonempty
      · apply hasThreeDetunedGoodBranch_of_overlapDebt δ₂ δ₃ δₓ g u hg
        unfold ThreeDetunedOverlapDebtPaid
        left
        have hinter := six_mul_detunedBadBranches_two_three_inter_card_eq
          δ₂ δ₃ g u hg hq2 hq3 hrow2 hrow3
        omega
      · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₃ δₓ g u
        have hzero : (detunedBadBranches δ₃ g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrow3
        omega
    · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂ δ₃ δₓ g u
      have hzero : (detunedBadBranches δ₂ g u).card = 0 := by
        rw [Finset.card_eq_zero]
        exact Finset.not_nonempty_iff_eq_empty.mp hrow2
      omega
  exact hgood.clearances

/-- Every third reduced denominator at least four satisfies the strict
`g/3` row-degree hypothesis. -/
theorem threeDetunedInstanceClearing_two_three_four_le
    (δ₂ δ₃ δₓ g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hq4 : 4 ≤ g / (Int.gcd δₓ g : ℤ)) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₃ δₓ g :=
  threeDetunedInstanceClearing_two_three_of_three_mul_badCount_lt
    δ₂ δ₃ δₓ g hg hq2 hq3 (by
      have hbound := seven_mul_badCount_le_two_mul δₓ g (by omega) hq4
      omega)

/-- General CRT overlap closure: a q-two row and a q-three row leave a
common good branch whenever the third row's universal bad degree is at most
`g/6`.  The exact `(2,3,6)` pattern is the first equality case. -/
theorem threeDetunedInstanceClearing_two_three_of_six_mul_badCount_le
    (δ₂ δ₃ δₓ g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hsmall : 6 * DetunedD3.badCount δₓ g ≤ g.toNat) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₃ δₓ g :=
  threeDetunedInstanceClearing_two_three_of_three_mul_badCount_lt
    δ₂ δ₃ δₓ g hg hq2 hq3 (by omega)

theorem threeDetunedInstanceClearing_two_three_six
    (δ₂ δ₃ δ₆ g : ℤ) (hg : 1 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hq6 : g / (Int.gcd δ₆ g : ℤ) = 6) :
    DetunedD3.ThreeDetunedInstanceClearing δ₂ δ₃ δ₆ g :=
  threeDetunedInstanceClearing_two_three_of_six_mul_badCount_le
    δ₂ δ₃ δ₆ g hg hq2 hq3 (by
      have hbad6 := six_mul_badCount_eq δ₆ g (by omega) hq6
      omega)

/-- Actual LRC consumer for the full q-two/q-three overlap credit. -/
theorem lonely14_of_three_detuned_two_three_of_three_mul_badCount_lt
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₃ iₓ : Fin 13) (h23 : i₂ ≠ i₃) (h2x : i₂ ≠ iₓ) (h3x : i₃ ≠ iₓ)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₃ → j ≠ iₓ → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq3 : g / (Int.gcd (v i₃) g : ℤ) = 3)
    (hsmall : 3 * DetunedD3.badCount (v iₓ) g < g.toNat) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ i₃ iₓ h23 h2x h3x hdvd
    (threeDetunedInstanceClearing_two_three_of_three_mul_badCount_lt
      (v i₂) (v i₃) (v iₓ) g (by omega) hq2 hq3 hsmall)

theorem lonely14_of_three_detuned_two_three_four_le
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₃ iₓ : Fin 13) (h23 : i₂ ≠ i₃) (h2x : i₂ ≠ iₓ) (h3x : i₃ ≠ iₓ)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₃ → j ≠ iₓ → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq3 : g / (Int.gcd (v i₃) g : ℤ) = 3)
    (hq4 : 4 ≤ g / (Int.gcd (v iₓ) g : ℤ)) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ i₃ iₓ h23 h2x h3x hdvd
    (threeDetunedInstanceClearing_two_three_four_le
      (v i₂) (v i₃) (v iₓ) g (by omega) hq2 hq3 hq4)

/-- Actual LRC consumer for the full q-two/q-three CRT overlap range. -/
theorem lonely14_of_three_detuned_two_three_of_six_mul_badCount_le
    (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₃ iₓ : Fin 13) (h23 : i₂ ≠ i₃) (h2x : i₂ ≠ iₓ) (h3x : i₃ ≠ iₓ)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₃ → j ≠ iₓ → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq3 : g / (Int.gcd (v i₃) g : ℤ) = 3)
    (hsmall : 6 * DetunedD3.badCount (v iₓ) g ≤ g.toNat) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₂ i₃ iₓ h23 h2x h3x hdvd
    (threeDetunedInstanceClearing_two_three_of_six_mul_badCount_le
      (v i₂) (v i₃) (v iₓ) g (by omega) hq2 hq3 hsmall)

theorem lonely14_of_three_detuned_two_three_six (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₃ i₆ : Fin 13) (h23 : i₂ ≠ i₃) (h26 : i₂ ≠ i₆) (h36 : i₃ ≠ i₆)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₃ → j ≠ i₆ → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq3 : g / (Int.gcd (v i₃) g : ℤ) = 3)
    (hq6 : g / (Int.gcd (v i₆) g : ℤ) = 6) :
    ∃ t : ℝ, Lonely 14 v t := by
  apply lonely14_of_three_detuned_two_three_of_six_mul_badCount_le
    cite v hv g hg i₂ i₃ i₆ h23 h26 h36 hdvd hq2 hq3
  have hbad6 := six_mul_badCount_eq (v i₆) g (by omega) hq6
  omega

/-- A selected q-two row and q-three row together with a third detuned row
closed either by the strict `g/3` bad-degree budget or by the exact q-three
double-overlap argument. -/
def HasTwoThreeCRTClosedPattern (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i₂ i₃ iₓ : Fin 13,
    i₂ ≠ i₃ ∧ i₂ ≠ iₓ ∧ i₃ ≠ iₓ ∧
    ¬ g ∣ v i₂ ∧ ¬ g ∣ v i₃ ∧ ¬ g ∣ v iₓ ∧
    g / (Int.gcd (v i₂) g : ℤ) = 2 ∧
    g / (Int.gcd (v i₃) g : ℤ) = 3 ∧
    (3 * DetunedD3.badCount (v iₓ) g < g.toNat ∨
      g / (Int.gcd (v iₓ) g : ℤ) = 3)

/-- The prior q-two/q-three dispatch with every CRT-small-third pattern
removed. -/
def DeepExceptionalDetunedDispatchAfterTwoThreeCRT : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    ((nonMultCard v g = 2 ∧ ¬ TwoAdicLiftTerminates v g) ∨
      (nonMultCard v g = 3 ∧
        (HasTwoWithSubNineCompanion v g ∨
          AllDetuningDenominatorsThree v g) ∧
        ¬ HasTwoThreeCRTClosedPattern v g)) →
    ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

/-- The CRT-small-third theorem supplies the previous exceptional dispatch;
all remaining cases are delegated to the strictly narrower interface. -/
theorem deepExceptionalDetunedDispatchTwoThree_of_afterCRT
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchAfterTwoThreeCRT) :
    DeepExceptionalDetunedDispatchTwoThree := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with hpair | ⟨hcard, hpattern⟩
  · exact hdeep v hv g hg (Or.inl hpair) hnongeneric
  · by_cases hcrt : HasTwoThreeCRTClosedPattern v g
    · obtain ⟨i₂, i₃, iₓ, h23, h2x, h3x, hδ2, hδ3, hδx,
        hq2, hq3, hclosed⟩ := hcrt
      have hselectedCard :
          ({i₂, i₃, iₓ} : Finset (Fin 13)).card = 3 :=
        Finset.card_eq_three.mpr ⟨i₂, i₃, iₓ, h23, h2x, h3x, rfl⟩
      have hsubset :
          ({i₂, i₃, iₓ} : Finset (Fin 13)) ⊆
            Finset.univ.filter (fun i => ¬ g ∣ v i) := by
        intro i hi
        simp only [Finset.mem_insert, Finset.mem_singleton] at hi
        simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        rcases hi with rfl | rfl | rfl
        · exact hδ2
        · exact hδ3
        · exact hδx
      have hfilterCard :
          (Finset.univ.filter (fun i => ¬ g ∣ v i)).card = 3 := by
        simpa [nonMultCard] using hcard
      have hselected :
          ({i₂, i₃, iₓ} : Finset (Fin 13)) =
            Finset.univ.filter (fun i => ¬ g ∣ v i) := by
        apply Finset.eq_of_subset_of_card_le hsubset
        rw [hfilterCard, hselectedCard]
      have hdvd : ∀ j, j ≠ i₂ → j ≠ i₃ → j ≠ iₓ → g ∣ v j := by
        intro j hj2 hj3 hjx
        by_contra hj
        have hjFilter : j ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          simp [hj]
        have hjSelected : j ∈ ({i₂, i₃, iₓ} : Finset (Fin 13)) := by
          rw [hselected]
          exact hjFilter
        simp only [Finset.mem_insert, Finset.mem_singleton] at hjSelected
        rcases hjSelected with h | h | h
        · exact hj2 h
        · exact hj3 h
        · exact hjx h
      rcases hclosed with hsmall | hq3x
      · exact lonely14_of_three_detuned_two_three_of_three_mul_badCount_lt
          cite v hv g hg i₂ i₃ iₓ h23 h2x h3x hdvd hq2 hq3 hsmall
      · exact lonely14_of_three_detuned_two_three_three
          cite v hv g hg i₂ i₃ iₓ h23 h2x h3x hdvd hq2 hq3 hq3x
    · exact hdeep v hv g hg
        (Or.inr ⟨hcard, hpattern, hcrt⟩) hnongeneric

/-- Sharpest relation-budget capstone after the first q-two phase closure. -/
theorem lrc14_from_afterTwoThreeCRT_and_relationBudget
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchAfterTwoThreeCRT)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_relationBudget cite
    (deepExceptionalDetunedDispatchTwoThree_of_afterCRT cite hdeep) hsupply

/-! ## Axiom audit -/

#print axioms detunedBadBranches_pair_modEq_of_q_le_seven
#print axioms detunedBadBranches_mem_of_modEq_reducedDenominator
#print axioms detunedBadBranches_eq_of_overlap_same_reducedDenominator
#print axioms two_mul_detunedBadBranches_card_eq_of_nonempty
#print axioms detunedBadBranches_two_three_overlap
#print axioms detunedBadBranches_two_three_inter_eq_modEq
#print axioms detunedBadBranches_two_three_inter_card_eq
#print axioms six_mul_detunedBadBranches_two_three_inter_card_eq
#print axioms not_hasThreeDetunedGoodBranch_two_two_three_iff_opposition
#print axioms no_three_detuned_goodBranch_three_twentySeven
#print axioms twoTwoThree_phaseOpposition_counterexample
#print axioms threeDetunedInstanceClearing_two_two_three_of_no_opposition
#print axioms lonely14_of_three_detuned_two_two_three_of_no_opposition
#print axioms threeDetunedInstanceClearing_two_three_three
#print axioms threeDetunedInstanceClearing_two_three_of_three_mul_badCount_lt
#print axioms threeDetunedInstanceClearing_two_three_four_le
#print axioms threeDetunedInstanceClearing_two_three_of_six_mul_badCount_le
#print axioms threeDetunedInstanceClearing_two_three_six
#print axioms lonely14_of_three_detuned_two_three_of_three_mul_badCount_lt
#print axioms lonely14_of_three_detuned_two_three_four_le
#print axioms lonely14_of_three_detuned_two_three_of_six_mul_badCount_le
#print axioms lonely14_of_three_detuned_two_three_six
#print axioms lonely14_of_three_detuned_two_three_three
#print axioms deepExceptionalDetunedDispatchTwoThree_of_afterCRT
#print axioms lrc14_from_afterTwoThreeCRT_and_relationBudget

end
end LRC14Grand
end LonelyRunner
