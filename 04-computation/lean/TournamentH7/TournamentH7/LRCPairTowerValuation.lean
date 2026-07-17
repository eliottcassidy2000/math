/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex (LRC multi-agent project, 2026-07-17-S46)
-/
import TournamentH7.LRCSelectedWitnessResidual
import TournamentH7.LRCPairTowerReduction

/-!
# Finite detuned unions and the first valuation-gap closure of the pair tower

The `d = 2` and `d = 3` detuned dispatches use the same finite-union
argument.  This module records that argument for an arbitrary finite set of
detuned coordinates.  It then combines the branch lemma with `LRCUpTo13` on
the complementary quotient family.  This is the reusable kernel-pure
consumer needed by repeated two-adic lifting; it does not freeze the number
of detuned coordinates at two or three.

The intended pair-tower application is the following exact debt calculation.
At the lift `G = 8g`, two harmonic multipliers of minimal two-adic valuation
have reduced denominator eight, while the original half-harmonic pair has
reduced denominator sixteen.  Their four bad-count fractions sum to

`2 * (2/8) + 2 * (3/16) = 7/8 < 1`.

The final theorem below is stated directly at the lifted modulus.  A separate
valuation-arithmetic adapter can supply its four reduced-denominator
hypotheses from a pair tower with exactly two minimal multipliers and a gap of
at least three.

Tournament-analysis audit: vertices are detuned proof obligations and branch
classes, with incidence `c` belongs to the bad row of `i`.  The relation is a
bipartite hypergraph, not a tournament.  Orienting obligation pairs would
discard the common branch which the proof must produce, so no arbitrary gauge
or Hamiltonian path is used.  The challenged assumption is that the detuned
set must have fixed cardinality; the finite-union lemma shows that only its
total exact bad-count debt matters.

No `sorry`; no `native_decide`.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- The union of the bad branch rows indexed by an arbitrary finite set of
detuned coordinates. -/
def detunedBadUnion (v : Fin 13 → ℤ) (g : ℤ) (u : ℝ)
    (D : Finset (Fin 13)) : Finset ℤ :=
  D.biUnion fun i => detunedBadBranches (v i) g u

/-- The arbitrary-finite-set form of the THM-678 branch union bound. -/
theorem detunedBadUnion_card_le (v : Fin 13 → ℤ) (g : ℤ) (u : ℝ)
    (D : Finset (Fin 13)) (hg : 1 ≤ g) :
    (detunedBadUnion v g u D).card ≤
      ∑ i ∈ D, DetunedD3.badCount (v i) g := by
  calc
    (detunedBadUnion v g u D).card ≤
        ∑ i ∈ D, (detunedBadBranches (v i) g u).card := by
          exact Finset.card_biUnion_le
    _ ≤ ∑ i ∈ D, DetunedD3.badCount (v i) g := by
      exact Finset.sum_le_sum fun i _ =>
        detunedBadBranches_card_le (v i) g u hg

/-- If the total bad-count debt of a finite detuned set is strictly below the
number of branch classes, one branch clears every coordinate in the set. -/
theorem exists_goodBranch_of_badCount_sum_lt
    (v : Fin 13 → ℤ) (g : ℤ) (u : ℝ) (D : Finset (Fin 13))
    (hg : 1 ≤ g)
    (hdebt : ∑ i ∈ D, DetunedD3.badCount (v i) g < g.toNat) :
    ∃ c ∈ Finset.Ico (0 : ℤ) g,
      ∀ i ∈ D, c ∉ detunedBadBranches (v i) g u := by
  have hbranches : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
    rw [Int.card_Ico]
    congr 1
    omega
  have hunion_lt : (detunedBadUnion v g u D).card <
      (Finset.Ico (0 : ℤ) g).card := by
    calc
      (detunedBadUnion v g u D).card ≤
          ∑ i ∈ D, DetunedD3.badCount (v i) g :=
        detunedBadUnion_card_le v g u D hg
      _ < g.toNat := hdebt
      _ = (Finset.Ico (0 : ℤ) g).card := hbranches.symm
  have hnotsub : ¬ Finset.Ico (0 : ℤ) g ⊆ detunedBadUnion v g u D := by
    intro hsub
    exact (Nat.not_le_of_lt hunion_lt) (Finset.card_le_card hsub)
  rw [Finset.not_subset] at hnotsub
  obtain ⟨c, hcIco, hcUnion⟩ := hnotsub
  refine ⟨c, hcIco, ?_⟩
  intro i hiD hbad
  apply hcUnion
  exact Finset.mem_biUnion.mpr ⟨i, hiD, hbad⟩

/-- Nonmembership in every bad row gives the simultaneous clearance
statements consumed by the LRC construction. -/
theorem exists_detunedClearances_of_badCount_sum_lt
    (v : Fin 13 → ℤ) (g : ℤ) (u : ℝ) (D : Finset (Fin 13))
    (hg : 1 ≤ g)
    (hdebt : ∑ i ∈ D, DetunedD3.badCount (v i) g < g.toNat) :
    ∃ c : ℤ, ∀ i ∈ D, ∀ n : ℤ,
      (1 : ℝ) / 14 ≤
        |(v i : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| := by
  obtain ⟨c, hcIco, hcgood⟩ :=
    exists_goodBranch_of_badCount_sum_lt v g u D hg hdebt
  refine ⟨c, ?_⟩
  intro i hiD n
  exact not_lt.mp fun hlt => hcgood i hiD (by
    rw [detunedBadBranches, Finset.mem_filter]
    exact ⟨hcIco, n, hlt⟩)

/-- General finite-detuned-set consumer.  The LRC citation is applied only to
the quotient coordinates outside `D`; the finite union bound clears all rows
inside `D` at one branch. -/
theorem lonely14_of_detunedFinset_badCount (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (D : Finset (Fin 13)) (hD : D.Nonempty)
    (hdvd : ∀ j, j ∉ D → g ∣ v j)
    (hdebt : ∑ i ∈ D, DetunedD3.badCount (v i) g < g.toNat) :
    ∃ t : ℝ, Lonely 14 v t := by
  let H : Finset (Fin 13) := Finset.univ \ D
  have hHcard : H.card = 13 - D.card := by
    simp [H, Finset.card_sdiff]
  have hDpos : 0 < D.card := hD.card_pos
  have hHle12 : H.card ≤ 12 := by
    rw [hHcard]
    omega
  let emb : Fin H.card ↪ Fin 13 := H.orderEmbOfFin rfl
  have hemb_mem : ∀ k, emb k ∈ H := fun k => H.orderEmbOfFin_mem rfl k
  let w : Fin H.card → ℤ := fun k => v (emb k) / g
  have hwne : ∀ k, w k ≠ 0 := by
    intro k hk
    have hkH := hemb_mem k
    have hkD : emb k ∉ D := by
      simpa [H] using hkH
    have hdiv := hdvd (emb k) hkD
    apply hv (emb k)
    have hfactor : v (emb k) = g * w k := by
      exact (Int.mul_ediv_cancel' hdiv).symm
    rw [hfactor, hk, mul_zero]
  obtain ⟨u, hu⟩ := cite H.card hHle12 w hwne
  obtain ⟨c, hc⟩ := exists_detunedClearances_of_badCount_sum_lt
    v g u D (by omega) hdebt
  refine ⟨(u + (c : ℝ)) / (g : ℝ), fun i n => ?_⟩
  by_cases hiD : i ∈ D
  · exact hc i hiD n
  · have hiH : i ∈ H := by
      simp [H, hiD]
    obtain ⟨k, hk⟩ : ∃ k, emb k = i := by
      have hrange : i ∈ Set.range emb := by
        rw [Finset.range_orderEmbOfFin]
        exact hiH
      exact hrange
    have hdiv := hdvd i hiD
    have hvi : (v i : ℝ) = (g : ℝ) * ((v i / g : ℤ) : ℝ) := by
      have : v i = g * (v i / g) := (Int.mul_ediv_cancel' hdiv).symm
      exact_mod_cast this
    have hval : (v i : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n =
        ((v i / g : ℤ) : ℝ) * u -
          (((n - (v i / g) * c : ℤ)) : ℝ) := by
      rw [hvi]
      push_cast
      field_simp
      ring
    rw [hval]
    have hwk : w k = v i / g := by
      show v (emb k) / g = v i / g
      rw [hk]
    have hcite := hu k (n - (v i / g) * c)
    rw [hwk] at hcite
    have hden : (H.card + 1 : ℕ) ≤ 14 := by omega
    have hpos : (0 : ℝ) < (H.card : ℝ) + 1 := by positivity
    have hthreshold : (1 : ℝ) / 14 ≤ 1 / ((H.card : ℝ) + 1) := by
      apply one_div_le_one_div_of_le hpos
      exact_mod_cast hden
    have hcast : ((H.card + 1 : ℕ) : ℝ) = (H.card : ℝ) + 1 := by
      push_cast
      ring
    rw [hcast] at hcite
    exact le_trans hthreshold hcite

/-! ## Exact valuation debt -/

/-- Universal bad-row fraction for a coordinate of dyadic reduced
denominator `2^r`.  Since powers of two are not divisible by seven, the
numerator is exactly `ceil(2^r/7)`. -/
def dyadicRowDebt (r : ℕ) : ℚ :=
  ((2 ^ r / 7 + 1 : ℕ) : ℚ) / (2 ^ r : ℕ)

/-- Exact union-bound debt at lift height `k`, expressed through the
histogram `h a` of harmonic multipliers having two-adic valuation `a`.

The original half-harmonic pair contributes two rows of denominator
`2^(k+1)`.  A harmonic multiplier of valuation `a < k` becomes a detuned row
of denominator `2^(k-a)`; multipliers of valuation at least `k` remain in the
harmonic quotient family. -/
def pairTowerValuationDebt (h : ℕ → ℕ) (k : ℕ) : ℚ :=
  2 * dyadicRowDebt (k + 1) +
    ∑ a ∈ Finset.range k, h a * dyadicRowDebt (k - a)

/-- A unique minimal multiplier reaches the selected `(2,4,4)` equality
wall at the first lift. -/
theorem unique_min_first_lift_debt :
    2 * dyadicRowDebt 2 + dyadicRowDebt 1 = 1 := by
  norm_num [dyadicRowDebt]

/-- With a valuation gap of two, the unique-minimum row closes by the plain
union bound at the second lift. -/
theorem unique_min_gap_two_debt :
    2 * dyadicRowDebt 3 + dyadicRowDebt 2 = 3 / 4 := by
  norm_num [dyadicRowDebt]

/-- Two minimal rows remain on a genuine equality wall at the second lift. -/
theorem two_min_second_lift_debt :
    2 * dyadicRowDebt 3 + 2 * dyadicRowDebt 2 = 1 := by
  norm_num [dyadicRowDebt]

/-- A gap of three moves the two-minimum profile below one at the third
lift; this is the `7/8` consumed by the theorem below. -/
theorem two_min_gap_three_debt :
    2 * dyadicRowDebt 4 + 2 * dyadicRowDebt 3 = 7 / 8 := by
  norm_num [dyadicRowDebt]

/-! ## The `(8,8,16,16)` valuation-gap leaf -/

/-- A reduced-denominator-sixteen row has exact THM-678 cost three
sixteenths of the branch classes. -/
theorem sixteen_mul_badCount_eq_three_mul_of_reducedDenominator_sixteen
    (x g : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd x g : ℤ) = 16) :
    16 * DetunedD3.badCount x g = 3 * g.toNat := by
  have hdvd : ((Int.gcd x g : ℤ)) ∣ g := Int.gcd_dvd_right x g
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq] at hfactor
  have htoNat : (g / (Int.gcd x g : ℤ)).toNat = 16 := by
    rw [hq]
    rfl
  rw [DetunedD3.badCount, htoNat]
  norm_num
  omega

/-- Exact valuation-debt arithmetic: `(8,8,16,16)` costs `7/8` of the
branch circle, hence is strictly generic. -/
theorem badCount_sum_eight_eight_sixteen_sixteen_lt
    (x₁ x₂ y₁ y₂ g : ℤ) (hg : 0 < g)
    (hx₁ : g / (Int.gcd x₁ g : ℤ) = 8)
    (hx₂ : g / (Int.gcd x₂ g : ℤ) = 8)
    (hy₁ : g / (Int.gcd y₁ g : ℤ) = 16)
    (hy₂ : g / (Int.gcd y₂ g : ℤ) = 16) :
    DetunedD3.badCount x₁ g + DetunedD3.badCount x₂ g +
      DetunedD3.badCount y₁ g + DetunedD3.badCount y₂ g < g.toNat := by
  have h8₁ := four_mul_badCount_eq_of_reducedDenominator_eight x₁ g hg hx₁
  have h8₂ := four_mul_badCount_eq_of_reducedDenominator_eight x₂ g hg hx₂
  have h16₁ := sixteen_mul_badCount_eq_three_mul_of_reducedDenominator_sixteen
    y₁ g hg hy₁
  have h16₂ := sixteen_mul_badCount_eq_three_mul_of_reducedDenominator_sixteen
    y₂ g hg hy₂
  omega

/-- Scaling two coprime integers by the same positive factor leaves the
reduced denominator unchanged.  This is the small gcd adapter used to turn
valuation-band factorizations into the denominator equalities consumed by
the finite detuned-set theorem. -/
theorem reducedDenominator_mul_of_gcd_eq_one
    (d a q : ℤ) (hd : 0 < d) (hq : 0 < q)
    (hcop : Int.gcd a q = 1) :
    (d * q) / (Int.gcd (d * a) (d * q) : ℤ) = q := by
  have hgcd : Int.gcd (d * a) (d * q) = d.natAbs * Int.gcd a q := by
    rw [Int.gcd_def, Int.natAbs_mul, Int.natAbs_mul, Nat.gcd_mul_left]
    rfl
  have hdabs : (d.natAbs : ℤ) = d := by
    rw [Int.natAbs_of_nonneg (le_of_lt hd)]
    exact Int.toNat_of_nonneg (le_of_lt hd)
  rw [hgcd, hcop, Nat.mul_one]
  norm_cast
  simpa [mul_comm] using Int.mul_ediv_cancel_left q (ne_of_gt hd)

/-- Kernel-pure `d = 4` leaf for the pair-tower valuation gap.  At the lifted
modulus, two rows have denominator eight, the original pair has denominator
sixteen, and every other speed is harmonic. -/
theorem lonely14_of_four_detuned_eight_eight_sixteen_sixteen
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ j₁ j₂ : Fin 13)
    (h12 : i₁ ≠ i₂) (h1j1 : i₁ ≠ j₁) (h1j2 : i₁ ≠ j₂)
    (h2j1 : i₂ ≠ j₁) (h2j2 : i₂ ≠ j₂) (hj12 : j₁ ≠ j₂)
    (hdvd : ∀ k, k ≠ i₁ → k ≠ i₂ → k ≠ j₁ → k ≠ j₂ → g ∣ v k)
    (hi₁ : g / (Int.gcd (v i₁) g : ℤ) = 8)
    (hi₂ : g / (Int.gcd (v i₂) g : ℤ) = 8)
    (hj₁ : g / (Int.gcd (v j₁) g : ℤ) = 16)
    (hj₂ : g / (Int.gcd (v j₂) g : ℤ) = 16) :
    ∃ t : ℝ, Lonely 14 v t := by
  let D : Finset (Fin 13) := {i₁, i₂, j₁, j₂}
  have hD : D.Nonempty := ⟨i₁, by simp [D]⟩
  have hDvd : ∀ k, k ∉ D → g ∣ v k := by
    intro k hk
    simp only [D, Finset.mem_insert, Finset.mem_singleton, not_or] at hk
    exact hdvd k hk.1 hk.2.1 hk.2.2.1 hk.2.2.2
  apply lonely14_of_detunedFinset_badCount cite v hv g hg D hD hDvd
  have hdebt := badCount_sum_eight_eight_sixteen_sixteen_lt
    (v i₁) (v i₂) (v j₁) (v j₂) g (by omega) hi₁ hi₂ hj₁ hj₂
  simpa [D, h12, h1j1, h1j2, h2j1, h2j2, hj12,
    Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hdebt

/-- The promised two-minimum valuation-gap closure in factorized form.

Write the original pair-tower modulus as `2d`.  The two original detuned
speeds are `d*a₁,d*a₂`, with their odd parts coprime to sixteen.  Exactly
two harmonic speeds have minimal valuation, `2d*m₁,2d*m₂`, with the
multipliers coprime to eight; every remaining speed is divisible by `16d`.
At the lift `G=16d=8(2d)` this is precisely the `(8,8,16,16)` leaf above.

The coprimality premises are the sign-insensitive Lean form of the valuation
statements `ν₂(aₛ)=ν₂(mₛ)=0`. -/
theorem lonely14_of_pairTower_two_min_gap_three
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (d : ℤ) (hd : 1 ≤ d)
    (i₁ i₂ j₁ j₂ : Fin 13)
    (h12 : i₁ ≠ i₂) (h1j1 : i₁ ≠ j₁) (h1j2 : i₁ ≠ j₂)
    (h2j1 : i₂ ≠ j₁) (h2j2 : i₂ ≠ j₂) (hj12 : j₁ ≠ j₂)
    (m₁ m₂ a₁ a₂ : ℤ)
    (hi₁ : v i₁ = (2 * d) * m₁) (hi₂ : v i₂ = (2 * d) * m₂)
    (hj₁ : v j₁ = d * a₁) (hj₂ : v j₂ = d * a₂)
    (hm₁ : Int.gcd m₁ 8 = 1) (hm₂ : Int.gcd m₂ 8 = 1)
    (ha₁ : Int.gcd a₁ 16 = 1) (ha₂ : Int.gcd a₂ 16 = 1)
    (hdvd : ∀ k, k ≠ i₁ → k ≠ i₂ → k ≠ j₁ → k ≠ j₂ →
      16 * d ∣ v k) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hd0 : 0 < d := by omega
  have hG : 2 ≤ 16 * d := by omega
  have hq8₁ : (16 * d) / (Int.gcd (v i₁) (16 * d) : ℤ) = 8 := by
    rw [hi₁]
    have h := reducedDenominator_mul_of_gcd_eq_one
      (2 * d) m₁ 8 (by omega) (by norm_num) hm₁
    convert h using 1 <;> ring
  have hq8₂ : (16 * d) / (Int.gcd (v i₂) (16 * d) : ℤ) = 8 := by
    rw [hi₂]
    have h := reducedDenominator_mul_of_gcd_eq_one
      (2 * d) m₂ 8 (by omega) (by norm_num) hm₂
    convert h using 1 <;> ring
  have hq16₁ : (16 * d) / (Int.gcd (v j₁) (16 * d) : ℤ) = 16 := by
    rw [hj₁]
    have h := reducedDenominator_mul_of_gcd_eq_one
      d a₁ 16 hd0 (by norm_num) ha₁
    convert h using 1 <;> ring
  have hq16₂ : (16 * d) / (Int.gcd (v j₂) (16 * d) : ℤ) = 16 := by
    rw [hj₂]
    have h := reducedDenominator_mul_of_gcd_eq_one
      d a₂ 16 hd0 (by norm_num) ha₂
    convert h using 1 <;> ring
  exact lonely14_of_four_detuned_eight_eight_sixteen_sixteen
    cite v hv (16 * d) hG i₁ i₂ j₁ j₂ h12 h1j1 h1j2 h2j1 h2j2 hj12
    hdvd hq8₁ hq8₂ hq16₁ hq16₂

/-- A named factorized leaf of the many-failure pair tower.  The original
modulus is `g = 2d`; two fresh harmonic rows have exact valuation one above
`d`, the original pair has valuation zero above `d`, and every other row is
quiet for three further lifts. -/
structure PairTowerTwoMinGapThreeLeaf (v : Fin 13 → ℤ) (g : ℤ) : Prop where
  d : ℤ
  hd : 1 ≤ d
  g_eq : g = 2 * d
  i₁ i₂ j₁ j₂ : Fin 13
  h12 : i₁ ≠ i₂
  h1j1 : i₁ ≠ j₁
  h1j2 : i₁ ≠ j₂
  h2j1 : i₂ ≠ j₁
  h2j2 : i₂ ≠ j₂
  hj12 : j₁ ≠ j₂
  m₁ m₂ a₁ a₂ : ℤ
  hi₁ : v i₁ = (2 * d) * m₁
  hi₂ : v i₂ = (2 * d) * m₂
  hj₁ : v j₁ = d * a₁
  hj₂ : v j₂ = d * a₂
  hm₁ : Int.gcd m₁ 8 = 1
  hm₂ : Int.gcd m₂ 8 = 1
  ha₁ : Int.gcd a₁ 16 = 1
  ha₂ : Int.gcd a₂ 16 = 1
  hdvd : ∀ k, k ≠ i₁ → k ≠ i₂ → k ≠ j₁ → k ≠ j₂ → 16 * d ∣ v k

/-- The named valuation-gap leaf is an unconditional theorem, not part of the
remaining pair supplier. -/
theorem lonely14_of_pairTowerTwoMinGapThreeLeaf
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ)
    (leaf : PairTowerTwoMinGapThreeLeaf v g) :
    ∃ t : ℝ, Lonely 14 v t :=
  lonely14_of_pairTower_two_min_gap_three cite v hv leaf.d leaf.hd
    leaf.i₁ leaf.i₂ leaf.j₁ leaf.j₂ leaf.h12 leaf.h1j1 leaf.h1j2
    leaf.h2j1 leaf.h2j2 leaf.hj12 leaf.m₁ leaf.m₂ leaf.a₁ leaf.a₂
    leaf.hi₁ leaf.hi₂ leaf.hj₁ leaf.hj₂ leaf.hm₁ leaf.hm₂ leaf.ha₁
    leaf.ha₂ leaf.hdvd

/-- Strictly reduced many-failure supplier: it is asked only after removing
the proved two-minimum, three-level valuation-gap leaf. -/
def ManyLiftFailureBeyondGapThreeSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    nonMultCard v g = 2 → 2 ≤ liftFailureCard v g →
    ¬ genericCount v g → ¬ PairTowerTwoMinGapThreeLeaf v g →
    ∃ t : ℝ, Lonely 14 v t

/-- Reinsert the proved valuation-gap leaf to recover the full many-failure
pair supplier. -/
theorem manyLiftFailurePairTowerSupply_of_beyondGapThree
    (cite : LRCUpTo13) (hremaining : ManyLiftFailureBeyondGapThreeSupply) :
    ManyLiftFailurePairTowerSupply := by
  intro v hv g hg hcard hmany hnongeneric
  by_cases hleaf : PairTowerTwoMinGapThreeLeaf v g
  · exact lonely14_of_pairTowerTwoMinGapThreeLeaf cite v hv g hleaf
  · exact hremaining v hv g hg hcard hmany hnongeneric hleaf

/-- The main capstone with the valuation-gap leaf removed from the pair-tower
hypothesis and THM-940's concrete signed-deviation socket on the dense core. -/
theorem lrc14_from_pairTowerBeyondGapThree_and_selectedWitnessSupplies_and_deviationBudget
    (cite : LRCUpTo13)
    (hremaining : ManyLiftFailureBeyondGapThreeSupply)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreDeviationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_manyLiftFailure_and_selectedWitnessSupplies_and_deviationBudget
    cite (manyLiftFailurePairTowerSupply_of_beyondGapThree cite hremaining)
    h22 h244 h333 hsupply

/-- The singleton-minimum part of the nonterminating pair tower reduces
exactly to the selected `(2,4,4)` supplier.

Write the original modulus as `2d`.  The unique minimal harmonic speed is
`2d*m`, the original half-harmonic pair is `d*a₁,d*a₂`, and all remaining
speeds are divisible by `4d`.  At the first lift `G=4d`, these are exactly one
q-two row and two q-four rows.  Their THM-678 debt is equality, so this theorem
uses the joint selected-witness supplier rather than claiming a phase-uniform
union bound. -/
theorem lonely14_of_pairTower_unique_min_selected
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (d : ℤ) (hd : 1 ≤ d) (i j₁ j₂ : Fin 13)
    (hij₁ : i ≠ j₁) (hij₂ : i ≠ j₂) (hj₁j₂ : j₁ ≠ j₂)
    (m a₁ a₂ : ℤ)
    (hi : v i = (2 * d) * m) (hj₁ : v j₁ = d * a₁)
    (hj₂ : v j₂ = d * a₂)
    (hm : Int.gcd m 2 = 1)
    (ha₁ : Int.gcd a₁ 4 = 1) (ha₂ : Int.gcd a₂ 4 = 1)
    (hdvd : ∀ k, k ≠ i → k ≠ j₁ → k ≠ j₂ → 4 * d ∣ v k) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hd0 : 0 < d := by omega
  have hG : 2 ≤ 4 * d := by omega
  have hq2 : (4 * d) / (Int.gcd (v i) (4 * d) : ℤ) = 2 := by
    rw [hi]
    have h := reducedDenominator_mul_of_gcd_eq_one
      (2 * d) m 2 (by omega) (by norm_num) hm
    convert h using 1 <;> ring
  have hq4₁ : (4 * d) / (Int.gcd (v j₁) (4 * d) : ℤ) = 4 := by
    rw [hj₁]
    have h := reducedDenominator_mul_of_gcd_eq_one
      d a₁ 4 hd0 (by norm_num) ha₁
    convert h using 1 <;> ring
  have hq4₂ : (4 * d) / (Int.gcd (v j₂) (4 * d) : ℤ) = 4 := by
    rw [hj₂]
    have h := reducedDenominator_mul_of_gcd_eq_one
      d a₂ 4 hd0 (by norm_num) ha₂
    convert h using 1 <;> ring
  have hδi : ¬ 4 * d ∣ v i :=
    not_dvd_of_one_lt_reducedDenominator (v i) (4 * d) 2 (by omega) hq2 (by omega)
  have hδ1 : ¬ 4 * d ∣ v j₁ :=
    not_dvd_of_one_lt_reducedDenominator (v j₁) (4 * d) 4 (by omega) hq4₁ (by omega)
  have hδ2 : ¬ 4 * d ∣ v j₂ :=
    not_dvd_of_one_lt_reducedDenominator (v j₂) (4 * d) 4 (by omega) hq4₂ (by omega)
  have hfilter :
      Finset.univ.filter (fun k : Fin 13 => ¬ 4 * d ∣ v k) = {i, j₁, j₂} := by
    ext k
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_insert, Finset.mem_singleton]
    constructor
    · intro hk
      by_cases hki : k = i
      · exact Or.inl hki
      by_cases hkj₁ : k = j₁
      · exact Or.inr (Or.inl hkj₁)
      by_cases hkj₂ : k = j₂
      · exact Or.inr (Or.inr hkj₂)
      exact False.elim (hk (hdvd k hki hkj₁ hkj₂))
    · rintro (rfl | rfl | rfl)
      · exact hδi
      · exact hδ1
      · exact hδ2
  have hnongeneric : ¬ genericCount v (4 * d) := by
    intro hgeneric
    rw [genericCount, hfilter,
      Finset.sum_insert (by simp [hij₁, hij₂]),
      Finset.sum_insert (by simp [hj₁j₂]), Finset.sum_singleton] at hgeneric
    have hb2 := two_mul_badCount_eq (v i) (4 * d) (by omega) hq2
    have hb4₁ := four_mul_badCount_eq_of_reducedDenominator_four
      (v j₁) (4 * d) (by omega) hq4₁
    have hb4₂ := four_mul_badCount_eq_of_reducedDenominator_four
      (v j₂) (4 * d) (by omega) hq4₂
    omega
  apply lonely14_of_three_detuned_selectedWitness
    v (4 * d) hG i j₁ j₂ hdvd
  exact h244 v hv (4 * d) hG hnongeneric i j₁ j₂
    hij₁ hij₂ hj₁j₂ hdvd hq2 hq4₁ hq4₂

/-! ## Axiom audit -/

#print axioms detunedBadUnion_card_le
#print axioms exists_goodBranch_of_badCount_sum_lt
#print axioms exists_detunedClearances_of_badCount_sum_lt
#print axioms lonely14_of_detunedFinset_badCount
#print axioms unique_min_first_lift_debt
#print axioms unique_min_gap_two_debt
#print axioms two_min_second_lift_debt
#print axioms two_min_gap_three_debt
#print axioms badCount_sum_eight_eight_sixteen_sixteen_lt
#print axioms reducedDenominator_mul_of_gcd_eq_one
#print axioms lonely14_of_four_detuned_eight_eight_sixteen_sixteen
#print axioms lonely14_of_pairTower_two_min_gap_three
#print axioms lonely14_of_pairTowerTwoMinGapThreeLeaf
#print axioms manyLiftFailurePairTowerSupply_of_beyondGapThree
#print axioms lrc14_from_pairTowerBeyondGapThree_and_selectedWitnessSupplies_and_deviationBudget
#print axioms lonely14_of_pairTower_unique_min_selected

end
end LRC14Grand
end LonelyRunner
