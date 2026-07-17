/-
  TournamentH7.LRCB5SubsetExpansion — THE DISCRETE HALF OF THE SINGULAR-SERIES
  IDENTIFICATION (death-star-2026-07-16-S35, HYP-7172; discharging the combinatorial
  face of codex-S28's named obligation in LRCB5RelationBudget).

  codex's budget module models the level-five singular series algebraically
  (equilibrium `2052/16807`, signed support-weights) but "intentionally does not
  assert the singular-series identity for the concrete discrete `B5`".  This module
  supplies the identity's DISCRETE half, exactly and kernel-pure:

  * `jointFail` / `momentS_eq_sum_jointFail` — the d-th factorial moment IS the sum
    of joint band-failure counts over d-subsets: `S_d = Σ_{|T|=d} N_T`
    (`choose` ↔ `powersetCard`, then a sum swap).
  * `B5_eq_subset_sum` — `B5 v q = Σ_{k ≤ 5} (−1)^k Σ_{|T|=k} N_T`: the quintic
    Bonferroni functional grouped by SUPPORT, the exact discrete analogue of the
    singular series' support decomposition.
  * `equilibrium_binomial` — `Σ_{d ≤ 5} (−1)^d C(13,d)/7^d = 2052/16807`: codex's
    `equilibrium` constant DERIVED (the truncated binomial at failure rate 1/7).
  * `deviation` / `B5_eq_equilibrium_add_deviation` — the deviation ledger:
    `D_T := N_T − (q−1)/7^|T|` (in ℚ) gives EXACTLY
        `(B5 v q : ℚ) = (q−1)·2052/16807 + Σ_{k ≤ 5} (−1)^k Σ_{|T|=k} D_T`.
    The relation masses of codex's `relationModel` are now DEFINED discrete
    quantities, not hypothesized reals.
  * `B5_pos_of_deviation_debt` — positivity from the deviation debt:
    `|Σ (−1)^k Σ D_T| < (q−1)·2052/16807 → 0 < B5 v q`.

  WHERE THE TRAP ENTERS (THM-939, stated here as the interface comment; the
  quantitative equidistribution of trapped D_T is the sharpened remaining
  obligation): for q beyond the family's magnitude scale, a mod-q coincidence
  pattern of support ≤ 5 whose top runner sits above the dense pair forces an
  integer relation of the same support — which `no_low_mass_relation_above_pair` /
  `no_unit_relation_high` forbid.  The analytic remainder is now: bound the
  CONCRETE ℚ-numbers D_T for trapped subsets.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDiscreteBonferroni

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **The joint band-failure count of a subset `T`**: multipliers `p ∈ (0, q)` at
which every runner in `T` misses the safe band. -/
def jointFail (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) : ℕ :=
  ((Finset.Ioo 0 q).filter fun p => ∀ i ∈ T, ¬ inBand v q p i).card

/-- The `d`-th factorial moment is the sum of joint-failure counts over
`d`-subsets: `S_d = Σ_{|T| = d} N_T`. -/
theorem momentS_eq_sum_jointFail (v : Fin 13 → ℤ) (q : ℕ) (d : ℕ) :
    momentS v q d
      = ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard d, jointFail v q T := by
  unfold momentS jointFail
  have hpoint : ∀ p : ℕ, (bandCount v q p).choose d
      = ((Finset.univ.filter fun i : Fin 13 => ¬ inBand v q p i).powersetCard d).card := by
    intro p
    rw [Finset.card_powersetCard]
    rfl
  have hswap : ∀ p : ℕ,
      ((Finset.univ.filter fun i : Fin 13 => ¬ inBand v q p i).powersetCard d).card
      = (((Finset.univ : Finset (Fin 13)).powersetCard d).filter
          (fun T => ∀ i ∈ T, ¬ inBand v q p i)).card := by
    intro p
    congr 1
    ext T
    constructor
    · intro hT
      obtain ⟨hsub, hcard⟩ := Finset.mem_powersetCard.mp hT
      rw [Finset.mem_filter, Finset.mem_powersetCard]
      exact ⟨⟨Finset.subset_univ T, hcard⟩,
        fun i hi => (Finset.mem_filter.mp (hsub hi)).2⟩
    · intro hT
      rw [Finset.mem_filter, Finset.mem_powersetCard] at hT
      obtain ⟨⟨_, hcard⟩, hfail⟩ := hT
      rw [Finset.mem_powersetCard]
      exact ⟨fun i hi => Finset.mem_filter.mpr ⟨Finset.mem_univ i, hfail i hi⟩, hcard⟩
  calc ∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d
      = ∑ p ∈ Finset.Ioo 0 q,
          (((Finset.univ : Finset (Fin 13)).powersetCard d).filter
            (fun T => ∀ i ∈ T, ¬ inBand v q p i)).card := by
        apply Finset.sum_congr rfl
        intro p _
        rw [hpoint p, hswap p]
    _ = ∑ p ∈ Finset.Ioo 0 q,
          ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard d,
            if (∀ i ∈ T, ¬ inBand v q p i) then 1 else 0 := by
        apply Finset.sum_congr rfl
        intro p _
        rw [Finset.card_filter]
    _ = ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard d,
          ∑ p ∈ Finset.Ioo 0 q,
            if (∀ i ∈ T, ¬ inBand v q p i) then 1 else 0 := Finset.sum_comm
    _ = ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard d,
          ((Finset.Ioo 0 q).filter fun p => ∀ i ∈ T, ¬ inBand v q p i).card := by
        apply Finset.sum_congr rfl
        intro T _
        rw [Finset.card_filter]

/-- **The subset expansion of the quintic Bonferroni functional**: `B5` grouped by
support — the exact discrete analogue of the singular series' support
decomposition. -/
theorem B5_eq_subset_sum (v : Fin 13 → ℤ) (q : ℕ) :
    B5 v q = ∑ k ∈ range 6, (-1 : ℤ) ^ k *
      ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k, (jointFail v q T : ℤ) := by
  unfold B5
  apply Finset.sum_congr rfl
  intro d _
  rw [momentS_eq_sum_jointFail]
  push_cast
  ring

/-- codex's `equilibrium` constant derived: the truncated binomial at failure
rate `1/7`. -/
theorem equilibrium_binomial :
    (∑ d ∈ range 6, (-1 : ℚ) ^ d * (Nat.choose 13 d : ℚ) / 7 ^ d)
      = 2052 / 16807 := by
  simp only [Finset.sum_range_succ, Finset.sum_range_zero,
    show Nat.choose 13 0 = 1 from rfl, show Nat.choose 13 1 = 13 from rfl,
    show Nat.choose 13 2 = 78 from rfl, show Nat.choose 13 3 = 286 from rfl,
    show Nat.choose 13 4 = 715 from rfl, show Nat.choose 13 5 = 1287 from rfl]
  norm_num

/-- **The deviation of a subset**: its joint-failure count against the equilibrium
share `(q−1)/7^|T|`, as an exact rational. -/
def deviation (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) : ℚ :=
  (jointFail v q T : ℚ) - ((q : ℚ) - 1) / 7 ^ T.card

/-- **THE DISCRETE IDENTIFICATION**: the quintic Bonferroni functional equals the
equilibrium share plus the signed deviation ledger, EXACTLY:
`B5 = (q−1)·2052/16807 + Σ_{k≤5} (−1)^k Σ_{|T|=k} D_T`. -/
theorem B5_eq_equilibrium_add_deviation (v : Fin 13 → ℤ) (q : ℕ) :
    (B5 v q : ℚ)
      = ((q : ℚ) - 1) * (2052 / 16807)
        + ∑ k ∈ range 6, (-1 : ℚ) ^ k *
            ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k, deviation v q T := by
  have hcard : ∀ k, ∀ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
      T.card = k := by
    intro k T hT
    exact (Finset.mem_powersetCard.mp hT).2
  have hcount : ∀ k : ℕ,
      ((Finset.univ : Finset (Fin 13)).powersetCard k).card = Nat.choose 13 k := by
    intro k
    rw [Finset.card_powersetCard]
    simp
  have hexp : (B5 v q : ℚ) = ∑ k ∈ range 6, (-1 : ℚ) ^ k *
      ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
        (jointFail v q T : ℚ) := by
    have := B5_eq_subset_sum v q
    have hcast : ((B5 v q : ℤ) : ℚ) = (B5 v q : ℚ) := rfl
    rw [this]
    push_cast
    ring
  rw [hexp]
  have hsplit : ∀ k ∈ range 6,
      (-1 : ℚ) ^ k * ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
        (jointFail v q T : ℚ)
      = (-1 : ℚ) ^ k * ((Nat.choose 13 k : ℚ) * (((q : ℚ) - 1) / 7 ^ k))
        + (-1 : ℚ) ^ k * ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
            deviation v q T := by
    intro k _
    have hterm : ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
        (jointFail v q T : ℚ)
        = ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
            (deviation v q T + ((q : ℚ) - 1) / 7 ^ k) := by
      apply Finset.sum_congr rfl
      intro T hT
      unfold deviation
      rw [hcard k T hT]
      ring
    rw [hterm, Finset.sum_add_distrib, Finset.sum_const, hcount, nsmul_eq_mul]
    ring
  rw [Finset.sum_congr rfl hsplit, Finset.sum_add_distrib]
  have hmain : ∑ k ∈ range 6,
      (-1 : ℚ) ^ k * ((Nat.choose 13 k : ℚ) * (((q : ℚ) - 1) / 7 ^ k))
      = ((q : ℚ) - 1) * (2052 / 16807) := by
    have hfactor : ∑ k ∈ range 6,
        (-1 : ℚ) ^ k * ((Nat.choose 13 k : ℚ) * (((q : ℚ) - 1) / 7 ^ k))
        = ((q : ℚ) - 1) * ∑ k ∈ range 6,
            (-1 : ℚ) ^ k * (Nat.choose 13 k : ℚ) / 7 ^ k := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro k _
      ring
    rw [hfactor, equilibrium_binomial]
  rw [hmain]

/-- **Positivity from the deviation debt**: if the signed deviation ledger stays
strictly inside the equilibrium share, the quintic certificate is positive.  This is
codex's `relationModel_pos_of_debt_lt` shape with the masses now DEFINED from the
discrete counts. -/
theorem B5_pos_of_deviation_debt (v : Fin 13 → ℤ) (q : ℕ)
    (hdebt : |∑ k ∈ range 6, (-1 : ℚ) ^ k *
        ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k, deviation v q T|
      < ((q : ℚ) - 1) * (2052 / 16807)) :
    0 < B5 v q := by
  have hid := B5_eq_equilibrium_add_deviation v q
  have habs := abs_lt.mp hdebt
  have hQ : (0 : ℚ) < (B5 v q : ℚ) := by
    rw [hid]
    linarith [habs.1]
  exact_mod_cast hQ

/-- The empty subset never deviates: `N_∅ = q − 1` exactly (for `1 ≤ q`). -/
theorem deviation_empty (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 ≤ q) :
    deviation v q ∅ = 0 := by
  unfold deviation jointFail
  have hall : (Finset.Ioo 0 q).filter
      (fun p => ∀ i ∈ (∅ : Finset (Fin 13)), ¬ inBand v q p i) = Finset.Ioo 0 q := by
    apply Finset.filter_true_of_mem
    intro p _ i hi
    exact absurd hi (by simp)
  rw [hall]
  simp only [Nat.card_Ioo, Nat.sub_zero, Finset.card_empty, pow_zero, div_one,
    Nat.cast_sub hq]
  simp

/-! ## Axiom audit -/
#print axioms momentS_eq_sum_jointFail
#print axioms B5_eq_subset_sum
#print axioms B5_eq_equilibrium_add_deviation
#print axioms B5_pos_of_deviation_debt

end LRC14Concrete
end LonelyRunner
