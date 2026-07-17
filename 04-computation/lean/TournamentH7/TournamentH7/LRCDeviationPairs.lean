/-
  TournamentH7.LRCDeviationPairs — THE PAIR RUNG of the deviation ledger
  (death-star-2026-07-17-S37, HYP-7181; the |T| = 2 floor of THM-940's discrete
  identification, extending THM-942A's unit-bijection technique).

  * `jointFail_anti` — N_T ≤ N_S for S ⊆ T (more conditions, fewer multipliers).
  * `jointFail_pair_lower` — the Bonferroni sandwich floor, subtraction-free:
        N_{i} + N_{j} ≤ N_{ij} + (q − 1).
  * `dilatePairCount_eq` — THE RESIDUE-LEVEL EXACT COUNT: at `q = 14Q`, the number of
    nonzero residues `s` with both `s` and `(2s) mod q` outside the safe band is
    EXACTLY `2·((Q−1)/2)` (ℕ-division).  Variable-modulus `%` is eliminated by hand
    (two branch rewrites); everything else is `omega`.
  * `jointFail_dilate_pair_eq` — THE TRANSFER: for `gcd(v i, q) = 1` and
    `v j ≡ 2·v i (mod q)` at `q = 14Q`:  `N_{ij} = 2·((Q−1)/2)`.  Hence
        D_{ij} = 2·((Q−1)/2) − (q−1)/49 = (5/7)·Q + O(1) — POSITIVE and Θ(q):
    the first formally exact quantification of the SYSTEMATIC dilate blocker (joint
    failures at rate ~1/14 against the independence rate 1/49).  kps's THM-934
    stratification ("only systematic structure kills B5") now has both poles formal:
    the trap (THM-939) forbids this shape above the dense pair; this theorem prices
    it exactly where it does occur.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDeviationSingles

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Joint failure is antitone in the subset: more runners to fail, fewer multipliers. -/
theorem jointFail_anti (v : Fin 13 → ℤ) (q : ℕ) {S T : Finset (Fin 13)}
    (hST : S ⊆ T) : jointFail v q T ≤ jointFail v q S := by
  unfold jointFail
  apply Finset.card_le_card
  intro p hp
  rw [Finset.mem_filter] at hp ⊢
  exact ⟨hp.1, fun i hi => hp.2 i (hST hi)⟩

/-- The Bonferroni pair floor, subtraction-free:
`N_{i} + N_{j} ≤ N_{ij} + (q − 1)`. -/
theorem jointFail_pair_lower (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13) :
    jointFail v q {i} + jointFail v q {j} ≤ jointFail v q {i, j} + (q - 1) := by
  unfold jointFail
  set A := (Finset.Ioo 0 q).filter (fun p => ∀ k ∈ ({i} : Finset (Fin 13)),
    ¬ inBand v q p k) with hA
  set B := (Finset.Ioo 0 q).filter (fun p => ∀ k ∈ ({j} : Finset (Fin 13)),
    ¬ inBand v q p k) with hB
  have hinter : A ∩ B ⊆ (Finset.Ioo 0 q).filter (fun p => ∀ k ∈ ({i, j} : Finset (Fin 13)),
      ¬ inBand v q p k) := by
    intro p hp
    rw [Finset.mem_inter, hA, hB, Finset.mem_filter, Finset.mem_filter] at hp
    rw [Finset.mem_filter]
    refine ⟨hp.1.1, fun k hk => ?_⟩
    rcases Finset.mem_insert.mp hk with rfl | hk
    · exact hp.1.2 k (Finset.mem_singleton_self k)
    · rcases Finset.mem_singleton.mp hk with rfl
      exact hp.2.2 k (Finset.mem_singleton_self k)
  have hunion : (A ∪ B).card ≤ q - 1 := by
    calc (A ∪ B).card ≤ (Finset.Ioo 0 q).card := by
          apply Finset.card_le_card
          intro p hp
          rcases Finset.mem_union.mp hp with h | h
          · exact Finset.mem_of_mem_filter p h
          · exact Finset.mem_of_mem_filter p h
      _ = q - 1 := by simp [Nat.card_Ioo]
  have hie := Finset.card_union_add_card_inter A B
  have hAB : (A ∩ B).card ≤ ((Finset.Ioo 0 q).filter (fun p => ∀ k ∈ ({i, j} : Finset (Fin 13)),
      ¬ inBand v q p k)).card := Finset.card_le_card hinter
  omega

/-- **The residue-level dilate-pair count**: at `q = 14Q`, the nonzero residues with
both `s` and `(2s) mod q` outside the safe band number exactly `2·((Q−1)/2)`. -/
theorem dilatePairCount_eq (Q : ℕ) (hQ : 1 ≤ Q) :
    ((Finset.Ioo 0 (14 * Q)).filter fun s =>
      ¬ (14 * Q ≤ 14 * s ∧ 14 * s ≤ 13 * (14 * Q)) ∧
      ¬ (14 * Q ≤ 14 * ((2 * s) % (14 * Q)) ∧
          14 * ((2 * s) % (14 * Q)) ≤ 13 * (14 * Q))).card
      = 2 * ((Q - 1) / 2) := by
  have hchar : (Finset.Ioo 0 (14 * Q)).filter (fun s =>
      ¬ (14 * Q ≤ 14 * s ∧ 14 * s ≤ 13 * (14 * Q)) ∧
      ¬ (14 * Q ≤ 14 * ((2 * s) % (14 * Q)) ∧
          14 * ((2 * s) % (14 * Q)) ≤ 13 * (14 * Q)))
      = Finset.Icc 1 ((Q - 1) / 2) ∪ Finset.Icc (27 * Q / 2 + 1) (14 * Q - 1) := by
    ext s
    simp only [Finset.mem_filter, Finset.mem_Ioo, Finset.mem_union, Finset.mem_Icc]
    constructor
    · rintro ⟨⟨hs0, hsq⟩, hout1, hout2⟩
      by_cases hlow : s < Q
      · left
        have hmod : (2 * s) % (14 * Q) = 2 * s := Nat.mod_eq_of_lt (by omega)
        rw [hmod] at hout2
        omega
      · right
        have hhigh : 13 * Q < s := by omega
        have h14 : 14 * Q ≤ 2 * s := by omega
        have hmod : (2 * s) % (14 * Q) = 2 * s - 14 * Q := by
          rw [Nat.mod_eq_sub_mod h14]
          exact Nat.mod_eq_of_lt (by omega)
        rw [hmod] at hout2
        omega
    · intro h
      rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · have hs0 : 0 < s := by omega
        have hlow : s < Q := by omega
        have hmod : (2 * s) % (14 * Q) = 2 * s := Nat.mod_eq_of_lt (by omega)
        refine ⟨⟨hs0, by omega⟩, by omega, ?_⟩
        rw [hmod]
        omega
      · have hhigh : 13 * Q < s := by omega
        have h14 : 14 * Q ≤ 2 * s := by omega
        have hmod : (2 * s) % (14 * Q) = 2 * s - 14 * Q := by
          rw [Nat.mod_eq_sub_mod h14]
          exact Nat.mod_eq_of_lt (by omega)
        refine ⟨⟨by omega, by omega⟩, by omega, ?_⟩
        rw [hmod]
        omega
  rw [hchar, Finset.card_union_of_disjoint]
  · rw [Nat.card_Icc, Nat.card_Icc]
    omega
  · rw [Finset.disjoint_left]
    intro s hs hs'
    rw [Finset.mem_Icc] at hs hs'
    omega

/-- **The dilate-pair transfer**: for `gcd(v i, q) = 1` and `v j ≡ 2·v i (mod q)` at
`q = 14Q`, the pair joint-failure count is exactly `2·((Q−1)/2)`.  With THM-940's
ledger this makes the dilate pair's deviation `(5/7)Q + O(1)` — positive and Θ(q),
the systematic blocker priced exactly. -/
theorem jointFail_dilate_pair_eq (v : Fin 13 → ℤ) (Q : ℕ) (i j : Fin 13)
    (hQ : 1 ≤ Q) (hij : i ≠ j)
    (hgcd : Int.gcd (v i) ((14 * Q : ℕ) : ℤ) = 1)
    (hdil : ((14 * Q : ℕ) : ℤ) ∣ (v j - 2 * v i)) :
    jointFail v (14 * Q) {i, j} = 2 * ((Q - 1) / 2) := by
  set q : ℕ := 14 * Q with hqdef
  have hq : 0 < q := by omega
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  set f : ℕ → ℕ := fun p => ((v i * (p : ℤ)) % (q : ℤ)).toNat with hf
  have hmod_nonneg : ∀ p : ℕ, (0 : ℤ) ≤ ((v i * (p : ℤ)) % (q : ℤ)) :=
    fun p => Int.emod_nonneg _ (by omega)
  have hmod_lt : ∀ p : ℕ, ((v i * (p : ℤ)) % (q : ℤ)) < (q : ℤ) :=
    fun p => Int.emod_lt_of_pos _ hqZ
  have hfval : ∀ p : ℕ, ((f p : ℕ) : ℤ) = (v i * (p : ℤ)) % (q : ℤ) := by
    intro p
    rw [hf]
    exact Int.toNat_of_nonneg (hmod_nonneg p)
  -- the i-condition transfers as in THM-942A
  have hband_i : ∀ p : ℕ, inBand v q p i ↔ (q ≤ 14 * f p ∧ 14 * f p ≤ 13 * q) := by
    intro p
    unfold inBand
    constructor
    · intro ⟨h1, h2⟩
      constructor
      · have h1' : (q : ℤ) ≤ 14 * ((f p : ℕ) : ℤ) := by rw [hfval]; exact h1
        exact_mod_cast h1'
      · have h2' : 14 * ((f p : ℕ) : ℤ) ≤ 13 * (q : ℤ) := by rw [hfval]; exact h2
        exact_mod_cast h2'
    · intro ⟨h1, h2⟩
      constructor
      · have h1' : (q : ℤ) ≤ 14 * ((f p : ℕ) : ℤ) := by exact_mod_cast h1
        rwa [hfval] at h1'
      · have h2' : 14 * ((f p : ℕ) : ℤ) ≤ 13 * (q : ℤ) := by exact_mod_cast h2
        rwa [hfval] at h2'
  -- the j-condition transfers through the dilate congruence
  have hjmod : ∀ p : ℕ, (v j * (p : ℤ)) % (q : ℤ) = (((2 * f p) % q : ℕ) : ℤ) := by
    intro p
    have hcong : v j * (p : ℤ) ≡ 2 * ((f p : ℕ) : ℤ) [ZMOD (q : ℤ)] := by
      calc v j * (p : ℤ)
          ≡ (2 * v i) * (p : ℤ) [ZMOD (q : ℤ)] := by
            apply Int.ModEq.mul_right
            exact (Int.modEq_iff_dvd.mpr (by simpa using hdil)).symm
        _ = 2 * (v i * (p : ℤ)) := by ring
        _ ≡ 2 * ((f p : ℕ) : ℤ) [ZMOD (q : ℤ)] := by
            apply Int.ModEq.mul_left
            show (v i * (p : ℤ)) % (q : ℤ) = ((f p : ℕ) : ℤ) % (q : ℤ)
            rw [hfval]
            exact (Int.emod_emod_of_dvd _ dvd_rfl).symm
    have hcast : (((2 * f p) % q : ℕ) : ℤ) = (2 * ((f p : ℕ) : ℤ)) % (q : ℤ) := by
      push_cast
      rfl
    rw [hcast]
    exact hcong
  have hband_j : ∀ p : ℕ, inBand v q p j ↔
      (q ≤ 14 * ((2 * f p) % q) ∧ 14 * ((2 * f p) % q) ≤ 13 * q) := by
    intro p
    unfold inBand
    rw [hjmod p]
    constructor
    · intro ⟨h1, h2⟩
      exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩
    · intro ⟨h1, h2⟩
      exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩
  -- rewrite the joint filter through the pair of residue conditions
  have hjf : jointFail v q {i, j}
      = ((Finset.Ioo 0 q).filter fun p =>
          ¬ (q ≤ 14 * f p ∧ 14 * f p ≤ 13 * q) ∧
          ¬ (q ≤ 14 * ((2 * f p) % q) ∧ 14 * ((2 * f p) % q) ≤ 13 * q)).card := by
    unfold jointFail
    congr 1
    apply Finset.filter_congr
    intro p _
    constructor
    · intro h
      exact ⟨fun hcon => (h i (Finset.mem_insert_self i {j})) ((hband_i p).mpr hcon),
        fun hcon => (h j (Finset.mem_insert_of_mem (Finset.mem_singleton_self j)))
          ((hband_j p).mpr hcon)⟩
    · intro ⟨h1, h2⟩ k hk
      rcases Finset.mem_insert.mp hk with rfl | hk
      · exact fun hcon => h1 ((hband_i p).mp hcon)
      · rcases Finset.mem_singleton.mp hk with rfl
        exact fun hcon => h2 ((hband_j p).mp hcon)
  rw [hjf, ← dilatePairCount_eq Q hQ]
  -- the bijection p ↦ f p between the multiplier filter and the residue filter
  obtain ⟨a, b, hab⟩ := Int.isCoprime_iff_gcd_eq_one.mpr hgcd
  have hva : v i * a ≡ 1 [ZMOD (q : ℤ)] :=
    Int.modEq_iff_dvd.mpr ⟨b, by linarith⟩
  apply Finset.card_bij (fun p _ => f p)
  · -- into
    intro p hp
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp
    obtain ⟨⟨hp0, hpq⟩, hcond⟩ := hp
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · rcases Nat.eq_zero_or_pos (f p) with h0 | hpos
      · exfalso
        have hzero : (v i * (p : ℤ)) % (q : ℤ) = 0 := by
          have hv := hfval p
          rw [h0] at hv
          simpa using hv.symm
        have hdvd : (q : ℤ) ∣ v i * (p : ℤ) := Int.dvd_of_emod_eq_zero hzero
        have hp_eq : (p : ℤ) = a * (v i * (p : ℤ)) + (b * (p : ℤ)) * (q : ℤ) := by
          linear_combination (-(p : ℤ)) * hab
        have hqp : (q : ℤ) ∣ (p : ℤ) := by
          rw [hp_eq]
          exact dvd_add (Dvd.dvd.mul_left hdvd a) (Dvd.intro_left _ rfl)
        have hppos : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp0
        have hple := Int.le_of_dvd hppos hqp
        have hplt : (p : ℤ) < (q : ℤ) := by exact_mod_cast hpq
        omega
      · exact hpos
    · have h1 : ((f p : ℕ) : ℤ) < (q : ℤ) := by rw [hfval]; exact hmod_lt p
      exact_mod_cast h1
    · have hqq : q = 14 * Q := hqdef
      constructor
      · intro hcon
        exact hcond.1 (by rw [hqq] at *; exact hcon)
      · intro hcon
        exact hcond.2 (by rw [hqq] at *; exact hcon)
  · -- injective
    intro p₁ hp₁ p₂ hp₂ hEq
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    have hmodeq : v i * (p₁ : ℤ) ≡ v i * (p₂ : ℤ) [ZMOD (q : ℤ)] := by
      show (v i * (p₁ : ℤ)) % (q : ℤ) = (v i * (p₂ : ℤ)) % (q : ℤ)
      have hc : ((f p₁ : ℕ) : ℤ) = ((f p₂ : ℕ) : ℤ) := by exact_mod_cast hEq
      rw [hfval, hfval] at hc
      exact hc
    have hpq : (p₁ : ℤ) ≡ (p₂ : ℤ) [ZMOD (q : ℤ)] := by
      calc (p₁ : ℤ) = 1 * (p₁ : ℤ) := (one_mul _).symm
        _ ≡ (v i * a) * (p₁ : ℤ) [ZMOD (q : ℤ)] := (hva.symm).mul_right _
        _ = a * (v i * (p₁ : ℤ)) := by ring
        _ ≡ a * (v i * (p₂ : ℤ)) [ZMOD (q : ℤ)] := hmodeq.mul_left _
        _ = (v i * a) * (p₂ : ℤ) := by ring
        _ ≡ 1 * (p₂ : ℤ) [ZMOD (q : ℤ)] := hva.mul_right _
        _ = (p₂ : ℤ) := one_mul _
    have h1 : (p₁ : ℤ) % (q : ℤ) = (p₁ : ℤ) :=
      Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hp₁.1.2)
    have h2 : (p₂ : ℤ) % (q : ℤ) = (p₂ : ℤ) :=
      Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hp₂.1.2)
    have hfin : (p₁ : ℤ) = (p₂ : ℤ) := by
      have hpq' : (p₁ : ℤ) % (q : ℤ) = (p₂ : ℤ) % (q : ℤ) := hpq
      rw [h1, h2] at hpq'
      exact hpq'
    exact_mod_cast hfin
  · -- surjective
    intro r hr
    rw [Finset.mem_filter, Finset.mem_Ioo] at hr
    obtain ⟨⟨hr0, hrq⟩, hrcond⟩ := hr
    set p : ℕ := ((a * (r : ℤ)) % (q : ℤ)).toNat with hpdef
    have hpmod_nonneg : (0 : ℤ) ≤ (a * (r : ℤ)) % (q : ℤ) := Int.emod_nonneg _ (by omega)
    have hpmod_lt : (a * (r : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
    have hpval : ((p : ℕ) : ℤ) = (a * (r : ℤ)) % (q : ℤ) := by
      rw [hpdef]
      exact Int.toNat_of_nonneg hpmod_nonneg
    have hp_ar : (p : ℤ) ≡ a * (r : ℤ) [ZMOD (q : ℤ)] := by
      show (p : ℤ) % (q : ℤ) = (a * (r : ℤ)) % (q : ℤ)
      rw [hpval]
      exact Int.emod_emod_of_dvd _ dvd_rfl
    have hvp : v i * (p : ℤ) ≡ (r : ℤ) [ZMOD (q : ℤ)] := by
      calc v i * (p : ℤ)
          ≡ v i * (a * (r : ℤ)) [ZMOD (q : ℤ)] := hp_ar.mul_left _
        _ = (v i * a) * (r : ℤ) := by ring
        _ ≡ 1 * (r : ℤ) [ZMOD (q : ℤ)] := hva.mul_right _
        _ = (r : ℤ) := one_mul _
    have hfp : f p = r := by
      have hfpZ : ((f p : ℕ) : ℤ) = (r : ℤ) := by
        rw [hfval]
        have hrsmall : (r : ℤ) % (q : ℤ) = (r : ℤ) :=
          Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hrq)
        have hvp' : (v i * (p : ℤ)) % (q : ℤ) = (r : ℤ) % (q : ℤ) := hvp
        rw [hrsmall] at hvp'
        exact hvp'
      exact_mod_cast hfpZ
    refine ⟨p, ?_, hfp⟩
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · rcases Nat.eq_zero_or_pos p with h0 | hpos
      · exfalso
        have hr' : f p = r := hfp
        rw [h0] at hr'
        have hf0 : f 0 = 0 := by
          rw [hf]
          simp
        omega
      · exact hpos
    · have hlt : ((p : ℕ) : ℤ) < (q : ℤ) := by rw [hpval]; exact hpmod_lt
      exact_mod_cast hlt
    · rw [hfp]
      exact hrcond

/-! ## Axiom audit -/
#print axioms jointFail_anti
#print axioms jointFail_pair_lower
#print axioms dilatePairCount_eq
#print axioms jointFail_dilate_pair_eq

end LRC14Concrete
end LonelyRunner
