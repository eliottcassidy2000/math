/-
klein-2026-07-01-S92 (HYP-3839) — the nest-chain telescoping identity.

Inclusion–exclusion over a NESTED family telescopes: when every d-fold intersection of a
subfamily T depends only on max T (the gcd-nest situation of THM-601 / HYP-3837/3838),

    ∑_{∅ ≠ T ⊆ Q} (-1)^(|T|+1) · g(max T)  =  g(min Q).

With g(v) = 2r/v this is the engine of the cap decomposition's nest chain. Companion to
PolygonPartitionDMNR.lean.
-/
import Mathlib

open Finset

namespace NestTelescope

variable {α : Type*} [LinearOrder α] {M : Type*} [AddCommGroup M]

/-- Total signed term: `(-1)^(|T|+1) • g(max T)`, reading `T.max : WithBot α` through
`recBotCoe` so the empty set contributes `0`. -/
noncomputable def term (g : α → M) (T : Finset α) : M :=
  ((-1 : ℤ) ^ (T.card + 1)) • (T.max.recBotCoe (0 : M) g)

lemma term_empty (g : α → M) : term g (∅ : Finset α) = 0 := by
  simp [term]

lemma term_nonempty (g : α → M) {T : Finset α} (h : T.Nonempty) :
    term g T = ((-1 : ℤ) ^ (T.card + 1)) • g (T.max' h) := by
  rw [term, ← Finset.coe_max' h]
  rfl

/-- **The telescoping identity.** -/
theorem sum_term_eq_min (Q : Finset α) (hQ : Q.Nonempty) (g : α → M) :
    ∑ T ∈ Q.powerset, term g T = g (Q.min' hQ) := by
  induction Q using Finset.induction_on_min with
  | empty => exact absurd hQ (by simp)
  | insert a s ha ih =>
    have has : a ∉ s := fun h => lt_irrefl a (ha a h)
    have hdisj : Disjoint s.powerset (s.powerset.image (insert a)) := by
      refine Finset.disjoint_left.mpr fun T hT hT2 => ?_
      obtain ⟨T', _, rfl⟩ := Finset.mem_image.mp hT2
      exact has (Finset.mem_powerset.mp hT (Finset.mem_insert_self a T'))
    have hinj : ∀ x ∈ s.powerset, ∀ y ∈ s.powerset, insert a x = insert a y → x = y := by
      intro x hx y hy hxy
      have hxs : a ∉ x := fun h => has (Finset.mem_powerset.mp hx h)
      have hys : a ∉ y := fun h => has (Finset.mem_powerset.mp hy h)
      have := congrArg (·.erase a) hxy
      simpa [Finset.erase_insert hxs, Finset.erase_insert hys] using this
    rw [Finset.powerset_insert, Finset.sum_union hdisj, Finset.sum_image hinj]
    have step : ∀ T ∈ s.powerset, term g (insert a T)
        = -term g T + (if T = ∅ then g a else 0) := by
      intro T hT
      have hTs : T ⊆ s := Finset.mem_powerset.mp hT
      have haT : a ∉ T := fun h => has (hTs h)
      have hcard : (insert a T).card = T.card + 1 := Finset.card_insert_of_notMem haT
      rcases T.eq_empty_or_nonempty with rfl | hTne
      · simp only [insert_empty, term_nonempty g (Finset.singleton_nonempty a),
          Finset.max'_singleton, Finset.card_singleton, term_empty, neg_zero, zero_add,
          if_pos rfl]
        norm_num
      · have hTne' : T ≠ ∅ := Finset.nonempty_iff_ne_empty.mp hTne
        have hins : (insert a T).Nonempty := ⟨a, Finset.mem_insert_self a T⟩
        have hmax : (insert a T).max = T.max := by
          rw [Finset.max_insert]
          obtain ⟨b, hb⟩ := hTne
          have hab : (a : WithBot α) ≤ T.max :=
            le_trans (WithBot.coe_le_coe.mpr (le_of_lt (ha b (hTs hb)))) (Finset.le_max hb)
          exact max_eq_right hab
        have hmax' : (insert a T).max' hins = T.max' hTne :=
          WithBot.coe_injective (by rw [Finset.coe_max', Finset.coe_max', hmax])
        rw [term_nonempty g hins, term_nonempty g hTne, hmax', hcard, if_neg hTne',
          add_zero, pow_succ, mul_comm, mul_smul, neg_one_zsmul]
    rw [Finset.sum_congr rfl step, Finset.sum_add_distrib]
    have h1 : ∑ T ∈ s.powerset, -term g T = -∑ T ∈ s.powerset, term g T := by
      simp
    have h2 : ∑ T ∈ s.powerset, (if T = ∅ then g a else 0) = g a := by
      rw [Finset.sum_ite_eq' s.powerset (∅ : Finset α) (fun _ => g a)]
      simp
    rw [h1, h2]
    have hmin : (insert a s).min' hQ = a := by
      apply le_antisymm (Finset.min'_le _ a (Finset.mem_insert_self a s))
      exact Finset.le_min' _ _ _ fun y hy => by
        rcases Finset.mem_insert.mp hy with rfl | hys
        · exact le_refl _
        · exact le_of_lt (ha y hys)
    rcases s.eq_empty_or_nonempty with rfl | hsne
    · simp [term, hmin]
    · rw [ih hsne, hmin]
      abel

end NestTelescope
