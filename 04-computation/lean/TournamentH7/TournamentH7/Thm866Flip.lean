/- Thm866Flip.lean — klein-2026-07-16-S317.
   THM-866 rung two, part I: the F3 arc-flip law (THM-855 F3).
   Flipping the arc between `w` and `l` sends `xLevel` to
   `xLevel + 4·(dev l − dev w) + 8` when `w` beat `l`; in particular a
   TIE-SPLITTING flip (`dev w = dev l`) is an exact +8 step.
   Pure finite algebra on the CompositionDefect definitions.  No sorries. -/
import TournamentH7.CompositionDefect

open Finset

namespace Tournament

variable {n : ℕ}

/-- The symmetric flip condition: the ordered pair `(i,j)` is one of the two
    slots between `w` and `l`. -/
def flipCond (w l i j : Fin n) : Prop := (i = w ∧ j = l) ∨ (i = l ∧ j = w)

instance (w l i j : Fin n) : Decidable (flipCond w l i j) := by
  unfold flipCond; infer_instance

theorem flipCond_symm {w l i j : Fin n} (h : flipCond w l i j) : flipCond w l j i := by
  rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · exact Or.inr ⟨h2, h1⟩
  · exact Or.inl ⟨h2, h1⟩

/-- The tournament obtained from `T` by reversing the arc between `w` and `l`. -/
def flip (T : Tournament n) (w l : Fin n) (hne : w ≠ l) : Tournament n where
  arc i j := if flipCond w l i j then !(T.arc i j) else T.arc i j
  irrefl := by
    intro i
    have hc : ¬ flipCond w l i i := by
      rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
      · exact hne (h1.symm.trans h2)
      · exact hne (h2.symm.trans h1)
    show (if flipCond w l i i then !(T.arc i i) else T.arc i i) = false
    rw [if_neg hc, T.irrefl i]
  total := by
    intro i j hij
    show (if flipCond w l i j then !(T.arc i j) else T.arc i j) = true ∨
         (if flipCond w l j i then !(T.arc j i) else T.arc j i) = true
    by_cases hc : flipCond w l i j
    · rw [if_pos hc, if_pos (flipCond_symm hc)]
      rcases T.total i j hij with h | h
      · right
        have hf : T.arc j i = false := by
          cases h2 : T.arc j i
          · rfl
          · exact absurd ⟨h, h2⟩ (T.asym i j)
        rw [hf]
        rfl
      · left
        have hf : T.arc i j = false := by
          cases h2 : T.arc i j
          · rfl
          · exact absurd ⟨h2, h⟩ (T.asym i j)
        rw [hf]
        rfl
    · rw [if_neg hc, if_neg (fun h => hc (flipCond_symm h))]
      exact T.total i j hij
  asym := by
    intro i j
    show ¬ ((if flipCond w l i j then !(T.arc i j) else T.arc i j) = true ∧
            (if flipCond w l j i then !(T.arc j i) else T.arc j i) = true)
    by_cases hc : flipCond w l i j
    · rw [if_pos hc, if_pos (flipCond_symm hc)]
      rintro ⟨h1, h2⟩
      have hij : i ≠ j := by
        rintro rfl
        rcases hc with ⟨a, b⟩ | ⟨a, b⟩
        · exact hne (a.symm.trans b)
        · exact hne (b.symm.trans a)
      have ha : T.arc i j = false := by
        cases h : T.arc i j
        · rfl
        · rw [h] at h1
          simp at h1
      have hb : T.arc j i = false := by
        cases h : T.arc j i
        · rfl
        · rw [h] at h2
          simp at h2
      rcases T.total i j hij with h | h
      · rw [ha] at h
        simp at h
      · rw [hb] at h
        simp at h
    · rw [if_neg hc, if_neg (fun h => hc (flipCond_symm h))]
      exact T.asym i j

variable (T : Tournament n) {w l : Fin n}

/-- Off the flipped pair, arcs are unchanged. -/
theorem flip_arc_other (hne : w ≠ l) {v j : Fin n} (hvw : v ≠ w) (hvl : v ≠ l) :
    (T.flip w l hne).arc v j = T.arc v j := by
  have hc : ¬ flipCond w l v j := by
    rintro (⟨h1, _⟩ | ⟨h1, _⟩)
    · exact hvw h1
    · exact hvl h1
  simp [flip, hc]

theorem flip_arc_w (hne : w ≠ l) {j : Fin n} (hj : j ≠ l) :
    (T.flip w l hne).arc w j = T.arc w j := by
  have hc : ¬ flipCond w l w j := by
    rintro (⟨_, h2⟩ | ⟨h1, _⟩)
    · exact hj h2
    · exact hne h1
  simp [flip, hc]

theorem flip_arc_wl (hne : w ≠ l) : (T.flip w l hne).arc w l = !(T.arc w l) := by
  have hc : flipCond w l w l := Or.inl ⟨rfl, rfl⟩
  simp [flip, hc]

theorem flip_arc_l (hne : w ≠ l) {j : Fin n} (hj : j ≠ w) :
    (T.flip w l hne).arc l j = T.arc l j := by
  have hc : ¬ flipCond w l l j := by
    rintro (⟨h1, _⟩ | ⟨_, h2⟩)
    · exact hne h1.symm
    · exact hj h2
  simp [flip, hc]

theorem flip_arc_lw (hne : w ≠ l) : (T.flip w l hne).arc l w = !(T.arc l w) := by
  have hc : flipCond w l l w := Or.inr ⟨rfl, rfl⟩
  simp [flip, hc]

/-- The loser's out-set loses exactly the flipped opponent. -/
theorem outDegree_flip_w (hne : w ≠ l) (hwl : T.arc w l = true) :
    (T.flip w l hne).outDegree w + 1 = T.outDegree w := by
  have hset : (Finset.univ.filter (fun j => (T.flip w l hne).arc w j = true))
      = (Finset.univ.filter (fun j => T.arc w j = true)).erase l := by
    ext j
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_erase]
    by_cases hj : j = l
    · subst hj
      rw [flip_arc_wl T hne, hwl]
      simp
    · rw [flip_arc_w T hne hj]
      constructor
      · exact fun h => ⟨hj, h⟩
      · exact fun h => h.2
  unfold outDegree
  rw [hset]
  have hmem : l ∈ Finset.univ.filter (fun j => T.arc w j = true) := by
    simp [hwl]
  rw [Finset.card_erase_of_mem hmem]
  have hpos : 0 < (Finset.univ.filter (fun j => T.arc w j = true)).card :=
    Finset.card_pos.mpr ⟨l, hmem⟩
  omega

/-- The winner's out-set gains exactly the flipped opponent. -/
theorem outDegree_flip_l (hne : w ≠ l) (hwl : T.arc w l = true) :
    (T.flip w l hne).outDegree l = T.outDegree l + 1 := by
  have hlw : T.arc l w = false := by
    cases h : T.arc l w with
    | true => exact absurd ⟨hwl, h⟩ (T.asym w l)
    | false => rfl
  have hset : (Finset.univ.filter (fun j => (T.flip w l hne).arc l j = true))
      = insert w (Finset.univ.filter (fun j => T.arc l j = true)) := by
    ext j
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_insert]
    by_cases hj : j = w
    · subst hj
      rw [flip_arc_lw T hne, hlw]
      simp
    · rw [flip_arc_l T hne hj]
      constructor
      · exact fun h => Or.inr h
      · rintro (h | h)
        · exact absurd h hj
        · exact h
  unfold outDegree
  rw [hset, Finset.card_insert_of_notMem (by simp [hlw])]

/-- Bystanders keep their out-degree. -/
theorem outDegree_flip_other (hne : w ≠ l) {v : Fin n} (hvw : v ≠ w) (hvl : v ≠ l) :
    (T.flip w l hne).outDegree v = T.outDegree v := by
  unfold outDegree
  congr 1
  apply Finset.filter_congr
  intro j _
  rw [flip_arc_other T hne hvw hvl]

/-- Centered-score changes: the loser drops by 2. -/
theorem dev_flip_w (hne : w ≠ l) (hwl : T.arc w l = true) :
    (T.flip w l hne).dev w = T.dev w - 2 := by
  have h := outDegree_flip_w T hne hwl
  unfold dev
  have : ((T.flip w l hne).outDegree w : ℤ) = (T.outDegree w : ℤ) - 1 := by
    omega
  rw [this]; ring

/-- The winner gains 2. -/
theorem dev_flip_l (hne : w ≠ l) (hwl : T.arc w l = true) :
    (T.flip w l hne).dev l = T.dev l + 2 := by
  have h := outDegree_flip_l T hne hwl
  unfold dev
  rw [h]; push_cast; ring

theorem dev_flip_other (hne : w ≠ l) {v : Fin n} (hvw : v ≠ w) (hvl : v ≠ l) :
    (T.flip w l hne).dev v = T.dev v := by
  unfold dev
  rw [outDegree_flip_other T hne hvw hvl]

/-- **THM-855 F3, the flip law**: reversing the arc `w → l` changes the axis
    level by exactly `4·(dev l − dev w) + 8`. -/
theorem xLevel_flip (hne : w ≠ l) (hwl : T.arc w l = true) :
    (T.flip w l hne).xLevel = T.xLevel + 4 * (T.dev l - T.dev w) + 8 := by
  unfold xLevel
  have hpt : ∀ v : Fin n, ((T.flip w l hne).dev v) ^ 2
      = (T.dev v) ^ 2 + ((if v = w then -4 * T.dev w + 4 else 0)
        + (if v = l then 4 * T.dev l + 4 else 0)) := by
    intro v
    by_cases hvw : v = w
    · subst hvw
      rw [if_pos rfl, if_neg hne, dev_flip_w T hne hwl]
      ring
    · by_cases hvl : v = l
      · subst hvl
        rw [if_neg hvw, if_pos rfl, dev_flip_l T hne hwl]
        ring
      · rw [if_neg hvw, if_neg hvl, dev_flip_other T hne hvw hvl]
        ring
  calc ∑ v : Fin n, ((T.flip w l hne).dev v) ^ 2
      = ∑ v : Fin n, ((T.dev v) ^ 2 + ((if v = w then -4 * T.dev w + 4 else 0)
          + (if v = l then 4 * T.dev l + 4 else 0))) :=
        Finset.sum_congr rfl fun v _ => hpt v
    _ = (∑ v : Fin n, (T.dev v) ^ 2)
        + ((∑ v : Fin n, if v = w then -4 * T.dev w + 4 else 0)
          + (∑ v : Fin n, if v = l then 4 * T.dev l + 4 else 0)) := by
        rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
    _ = (∑ v : Fin n, (T.dev v) ^ 2) + ((-4 * T.dev w + 4) + (4 * T.dev l + 4)) := by
        rw [Finset.sum_ite_eq' Finset.univ w (fun _ => -4 * T.dev w + 4),
            Finset.sum_ite_eq' Finset.univ l (fun _ => 4 * T.dev l + 4)]
        simp
    _ = (∑ v : Fin n, (T.dev v) ^ 2) + 4 * (T.dev l - T.dev w) + 8 := by ring

/-- **THM-866 step (3)**: a TIE-SPLITTING flip is an exact `+8` step. -/
theorem xLevel_flip_tie (hne : w ≠ l) (hwl : T.arc w l = true)
    (htie : T.dev w = T.dev l) :
    (T.flip w l hne).xLevel = T.xLevel + 8 := by
  rw [xLevel_flip T hne hwl, htie]; ring

/-- Ties in `dev` are exactly ties in score. -/
theorem dev_eq_iff_outDegree_eq {u v : Fin n} :
    T.dev u = T.dev v ↔ T.outDegree u = T.outDegree v := by
  unfold dev
  constructor
  · intro h; omega
  · intro h; rw [h]

end Tournament
