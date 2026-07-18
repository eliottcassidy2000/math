import TournamentH7.LRCQ25Obstruction

namespace LonelyRunner
namespace LRC14Concrete

theorem q_pos_of_hasQWitness {v : Fin 13 → ℤ} {q : ℕ}
    (h : HasQWitness v q) : 0 < q := by
  obtain ⟨p, hp, _⟩ := h
  exact lt_trans (Finset.mem_Ioo.mp hp).1 (Finset.mem_Ioo.mp hp).2

theorem speeds_ne_zero_of_hasQWitness {v : Fin 13 → ℤ} {q : ℕ}
    (h : HasQWitness v q) : ∀ i, v i ≠ 0 := by
  have hq := q_pos_of_hasQWitness h
  obtain ⟨p, _hp, hband⟩ := h
  intro i hzero
  have hlo := (hband i).1
  simp [hzero] at hlo
  omega

theorem mreach_ge_of_hasQWitness {v : Fin 13 → ℤ} {q : ℕ}
    (h : HasQWitness v q) : (1 : ℝ) / 14 ≤ Mreach v := by
  have hq := q_pos_of_hasQWitness h
  obtain ⟨p, _hp, hband⟩ := h
  exact mreach_ge_of_pairsum_band v (p : ℤ) (q : ℤ)
    (by exact_mod_cast hq) hband

theorem lonely_of_hasQWitness {v : Fin 13 → ℤ} {q : ℕ}
    (h : HasQWitness v q) : ∃ t : ℝ, Lonely 14 v t :=
  lonely_of_Mreach_ge v (speeds_ne_zero_of_hasQWitness h)
    (mreach_ge_of_hasQWitness h)

theorem lonely_of_liveCount_pos {v : Fin 13 → ℤ} {q : ℕ}
    (h : 0 < liveCount v q) : ∃ t : ℝ, Lonely 14 v t :=
  lonely_of_hasQWitness ((hasQWitness_iff_liveCount_pos v q).2 h)

#print axioms q_pos_of_hasQWitness
#print axioms speeds_ne_zero_of_hasQWitness
#print axioms mreach_ge_of_hasQWitness
#print axioms lonely_of_hasQWitness
#print axioms lonely_of_liveCount_pos

end LRC14Concrete
end LonelyRunner
