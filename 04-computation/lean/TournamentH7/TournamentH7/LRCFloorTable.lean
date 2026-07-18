/- LRCFloorTable.lean -- opus-2026-07-17-S347 (HYP-7300 / THM-964 (S)+(P)).
   THE PER-CLASS FLOOR, from the kernel-proved folded identity (muNum_folded,
   THM-965).  fold(r) = r(14−r) sits in [0, 49] (max (r−7)² ≥ 0), so

     14 · muNum a b  ≥  4ab − 49,

   i.e. the pair overlap is ≥ (4ab − 49)/(196 ab) = 1/49 − 1/(4ab), monotone
   increasing in ab — the whole THM-964 sawtooth floor table follows by
   evaluating this bound (no case sweep).  Kernel-pure: no sorry, no
   native_decide. -/
import Mathlib
import TournamentH7.LRCFoldedIdentity

namespace LonelyRunner.LRC14.Hunter

/-- `fold r = r(14−r)` on a residue `n % 14` lies in `[0, 49]`. -/
theorem fold_mem (n : ℕ) :
    0 ≤ ((n % 14 : ℕ) : ℤ) * (14 - ((n % 14 : ℕ) : ℤ))
      ∧ ((n % 14 : ℕ) : ℤ) * (14 - ((n % 14 : ℕ) : ℤ)) ≤ 49 := by
  have hlt : (n % 14 : ℕ) < 14 := Nat.mod_lt _ (by norm_num)
  have h13 : ((n % 14 : ℕ) : ℤ) ≤ 13 := by exact_mod_cast Nat.lt_succ_iff.mp hlt
  have h0 : (0 : ℤ) ≤ ((n % 14 : ℕ) : ℤ) := Int.natCast_nonneg _
  set r : ℤ := ((n % 14 : ℕ) : ℤ)
  refine ⟨?_, ?_⟩
  · nlinarith
  · nlinarith [sq_nonneg (r - 7)]

/-- **THE PER-CLASS FLOOR (muNum units):** `4ab − 49 ≤ 14·muNum a b`.
The pair-overlap sum is within 49 of `4ab` — the exact deviation is two
folds (THM-965), each in `[0, 49]`. -/
theorem muNum_lower (a b : ℕ) (ha : 1 ≤ a) (hab : a ≤ b) :
    (4 * (a : ℤ) * b : ℤ) - 49 ≤ 14 * muNum a b := by
  rw [muNum_folded a b ha hab]
  have hS := fold_mem (a + b)
  have hD := fold_mem (b - a)
  linarith [hS.1, hD.2]

/-- **THE MONOTONE FLOOR:** the muNum lower bound is increasing in `ab`, so a
uniform product lower bound `4ab ≥ Q` gives `14·muNum a b ≥ Q − 49`.  This is
the shape THM-964 (P) consumes: each Hunter path-tree edge of a comparable
7-block carries a size-dependent floor whose sum is the global μ(U₇) bound. -/
theorem muNum_lower_of_prod (a b : ℕ) (ha : 1 ≤ a) (hab : a ≤ b)
    (Q : ℤ) (hQ : Q ≤ 4 * (a : ℤ) * b) :
    Q - 49 ≤ 14 * muNum a b :=
  le_trans (by linarith) (muNum_lower a b ha hab)

/-- The rational overlap floor for coprime pairs: `μ = muNum/(14ab)` gives
`μ ≥ 1/49 − 1/(4ab)`.  (Stated as the muNum inequality scaled; the measure
identification `μ·14ab = muNum` is boxeph's LEM-042 defining relation.) -/
theorem overlap_floor_rat (a b : ℕ) (ha : 1 ≤ a) (hab : a ≤ b) :
    (1 : ℚ) / 49 - 1 / (4 * a * b) ≤ (muNum a b : ℚ) / (14 * a * b) := by
  have ha0 : (0 : ℚ) < a := by exact_mod_cast ha
  have hb0 : (0 : ℚ) < b := by exact_mod_cast (le_trans ha hab)
  have hlowQ : (4 * (a : ℚ) * b - 49 : ℚ) ≤ 14 * (muNum a b : ℚ) := by
    have := muNum_lower a b ha hab
    exact_mod_cast this
  have key : (1 : ℚ) / 49 - 1 / (4 * a * b) = (4 * a * b - 49) / (196 * a * b) := by
    field_simp
    ring
  have hrhs : (muNum a b : ℚ) / (14 * a * b) = (14 * muNum a b) / (196 * a * b) := by
    rw [div_eq_div_iff (by positivity) (by positivity)]
    ring
  rw [key, hrhs]
  gcongr

/-! ## Axiom audit -/
#print axioms fold_mem
#print axioms muNum_lower
#print axioms muNum_lower_of_prod
#print axioms overlap_floor_rat

end LonelyRunner.LRC14.Hunter
