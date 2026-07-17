/- Thm892Shadows.lean — klein-2026-07-16-S318.
   THM-892 (K), the rational heart (boxeph manifest item 6, first lemma):
   the tent K(j/P) = (j/P)(1−j/P) on Z_P has discrete second difference
   −2/P² away from 0 and −2/P² + 2/P at 0 — for EVERY P ≥ 2, general proof,
   no decide.  The csc² kernel identity k̂_P(n) = −csc²(πn/P)/(2P²) of THM-892
   is the DFT diagonalization of exactly this statement (−4sin²(πn/P)·k̂(n)
   = 2/P²), so this module is the kernel identity's kernel-checkable core.
   No sorries, no native_decide. -/
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

namespace LonelyRunner
namespace LRC14
namespace Thm892

/-- The tent value at a residue: `K(j/P) = (j/P)·(1 − j/P)`, cleared to ℚ. -/
def tent (P : ℕ) (j : ZMod P) : ℚ :=
  (j.val : ℚ) * ((P : ℚ) - (j.val : ℚ)) / (P : ℚ) ^ 2

/-- The ℕ-representative form used by the wrap-around computation. -/
def tentN (P v : ℕ) : ℚ := (v : ℚ) * ((P : ℚ) - (v : ℚ)) / (P : ℚ) ^ 2

theorem tent_eq_tentN (P : ℕ) (j : ZMod P) : tent P j = tentN P j.val := rfl

/-- The representative-level second difference: for `v < P`,
    `K((v + (P−1)) mod P) + K((v+1) mod P) − 2·K(v)` equals `−2/P²`,
    plus the extra `2/P` exactly at `v = 0`. -/
theorem tentN_second_difference (P v : ℕ) (h2 : 2 ≤ P) (hv : v < P) :
    tentN P ((v + (P - 1)) % P) + tentN P ((v + 1) % P) - 2 * tentN P v
      = -2 / (P : ℚ) ^ 2 + (if v = 0 then 2 / (P : ℚ) else 0) := by
  have hP0 : (P : ℚ) ≠ 0 := by
    have hp : 0 < P := by omega
    exact_mod_cast hp.ne'
  rcases Nat.eq_zero_or_pos v with rfl | hvpos
  · -- v = 0: left neighbour is P−1, right neighbour is 1
    have e1 : (0 + (P - 1)) % P = P - 1 := by
      have h' : 0 + (P - 1) = P - 1 := by omega
      rw [h', Nat.mod_eq_of_lt (by omega)]
    have e2 : (0 + 1) % P = 1 := by
      have h' : 0 + 1 = 1 := by omega
      rw [h', Nat.mod_eq_of_lt (by omega)]
    rw [e1, e2, if_pos rfl]
    unfold tentN
    have hc : ((P - 1 : ℕ) : ℚ) = (P : ℚ) - 1 := by
      have h1 : (1:ℕ) ≤ P := by omega
      push_cast [Nat.cast_sub h1]
      ring
    rw [hc]
    field_simp
    ring
  · rcases Nat.lt_or_ge v (P - 1) with hmid | htop
    · -- 1 ≤ v ≤ P−2: no wrap on either side
      have e1 : (v + (P - 1)) % P = v - 1 := by
        have h' : v + (P - 1) = (v - 1) + P := by omega
        rw [h', Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
      have e2 : (v + 1) % P = v + 1 := Nat.mod_eq_of_lt (by omega)
      rw [e1, e2, if_neg (by omega)]
      unfold tentN
      have hc1 : ((v - 1 : ℕ) : ℚ) = (v : ℚ) - 1 := by
        have h1 : (1:ℕ) ≤ v := hvpos
        push_cast [Nat.cast_sub h1]
        ring
      have hc2 : ((v + 1 : ℕ) : ℚ) = (v : ℚ) + 1 := by push_cast; ring
      rw [hc1, hc2]
      field_simp
      ring
    · -- v = P−1: right neighbour wraps to 0
      have hvtop : v = P - 1 := by omega
      subst hvtop
      have e1 : ((P - 1) + (P - 1)) % P = P - 2 := by
        have h' : (P - 1) + (P - 1) = (P - 2) + P := by omega
        rw [h', Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
      have e2 : ((P - 1) + 1) % P = 0 := by
        have h' : (P - 1) + 1 = P := by omega
        rw [h', Nat.mod_self]
      rw [e1, e2, if_neg (by omega)]
      unfold tentN
      have hc1 : ((P - 2 : ℕ) : ℚ) = (P : ℚ) - 2 := by
        have h1 : (2:ℕ) ≤ P := h2
        push_cast [Nat.cast_sub h1]
        ring
      have hc2 : ((P - 1 : ℕ) : ℚ) = (P : ℚ) - 1 := by
        have h1 : (1:ℕ) ≤ P := by omega
        push_cast [Nat.cast_sub h1]
        ring
      rw [hc1, hc2]
      field_simp
      ring

/-- **THM-892 (K), the rational heart on `ZMod P`**: for every `P ≥ 2` and
    `j : ZMod P`, `K(j−1) + K(j+1) − 2·K(j) = −2/P² + (2/P)·[j = 0]`. -/
theorem tent_second_difference (P : ℕ) (h2 : 2 ≤ P) (j : ZMod P) :
    tent P (j - 1) + tent P (j + 1) - 2 * tent P j
      = -2 / (P : ℚ) ^ 2 + (if j = 0 then 2 / (P : ℚ) else 0) := by
  haveI : NeZero P := ⟨by omega⟩
  obtain ⟨n, rfl⟩ : ∃ n, P = n + 1 := ⟨P - 1, by omega⟩
  have hval_add : (j + 1).val = (j.val + 1) % (n + 1) := by
    rw [ZMod.val_add, ZMod.val_one_eq_one_mod]
    congr 1
    have h1 : 1 % (n + 1) = 1 := Nat.mod_eq_of_lt (by omega)
    rw [h1]
  have hval_sub : (j - 1).val = (j.val + n) % (n + 1) := by
    rw [sub_eq_add_neg, ZMod.val_add, ZMod.val_neg_one]
  have hkey := tentN_second_difference (n + 1) j.val h2 (ZMod.val_lt j)
  have hn1 : (n + 1) - 1 = n := by omega
  rw [hn1] at hkey
  rw [tent_eq_tentN, tent_eq_tentN, tent_eq_tentN, hval_add, hval_sub, hkey]
  congr 1
  simp only [ZMod.val_eq_zero]

end Thm892
end LRC14
end LonelyRunner
