/-
  TournamentH7.RedeiFromOCF — Rédei's H ≥ 1 derived from OCF

  Project canon Rédei's theorem has two parts:
    (R1) Every tournament has at least one Hamilton path: H(T) ≥ 1.
    (R2) The number of Hamilton paths is always odd.

  Currently both are axiomatised in `Redei.lean`.  But (R1) actually
  follows from OCF (THM-326): H(T) = 1 + 2 α_1(T) + 4 α_2(T) + …
  Since each α_k(T) ≥ 0, this is at least 1.

  This module formalises this observation.
-/

import TournamentH7.OCF
import TournamentH7.SCC
import Mathlib.Algebra.Ring.Parity

namespace Tournament

variable {n : ℕ}

/-! ### Rédei's existence (H ≥ 1) — derived from OCF, no axiom needed -/

/-- **Theorem (Rédei's existence, oracle-S10).**

    For every tournament T, `H(T) ≥ 1`. Derived from OCF
    (without needing the separate `redei_existence` axiom). -/
theorem H_ge_one (T : Tournament n) : 1 ≤ H T := by
  have hocf := ocf T
  -- hocf : H T = 1 + 2 α_1 + 4 α_2 + 8 α_3 + 16 α_4
  -- Each α_k ≥ 0 (natural numbers), so RHS ≥ 1.
  omega

/-! ### Corollaries -/

/-- H(T) ≠ 0 (Rédei's existence). -/
theorem H_ne_zero (T : Tournament n) : H T ≠ 0 := by
  have := H_ge_one T
  omega

/-! ### Rédei's parity (H is odd) — also derived from OCF -/

/-- **Theorem (Rédei's parity, oracle-S11).**

    For every tournament T, `H(T)` is odd. Derived from OCF:
    `H = 1 + 2 α_1 + 4 α_2 + 8 α_3 + 16 α_4 = 1 + 2 * (α_1 + 2 α_2 + …)`,
    which is `1 + even`, hence odd.

    NO `redei_parity` axiom is needed! -/
theorem H_odd_from_ocf (T : Tournament n) : Odd (H T) := by
  have hocf := ocf T
  -- H T = 1 + 2 α_1 + 4 α_2 + 8 α_3 + 16 α_4
  -- = 1 + 2(α_1 + 2 α_2 + 4 α_3 + 8 α_4)
  refine ⟨alphaCount 1 T + 2 * alphaCount 2 T
            + 4 * alphaCount 3 T + 8 * alphaCount 4 T, ?_⟩
  omega

/-- OCF leaves the explicit mod-2 residue `1` for the Hamiltonian-path count. -/
theorem H_mod_two_eq_one_from_ocf (T : Tournament n) :
    H T % 2 = 1 := by
  rcases H_odd_from_ocf T with ⟨k, hk⟩
  rw [hk]
  omega

/-- H(T) is not even. -/
theorem H_not_even (T : Tournament n) : ¬ Even (H T) := by
  rw [Nat.not_even_iff_odd]
  exact H_odd_from_ocf T

/-- H(T) ≠ 2 (since H is odd). -/
theorem H_ne_two_from_ocf (T : Tournament n) : H T ≠ 2 := by
  intro h
  have := H_odd_from_ocf T
  rw [h] at this
  exact (by decide : ¬ Odd 2) this

/-- H(T) ≠ any even number. -/
theorem H_ne_even_from_ocf (T : Tournament n) (m : ℕ) (hm : Even m) : H T ≠ m := by
  intro h
  have hodd := H_odd_from_ocf T
  rw [h] at hodd
  exact (Nat.not_odd_iff_even.mpr hm) hodd

end Tournament
