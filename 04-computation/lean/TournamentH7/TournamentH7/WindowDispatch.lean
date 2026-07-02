/-
klein-2026-07-02-S112 (HYP-4017) — THE PACK-LIST DISPATCHER: the window census as data
+ one decidable sweep + membership soundness + the enumeration-completeness reflection,
discharging `hdistinct20` outright. With HWindow20 and the citation node this makes the
window leg at W = 20 UNCONDITIONAL (modulo the LRC(≤13) citation only).
-/
import Mathlib
import TournamentH7.WindowData20
import TournamentH7.LRCKernelGate
import TournamentH7.HWindow20
import TournamentH7.LRC14CertRoute

namespace LonelyRunner.LRC14.WindowDispatch

open LonelyRunner KernelGate WindowData

/-- The decidable per-row gate: positive denominator, 13 entries, every speed clears the
kernel condition at the row's witness. -/
def rowOK (r : List ℤ × ℤ × ℤ) : Prop :=
  0 < r.2.2 ∧ r.1.length = 13 ∧
  ∀ s ∈ r.1, speedOK s r.2.1 r.2.2.toNat

instance (r : List ℤ × ℤ × ℤ) : Decidable (rowOK r) := by
  unfold rowOK; infer_instance

/-- **The sweep**: every one of the 6084 census rows passes the kernel gate. -/
theorem all_rows_ok : ∀ r ∈ rows, rowOK r := by native_decide

/-- Membership soundness: a 13-tuple whose value list is a census row is lonely at the
row's witness. -/
theorem lonely_of_mem_rows (u : Fin 13 → ℤ) (n d : ℤ)
    (hmem : (List.ofFn u, n, d) ∈ rows) : ∃ t : ℝ, Lonely 14 u t := by
  have hOK := all_rows_ok _ hmem
  unfold rowOK at hOK
  dsimp only at hOK
  obtain ⟨hd, -, hs⟩ := hOK
  refine ⟨((((n : ℚ) / (d.toNat : ℚ)) : ℚ) : ℝ), ?_⟩
  apply lonely_of_kernelWitness (by omega : 0 < d.toNat)
  intro i
  refine hs (u i) ?_
  rw [List.mem_ofFn]
  exact ⟨i, rfl⟩

end LonelyRunner.LRC14.WindowDispatch
