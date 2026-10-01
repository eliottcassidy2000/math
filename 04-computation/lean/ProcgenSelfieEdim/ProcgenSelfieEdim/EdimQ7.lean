import ProcgenSelfieEdim.EdimQ7P0
import ProcgenSelfieEdim.EdimQ7P1

set_option autoImplicit false

/-!
# THM-4525: `edim_m(Q_7) ≤ 19`

The 448 edge keys of the explicit 19-set (`q7set`) equal the precomputed numbers
(`q7_block_0` … `q7_block_6`), which are pairwise distinct (verified merge sort,
9 rounds, then a strict-increase check). Hence the set is edge-multiset resolving.
-/

namespace ProcgenSelfieEdim

theorem q7_distinct : strictIncR (msortR 9 q7keys) = true := by
  decide +kernel

theorem q7_keys_eq : edgeKeys 32 7 q7set = q7keys := by
  apply eq_of_blocks 64 7 _ _ (by rw [q7_lengths.1]; decide) (by rw [q7_lengths.2]; decide)
  intro i hi
  match i, hi with
  | 0, _ => exact eq_of_eqListR _ _ q7_block_0
  | 1, _ => exact eq_of_eqListR _ _ q7_block_1
  | 2, _ => exact eq_of_eqListR _ _ q7_block_2
  | 3, _ => exact eq_of_eqListR _ _ q7_block_3
  | 4, _ => exact eq_of_eqListR _ _ q7_block_4
  | 5, _ => exact eq_of_eqListR _ _ q7_block_5
  | 6, _ => exact eq_of_eqListR _ _ q7_block_6
  | _ + 7, h => exact absurd h (by omega)

/-- **THM-4525:** the explicit 19-set is edge-multiset resolving in `Q_7`, so
`edim_m(Q_7) ≤ 19` (the paper's certificate had 63 vertices). -/
theorem q7set_resolving : Resolving 7 q7set := by
  apply resolving_of_edgeKeys 32 7 q7set
  rw [q7_keys_eq]
  exact nodup_of_msortR 9 q7keys q7_distinct

end ProcgenSelfieEdim
