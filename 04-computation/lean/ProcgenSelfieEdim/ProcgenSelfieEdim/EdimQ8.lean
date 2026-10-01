import ProcgenSelfieEdim.EdimQ8P0
import ProcgenSelfieEdim.EdimQ8P1
import ProcgenSelfieEdim.EdimQ8P2
import ProcgenSelfieEdim.EdimQ8P3

set_option autoImplicit false

/-!
# THM-4525: `edim_m(Q_8) ≤ 26`

The 1024 edge keys of the explicit 26-set (`q8set`) equal the precomputed numbers
(`q8_block_0` … `q8_block_15`), which are pairwise distinct (verified merge sort,
10 rounds, then a strict-increase check). Hence the set is edge-multiset resolving.
-/

namespace ProcgenSelfieEdim

theorem q8_distinct : strictIncR (msortR 10 q8keys) = true := by
  decide +kernel

theorem q8_keys_eq : edgeKeys 32 8 q8set = q8keys := by
  apply eq_of_blocks 64 16 _ _ (by rw [q8_lengths.1]; decide) (by rw [q8_lengths.2]; decide)
  intro i hi
  match i, hi with
  | 0, _ => exact eq_of_eqListR _ _ q8_block_0
  | 1, _ => exact eq_of_eqListR _ _ q8_block_1
  | 2, _ => exact eq_of_eqListR _ _ q8_block_2
  | 3, _ => exact eq_of_eqListR _ _ q8_block_3
  | 4, _ => exact eq_of_eqListR _ _ q8_block_4
  | 5, _ => exact eq_of_eqListR _ _ q8_block_5
  | 6, _ => exact eq_of_eqListR _ _ q8_block_6
  | 7, _ => exact eq_of_eqListR _ _ q8_block_7
  | 8, _ => exact eq_of_eqListR _ _ q8_block_8
  | 9, _ => exact eq_of_eqListR _ _ q8_block_9
  | 10, _ => exact eq_of_eqListR _ _ q8_block_10
  | 11, _ => exact eq_of_eqListR _ _ q8_block_11
  | 12, _ => exact eq_of_eqListR _ _ q8_block_12
  | 13, _ => exact eq_of_eqListR _ _ q8_block_13
  | 14, _ => exact eq_of_eqListR _ _ q8_block_14
  | 15, _ => exact eq_of_eqListR _ _ q8_block_15
  | _ + 16, h => exact absurd h (by omega)

/-- **THM-4525:** the explicit 26-set is edge-multiset resolving in `Q_8`, so
`edim_m(Q_8) ≤ 26` (the paper's certificate had 115 vertices). -/
theorem q8set_resolving : Resolving 8 q8set := by
  apply resolving_of_edgeKeys 32 8 q8set
  rw [q8_keys_eq]
  exact nodup_of_msortR 10 q8keys q8_distinct

end ProcgenSelfieEdim
