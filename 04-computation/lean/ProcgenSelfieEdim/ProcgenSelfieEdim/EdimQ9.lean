import ProcgenSelfieEdim.EdimQ9P0
import ProcgenSelfieEdim.EdimQ9P1
import ProcgenSelfieEdim.EdimQ9P2
import ProcgenSelfieEdim.EdimQ9P3
import ProcgenSelfieEdim.EdimQ9P4
import ProcgenSelfieEdim.EdimQ9P5
import ProcgenSelfieEdim.EdimQ9P6
import ProcgenSelfieEdim.EdimQ9P7
import ProcgenSelfieEdim.EdimQ9P8
import ProcgenSelfieEdim.EdimQ9P9
import ProcgenSelfieEdim.EdimQ9P10
import ProcgenSelfieEdim.EdimQ9P11

set_option autoImplicit false

/-!
# THM-4525: `edim_m(Q_9) ≤ 38`

The 2304 edge keys of the explicit 38-set (`q9set`) equal the precomputed numbers
(`q9_block_0` … `q9_block_35`), which are pairwise distinct (verified merge sort,
12 rounds, then a strict-increase check). Hence the set is edge-multiset resolving.
-/

namespace ProcgenSelfieEdim

theorem q9_distinct : strictIncR (msortR 12 q9keys) = true := by
  decide +kernel

theorem q9_keys_eq : edgeKeys 64 9 q9set = q9keys := by
  apply eq_of_blocks 64 36 _ _ (by rw [q9_lengths.1]; decide) (by rw [q9_lengths.2]; decide)
  intro i hi
  match i, hi with
  | 0, _ => exact eq_of_eqListR _ _ q9_block_0
  | 1, _ => exact eq_of_eqListR _ _ q9_block_1
  | 2, _ => exact eq_of_eqListR _ _ q9_block_2
  | 3, _ => exact eq_of_eqListR _ _ q9_block_3
  | 4, _ => exact eq_of_eqListR _ _ q9_block_4
  | 5, _ => exact eq_of_eqListR _ _ q9_block_5
  | 6, _ => exact eq_of_eqListR _ _ q9_block_6
  | 7, _ => exact eq_of_eqListR _ _ q9_block_7
  | 8, _ => exact eq_of_eqListR _ _ q9_block_8
  | 9, _ => exact eq_of_eqListR _ _ q9_block_9
  | 10, _ => exact eq_of_eqListR _ _ q9_block_10
  | 11, _ => exact eq_of_eqListR _ _ q9_block_11
  | 12, _ => exact eq_of_eqListR _ _ q9_block_12
  | 13, _ => exact eq_of_eqListR _ _ q9_block_13
  | 14, _ => exact eq_of_eqListR _ _ q9_block_14
  | 15, _ => exact eq_of_eqListR _ _ q9_block_15
  | 16, _ => exact eq_of_eqListR _ _ q9_block_16
  | 17, _ => exact eq_of_eqListR _ _ q9_block_17
  | 18, _ => exact eq_of_eqListR _ _ q9_block_18
  | 19, _ => exact eq_of_eqListR _ _ q9_block_19
  | 20, _ => exact eq_of_eqListR _ _ q9_block_20
  | 21, _ => exact eq_of_eqListR _ _ q9_block_21
  | 22, _ => exact eq_of_eqListR _ _ q9_block_22
  | 23, _ => exact eq_of_eqListR _ _ q9_block_23
  | 24, _ => exact eq_of_eqListR _ _ q9_block_24
  | 25, _ => exact eq_of_eqListR _ _ q9_block_25
  | 26, _ => exact eq_of_eqListR _ _ q9_block_26
  | 27, _ => exact eq_of_eqListR _ _ q9_block_27
  | 28, _ => exact eq_of_eqListR _ _ q9_block_28
  | 29, _ => exact eq_of_eqListR _ _ q9_block_29
  | 30, _ => exact eq_of_eqListR _ _ q9_block_30
  | 31, _ => exact eq_of_eqListR _ _ q9_block_31
  | 32, _ => exact eq_of_eqListR _ _ q9_block_32
  | 33, _ => exact eq_of_eqListR _ _ q9_block_33
  | 34, _ => exact eq_of_eqListR _ _ q9_block_34
  | 35, _ => exact eq_of_eqListR _ _ q9_block_35
  | _ + 36, h => exact absurd h (by omega)

/-- **THM-4525:** the explicit 38-set is edge-multiset resolving in `Q_9`, so
`edim_m(Q_9) ≤ 38` (the paper's certificate had 246 vertices). -/
theorem q9set_resolving : Resolving 9 q9set := by
  apply resolving_of_edgeKeys 64 9 q9set
  rw [q9_keys_eq]
  exact nodup_of_msortR 12 q9keys q9_distinct

end ProcgenSelfieEdim
