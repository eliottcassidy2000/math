import ProcgenSelfieEdim.EdimQ8Data

set_option autoImplicit false

/-! `Q_8` certificate, key blocks 0-3 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q8_block_0 : eqListR (blk 64 0 (edgeKeys 32 8 q8set)) (blk 64 0 q8keys) = true := by
  decide +kernel

theorem q8_block_1 : eqListR (blk 64 1 (edgeKeys 32 8 q8set)) (blk 64 1 q8keys) = true := by
  decide +kernel

theorem q8_block_2 : eqListR (blk 64 2 (edgeKeys 32 8 q8set)) (blk 64 2 q8keys) = true := by
  decide +kernel

theorem q8_block_3 : eqListR (blk 64 3 (edgeKeys 32 8 q8set)) (blk 64 3 q8keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
