import ProcgenSelfieEdim.EdimQ8Data

set_option autoImplicit false

/-! `Q_8` certificate, key blocks 12-15 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q8_block_12 : eqListR (blk 64 12 (edgeKeys 32 8 q8set)) (blk 64 12 q8keys) = true := by
  decide +kernel

theorem q8_block_13 : eqListR (blk 64 13 (edgeKeys 32 8 q8set)) (blk 64 13 q8keys) = true := by
  decide +kernel

theorem q8_block_14 : eqListR (blk 64 14 (edgeKeys 32 8 q8set)) (blk 64 14 q8keys) = true := by
  decide +kernel

theorem q8_block_15 : eqListR (blk 64 15 (edgeKeys 32 8 q8set)) (blk 64 15 q8keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
