import ProcgenSelfieEdim.EdimQ8Data

set_option autoImplicit false

/-! `Q_8` certificate, key blocks 8-11 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q8_block_8 : eqListR (blk 64 8 (edgeKeys 32 8 q8set)) (blk 64 8 q8keys) = true := by
  decide +kernel

theorem q8_block_9 : eqListR (blk 64 9 (edgeKeys 32 8 q8set)) (blk 64 9 q8keys) = true := by
  decide +kernel

theorem q8_block_10 : eqListR (blk 64 10 (edgeKeys 32 8 q8set)) (blk 64 10 q8keys) = true := by
  decide +kernel

theorem q8_block_11 : eqListR (blk 64 11 (edgeKeys 32 8 q8set)) (blk 64 11 q8keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
