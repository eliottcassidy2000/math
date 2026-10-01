import ProcgenSelfieEdim.EdimQ8Data

set_option autoImplicit false

/-! `Q_8` certificate, key blocks 4-7 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q8_block_4 : eqListR (blk 64 4 (edgeKeys 32 8 q8set)) (blk 64 4 q8keys) = true := by
  decide +kernel

theorem q8_block_5 : eqListR (blk 64 5 (edgeKeys 32 8 q8set)) (blk 64 5 q8keys) = true := by
  decide +kernel

theorem q8_block_6 : eqListR (blk 64 6 (edgeKeys 32 8 q8set)) (blk 64 6 q8keys) = true := by
  decide +kernel

theorem q8_block_7 : eqListR (blk 64 7 (edgeKeys 32 8 q8set)) (blk 64 7 q8keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
