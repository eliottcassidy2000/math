import ProcgenSelfieEdim.EdimQ7Data

set_option autoImplicit false

/-! `Q_7` certificate, key blocks 0-3 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q7_block_0 : eqListR (blk 64 0 (edgeKeys 32 7 q7set)) (blk 64 0 q7keys) = true := by
  decide +kernel

theorem q7_block_1 : eqListR (blk 64 1 (edgeKeys 32 7 q7set)) (blk 64 1 q7keys) = true := by
  decide +kernel

theorem q7_block_2 : eqListR (blk 64 2 (edgeKeys 32 7 q7set)) (blk 64 2 q7keys) = true := by
  decide +kernel

theorem q7_block_3 : eqListR (blk 64 3 (edgeKeys 32 7 q7set)) (blk 64 3 q7keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
