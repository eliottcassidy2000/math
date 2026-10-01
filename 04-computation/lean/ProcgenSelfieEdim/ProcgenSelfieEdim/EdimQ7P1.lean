import ProcgenSelfieEdim.EdimQ7Data

set_option autoImplicit false

/-! `Q_7` certificate, key blocks 4-6 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q7_block_4 : eqListR (blk 64 4 (edgeKeys 32 7 q7set)) (blk 64 4 q7keys) = true := by
  decide +kernel

theorem q7_block_5 : eqListR (blk 64 5 (edgeKeys 32 7 q7set)) (blk 64 5 q7keys) = true := by
  decide +kernel

theorem q7_block_6 : eqListR (blk 64 6 (edgeKeys 32 7 q7set)) (blk 64 6 q7keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
