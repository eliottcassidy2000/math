import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 21-23 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_21 : eqListR (blk 64 21 (edgeKeys 64 9 q9set)) (blk 64 21 q9keys) = true := by
  decide +kernel

theorem q9_block_22 : eqListR (blk 64 22 (edgeKeys 64 9 q9set)) (blk 64 22 q9keys) = true := by
  decide +kernel

theorem q9_block_23 : eqListR (blk 64 23 (edgeKeys 64 9 q9set)) (blk 64 23 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
