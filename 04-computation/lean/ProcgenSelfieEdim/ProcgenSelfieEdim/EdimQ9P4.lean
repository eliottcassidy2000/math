import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 12-14 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_12 : eqListR (blk 64 12 (edgeKeys 64 9 q9set)) (blk 64 12 q9keys) = true := by
  decide +kernel

theorem q9_block_13 : eqListR (blk 64 13 (edgeKeys 64 9 q9set)) (blk 64 13 q9keys) = true := by
  decide +kernel

theorem q9_block_14 : eqListR (blk 64 14 (edgeKeys 64 9 q9set)) (blk 64 14 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
