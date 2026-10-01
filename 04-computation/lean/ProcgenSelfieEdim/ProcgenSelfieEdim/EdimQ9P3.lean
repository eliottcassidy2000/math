import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 9-11 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_9 : eqListR (blk 64 9 (edgeKeys 64 9 q9set)) (blk 64 9 q9keys) = true := by
  decide +kernel

theorem q9_block_10 : eqListR (blk 64 10 (edgeKeys 64 9 q9set)) (blk 64 10 q9keys) = true := by
  decide +kernel

theorem q9_block_11 : eqListR (blk 64 11 (edgeKeys 64 9 q9set)) (blk 64 11 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
