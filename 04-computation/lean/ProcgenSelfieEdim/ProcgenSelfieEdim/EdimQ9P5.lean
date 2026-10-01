import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 15-17 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_15 : eqListR (blk 64 15 (edgeKeys 64 9 q9set)) (blk 64 15 q9keys) = true := by
  decide +kernel

theorem q9_block_16 : eqListR (blk 64 16 (edgeKeys 64 9 q9set)) (blk 64 16 q9keys) = true := by
  decide +kernel

theorem q9_block_17 : eqListR (blk 64 17 (edgeKeys 64 9 q9set)) (blk 64 17 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
