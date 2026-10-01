import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 0-2 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_0 : eqListR (blk 64 0 (edgeKeys 64 9 q9set)) (blk 64 0 q9keys) = true := by
  decide +kernel

theorem q9_block_1 : eqListR (blk 64 1 (edgeKeys 64 9 q9set)) (blk 64 1 q9keys) = true := by
  decide +kernel

theorem q9_block_2 : eqListR (blk 64 2 (edgeKeys 64 9 q9set)) (blk 64 2 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
