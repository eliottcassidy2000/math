import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 6-8 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_6 : eqListR (blk 64 6 (edgeKeys 64 9 q9set)) (blk 64 6 q9keys) = true := by
  decide +kernel

theorem q9_block_7 : eqListR (blk 64 7 (edgeKeys 64 9 q9set)) (blk 64 7 q9keys) = true := by
  decide +kernel

theorem q9_block_8 : eqListR (blk 64 8 (edgeKeys 64 9 q9set)) (blk 64 8 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
