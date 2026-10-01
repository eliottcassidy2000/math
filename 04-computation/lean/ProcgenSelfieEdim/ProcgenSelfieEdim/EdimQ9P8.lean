import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 24-26 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_24 : eqListR (blk 64 24 (edgeKeys 64 9 q9set)) (blk 64 24 q9keys) = true := by
  decide +kernel

theorem q9_block_25 : eqListR (blk 64 25 (edgeKeys 64 9 q9set)) (blk 64 25 q9keys) = true := by
  decide +kernel

theorem q9_block_26 : eqListR (blk 64 26 (edgeKeys 64 9 q9set)) (blk 64 26 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
