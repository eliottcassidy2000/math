import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 30-32 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_30 : eqListR (blk 64 30 (edgeKeys 64 9 q9set)) (blk 64 30 q9keys) = true := by
  decide +kernel

theorem q9_block_31 : eqListR (blk 64 31 (edgeKeys 64 9 q9set)) (blk 64 31 q9keys) = true := by
  decide +kernel

theorem q9_block_32 : eqListR (blk 64 32 (edgeKeys 64 9 q9set)) (blk 64 32 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
