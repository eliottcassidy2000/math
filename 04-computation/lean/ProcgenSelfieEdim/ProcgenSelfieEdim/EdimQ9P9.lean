import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 27-29 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_27 : eqListR (blk 64 27 (edgeKeys 64 9 q9set)) (blk 64 27 q9keys) = true := by
  decide +kernel

theorem q9_block_28 : eqListR (blk 64 28 (edgeKeys 64 9 q9set)) (blk 64 28 q9keys) = true := by
  decide +kernel

theorem q9_block_29 : eqListR (blk 64 29 (edgeKeys 64 9 q9set)) (blk 64 29 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
