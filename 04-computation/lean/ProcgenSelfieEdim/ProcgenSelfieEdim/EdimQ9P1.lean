import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 3-5 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_3 : eqListR (blk 64 3 (edgeKeys 64 9 q9set)) (blk 64 3 q9keys) = true := by
  decide +kernel

theorem q9_block_4 : eqListR (blk 64 4 (edgeKeys 64 9 q9set)) (blk 64 4 q9keys) = true := by
  decide +kernel

theorem q9_block_5 : eqListR (blk 64 5 (edgeKeys 64 9 q9set)) (blk 64 5 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
