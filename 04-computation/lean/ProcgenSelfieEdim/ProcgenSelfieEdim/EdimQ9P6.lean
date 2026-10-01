import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 18-20 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_18 : eqListR (blk 64 18 (edgeKeys 64 9 q9set)) (blk 64 18 q9keys) = true := by
  decide +kernel

theorem q9_block_19 : eqListR (blk 64 19 (edgeKeys 64 9 q9set)) (blk 64 19 q9keys) = true := by
  decide +kernel

theorem q9_block_20 : eqListR (blk 64 20 (edgeKeys 64 9 q9set)) (blk 64 20 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
