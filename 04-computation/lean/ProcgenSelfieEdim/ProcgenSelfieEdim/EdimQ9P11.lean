import ProcgenSelfieEdim.EdimQ9Data

set_option autoImplicit false

/-! `Q_9` certificate, key blocks 33-35 (one kernel check per block). -/

namespace ProcgenSelfieEdim

theorem q9_block_33 : eqListR (blk 64 33 (edgeKeys 64 9 q9set)) (blk 64 33 q9keys) = true := by
  decide +kernel

theorem q9_block_34 : eqListR (blk 64 34 (edgeKeys 64 9 q9set)) (blk 64 34 q9keys) = true := by
  decide +kernel

theorem q9_block_35 : eqListR (blk 64 35 (edgeKeys 64 9 q9set)) (blk 64 35 q9keys) = true := by
  decide +kernel

end ProcgenSelfieEdim
