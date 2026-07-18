import TournamentH7.LRCPairRatioTau5Replay3

namespace LonelyRunner
namespace LRCPairRatioTau5

set_option maxRecDepth 100000

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_12_true :
    edgeTableCheckBlock 12 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_12_true :
    localTableCheckBlock 12 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_13_true :
    edgeTableCheckBlock 13 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_13_true :
    localTableCheckBlock 13 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_14_true :
    edgeTableCheckBlock 14 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_14_true :
    localTableCheckBlock 14 = true := by
  decide +kernel

end LRCPairRatioTau5
end LonelyRunner
