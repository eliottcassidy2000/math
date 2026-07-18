import TournamentH7.LRCPairRatioTau5Replay0

namespace LonelyRunner
namespace LRCPairRatioTau5

set_option maxRecDepth 100000

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_3_true :
    edgeTableCheckBlock 3 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_3_true :
    localTableCheckBlock 3 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_4_true :
    edgeTableCheckBlock 4 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_4_true :
    localTableCheckBlock 4 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_5_true :
    edgeTableCheckBlock 5 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_5_true :
    localTableCheckBlock 5 = true := by
  decide +kernel

end LRCPairRatioTau5
end LonelyRunner
