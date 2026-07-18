import TournamentH7.LRCPairRatioTau5Replay4

namespace LonelyRunner
namespace LRCPairRatioTau5

set_option maxRecDepth 100000

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_15_true :
    edgeTableCheckBlock 15 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_15_true :
    localTableCheckBlock 15 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_16_true :
    edgeTableCheckBlock 16 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_16_true :
    localTableCheckBlock 16 = true := by
  decide +kernel

end LRCPairRatioTau5
end LonelyRunner
