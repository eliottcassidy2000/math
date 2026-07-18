import TournamentH7.LRCPairRatioTau5Replay2

namespace LonelyRunner
namespace LRCPairRatioTau5

set_option maxRecDepth 100000

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_9_true :
    edgeTableCheckBlock 9 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_9_true :
    localTableCheckBlock 9 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_10_true :
    edgeTableCheckBlock 10 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_10_true :
    localTableCheckBlock 10 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_11_true :
    edgeTableCheckBlock 11 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_11_true :
    localTableCheckBlock 11 = true := by
  decide +kernel

end LRCPairRatioTau5
end LonelyRunner
