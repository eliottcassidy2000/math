import TournamentH7.LRCPairRatioTau5Certificate

namespace LonelyRunner
namespace LRCPairRatioTau5

set_option maxRecDepth 100000

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_zero_true :
    edgeTableCheckBlock 0 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_zero_true :
    localTableCheckBlock 0 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_1_true :
    edgeTableCheckBlock 1 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_1_true :
    localTableCheckBlock 1 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_2_true :
    edgeTableCheckBlock 2 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_2_true :
    localTableCheckBlock 2 = true := by
  decide +kernel

end LRCPairRatioTau5
end LonelyRunner
