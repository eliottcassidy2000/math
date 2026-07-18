import TournamentH7.LRCPairRatioTau5Replay1

namespace LonelyRunner
namespace LRCPairRatioTau5

set_option maxRecDepth 100000

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_6_true :
    edgeTableCheckBlock 6 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_6_true :
    localTableCheckBlock 6 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_7_true :
    edgeTableCheckBlock 7 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_7_true :
    localTableCheckBlock 7 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem edgeTableCheckBlock_8_true :
    edgeTableCheckBlock 8 = true := by
  decide +kernel

set_option maxHeartbeats 50000000 in
theorem localTableCheckBlock_8_true :
    localTableCheckBlock 8 = true := by
  decide +kernel

end LRCPairRatioTau5
end LonelyRunner
