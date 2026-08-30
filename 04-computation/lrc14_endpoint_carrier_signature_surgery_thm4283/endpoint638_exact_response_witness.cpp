// Serialize the exact nine-mask response witness at the first THM-4283
// carrier boundary.  This is a local endpoint certificate only; no lower
// layer or proof-graph consequence is inferred here.

#include "carrier_scan_support.cpp"

#include <set>

namespace {

constexpr u32 kPriorRepair = UINT32_C(0x014c9084);
constexpr u64 kNested8997Fnv = UINT64_C(0x8e1860a25d0fcf87);
constexpr u64 kFailureFnv638 = UINT64_C(0x917d107c4536efc9);

u64 local_mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 5,
                "usage: endpoint638-witness JOINT421 BASE8951 ADDITIONS45 "
                "WITNESS_OUT");
        init_choose8_local();
        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        std::vector<u32> carrier =
            read_masks(argv[2], kCarrierCount, kCarrierFnv);
        const std::vector<u32> additions = read_additions(argv[3]);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second,
                    "carrier addition overlaps inherited carrier");
            carrier.push_back(mask);
        }
        require(carrier.size() == kCarrierCount + kAdditionCount &&
                    local_mask_fnv(carrier) == kAugmentedCarrierFnv,
                "base8996 identity changed");
        require(distinct.insert(kPriorRepair).second,
                "prior endpoint-644 repair overlaps carrier");
        carrier.push_back(kPriorRepair);
        require(carrier.size() == 8997 &&
                    local_mask_fnv(carrier) == kNested8997Fnv,
                "nested8997 identity changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const Pair pair{256, 638};
        const PairAudit boundary =
            audit_pair(cells, joint, carrier, joint_set, pair);
        require(boundary.failures == 40 &&
                    boundary.failure_fnv == kFailureFnv638 &&
                    boundary.complete_response_classes == 315 &&
                    boundary.full_response_masks == 0 &&
                    boundary.exact_minimum_replacements == 9 &&
                    boundary.exact_witness_masks.size() == 9,
                "endpoint-638 exact response boundary changed");
        for (u32 mask : boundary.exact_witness_masks)
            require(distinct.insert(mask).second,
                    "exact response witness overlaps nested carrier");

        std::ofstream witness_output(argv[4]);
        require(static_cast<bool>(witness_output),
                "cannot create endpoint-638 witness ledger");
        for (u32 mask : boundary.exact_witness_masks)
            witness_output << std::hex << std::setw(8) << std::setfill('0')
                           << mask << std::dec << std::setfill(' ') << '\n';
        require(witness_output.good(), "endpoint-638 witness write failed");

        std::vector<u32> repaired = carrier;
        repaired.insert(repaired.end(), boundary.exact_witness_masks.begin(),
                        boundary.exact_witness_masks.end());
        const PairAudit replay =
            audit_pair(cells, joint, repaired, joint_set, pair);
        require(replay.failures == 0,
                "nine-mask endpoint-638 response does not close pair");

        std::cout << "THM4283_ENDPOINT638_EXACT_RESPONSE_WITNESS_V1\n"
                  << "PAIR 256,638 BASE8997 " << carrier.size() << " FNV "
                  << std::hex << local_mask_fnv(carrier) << std::dec
                  << " FAILURES " << boundary.failures << " FAILURE_FNV "
                  << std::hex << boundary.failure_fnv << std::dec
                  << " RESPONSE_CLASSES "
                  << boundary.complete_response_classes
                  << " FULL_RESPONDERS " << boundary.full_response_masks
                  << " EXACT_MINIMUM "
                  << boundary.exact_minimum_replacements << '\n'
                  << "WITNESSES " << boundary.exact_witness_masks.size()
                  << " FNV " << std::hex
                  << local_mask_fnv(boundary.exact_witness_masks)
                  << " MASKS";
        for (u32 mask : boundary.exact_witness_masks)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << mask;
        std::cout << std::dec << std::setfill(' ') << '\n'
                  << "REPAIRED9006_FNV " << std::hex
                  << local_mask_fnv(repaired) << std::dec
                  << " REPLAY_FAILURES " << replay.failures << '\n'
                  << "SCOPE LOCAL_PAIR_RESPONSE_ONLY_NO_LOWER_LAYER_CLAIM\n"
                  << "VERDICT PASS EXACT_MINIMUM_NINE_AND_EXPLICIT_WITNESS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_ENDPOINT638_WITNESS_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
