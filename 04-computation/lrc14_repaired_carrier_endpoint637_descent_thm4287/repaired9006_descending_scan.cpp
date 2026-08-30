// Exact descending scan of THM-4283's repaired 9,006-mask carrier against
// the post-THM-4283 residual.  This primary path inherits the independently
// audited joint-exposure implementation pinned by THM-4283, but changes the
// consumer: every row is now drawn from THM-4283's typed final residual.
//
// The program completes a whole endpoint layer before stopping.  A failure
// therefore establishes a layer boundary, not merely a first-row boundary.

#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

#include <algorithm>
#include <set>
#include <tuple>

namespace {

constexpr std::size_t kCurrentResidualCount = 22682;
constexpr u64 kCurrentResidualFnv = UINT64_C(0xf7563445f15efebf);
constexpr u32 kPriorRepair = UINT32_C(0x014c9084);
constexpr std::size_t kBoundaryWitnessCount = 9;
constexpr u64 kBoundaryWitnessFnv = UINT64_C(0x02b936529030e4bc);
constexpr u64 kPriorCarrierFnv = UINT64_C(0x8e1860a25d0fcf87);
constexpr u64 kRepairedCarrierFnv = UINT64_C(0xfdc1c57ae4dc1bb6);

std::vector<Pair> read_current_band(const std::filesystem::path& path,
                                    int lower_endpoint) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open current residual");
    std::vector<Pair> all;
    std::vector<Pair> band;
    FnvLocal all_ledger;
    int maximum_endpoint = 0;
    std::vector<Pair> top_rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed current-residual row");
        std::size_t used_q = 0;
        std::size_t used_r = 0;
        Pair pair;
        pair.q = std::stoi(line.substr(0, comma), &used_q);
        pair.r = std::stoi(line.substr(comma + 1), &used_r);
        require(used_q == comma && used_r == line.size() - comma - 1 &&
                    pair.q > 0 && pair.q < pair.r,
                "invalid current-residual pair");
        if (!all.empty())
            require(std::tie(all.back().q, all.back().r) <
                        std::tie(pair.q, pair.r),
                    "current-residual order changed");
        all.push_back(pair);
        all_ledger.add(pair.q);
        all_ledger.add(pair.r);
        if (pair.r > maximum_endpoint) {
            maximum_endpoint = pair.r;
            top_rows.assign(1, pair);
        } else if (pair.r == maximum_endpoint) {
            top_rows.push_back(pair);
        }
        if (pair.r >= lower_endpoint) band.push_back(pair);
    }
    require(all.size() == kCurrentResidualCount &&
                all_ledger.state == kCurrentResidualFnv,
            "post-THM4283 residual identity changed");
    require(maximum_endpoint == 637 && top_rows.size() == 3 &&
                top_rows[0].q == 100 && top_rows[0].r == 637 &&
                top_rows[1].q == 294 && top_rows[1].r == 637 &&
                top_rows[2].q == 520 && top_rows[2].r == 637,
            "post-THM4283 residual boundary changed");
    std::sort(band.begin(), band.end(), [](const Pair& left, const Pair& right) {
        if (left.r != right.r) return left.r > right.r;
        return left.q < right.q;
    });
    require(!band.empty(), "selected current-residual band is empty");
    require(band.front().q == 100 && band.front().r == 637,
            "top current-residual row changed");
    return band;
}

u64 pair_list_fnv(const std::vector<Pair>& pairs) {
    FnvLocal ledger;
    for (const Pair pair : pairs) {
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    return ledger.state;
}

u64 body_fnv(const std::vector<u32>& bodies) {
    FnvLocal ledger;
    for (u32 body : bodies) ledger.add(body);
    return ledger.state;
}

u64 local_mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

PairAudit audit_pair_without_response_search(
    const std::vector<Cell>& cells,
    const std::vector<u32>& joint,
    const std::vector<u32>& carrier,
    const std::unordered_set<u32>& joint_set,
    const Pair pair) {
    const ActiveUniverse active = build_active_universe(cells, pair.q, pair.r);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    active_joint.reserve(joint.size());
    active_nonjoint.reserve(carrier.size() - joint.size());
    PairAudit audit;
    FnvLocal active_carrier_ledger;
    FnvLocal inactive_joint_ledger;
    for (std::size_t index = 0; index < joint.size(); ++index) {
        const u32 mask = joint[index];
        if (active.active[colex_rank8_local(mask)]) {
            active_joint.push_back(mask);
        } else {
            ++audit.inactive_joint;
            inactive_joint_ledger.add(index);
            inactive_joint_ledger.add(mask);
        }
    }
    for (u32 mask : carrier) {
        if (!active.active[colex_rank8_local(mask)]) continue;
        ++audit.active_carrier;
        active_carrier_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    audit.active_joint = active_joint.size();
    audit.active_nonjoint = active_nonjoint.size();
    require(audit.active_carrier ==
                audit.active_joint + audit.active_nonjoint,
            "carrier partition changed");

    FnvLocal exposed_ledger;
    FnvLocal response_ledger;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        bool joint_hit = false;
        for (u32 mask : active_joint) {
            ++audit.joint_checks;
            if ((body & mask) == 0) {
                joint_hit = true;
                break;
            }
        }
        if (!joint_hit) {
            ++audit.exposed;
            exposed_ledger.add(body);
            u64 hits = 0;
            u32 least = 0;
            for (u32 mask : active_nonjoint) {
                ++audit.nonjoint_checks;
                if ((body & mask) != 0) continue;
                if (hits == 0 || mask < least) least = mask;
                ++hits;
            }
            response_ledger.add(body);
            response_ledger.add(hits);
            response_ledger.add(least);
            audit.hit_incidences += hits;
            if (hits == 0) {
                ++audit.failures;
                audit.failure_bodies.push_back(body);
            } else {
                if (hits < audit.minimum_hits) {
                    audit.minimum_hits = hits;
                    audit.minimum_hit_body = body;
                }
                if (hits > audit.maximum_hits) {
                    audit.maximum_hits = hits;
                    audit.maximum_hit_body = body;
                }
            }
        }
        if (ordinal + 1 < kBodyCount) body = next_same_popcount(body);
    }
    require(body == UINT32_C(0x3fe00000), "body enumeration changed");
    audit.exposed_fnv = exposed_ledger.state;
    audit.response_fnv = response_ledger.state;
    audit.active_carrier_fnv = active_carrier_ledger.state;
    audit.inactive_joint_fnv = inactive_joint_ledger.state;
    if (audit.exposed == 0) audit.minimum_hits = 0;
    FnvLocal failure_ledger;
    for (u32 failure : audit.failure_bodies) failure_ledger.add(failure);
    audit.failure_fnv = failure_ledger.state;
    return audit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 8,
                "usage: repaired9006-scan JOINT421 BASE8951 ADDITIONS45 "
                "WITNESS9 CURRENT22682 LOWER_ENDPOINT FAILURE_LEDGER");
        init_choose8_local();
        const int lower_endpoint = std::stoi(argv[6]);
        require(lower_endpoint > 0 && lower_endpoint <= 637,
                "lower endpoint must lie in 1..637");

        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        std::vector<u32> carrier =
            read_masks(argv[2], kCarrierCount, kCarrierFnv);
        const std::vector<u32> additions = read_additions(argv[3]);
        const std::vector<u32> boundary_witness =
            read_masks(argv[4], kBoundaryWitnessCount, kBoundaryWitnessFnv);

        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second,
                    "45-mask addition overlaps inherited carrier");
            carrier.push_back(mask);
        }
        require(carrier.size() == kCarrierCount + kAdditionCount &&
                    local_mask_fnv(carrier) == kAugmentedCarrierFnv,
                "8,996-mask carrier identity changed");
        const std::vector<u32> base_carrier = carrier;
        require(std::popcount(kPriorRepair) == 8 &&
                    distinct.insert(kPriorRepair).second,
                "endpoint-644 repair rank or novelty changed");
        carrier.push_back(kPriorRepair);
        require(carrier.size() == 8997 &&
                    local_mask_fnv(carrier) == kPriorCarrierFnv,
                "prior 8,997-mask carrier identity changed");
        const std::vector<u32> prior_carrier = carrier;
        for (u32 mask : boundary_witness) {
            require(distinct.insert(mask).second,
                    "endpoint-638 witness overlaps repaired carrier");
            carrier.push_back(mask);
        }
        require(carrier.size() == 9006 &&
                    local_mask_fnv(carrier) == kRepairedCarrierFnv,
                "repaired 9,006-mask carrier identity changed");

        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            require(distinct.contains(mask),
                    "joint mask absent from repaired carrier");
        const std::vector<Pair> band =
            read_current_band(argv[5], lower_endpoint);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");

        std::ofstream failure_output(argv[7]);
        require(static_cast<bool>(failure_output),
                "cannot create failure ledger");
        failure_output << "q,r,body_hex\n";

        std::cout << "THM4287_REPAIRED9006_DESCENDING_SCAN_V1\n"
                  << "CURRENT_RESIDUAL " << kCurrentResidualCount << " FNV "
                  << std::hex << kCurrentResidualFnv << std::dec
                  << " REQUESTED_LOWER " << lower_endpoint << " SELECTED "
                  << band.size() << " PAIR_FNV " << std::hex
                  << pair_list_fnv(band) << " BASE8996_FNV "
                  << kAugmentedCarrierFnv << " PRIOR8997_FNV "
                  << kPriorCarrierFnv << " CARRIER9006_FNV "
                  << kRepairedCarrierFnv << std::dec << '\n'
                  << "POOL_CELLS " << cells.size() << " BODY_UNIVERSE "
                  << kBodyCount << " REPAIR_UNIVERSE " << EXPECTED_REPAIRS
                  << '\n';

        FnvLocal pair_ledger;
        u64 total_failures = 0;
        std::size_t completed_rows = 0;
        bool stopped = false;
        std::size_t index = 0;
        while (index < band.size() && !stopped) {
            const int endpoint = band[index].r;
            u64 layer_failures = 0;
            std::size_t layer_rows = 0;
            FnvLocal layer_ledger;
            while (index < band.size() && band[index].r == endpoint) {
                const Pair pair = band[index++];
                const PairAudit audit = audit_pair_without_response_search(
                    cells, joint, carrier, joint_set, pair);
                PairAudit base_audit;
                PairAudit prior_audit;
                if (pair.r == 637) {
                    base_audit = audit_pair_without_response_search(
                        cells, joint, base_carrier, joint_set, pair);
                    prior_audit = audit_pair_without_response_search(
                        cells, joint, prior_carrier, joint_set, pair);
                }
                for (u32 body : audit.failure_bodies)
                    failure_output << pair.q << ',' << pair.r << ','
                                   << std::hex << std::setw(8)
                                   << std::setfill('0') << body << std::dec
                                   << std::setfill(' ') << '\n';
                layer_failures += audit.failures;
                total_failures += audit.failures;
                ++layer_rows;
                ++completed_rows;
                layer_ledger.add(pair.q);
                layer_ledger.add(pair.r);
                pair_ledger.add(pair.q);
                pair_ledger.add(pair.r);
                pair_ledger.add(audit.active_carrier);
                pair_ledger.add(audit.active_carrier_fnv);
                pair_ledger.add(audit.active_joint);
                pair_ledger.add(audit.active_nonjoint);
                pair_ledger.add(audit.inactive_joint);
                pair_ledger.add(audit.inactive_joint_fnv);
                pair_ledger.add(audit.exposed);
                pair_ledger.add(audit.exposed_fnv);
                pair_ledger.add(audit.failures);
                pair_ledger.add(audit.failure_fnv);
                pair_ledger.add(audit.hit_incidences);
                pair_ledger.add(audit.minimum_hits);
                pair_ledger.add(audit.minimum_hit_body);
                pair_ledger.add(audit.maximum_hits);
                pair_ledger.add(audit.maximum_hit_body);
                pair_ledger.add(audit.response_fnv);
                pair_ledger.add(pair.r == 637);
                if (pair.r == 637) {
                    pair_ledger.add(base_audit.failures);
                    pair_ledger.add(base_audit.failure_fnv);
                    pair_ledger.add(prior_audit.failures);
                    pair_ledger.add(prior_audit.failure_fnv);
                }
                std::cout << "PAIR " << pair.q << ',' << pair.r
                          << " ACTIVE_CARRIER " << audit.active_carrier
                          << " ACTIVE_CARRIER_FNV " << std::hex
                          << audit.active_carrier_fnv << std::dec
                          << " ACTIVE_JOINT " << audit.active_joint
                          << " ACTIVE_NONJOINT " << audit.active_nonjoint
                          << " INACTIVE_JOINT " << audit.inactive_joint
                          << " INACTIVE_JOINT_FNV " << std::hex
                          << audit.inactive_joint_fnv << std::dec
                          << " EXPOSED " << audit.exposed << " EXPOSED_FNV "
                          << std::hex << audit.exposed_fnv << std::dec
                          << " FAILURES " << audit.failures
                          << " FAILURE_FNV " << std::hex << audit.failure_fnv
                          << std::dec << " HIT_INCIDENCES "
                          << audit.hit_incidences << " POSITIVE_HIT_RANGE "
                          << audit.minimum_hits << ".." << audit.maximum_hits
                          << " MIN_BODY " << std::hex << std::setw(8)
                          << std::setfill('0') << audit.minimum_hit_body
                          << " MAX_BODY " << std::setw(8)
                          << audit.maximum_hit_body << " RESPONSE_FNV "
                          << audit.response_fnv << std::dec
                          << std::setfill(' ');
                if (pair.r == 637)
                    std::cout << " BASE8996_FAILURES " << base_audit.failures
                              << " BASE8996_FAILURE_FNV " << std::hex
                              << base_audit.failure_fnv << std::dec
                              << " PRIOR8997_FAILURES "
                              << prior_audit.failures
                              << " PRIOR8997_FAILURE_FNV " << std::hex
                              << prior_audit.failure_fnv << std::dec;
                std::cout << '\n';
            }
            std::cout << "LAYER " << endpoint << " ROWS " << layer_rows
                      << " FNV " << std::hex << layer_ledger.state << std::dec
                      << " FAILURES " << layer_failures << '\n';
            stopped = layer_failures != 0;
        }
        require(failure_output.good(), "failure-ledger write failed");
        std::cout << "PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << " COMPLETED_ROWS " << completed_rows
                  << " TOTAL_FAILURES " << total_failures
                  << " STOPPED_AT_FAILURE " << stopped << '\n'
                  << "SCOPE FIXED_POOL_TYPED_RESIDUAL_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_DESCENDING_SCAN\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4287_REPAIRED9006_SCAN_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
