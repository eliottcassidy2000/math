// Exact carrier-via-joint-exposure scan for THM-4282.
//
// The 421-mask THM-4281 joint deck is placed first as a fast body-cover
// sidecar.  At each requested pair we enumerate every labelled nine-body.
// A body missed by the active part of the joint deck is an exact exposed
// obligation; we then test every active nonjoint mask of the frozen final
// 8,951-mask THM-4281 carrier on that obligation.  Reordering is used only
// for the body-cover test and cannot change its truth value.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <bit>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <unordered_set>

namespace {

constexpr std::size_t kJointCount = 421;
constexpr std::size_t kCarrierCount = 8951;
constexpr std::size_t kAdditionCount = 45;
constexpr std::size_t kBandPairCount = 90;
constexpr std::size_t kBodyCount = 14307150;
constexpr u64 kJointFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr u64 kCarrierFnv = UINT64_C(0x188f82ab9dd1695a);
constexpr u64 kAdditionFnv = UINT64_C(0xec083b65cc8c34e3);
constexpr u64 kAugmentedCarrierFnv = UINT64_C(0xfd899660f14b311c);
constexpr u64 kBandPairFnv = UINT64_C(0x942995bee7469430);

struct Pair {
    int q = 0;
    int r = 0;
};

std::vector<u32> read_masks(const std::filesystem::path& path,
                            std::size_t expected_count,
                            u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    FnvLocal ledger;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "mask rank/distinctness changed");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(masks.size() == expected_count && ledger.state == expected_fnv,
            "mask ledger identity changed");
    return masks;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_additions(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open carrier additions");
    std::vector<u32> masks;
    std::set<u32> distinct;
    FnvLocal ledger;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid addition token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "addition rank/distinctness changed");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(masks.size() == kAdditionCount && ledger.state == kAdditionFnv,
            "canonical 45-mask addition ledger identity changed");
    return masks;
}

std::vector<Pair> read_pairs(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair ledger");
    std::vector<Pair> pairs;
    std::set<std::pair<int, int>> distinct;
    FnvLocal ledger;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        Pair pair;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "invalid pair delimiter");
        std::size_t used_q = 0;
        std::size_t used_r = 0;
        pair.q = std::stoi(line.substr(0, comma), &used_q);
        pair.r = std::stoi(line.substr(comma + 1), &used_r);
        require(used_q == comma && used_r == line.size() - comma - 1 &&
                    pair.q > 0 && pair.q < pair.r &&
                    distinct.insert({pair.q, pair.r}).second,
                "invalid or duplicate pair row");
        pairs.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(pairs.size() == kBandPairCount && ledger.state == kBandPairFnv,
            "canonical 90-row pair ledger identity changed");
    return pairs;
}

u32 next_same_popcount(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

struct PairAudit {
    u64 active_carrier = 0;
    u64 active_joint = 0;
    u64 active_nonjoint = 0;
    u64 inactive_joint = 0;
    u64 joint_checks = 0;
    u64 nonjoint_checks = 0;
    u64 exposed = 0;
    u64 exposed_fnv = 0;
    u64 failures = 0;
    u64 hit_incidences = 0;
    u64 minimum_hits = UINT64_MAX;
    u64 maximum_hits = 0;
    u32 minimum_hit_body = 0;
    u32 maximum_hit_body = 0;
    u64 response_fnv = 0;
    u64 active_carrier_fnv = 0;
    u64 inactive_joint_fnv = 0;
    std::vector<u32> failure_bodies;
    u64 failure_fnv = 0;
    u64 complete_response_classes = 0;
    u64 complete_response_fnv = 0;
    u64 complete_active_response_fnv = 0;
    u64 full_response_masks = 0;
    u32 least_full_response = 0;
    u64 exact_minimum_replacements = 0;
    std::vector<u32> exact_witness_masks;
};

PairAudit audit_pair(const std::vector<Cell>& cells,
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
    if (!audit.failure_bodies.empty() && audit.failure_bodies.size() < 64) {
        struct ResponseClass {
            u64 count = 0;
            u32 least = 0;
        };
        std::map<u64, ResponseClass> classes;
        FnvLocal active_response_ledger;
        const u64 full = (UINT64_C(1) << audit.failure_bodies.size()) - 1;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if (!active.active[rank]) continue;
            const u32 mask = unrank_colex8(rank);
            u64 response = 0;
            for (std::size_t index = 0;
                 index < audit.failure_bodies.size(); ++index)
                if ((mask & audit.failure_bodies[index]) == 0)
                    response |= UINT64_C(1) << index;
            active_response_ledger.add(mask);
            active_response_ledger.add(response);
            ResponseClass& response_class = classes[response];
            ++response_class.count;
            if (response_class.count == 1 || mask < response_class.least)
                response_class.least = mask;
        }
        FnvLocal class_ledger;
        for (const auto& [response, response_class] : classes) {
            class_ledger.add(response);
            class_ledger.add(response_class.count);
            class_ledger.add(response_class.least);
        }
        audit.complete_response_classes = classes.size();
        audit.complete_response_fnv = class_ledger.state;
        audit.complete_active_response_fnv = active_response_ledger.state;
        const auto full_class = classes.find(full);
        if (full_class != classes.end()) {
            audit.full_response_masks = full_class->second.count;
            audit.least_full_response = full_class->second.least;
            require(std::find(carrier.begin(), carrier.end(),
                              audit.least_full_response) == carrier.end(),
                    "full responder unexpectedly already in failed carrier");
        }

        struct Predecessor {
            u64 previous = 0;
            u64 pattern = 0;
        };
        std::vector<u64> patterns;
        for (const auto& [pattern, response_class] : classes)
            if (pattern != 0) patterns.push_back(pattern);
        std::sort(patterns.begin(), patterns.end(), [&](u64 left, u64 right) {
            const unsigned left_cover = std::popcount(left);
            const unsigned right_cover = std::popcount(right);
            if (left_cover != right_cover) return left_cover > right_cover;
            return classes.at(left).least < classes.at(right).least;
        });
        std::vector<u64> maximal;
        for (u64 pattern : patterns) {
            bool dominated = false;
            for (u64 prior : maximal)
                if ((pattern & ~prior) == 0) {
                    dominated = true;
                    break;
                }
            if (!dominated) maximal.push_back(pattern);
        }
        std::unordered_map<u64, Predecessor> predecessor;
        predecessor.emplace(0, Predecessor{});
        std::vector<u64> frontier{0};
        for (u64 depth = 1; !frontier.empty() && !predecessor.contains(full);
             ++depth) {
            std::vector<u64> next;
            for (u64 covered : frontier) {
                for (u64 pattern : maximal) {
                    const u64 joined = covered | pattern;
                    if (joined == covered || predecessor.contains(joined))
                        continue;
                    predecessor.emplace(joined,
                                        Predecessor{covered, pattern});
                    next.push_back(joined);
                }
            }
            frontier.swap(next);
            if (predecessor.contains(full))
                audit.exact_minimum_replacements = depth;
        }
        if (audit.exact_minimum_replacements != 0) {
            u64 state = full;
            while (state != 0) {
                const Predecessor step = predecessor.at(state);
                audit.exact_witness_masks.push_back(
                    classes.at(step.pattern).least);
                state = step.previous;
            }
            std::sort(audit.exact_witness_masks.begin(),
                      audit.exact_witness_masks.end());
        }
    }
    return audit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 5 || argc == 6,
                "usage: final8951-joint-exposure JOINT_DECK FINAL_CARRIER "
                "CANONICAL_90_PAIRS CANONICAL_45_ADDITIONS [FAILURES_OUT]");
        init_choose8_local();
        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        std::vector<u32> carrier =
            read_masks(argv[2], kCarrierCount, kCarrierFnv);
        std::set<u32> carrier_set(carrier.begin(), carrier.end());
        const std::vector<u32> additions = read_additions(argv[4]);
        for (u32 mask : additions) {
            require(carrier_set.insert(mask).second,
                    "carrier addition overlaps prior carrier");
            carrier.push_back(mask);
        }
        require(carrier.size() == kCarrierCount + kAdditionCount &&
                    mask_fnv(carrier) == kAugmentedCarrierFnv,
                "canonical augmented carrier identity changed");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        require(joint_set.size() == joint.size(), "joint set changed");
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(),
                    "joint mask absent from final carrier");
        const std::vector<Pair> pairs = read_pairs(argv[3]);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");

        std::cout << "THM4282_FINAL8996_JOINT_EXPOSURE_SCAN_V2\n"
                  << "JOINT " << joint.size() << " FNV " << std::hex
                  << kJointFnv << std::dec << " BASE_CARRIER "
                  << kCarrierCount << " FNV " << std::hex << kCarrierFnv
                  << std::dec << " ADDITIONS " << additions.size()
                  << " AUGMENTED_CARRIER " << carrier.size() << " FNV "
                  << std::hex << mask_fnv(carrier) << std::dec << " PAIRS "
                  << pairs.size() << '\n';
        FnvLocal pair_ledger;
        u64 total_failures = 0;
        std::ofstream failure_output;
        if (argc == 6) {
            failure_output.open(argv[5]);
            require(static_cast<bool>(failure_output),
                    "cannot create failure ledger");
            failure_output << "q,r,body_hex\n";
        }
        for (const Pair pair : pairs) {
            const PairAudit audit =
                audit_pair(cells, joint, carrier, joint_set, pair);
            total_failures += audit.failures;
            for (u32 body : audit.failure_bodies)
                if (failure_output.is_open())
                    failure_output << pair.q << ',' << pair.r << ','
                                   << std::hex << std::setw(8)
                                   << std::setfill('0') << body << std::dec
                                   << std::setfill(' ') << '\n';
            pair_ledger.add(pair.q);
            pair_ledger.add(pair.r);
            pair_ledger.add(audit.active_carrier);
            pair_ledger.add(audit.active_carrier_fnv);
            pair_ledger.add(audit.active_joint);
            pair_ledger.add(audit.inactive_joint);
            pair_ledger.add(audit.inactive_joint_fnv);
            pair_ledger.add(audit.exposed);
            pair_ledger.add(audit.exposed_fnv);
            pair_ledger.add(audit.failures);
            pair_ledger.add(audit.hit_incidences);
            pair_ledger.add(audit.minimum_hits);
            pair_ledger.add(audit.minimum_hit_body);
            pair_ledger.add(audit.maximum_hits);
            pair_ledger.add(audit.maximum_hit_body);
            pair_ledger.add(audit.response_fnv);
            pair_ledger.add(audit.failure_fnv);
            pair_ledger.add(audit.complete_response_classes);
            pair_ledger.add(audit.complete_response_fnv);
            pair_ledger.add(audit.complete_active_response_fnv);
            pair_ledger.add(audit.full_response_masks);
            pair_ledger.add(audit.least_full_response);
            pair_ledger.add(audit.exact_minimum_replacements);
            for (u32 mask : audit.exact_witness_masks) pair_ledger.add(mask);
            std::cout << "PAIR " << pair.q << ',' << pair.r
                      << " ACTIVE_CARRIER " << audit.active_carrier
                      << " ACTIVE_CARRIER_FNV " << std::hex
                      << audit.active_carrier_fnv << std::dec
                      << " ACTIVE_JOINT " << audit.active_joint
                      << " ACTIVE_NONJOINT " << audit.active_nonjoint
                      << " INACTIVE_JOINT " << audit.inactive_joint
                      << " INACTIVE_JOINT_FNV " << std::hex
                      << audit.inactive_joint_fnv << std::dec << '\n'
                      << "EXPOSED " << audit.exposed << " EXPOSED_FNV "
                      << std::hex << audit.exposed_fnv << std::dec
                      << " JOINT_CHECKS " << audit.joint_checks
                      << " NONJOINT_CHECKS " << audit.nonjoint_checks
                      << " FAILURES " << audit.failures
                      << " HIT_INCIDENCES " << audit.hit_incidences
                      << " HIT_RANGE " << audit.minimum_hits << ".."
                      << audit.maximum_hits << " MIN_BODY " << std::hex
                      << std::setw(8) << std::setfill('0')
                      << audit.minimum_hit_body << " MAX_BODY "
                      << std::setw(8) << audit.maximum_hit_body
                      << " RESPONSE_FNV " << audit.response_fnv << std::dec
                      << std::setfill(' ') << '\n';
            if (!audit.failure_bodies.empty()) {
                std::cout << "FAILURE_BODIES " << audit.failure_bodies.size()
                          << " FNV " << std::hex << audit.failure_fnv
                          << " MASKS";
                if (audit.failure_bodies.size() < 64) {
                    for (u32 body : audit.failure_bodies)
                        std::cout << ' ' << std::setw(8) << std::setfill('0')
                                  << body;
                } else {
                    std::cout << " OMITTED_GT63 FIRST " << std::setw(8)
                              << std::setfill('0')
                              << audit.failure_bodies.front() << " LAST "
                              << std::setw(8) << audit.failure_bodies.back();
                }
                std::cout << std::dec << std::setfill(' ') << '\n'
                          << "COMPLETE_ACTIVE_RESPONSE_CLASSES ";
                if (audit.failure_bodies.size() < 64) {
                    std::cout << audit.complete_response_classes
                              << " CLASS_FNV " << std::hex
                              << audit.complete_response_fnv
                              << " ACTIVE_RESPONSE_FNV "
                              << audit.complete_active_response_fnv << std::dec
                              << " FULL_RESPONSE_MASKS "
                              << audit.full_response_masks << " LEAST_FULL "
                              << std::hex << std::setw(8) << std::setfill('0')
                              << audit.least_full_response << std::dec
                              << std::setfill(' ')
                              << " EXACT_MINIMUM "
                              << audit.exact_minimum_replacements
                              << " WITNESS";
                    for (u32 mask : audit.exact_witness_masks)
                        std::cout << ' ' << std::hex << std::setw(8)
                                  << std::setfill('0') << mask << std::dec
                                  << std::setfill(' ');
                } else {
                    std::cout << "SKIPPED_FAILURE_COUNT_GT63";
                }
                std::cout << '\n';
            }
        }
        require(!failure_output.is_open() || failure_output.good(),
                "failed writing failure ledger");
        std::cout << "PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << " TOTAL_FAILURES " << total_failures << '\n'
                  << "VERDICT "
                  << (total_failures == 0 ? "PASS ALL_PAIRS_CLOSED"
                                          : "FAIL RESISTANT_BODY_FOUND")
                  << '\n';
        return total_failures == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "FINAL8951_EXPOSURE_ERROR " << error.what() << '\n';
        return 1;
    }
}
