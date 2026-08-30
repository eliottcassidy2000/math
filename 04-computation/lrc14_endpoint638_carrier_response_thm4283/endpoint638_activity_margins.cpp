#include <filesystem>
// Exact activity-margin audit for THM-4283.
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
        require(argc == 4, "usage: margins JOINT BASE ADDITIONS");
        init_choose8_local();
        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        std::vector<u32> carrier =
            read_masks(argv[2], kCarrierCount, kCarrierFnv);
        std::set<u32> carrier_set(carrier.begin(), carrier.end());
        for (u32 mask : read_additions(argv[3])) {
            require(carrier_set.insert(mask).second, "addition overlap");
            carrier.push_back(mask);
        }
        constexpr u32 repair = UINT32_C(0x014c9084);
        require(carrier_set.insert(repair).second && std::popcount(repair) == 8,
                "repair invalid/duplicate");
        carrier.push_back(repair);
        const std::vector<Cell> cells = build_pool_cells();
        const std::array<Pair, 8> pairs{{
            {220,644}, {256,644}, {258,644}, {294,644},
            {366,644}, {416,644}, {512,644}, {256,638}}};
        FnvLocal global;
        for (const Pair pair : pairs) {
            const i64 g = std::gcd(pair.q, pair.r);
            const PrimitivePair primitive =
                build_primitive(pair.q / g, pair.r / g);
            const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
            std::vector<i128> masses(EXPECTED_REPAIRS, 0);
            u64 operations = 0;
            for (const auto& [mask, value] : atoms.mass)
                add_supersets_pair(mask, 8 - std::popcount(mask), 0, 0,
                                   value, masses, operations);
            require(operations == UINT64_C(152170690),
                    "zeta operation count changed");
            const i128 denominator =
                static_cast<i128>(primitive.grid) * g * COMMON;
            auto margin = [&](u32 mask) {
                return static_cast<i128>(63) *
                           masses[colex_rank8_local(mask)] -
                       static_cast<i128>(4) * denominator;
            };
            u64 active_all = 0;
            FnvLocal active_all_ledger;
            bool min_all_set = false, max_inactive_all_set = false;
            i128 min_all = 0, max_inactive_all = 0;
            u32 min_all_mask = 0, max_inactive_all_mask = 0;
            u32 mask = (u32{1} << 8) - 1;
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
                const i128 m = static_cast<i128>(63) * masses[rank] -
                               static_cast<i128>(4) * denominator;
                if (m >= 0) {
                    ++active_all;
                    active_all_ledger.add(mask);
                    if (!min_all_set || m < min_all ||
                        (m == min_all && mask < min_all_mask)) {
                        min_all_set = true; min_all = m; min_all_mask = mask;
                    }
                } else if (!max_inactive_all_set || m > max_inactive_all ||
                           (m == max_inactive_all &&
                            mask < max_inactive_all_mask)) {
                    max_inactive_all_set = true;
                    max_inactive_all = m;
                    max_inactive_all_mask = mask;
                }
                if (rank + 1 < EXPECTED_REPAIRS)
                    mask = next_same_popcount(mask);
            }
            u64 active_carrier = 0, inactive_carrier = 0;
            FnvLocal active_carrier_ledger;
            bool min_carrier_set = false, max_inactive_carrier_set = false;
            i128 min_carrier = 0, max_inactive_carrier = 0;
            u32 min_carrier_mask = 0, max_inactive_carrier_mask = 0;
            for (u32 carrier_mask : carrier) {
                const i128 m = margin(carrier_mask);
                if (m >= 0) {
                    ++active_carrier;
                    active_carrier_ledger.add(carrier_mask);
                    if (!min_carrier_set || m < min_carrier ||
                        (m == min_carrier &&
                         carrier_mask < min_carrier_mask)) {
                        min_carrier_set = true;
                        min_carrier = m;
                        min_carrier_mask = carrier_mask;
                    }
                } else {
                    ++inactive_carrier;
                    if (!max_inactive_carrier_set ||
                        m > max_inactive_carrier ||
                        (m == max_inactive_carrier &&
                         carrier_mask < max_inactive_carrier_mask)) {
                        max_inactive_carrier_set = true;
                        max_inactive_carrier = m;
                        max_inactive_carrier_mask = carrier_mask;
                    }
                }
            }
            std::vector<std::pair<std::size_t,u32>> inactive_joint;
            for (std::size_t index = 0; index < joint.size(); ++index)
                if (margin(joint[index]) < 0)
                    inactive_joint.push_back({index,joint[index]});
            global.add(pair.q); global.add(pair.r);
            global.add(active_all); global.add(active_all_ledger.state);
            global.add(active_carrier); global.add(active_carrier_ledger.state);
            global.add(inactive_carrier);
            global.add(static_cast<u64>(min_all));
            global.add(static_cast<u64>(max_inactive_all));
            global.add(static_cast<u64>(min_carrier));
            global.add(static_cast<u64>(max_inactive_carrier));
            global.add(static_cast<u64>(margin(repair)));
            for (const auto& [index,joint_mask] : inactive_joint) {
                global.add(index); global.add(joint_mask);
                global.add(static_cast<u64>(margin(joint_mask)));
            }
            std::cout << "PAIR " << pair.q << ',' << pair.r
                      << " DEN " << decimal(denominator)
                      << " ACTIVE_ALL " << active_all
                      << " ACTIVE_ALL_FNV " << std::hex
                      << active_all_ledger.state << std::dec
                      << " ALL_MIN_ACTIVE " << std::hex << std::setw(8)
                      << std::setfill('0') << min_all_mask << std::dec
                      << std::setfill(' ') << ':' << decimal(min_all)
                      << " ALL_CLOSEST_INACTIVE " << std::hex << std::setw(8)
                      << std::setfill('0') << max_inactive_all_mask << std::dec
                      << std::setfill(' ') << ':'
                      << decimal(max_inactive_all) << '\n'
                      << "CARRIER ACTIVE " << active_carrier
                      << " INACTIVE " << inactive_carrier
                      << " ACTIVE_FNV " << std::hex
                      << active_carrier_ledger.state << std::dec
                      << " MIN_ACTIVE " << std::hex << std::setw(8)
                      << std::setfill('0') << min_carrier_mask << std::dec
                      << std::setfill(' ') << ':' << decimal(min_carrier)
                      << " CLOSEST_INACTIVE " << std::hex << std::setw(8)
                      << std::setfill('0') << max_inactive_carrier_mask
                      << std::dec << std::setfill(' ') << ':'
                      << decimal(max_inactive_carrier)
                      << " REPAIR_014c9084 " << decimal(margin(repair))
                      << '\n'
                      << "INACTIVE_JOINT";
            for (const auto& [index,joint_mask] : inactive_joint)
                std::cout << ' ' << index << ':' << std::hex
                          << std::setw(8) << std::setfill('0') << joint_mask
                          << std::dec << std::setfill(' ') << ':'
                          << decimal(margin(joint_mask));
            std::cout << '\n';
            if (pair.q == 256 && pair.r == 638) {
                std::cout << "NEXT_WITNESS";
                for (u32 witness : std::array<u32,9>{{
                     UINT32_C(0x02203226), UINT32_C(0x081e1084),
                     UINT32_C(0x08a89440), UINT32_C(0x180a8281),
                     UINT32_C(0x18261042), UINT32_C(0x18a0d040),
                     UINT32_C(0x1a82a200), UINT32_C(0x202a9440),
                     UINT32_C(0x280a0a88)}})
                    std::cout << ' ' << std::hex << std::setw(8)
                              << std::setfill('0') << witness << std::dec
                              << std::setfill(' ') << ':'
                              << decimal(margin(witness));
                std::cout << '\n';
            }
        }
        std::cout << "LEDGER " << std::hex << global.state << std::dec << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << '\n';
        return 1;
    }
}
