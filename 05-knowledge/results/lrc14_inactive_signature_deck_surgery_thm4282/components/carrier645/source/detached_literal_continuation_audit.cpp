// Detached literal/body-cover audit of the final8951 continuation.
//
// Activity is reconstructed from the actual q/r/fixed-pool joint wall
// arrangement.  This translation unit intentionally does not include or call
// response_pattern_atlas.cpp, PrimitivePair, or build_active_universe().

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define main detached_wall_dependency_unused_main
#include "04-computation/lrc14_endpoint_cascade_direct_wall_body_audit.cpp"
#undef main
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <atomic>
#include <bit>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <set>
#include <unordered_set>

namespace {

constexpr u64 kRank8Count = UINT64_C(5852925);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr std::size_t kCarrierCount = 8951;
constexpr std::size_t kJointCount = 421;
constexpr u64 kCarrierFnv = UINT64_C(0x188f82ab9dd1695a);
constexpr u64 kJointFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kResidualCount = 23637;
constexpr u64 kResidualFnv = UINT64_C(0xe8b363d2b3d9ba6a);
constexpr std::size_t kTargetPairCount = 90;
constexpr u64 kTargetPairFnv = UINT64_C(0x942995bee7469430);
constexpr std::size_t kAdditionCount = 45;
constexpr u64 kAdditionFnv = UINT64_C(0xec083b65cc8c34e3);
constexpr u32 kFirstAddition = UINT32_C(0x084a0a81);
constexpr u32 kSecondAddition = UINT32_C(0x0016580c);
constexpr std::array<std::size_t, 6> kStageCuts{0, 1, 2, 11, 42, 45};

std::array<std::array<u64, 9>, 31> local_choose{};

void init_local_choose() {
    for (int n = 0; n <= 30; ++n) {
        local_choose[n][0] = 1;
        for (int k = 1; k <= 8; ++k)
            local_choose[n][k] =
                n == 0 ? 0 : local_choose[n - 1][k] +
                                   local_choose[n - 1][k - 1];
    }
    require(local_choose[30][8] == kRank8Count,
            "detached rank-eight universe changed");
}

u64 local_rank8(u32 mask) {
    require(std::popcount(mask) == 8, "colex rank off rank eight");
    u64 rank = 0;
    int ordinal = 1;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += local_choose[bit][ordinal++];
    }
    require(ordinal == 9 && rank < kRank8Count, "colex rank escaped");
    return rank;
}

u32 local_unrank8(u64 rank) {
    require(rank < kRank8Count, "colex unrank escaped");
    const u64 original = rank;
    u32 mask = 0;
    int ceiling = 29;
    for (int ordinal = 8; ordinal >= 1; --ordinal) {
        while (ceiling >= 0 && local_choose[ceiling][ordinal] > rank)
            --ceiling;
        require(ceiling >= ordinal - 1, "colex unrank failed");
        mask |= u32{1} << ceiling;
        rank -= local_choose[ceiling][ordinal];
        --ceiling;
    }
    require(rank == 0 && local_rank8(mask) == original,
            "colex round trip failed");
    return mask;
}

std::vector<u32> read_masks_detached(const std::filesystem::path& path,
                                     std::size_t expected_count,
                                     u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    Fnv ledger;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "bad mask token");
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

std::vector<u32> read_additions(const std::filesystem::path& path) {
    std::vector<u32> masks =
        read_masks_detached(path, kAdditionCount, kAdditionFnv);
    require(masks[0] == kFirstAddition && masks[1] == kSecondAddition,
            "seed-addition identity/order changed");
    return masks;
}

struct PairKey {
    int q = 0;
    int r = 0;
    auto operator<=>(const PairKey&) const = default;
};

std::vector<PairKey> read_target_pairs(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open combined residual");
    std::vector<PairKey> all;
    Fnv all_ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "empty residual row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed residual row");
        PairKey pair{std::stoi(line.substr(0, comma)),
                     std::stoi(line.substr(comma + 1))};
        require(pair.q > 0 && pair.q < pair.r &&
                    (all.empty() || all.back() < pair),
                "residual order/coordinates changed");
        all.push_back(pair);
        all_ledger.add(pair.q);
        all_ledger.add(pair.r);
    }
    require(all.size() == kResidualCount &&
                all_ledger.state == kResidualFnv,
            "combined residual identity changed");
    std::vector<PairKey> target;
    Fnv target_ledger;
    for (const PairKey pair : all) {
        if (pair.r < 645 || pair.r > 663) continue;
        target.push_back(pair);
        target_ledger.add(pair.q);
        target_ledger.add(pair.r);
    }
    require(target.size() == kTargetPairCount &&
                target_ledger.state == kTargetPairFnv,
            "endpoint 663..645 pair slice changed");
    return target;
}

struct LiteralAtoms {
    i64 grid = 0;
    i64 pair_ticks = 0;
    u64 cells = 0;
    std::vector<std::pair<u32, i64>> masses;
};

LiteralAtoms build_literal_atoms(const PairKey pair) {
    const Geometry geometry = build_joint_geometry(pair.q, pair.r);
    std::map<u32, i64> grouped;
    i64 safe_sum = 0;
    for (const Cell& cell : geometry.cells) {
        if (!cell.pair_safe) continue;
        safe_sum += cell.width;
        if (std::popcount(cell.failed_pool) <= 8)
            grouped[cell.failed_pool] += cell.width;
    }
    require(safe_sum == geometry.pair_ticks,
            "literal safe cells do not telescope");
    LiteralAtoms answer;
    answer.grid = geometry.grid;
    answer.pair_ticks = geometry.pair_ticks;
    answer.cells = geometry.cells.size();
    answer.masses.assign(grouped.begin(), grouped.end());
    return answer;
}

i64 literal_mass(const LiteralAtoms& atoms, u32 mask) {
    i64 mass = 0;
    for (const auto [failure, value] : atoms.masses)
        if ((failure & ~mask) == 0) mass += value;
    return mass;
}

i128 literal_margin(const LiteralAtoms& atoms, u32 mask) {
    return static_cast<i128>(63) * literal_mass(atoms, mask) -
           static_cast<i128>(4) * atoms.grid;
}

struct StageAudit {
    u64 failures = 0;
    u64 failure_fnv = 0;
    u64 hit_incidences = 0;
    u64 minimum_hits = 0;
    u64 maximum_hits = 0;
    u32 minimum_body = 0;
    u32 maximum_body = 0;
    u64 response_fnv = 0;
    std::vector<u32> failure_bodies;
};

struct ExposureDatum {
    u32 body = 0;
    std::array<u64, 6> hits{};
    std::array<u32, 6> least{};
};

struct PairAuditDetached {
    PairKey pair;
    i64 grid = 0;
    i64 pair_ticks = 0;
    u64 cells = 0;
    u64 atom_classes = 0;
    u64 active_base = 0;
    u64 active_base_fnv = 0;
    u64 active_joint = 0;
    u64 active_nonjoint = 0;
    u64 inactive_joint = 0;
    u64 inactive_joint_fnv = 0;
    u64 equalities = 0;
    std::vector<i128> addition_margin;
    std::vector<bool> addition_active;
    u64 joint_checks = 0;
    u64 nonjoint_checks = 0;
    u64 exposed = 0;
    u64 exposed_fnv = 0;
    std::array<StageAudit, 6> stages;
};

u32 next_same_popcount_detached(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

StageAudit summarize_stage(const std::vector<ExposureDatum>& exposed,
                           std::size_t stage) {
    StageAudit answer;
    Fnv response_ledger;
    Fnv failure_ledger;
    bool positive_seen = false;
    for (const ExposureDatum& row : exposed) {
        const u64 hits = row.hits[stage];
        response_ledger.add(row.body);
        response_ledger.add(hits);
        response_ledger.add(row.least[stage]);
        answer.hit_incidences += hits;
        if (hits == 0) {
            answer.failure_bodies.push_back(row.body);
            failure_ledger.add(row.body);
            continue;
        }
        if (!positive_seen || hits < answer.minimum_hits) {
            positive_seen = true;
            answer.minimum_hits = hits;
            answer.minimum_body = row.body;
        }
        if (hits > answer.maximum_hits) {
            answer.maximum_hits = hits;
            answer.maximum_body = row.body;
        }
    }
    answer.failures = answer.failure_bodies.size();
    answer.failure_fnv = failure_ledger.state;
    answer.response_fnv = response_ledger.state;
    return answer;
}

PairAuditDetached audit_pair_detached(
    const PairKey pair, const std::vector<u32>& carrier,
    const std::vector<u32>& joint, const std::unordered_set<u32>& joint_set,
    const std::vector<u32>& additions) {
    const LiteralAtoms atoms = build_literal_atoms(pair);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    PairAuditDetached answer;
    answer.pair = pair;
    answer.grid = atoms.grid;
    answer.pair_ticks = atoms.pair_ticks;
    answer.cells = atoms.cells;
    answer.atom_classes = atoms.masses.size();
    Fnv active_ledger;
    Fnv inactive_joint_ledger;
    for (u32 mask : carrier) {
        const i128 margin = literal_margin(atoms, mask);
        answer.equalities += margin == 0;
        if (margin < 0) continue;
        ++answer.active_base;
        active_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    for (std::size_t index = 0; index < joint.size(); ++index) {
        const i128 margin = literal_margin(atoms, joint[index]);
        if (margin >= 0) {
            active_joint.push_back(joint[index]);
        } else {
            ++answer.inactive_joint;
            inactive_joint_ledger.add(index);
            inactive_joint_ledger.add(joint[index]);
        }
    }
    answer.active_base_fnv = active_ledger.state;
    answer.active_joint = active_joint.size();
    answer.active_nonjoint = active_nonjoint.size();
    answer.inactive_joint_fnv = inactive_joint_ledger.state;
    require(answer.active_base == answer.active_joint +
                                      answer.active_nonjoint,
            "detached carrier partition changed");
    answer.addition_margin.resize(additions.size());
    answer.addition_active.resize(additions.size());
    for (std::size_t index = 0; index < additions.size(); ++index) {
        answer.addition_margin[index] = literal_margin(atoms, additions[index]);
        answer.addition_active[index] = answer.addition_margin[index] >= 0;
    }

    std::vector<ExposureDatum> exposed;
    u32 body = (u32{1} << 9) - 1;
    for (u64 ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        bool joint_hit = false;
        for (u32 mask : active_joint) {
            ++answer.joint_checks;
            if ((body & mask) == 0) {
                joint_hit = true;
                break;
            }
        }
        if (!joint_hit) {
            ExposureDatum row;
            row.body = body;
            for (u32 mask : active_nonjoint) {
                ++answer.nonjoint_checks;
                if ((body & mask) != 0) continue;
                ++row.hits[0];
                if (row.hits[0] == 1 || mask < row.least[0])
                    row.least[0] = mask;
            }
            std::size_t previous_cut = 0;
            for (std::size_t stage = 1; stage < kStageCuts.size(); ++stage) {
                row.hits[stage] = row.hits[stage - 1];
                row.least[stage] = row.least[stage - 1];
                for (std::size_t index = previous_cut;
                     index < kStageCuts[stage]; ++index) {
                    if (!answer.addition_active[index] ||
                        (body & additions[index]) != 0)
                        continue;
                    ++row.hits[stage];
                    if (row.least[stage] == 0 ||
                        additions[index] < row.least[stage])
                        row.least[stage] = additions[index];
                }
                previous_cut = kStageCuts[stage];
            }
            exposed.push_back(row);
        }
        if (ordinal + 1 < kBodyCount)
            body = next_same_popcount_detached(body);
    }
    require(body == UINT32_C(0x3fe00000),
            "detached body endpoint changed");
    answer.exposed = exposed.size();
    Fnv exposed_ledger;
    for (const ExposureDatum& row : exposed) exposed_ledger.add(row.body);
    answer.exposed_fnv = exposed_ledger.state;
    for (std::size_t stage = 0; stage < answer.stages.size(); ++stage)
        answer.stages[stage] = summarize_stage(exposed, stage);
    return answer;
}

template <class Callback>
void choose_extra(const std::vector<unsigned char>& positions,
                  int start, int need, u32 mask, Callback& callback) {
    if (need == 0) {
        callback(mask);
        return;
    }
    for (int index = start;
         index <= static_cast<int>(positions.size()) - need; ++index)
        choose_extra(positions, index + 1, need - 1,
                     mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
u64 enumerate_supersets(u32 atom, Callback callback) {
    const int arity = std::popcount(atom);
    require(arity <= 8, "literal atom exceeds rank eight");
    std::vector<unsigned char> complement;
    for (int bit = 0; bit < 30; ++bit)
        if ((atom & (u32{1} << bit)) == 0)
            complement.push_back(static_cast<unsigned char>(bit));
    u64 emitted = 0;
    auto checked = [&](u32 extra) {
        const u32 mask = atom | extra;
        callback(mask, local_rank8(mask));
        ++emitted;
    };
    choose_extra(complement, 0, 8 - arity, 0, checked);
    return emitted;
}

struct CompleteResponseAudit {
    u64 active = 0;
    u64 active_fnv = 0;
    u64 equalities = 0;
    u64 zeta_operations = 0;
    u64 classes = 0;
    u64 class_fnv = 0;
    u64 active_response_fnv = 0;
    u64 response_union = 0;
    u64 full_response_masks = 0;
    u32 least_full = 0;
    u64 maximal_classes = 0;
    u64 maximal_class_fnv = 0;
    u64 minimum_cover = 0;
    std::vector<u32> minimum_witness;
};

CompleteResponseAudit complete_response_detached(
    const PairKey pair, const std::vector<u32>& failures) {
    require(!failures.empty() && failures.size() < 64,
            "boundary failure quotient outside one word");
    const LiteralAtoms atoms = build_literal_atoms(pair);
    std::vector<i64> masses(kRank8Count, 0);
    u64 operations = 0;
    for (const auto [atom, value] : atoms.masses)
        operations += enumerate_supersets(
            atom, [&](u32, u64 rank) { masses[rank] += value; });
    struct ResponseClass {
        u64 count = 0;
        u32 least = 0;
    };
    std::map<u64, ResponseClass> classes;
    CompleteResponseAudit answer;
    answer.zeta_operations = operations;
    Fnv active_ledger;
    Fnv active_response_ledger;
    const u64 full = (UINT64_C(1) << failures.size()) - 1;
    for (u64 rank = 0; rank < kRank8Count; ++rank) {
        const i128 margin = static_cast<i128>(63) * masses[rank] -
                            static_cast<i128>(4) * atoms.grid;
        answer.equalities += margin == 0;
        if (margin < 0) continue;
        const u32 mask = local_unrank8(rank);
        ++answer.active;
        active_ledger.add(mask);
        u64 response = 0;
        for (std::size_t bit = 0; bit < failures.size(); ++bit)
            if ((mask & failures[bit]) == 0)
                response |= UINT64_C(1) << bit;
        answer.response_union |= response;
        active_response_ledger.add(mask);
        active_response_ledger.add(response);
        ResponseClass& response_class = classes[response];
        ++response_class.count;
        if (response_class.count == 1 || mask < response_class.least)
            response_class.least = mask;
    }
    answer.active_fnv = active_ledger.state;
    answer.active_response_fnv = active_response_ledger.state;
    answer.classes = classes.size();
    Fnv class_ledger;
    for (const auto& [response, response_class] : classes) {
        class_ledger.add(response);
        class_ledger.add(response_class.count);
        class_ledger.add(response_class.least);
    }
    answer.class_fnv = class_ledger.state;
    const auto full_class = classes.find(full);
    if (full_class != classes.end()) {
        answer.full_response_masks = full_class->second.count;
        answer.least_full = full_class->second.least;
    }
    require(answer.response_union == full,
            "complete boundary response union misses an obligation");

    std::vector<std::pair<u64, u32>> maximal;
    for (const auto& [response, response_class] : classes) {
        if (response == 0) continue;
        bool dominated = false;
        for (const auto& [other, ignored] : classes) {
            (void)ignored;
            if (response != other && (response | other) == other) {
                dominated = true;
                break;
            }
        }
        if (!dominated) maximal.emplace_back(response, response_class.least);
    }
    answer.maximal_classes = maximal.size();
    Fnv maximal_ledger;
    for (const auto [response, least] : maximal) {
        maximal_ledger.add(response);
        maximal_ledger.add(least);
    }
    answer.maximal_class_fnv = maximal_ledger.state;

    std::sort(maximal.begin(), maximal.end(), [](const auto& left,
                                                  const auto& right) {
        const int left_size = std::popcount(left.first);
        const int right_size = std::popcount(right.first);
        if (left_size != right_size) return left_size > right_size;
        if (left.first != right.first) return left.first < right.first;
        return left.second < right.second;
    });
    std::vector<std::vector<std::size_t>> containing(failures.size());
    for (std::size_t index = 0; index < maximal.size(); ++index)
        for (std::size_t bit = 0; bit < failures.size(); ++bit)
            if ((maximal[index].first & (UINT64_C(1) << bit)) != 0)
                containing[bit].push_back(index);

    std::vector<std::size_t> path;
    std::vector<std::size_t> winning_path;
    for (std::size_t depth = 1; depth <= failures.size(); ++depth) {
        std::unordered_set<u64> dead;
        std::function<bool(u64, std::size_t)> search =
            [&](u64 state, std::size_t remaining) -> bool {
            if (state == full) {
                winning_path = path;
                return true;
            }
            if (remaining == 0) return false;
            const u64 uncovered = full & ~state;
            int maximum_gain = 0;
            for (const auto [response, ignored] : maximal) {
                (void)ignored;
                maximum_gain = std::max(
                    maximum_gain, std::popcount(response & uncovered));
            }
            if (maximum_gain == 0 ||
                (std::popcount(uncovered) + maximum_gain - 1) /
                        maximum_gain >
                    static_cast<int>(remaining))
                return false;
            const u64 key = state |
                            (static_cast<u64>(remaining) << failures.size());
            if (dead.contains(key)) return false;
            std::size_t chosen_bit = failures.size();
            std::size_t chosen_count = std::numeric_limits<std::size_t>::max();
            for (std::size_t bit = 0; bit < failures.size(); ++bit) {
                if ((uncovered & (UINT64_C(1) << bit)) == 0) continue;
                std::size_t count = 0;
                for (std::size_t index : containing[bit])
                    count += (maximal[index].first & uncovered) != 0;
                if (count < chosen_count) {
                    chosen_bit = bit;
                    chosen_count = count;
                }
            }
            require(chosen_bit < failures.size() && chosen_count > 0,
                    "OR-cover branching found an unreachable obligation");
            for (std::size_t index : containing[chosen_bit]) {
                const u64 next = state | maximal[index].first;
                if (next == state) continue;
                path.push_back(index);
                if (search(next, remaining - 1)) return true;
                path.pop_back();
            }
            dead.insert(key);
            return false;
        };
        path.clear();
        if (!search(0, depth)) continue;
        answer.minimum_cover = depth;
        for (std::size_t index : winning_path)
            answer.minimum_witness.push_back(maximal[index].second);
        break;
    }
    require(answer.minimum_cover > 0,
            "exact branch-and-bound found no OR cover");
    return answer;
}

void require_witness_cover(const PairKey pair,
                           const std::vector<u32>& failures,
                           const std::vector<u32>& witnesses) {
    const LiteralAtoms atoms = build_literal_atoms(pair);
    for (u32 witness : witnesses)
        require(std::popcount(witness) == 8 &&
                    literal_margin(atoms, witness) >= 0,
                "claimed witness is not literally active");
    for (u32 body : failures) {
        bool hit = false;
        for (u32 witness : witnesses)
            hit = hit || ((body & witness) == 0);
        require(hit, "claimed witness family misses a failure body");
    }
}

void add_i128_ledger(Fnv& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

void write_pair_rows(const std::filesystem::path& path,
                     const std::vector<PairAuditDetached>& rows) {
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create pair-row TSV");
    output << "q\tr\tgrid\tcells\tatoms\tactive_base\tactive_base_fnv"
              "\tactive_joint\tactive_nonjoint\tinactive_joint"
              "\tinactive_joint_fnv\tactive_additions\taddition_activity_fnv"
              "\texposed\texposed_fnv";
    for (std::size_t stage = 0; stage < kStageCuts.size(); ++stage)
        output << "\tstage" << stage << "_cut" << kStageCuts[stage]
               << "_fail\tstage" << stage << "_failure_fnv";
    output << '\n';
    for (const PairAuditDetached& row : rows) {
        output << row.pair.q << '\t' << row.pair.r << '\t' << row.grid
               << '\t' << row.cells << '\t' << row.atom_classes << '\t'
               << row.active_base << '\t' << std::hex
               << row.active_base_fnv << std::dec << '\t'
               << row.active_joint << '\t' << row.active_nonjoint << '\t'
               << row.inactive_joint << '\t' << std::hex
               << row.inactive_joint_fnv << std::dec;
        Fnv addition_activity_ledger;
        u64 active_additions = 0;
        for (std::size_t index = 0; index < row.addition_active.size(); ++index) {
            addition_activity_ledger.add(index);
            add_i128_ledger(addition_activity_ledger,
                            row.addition_margin[index]);
            addition_activity_ledger.add(row.addition_active[index]);
            active_additions += row.addition_active[index];
        }
        output << '\t' << active_additions << '\t' << std::hex
               << addition_activity_ledger.state << std::dec << '\t'
               << row.exposed << '\t' << std::hex << row.exposed_fnv
               << std::dec;
        for (const StageAudit& stage : row.stages)
            output << '\t' << stage.failures << '\t' << std::hex
                   << stage.failure_fnv << std::dec;
        output << '\n';
    }
    require(output.good(), "failed writing pair-row TSV");
}

std::string masks_text(const std::vector<u32>& masks) {
    std::ostringstream output;
    for (u32 mask : masks)
        output << ' ' << std::hex << std::setw(8) << std::setfill('0')
               << mask;
    return output.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 6,
                "usage: detached-audit CARRIER JOINT RESIDUAL ADDITIONS45 ROWS_TSV");
        init_local_choose();
        const std::vector<u32> carrier =
            read_masks_detached(argv[1], kCarrierCount, kCarrierFnv);
        const std::vector<u32> joint =
            read_masks_detached(argv[2], kJointCount, kJointFnv);
        const std::vector<u32> additions = read_additions(argv[4]);
        const std::vector<PairKey> pairs = read_target_pairs(argv[3]);
        const std::set<u32> carrier_set(carrier.begin(), carrier.end());
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            require(carrier_set.contains(mask),
                    "joint mask absent from reconstructed carrier");
        for (u32 mask : additions)
            require(!carrier_set.contains(mask) && !joint_set.contains(mask),
                    "continuation mask is not novel nonjoint");
        Fnv final_carrier_ledger;
        for (u32 mask : carrier) final_carrier_ledger.add(mask);
        for (u32 mask : additions) final_carrier_ledger.add(mask);
        require(final_carrier_ledger.state == UINT64_C(0xfd899660f14b311c),
                "final8996 ordered carrier identity changed");

        std::vector<PairAuditDetached> audits(pairs.size());
        std::atomic<std::size_t> next{0};
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned threads =
            std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
        std::vector<std::thread> workers;
        std::mutex progress_mutex;
        for (unsigned lane = 0; lane < threads; ++lane) {
            workers.emplace_back([&]() {
                for (;;) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= pairs.size()) break;
                    audits[index] = audit_pair_detached(
                        pairs[index], carrier, joint, joint_set, additions);
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    std::cerr << "PROGRESS " << index + 1 << '/' << pairs.size()
                              << " PAIR " << pairs[index].q << ','
                              << pairs[index].r << '\n';
                }
            });
        }
        for (std::thread& worker : workers) worker.join();

        u64 final_failures = 0;
        Fnv pair_ledger;
        for (const PairAuditDetached& row : audits) {
            final_failures += row.stages.back().failures;
            pair_ledger.add(row.pair.q);
            pair_ledger.add(row.pair.r);
            pair_ledger.add(row.active_base);
            pair_ledger.add(row.active_base_fnv);
            pair_ledger.add(row.active_joint);
            pair_ledger.add(row.active_nonjoint);
            pair_ledger.add(row.inactive_joint);
            pair_ledger.add(row.inactive_joint_fnv);
            for (std::size_t index = 0; index < additions.size(); ++index) {
                pair_ledger.add(index);
                add_i128_ledger(pair_ledger, row.addition_margin[index]);
                pair_ledger.add(row.addition_active[index]);
            }
            pair_ledger.add(row.exposed);
            pair_ledger.add(row.exposed_fnv);
            for (const StageAudit& stage : row.stages) {
                pair_ledger.add(stage.failures);
                pair_ledger.add(stage.failure_fnv);
                pair_ledger.add(stage.hit_incidences);
                pair_ledger.add(stage.minimum_hits);
                pair_ledger.add(stage.minimum_body);
                pair_ledger.add(stage.maximum_hits);
                pair_ledger.add(stage.maximum_body);
                pair_ledger.add(stage.response_fnv);
            }
        }
        write_pair_rows(argv[5], audits);

        const auto find_pair = [&](int q, int r) -> const PairAuditDetached& {
            const auto found = std::find_if(
                audits.begin(), audits.end(), [&](const auto& row) {
                    return row.pair.q == q && row.pair.r == r;
                });
            require(found != audits.end(), "boundary pair missing");
            return *found;
        };
        const PairAuditDetached& boundary660 = find_pair(256, 660);
        const PairAuditDetached& boundary657 = find_pair(256, 657);
        const PairAuditDetached& boundary650 = find_pair(256, 650);
        const PairAuditDetached& boundary648 = find_pair(256, 648);
        const PairAuditDetached& small294 = find_pair(294, 650);
        const PairAuditDetached& small366 = find_pair(366, 650);
        const PairAuditDetached& small416 = find_pair(416, 650);
        const PairAuditDetached& small512 = find_pair(512, 650);
        const std::vector<u32> expected660{
            UINT32_C(0x05346408), UINT32_C(0x07146048),
            UINT32_C(0x07147008), UINT32_C(0x0714e008),
            UINT32_C(0x07346008)};
        require(boundary660.stages[0].failure_bodies == expected660 &&
                    boundary660.stages[1].failures == 0,
                "(256,660) staged failure fibre changed");
        require(boundary657.stages[1].failures == 1 &&
                    boundary657.stages[2].failures == 0,
                "(256,657) staged failure fibre changed");
        const std::map<int, u64> expected650{{256, 368}, {294, 2},
                                             {366, 18}, {416, 14},
                                             {512, 1}};
        std::map<int, u64> observed650;
        for (const PairAuditDetached& row : audits)
            if (row.pair.r == 650 && row.stages[2].failures != 0)
                observed650.emplace(row.pair.q, row.stages[2].failures);
        require(observed650 == expected650,
                "endpoint650 two-seed resistant classification changed");
        require(boundary650.stages[3].failures == 358 &&
                    boundary650.stages[4].failures == 0 &&
                    small294.stages[3].failures == 0 &&
                    small366.stages[3].failures == 0 &&
                    small416.stages[3].failures == 0 &&
                    small512.stages[3].failures == 0,
                "endpoint650 nine-mask/greedy continuation changed");
        std::vector<PairKey> prefinal_resistant;
        for (const PairAuditDetached& row : audits)
            if (row.stages[4].failures != 0)
                prefinal_resistant.push_back(row.pair);
        require(prefinal_resistant == std::vector<PairKey>{{256, 648}} &&
                    boundary648.stages[4].failures == 6 &&
                    boundary648.stages[5].failures == 0,
                "pre-final endpoint648 boundary classification changed");
        require(final_failures == 0,
                "45-mask continuation misses a 663..645 body");

        const CompleteResponseAudit quotient660 = complete_response_detached(
            {256,660}, boundary660.stages[0].failure_bodies);
        const CompleteResponseAudit quotient657 = complete_response_detached(
            {256,657}, boundary657.stages[1].failure_bodies);
        const CompleteResponseAudit quotient294 = complete_response_detached(
            {294,650}, small294.stages[2].failure_bodies);
        const CompleteResponseAudit quotient366 = complete_response_detached(
            {366,650}, small366.stages[2].failure_bodies);
        const CompleteResponseAudit quotient416 = complete_response_detached(
            {416,650}, small416.stages[2].failure_bodies);
        const CompleteResponseAudit quotient512 = complete_response_detached(
            {512,650}, small512.stages[2].failure_bodies);
        const CompleteResponseAudit quotient648 = complete_response_detached(
            {256,648}, boundary648.stages[4].failure_bodies);
        require(quotient660.full_response_masks > 0 &&
                    quotient660.least_full == kFirstAddition,
                "(256,660) least full literal responder changed");
        require(quotient657.full_response_masks > 0 &&
                    quotient657.least_full == kSecondAddition,
                "(256,657) least full literal responder changed");
        require(quotient294.minimum_cover == 1 &&
                    quotient366.minimum_cover == 4 &&
                    quotient416.minimum_cover == 3 &&
                    quotient512.minimum_cover == 1 &&
                    quotient648.minimum_cover == 3,
                "detached exact minimum response cover changed");
        const std::vector<u32> witness294(additions.begin() + 2,
                                          additions.begin() + 3);
        const std::vector<u32> witness366(additions.begin() + 3,
                                          additions.begin() + 7);
        const std::vector<u32> witness416(additions.begin() + 7,
                                          additions.begin() + 10);
        const std::vector<u32> witness512(additions.begin() + 10,
                                          additions.begin() + 11);
        const std::vector<u32> witness648(additions.begin() + 42,
                                          additions.begin() + 45);
        require_witness_cover({294,650}, small294.stages[2].failure_bodies,
                              witness294);
        require_witness_cover({366,650}, small366.stages[2].failure_bodies,
                              witness366);
        require_witness_cover({416,650}, small416.stages[2].failure_bodies,
                              witness416);
        require_witness_cover({512,650}, small512.stages[2].failure_bodies,
                              witness512);
        require_witness_cover({256,648}, boundary648.stages[4].failure_bodies,
                              witness648);

        std::cout << "THM4282_FINAL8951_CONTINUATION_DETACHED_LITERAL_AUDIT_V2\n"
                  << "BASE_CARRIER " << carrier.size() << " FNV " << std::hex
                  << kCarrierFnv << " JOINT " << std::dec << joint.size()
                  << " FNV " << std::hex << kJointFnv << std::dec
                  << " ADDITIONS " << std::hex << std::setw(8)
                  << std::setfill('0') << additions[0] << ',' << std::setw(8)
                  << additions[1] << std::dec << std::setfill(' ')
                  << " TOTAL " << additions.size() << " ADDITION_FNV "
                  << std::hex << kAdditionFnv << " FINAL_CARRIER 8996 FNV "
                  << final_carrier_ledger.state << std::dec << '\n'
                  << "COMBINED_RESIDUAL " << kResidualCount << " FNV "
                  << std::hex << kResidualFnv << std::dec
                  << " TARGET_ENDPOINTS 663..645 PAIRS " << pairs.size()
                  << " FNV " << std::hex << kTargetPairFnv << std::dec
                  << " THREADS " << threads << '\n'
                  << "PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << " FINAL_FAILURES " << final_failures << '\n';
        for (const PairAuditDetached* boundary :
             {&boundary660, &boundary657, &boundary650, &boundary648}) {
            std::cout << "BOUNDARY " << boundary->pair.q << ','
                      << boundary->pair.r << " GRID " << boundary->grid
                      << " CELLS " << boundary->cells << " ATOMS "
                      << boundary->atom_classes << " ACTIVE_BASE "
                      << boundary->active_base << " ACTIVE_BASE_FNV "
                      << std::hex << boundary->active_base_fnv << std::dec
                      << " ACTIVE_JOINT " << boundary->active_joint
                      << " EXPOSED " << boundary->exposed;
            for (std::size_t stage = 0; stage < kStageCuts.size(); ++stage)
                std::cout << " CUT" << kStageCuts[stage] << "_FAIL "
                          << boundary->stages[stage].failures;
            std::cout << '\n';
        }
        std::cout << "BOUNDARY_FAILURE_BODIES PAIR 256,660 CUT 0 COUNT "
                  << boundary660.stages[0].failure_bodies.size() << " FNV "
                  << std::hex << boundary660.stages[0].failure_fnv << " MASKS"
                  << masks_text(boundary660.stages[0].failure_bodies)
                  << std::dec << '\n'
                  << "BOUNDARY_FAILURE_BODIES PAIR 256,657 CUT 1 COUNT "
                  << boundary657.stages[1].failure_bodies.size() << " FNV "
                  << std::hex << boundary657.stages[1].failure_fnv << " MASKS"
                  << masks_text(boundary657.stages[1].failure_bodies)
                  << std::dec << '\n';
        const auto print_quotient = [&](const PairKey pair,
                                        const CompleteResponseAudit& q) {
            std::cout << "COMPLETE_LITERAL_RESPONSE PAIR " << pair.q << ','
                      << pair.r << " RANK8 " << kRank8Count << " ACTIVE "
                      << q.active << " ACTIVE_FNV " << std::hex << q.active_fnv
                      << std::dec << " EQUALITIES " << q.equalities
                      << " ZETA_OPERATIONS " << q.zeta_operations
                      << " CLASSES " << q.classes << " CLASS_FNV " << std::hex
                      << q.class_fnv << " ACTIVE_RESPONSE_FNV "
                      << q.active_response_fnv << " UNION " << q.response_union
                      << std::dec << " FULL_RESPONSE_MASKS "
                      << q.full_response_masks << " LEAST_FULL " << std::hex
                      << std::setw(8) << std::setfill('0') << q.least_full
                      << std::dec << std::setfill(' ') << " MAXIMAL_CLASSES "
                      << q.maximal_classes << " MAXIMAL_FNV " << std::hex
                      << q.maximal_class_fnv << std::dec << " MINIMUM_COVER "
                      << q.minimum_cover << " BFS_WITNESS"
                      << masks_text(q.minimum_witness) << std::dec << '\n';
        };
        print_quotient({256,660}, quotient660);
        print_quotient({256,657}, quotient657);
        print_quotient({294,650}, quotient294);
        print_quotient({366,650}, quotient366);
        print_quotient({416,650}, quotient416);
        print_quotient({512,650}, quotient512);
        print_quotient({256,648}, quotient648);
        std::cout << "CLAIMED_MINIMUM_WITNESSES PAIR 294,650 MASKS"
                  << masks_text(witness294) << '\n'
                  << "CLAIMED_MINIMUM_WITNESSES PAIR 366,650 MASKS"
                  << masks_text(witness366) << '\n'
                  << "CLAIMED_MINIMUM_WITNESSES PAIR 416,650 MASKS"
                  << masks_text(witness416) << '\n'
                  << "CLAIMED_MINIMUM_WITNESSES PAIR 512,650 MASKS"
                  << masks_text(witness512) << '\n'
                  << "CLAIMED_MINIMUM_WITNESSES PAIR 256,648 MASKS"
                  << masks_text(witness648) << std::dec << '\n';
        std::cout << "VERDICT PASS DETACHED_LITERAL_45_MASK_CONTINUATION_663_TO_645\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DETACHED_CONTINUATION_ERROR " << error.what() << '\n';
        return 1;
    }
}
