// Exact generic deletion/repair surgery for a frozen rank-eight deck.
// This instance is locked to THM-4281 target (520,663) and the five inactive
// deck positions read independently from the frozen full-signature CSV.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace {

constexpr int kTargetQ = 520;
constexpr int kTargetR = 663;
constexpr std::size_t kSurgeryDeckSize = 421;
constexpr std::size_t kSurgeryBodyCount = 14307150;
constexpr u64 kSurgeryDeckFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr u64 kSurgeryFinalResidualFnv = UINT64_C(0x80ec0687d8c7dba7);
constexpr std::size_t kSurgerySignatureRows = 24223;
constexpr std::array<std::size_t, 5> kExpectedDeleted{
    57, 107, 222, 275, 345};
constexpr int kDualDenominator = 11;
constexpr std::array<int, 53> kDualNumerators{
    2,0,0,0,0,0,0,1,3,3,6,1,2,1,2,0,3,3,3,0,1,0,0,0,1,2,0,
    0,0,2,4,0,2,0,1,0,2,1,3,0,0,0,1,2,0,0,0,0,0,4,1,0,0};

struct SurgerySignature {
    std::array<u64, 7> words{};
    u64 inactive_count = 0;
};

std::vector<std::string> surgery_split(const std::string& value,
                                       char delimiter) {
    std::vector<std::string> pieces;
    std::stringstream stream(value);
    std::string piece;
    while (std::getline(stream, piece, delimiter)) pieces.push_back(piece);
    return pieces;
}

std::vector<u32> surgery_read_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open surgery deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid surgery deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "surgery deck rank/distinctness changed");
        deck.push_back(mask);
    }
    FnvLocal ledger;
    for (u32 mask : deck) ledger.add(mask);
    require(deck.size() == kSurgeryDeckSize &&
                ledger.state == kSurgeryDeckFnv,
            "surgery deck identity changed");
    return deck;
}

SurgerySignature surgery_read_signature(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open full-signature CSV");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "full-signature header changed");
    std::size_t rows = 0;
    bool found = false;
    SurgerySignature target;
    FnvLocal edge_ledger;
    int previous_q = 0;
    int previous_r = 0;
    while (std::getline(input, line)) {
        require(!line.empty(), "empty full-signature row");
        const std::vector<std::string> fields = surgery_split(line, ',');
        require(fields.size() == 10, "malformed full-signature row");
        const int q = std::stoi(fields[0]);
        const int r = std::stoi(fields[1]);
        const u64 count = std::stoull(fields[2]);
        require(q > 0 && q < r && count > 0,
                "invalid full-signature coordinates/count");
        require(rows == 0 || std::pair{previous_q, previous_r} <
                                 std::pair{q, r},
                "full-signature rows not strictly ordered");
        previous_q = q;
        previous_r = r;
        edge_ledger.add(static_cast<u64>(q));
        edge_ledger.add(static_cast<u64>(r));
        std::array<u64, 7> words{};
        u64 reconstructed_count = 0;
        for (std::size_t word = 0; word < words.size(); ++word) {
            std::size_t used = 0;
            words[word] = std::stoull(fields[3 + word], &used, 16);
            require(used == fields[3 + word].size(),
                    "bad signature word token");
            reconstructed_count += std::popcount(words[word]);
        }
        require(reconstructed_count == count && (words[6] >> 37) == 0,
                "signature count/padding changed");
        if (q == kTargetQ && r == kTargetR) {
            require(!found, "duplicate target signature");
            found = true;
            target.words = words;
            target.inactive_count = count;
        }
        ++rows;
    }
    require(rows == kSurgerySignatureRows &&
                edge_ledger.state == kSurgeryFinalResidualFnv && found,
            "full-signature universe/target changed");
    return target;
}

u32 surgery_next_same_popcount(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

struct LostBody {
    u32 body = 0;
    unsigned deleted_response = 0;
};

std::vector<LostBody> surgery_enumerate_lost_bodies(
    const std::vector<u32>& deck,
    const std::array<bool, kSurgeryDeckSize>& deleted,
    const std::filesystem::path& output_path,
    u64& retained_checks,
    u64& obligation_fnv,
    std::array<u64, 32>& response_histogram) {
    std::ofstream output(output_path);
    require(static_cast<bool>(output), "cannot create lost-body CSV");
    output << "ordinal,body_hex,deleted_response_hex\n";
    std::vector<LostBody> lost;
    FnvLocal ledger;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < kSurgeryBodyCount; ++ordinal) {
        bool retained_hit = false;
        for (std::size_t index = 0; index < deck.size(); ++index) {
            if (deleted[index]) continue;
            ++retained_checks;
            if ((body & deck[index]) == 0) {
                retained_hit = true;
                break;
            }
        }
        if (!retained_hit) {
            unsigned response = 0;
            unsigned local = 0;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (!deleted[index]) continue;
                if ((body & deck[index]) == 0) response |= 1u << local;
                ++local;
            }
            require(local == kExpectedDeleted.size() && response != 0,
                    "original deck failed body coverage");
            ledger.add(body);
            ledger.add(response);
            ++response_histogram[response];
            output << lost.size() << ',' << std::hex << std::setw(8)
                   << std::setfill('0') << body << ',' << std::setw(2)
                   << response << std::dec << std::setfill(' ') << '\n';
            lost.push_back({body, response});
        }
        if (ordinal + 1 < kSurgeryBodyCount)
            body = surgery_next_same_popcount(body);
    }
    require(body == UINT32_C(0x3fe00000) && output.good(),
            "body enumeration/output changed");
    obligation_fnv = ledger.state;
    return lost;
}

struct SurgeryClass {
    u64 pattern = 0;
    u64 multiplicity = 0;
    u32 least_mask = 0;
    unsigned cover = 0;
    bool maximal = false;
};

struct CoverSearch {
    const std::vector<u64>& patterns;
    u64 nodes = 0;
    u64 bound_prunes = 0;
    u64 memo_prunes = 0;
    std::unordered_map<u64, int> failed;
    std::vector<std::size_t> chosen;

    bool solve(u64 uncovered, int remaining) {
        ++nodes;
        if (uncovered == 0) return true;
        if (remaining == 0) return false;
        const auto memo = failed.find(uncovered);
        if (memo != failed.end() && memo->second >= remaining) {
            ++memo_prunes;
            return false;
        }
        unsigned maximum_gain = 0;
        for (u64 pattern : patterns)
            maximum_gain = std::max(maximum_gain,
                                    static_cast<unsigned>(
                                        std::popcount(pattern & uncovered)));
        if (maximum_gain == 0 ||
            (std::popcount(uncovered) + maximum_gain - 1) / maximum_gain >
                static_cast<unsigned>(remaining)) {
            ++bound_prunes;
            failed[uncovered] = std::max(failed[uncovered], remaining);
            return false;
        }

        unsigned selected_bit = 64;
        std::size_t selected_options = patterns.size() + 1;
        u64 scan = uncovered;
        while (scan) {
            const unsigned bit = std::countr_zero(scan);
            scan &= scan - 1;
            std::size_t options = 0;
            for (u64 pattern : patterns)
                options += ((pattern >> bit) & 1) != 0;
            if (options < selected_options) {
                selected_options = options;
                selected_bit = bit;
            }
        }
        require(selected_bit < 64 && selected_options > 0,
                "uncovered obligation has no response");

        std::vector<std::size_t> branches;
        for (std::size_t index = 0; index < patterns.size(); ++index)
            if ((patterns[index] >> selected_bit) & 1) branches.push_back(index);
        std::sort(branches.begin(), branches.end(), [&](std::size_t left,
                                                        std::size_t right) {
            const unsigned left_gain = std::popcount(patterns[left] & uncovered);
            const unsigned right_gain = std::popcount(patterns[right] & uncovered);
            if (left_gain != right_gain) return left_gain > right_gain;
            return patterns[left] < patterns[right];
        });
        std::unordered_set<u64> seen_restrictions;
        for (std::size_t index : branches) {
            const u64 restriction = patterns[index] & uncovered;
            if (!seen_restrictions.insert(restriction).second) continue;
            chosen.push_back(index);
            if (solve(uncovered & ~patterns[index], remaining - 1)) return true;
            chosen.pop_back();
        }
        failed[uncovered] = std::max(failed[uncovered], remaining);
        return false;
    }
};

bool surgery_find_packing(const std::array<u64, 64>& incompatible,
                          u64 candidates, unsigned need,
                          std::vector<unsigned>& chosen) {
    if (need == 0) return true;
    while (std::popcount(candidates) >= need) {
        const unsigned vertex = std::countr_zero(candidates);
        candidates &= candidates - 1;
        chosen.push_back(vertex);
        if (surgery_find_packing(incompatible,
                                 candidates & incompatible[vertex],
                                 need - 1, chosen))
            return true;
        chosen.pop_back();
    }
    return false;
}

struct CoverageAudit {
    u64 checks = 0;
    u64 maximum_prefix = 0;
    u64 failures = 0;
    u64 row_fnv = 0;
};

CoverageAudit surgery_rescan_body_coverage(const std::vector<u32>& deck) {
    CoverageAudit audit;
    FnvLocal row_ledger;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < kSurgeryBodyCount; ++ordinal) {
        u64 prefix = 0;
        for (u32 repair : deck) {
            ++prefix;
            ++audit.checks;
            if ((body & repair) == 0) break;
        }
        if (prefix == deck.size() && (body & deck.back()) != 0) {
            ++audit.failures;
            prefix = 0;
        } else {
            audit.maximum_prefix = std::max(audit.maximum_prefix, prefix);
        }
        row_ledger.add(body);
        row_ledger.add(prefix);
        if (ordinal + 1 < kSurgeryBodyCount)
            body = surgery_next_same_popcount(body);
    }
    require(body == UINT32_C(0x3fe00000), "rebuilt body endpoint changed");
    audit.row_fnv = row_ledger.state;
    return audit;
}

void surgery_write_deck(const std::filesystem::path& path,
                        const std::vector<u32>& deck) {
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create rebuilt deck");
    for (u32 mask : deck)
        output << std::hex << std::setw(8) << std::setfill('0') << mask
               << std::dec << std::setfill(' ') << '\n';
    require(output.good(), "failed writing rebuilt deck");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 7,
                "usage: signature-surgery SIGNATURE_CSV DECK LOST_CSV "
                "CLASS_TSV ACTIVE_MASK_TSV REBUILT_DECK");
        init_choose8_local();
        const SurgerySignature signature = surgery_read_signature(argv[1]);
        const std::vector<u32> deck = surgery_read_deck(argv[2]);

        std::array<bool, kSurgeryDeckSize> deleted{};
        std::vector<std::size_t> deleted_indices;
        for (std::size_t word = 0; word < signature.words.size(); ++word) {
            u64 bits = signature.words[word];
            while (bits) {
                const unsigned bit = std::countr_zero(bits);
                const std::size_t index = 64 * word + bit;
                require(index < deck.size(), "signature bit outside deck");
                deleted[index] = true;
                deleted_indices.push_back(index);
                bits &= bits - 1;
            }
        }
        require(signature.inactive_count == kExpectedDeleted.size() &&
                    deleted_indices == std::vector<std::size_t>(
                        kExpectedDeleted.begin(), kExpectedDeleted.end()),
                "target inactive signature changed");

        u64 retained_checks = 0;
        u64 obligation_fnv = 0;
        std::array<u64, 32> response_histogram{};
        const std::vector<LostBody> lost = surgery_enumerate_lost_bodies(
            deck, deleted, argv[3], retained_checks, obligation_fnv,
            response_histogram);
        require(!lost.empty() && lost.size() < 64,
                "lost-body universe outside one-word exact solver");
        const u64 full_universe = (UINT64_C(1) << lost.size()) - 1;

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cell count changed");
        const ActiveUniverse active =
            build_active_universe(cells, kTargetQ, kTargetR);
        u64 inactive_deck_count = 0;
        for (std::size_t index = 0; index < deck.size(); ++index) {
            const u64 rank = colex_rank8_local(deck[index]);
            require(rank < active.active.size(), "deck colex rank escaped");
            const bool inactive = active.active[rank] == 0;
            require(inactive == deleted[index],
                    "signature disagrees with complete active universe");
            inactive_deck_count += inactive;
        }
        require(inactive_deck_count == kExpectedDeleted.size(),
                "target inactive deck count changed");

        std::vector<u64> response_by_rank(EXPECTED_REPAIRS, 0);
        u64 complement_checks = 0;
        u64 active_incidences = 0;
        for (std::size_t obligation = 0; obligation < lost.size(); ++obligation) {
            enumerate_disjoint_repairs(lost[obligation].body,
                [&](u32, u64 rank) {
                    ++complement_checks;
                    if (!active.active[rank]) return;
                    response_by_rank[rank] |= UINT64_C(1) << obligation;
                    ++active_incidences;
                });
        }
        require(complement_checks == lost.size() * DISJOINT_REPAIRS_PER_BODY,
                "complement enumeration count changed");

        std::map<u64, SurgeryClass> class_map;
        FnvLocal active_response_ledger;
        u64 active_replayed = 0;
        u64 response_union = 0;
        u64 response_incidences_replayed = 0;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if (!active.active[rank]) continue;
            ++active_replayed;
            const u32 mask = unrank_colex8(rank);
            const u64 pattern = response_by_rank[rank];
            response_union |= pattern;
            response_incidences_replayed += std::popcount(pattern);
            active_response_ledger.add(mask);
            active_response_ledger.add(pattern);
            SurgeryClass& response_class = class_map[pattern];
            response_class.pattern = pattern;
            ++response_class.multiplicity;
            if (response_class.multiplicity == 1 || mask < response_class.least_mask)
                response_class.least_mask = mask;
        }
        require(active_replayed == active.count &&
                    response_incidences_replayed == active_incidences &&
                    response_union == full_universe,
                "active response quotient lost accounting/coverage");

        std::vector<u64> ordered_patterns;
        ordered_patterns.reserve(class_map.size());
        for (auto& [pattern, response_class] : class_map) {
            response_class.cover = std::popcount(pattern);
            ordered_patterns.push_back(pattern);
        }
        std::vector<u64> dominance_order;
        for (u64 pattern : ordered_patterns)
            if (pattern != 0) dominance_order.push_back(pattern);
        std::sort(dominance_order.begin(), dominance_order.end(),
                  [](u64 left, u64 right) {
            const unsigned lc = std::popcount(left);
            const unsigned rc = std::popcount(right);
            if (lc != rc) return lc > rc;
            return left < right;
        });
        std::vector<u64> maximal_patterns;
        for (u64 pattern : dominance_order) {
            bool dominated = false;
            for (u64 maximal : maximal_patterns) {
                if ((pattern & ~maximal) == 0) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) maximal_patterns.push_back(pattern);
        }
        for (u64 pattern : maximal_patterns) class_map[pattern].maximal = true;

        std::ofstream class_output(argv[4]);
        require(static_cast<bool>(class_output), "cannot create class TSV");
        class_output << "class\tpattern_hex\tleast_mask_hex\tmultiplicity"
                        "\tcover\tmaximal\n";
        std::unordered_map<u64, std::size_t> class_id;
        FnvLocal class_ledger;
        std::size_t class_index = 0;
        for (const auto& [pattern, response_class] : class_map) {
            class_id[pattern] = class_index;
            class_output << class_index << '\t' << std::hex << std::setw(14)
                         << std::setfill('0') << pattern << '\t'
                         << std::setw(8) << response_class.least_mask << std::dec
                         << std::setfill(' ') << '\t'
                         << response_class.multiplicity << '\t'
                         << response_class.cover << '\t'
                         << response_class.maximal << '\n';
            class_ledger.add(pattern);
            class_ledger.add(response_class.least_mask);
            class_ledger.add(response_class.multiplicity);
            class_ledger.add(response_class.cover);
            class_ledger.add(response_class.maximal);
            ++class_index;
        }
        require(class_output.good(), "failed writing class TSV");

        std::ofstream active_output(argv[5]);
        require(static_cast<bool>(active_output), "cannot create active-mask TSV");
        active_output << "colex_rank\tmask_hex\tclass\tpattern_hex\n";
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if (!active.active[rank]) continue;
            const u32 mask = unrank_colex8(rank);
            const u64 pattern = response_by_rank[rank];
            active_output << rank << '\t' << std::hex << std::setw(8)
                          << std::setfill('0') << mask << std::dec
                          << std::setfill(' ') << '\t' << class_id.at(pattern)
                          << '\t' << std::hex << std::setw(14)
                          << std::setfill('0') << pattern << std::dec
                          << std::setfill(' ') << '\n';
        }
        require(active_output.good(), "failed writing active-mask TSV");

        u64 uncovered = full_universe;
        std::vector<std::size_t> greedy;
        while (uncovered) {
            unsigned best_gain = 0;
            std::size_t best = maximal_patterns.size();
            for (std::size_t index = 0; index < maximal_patterns.size(); ++index) {
                const unsigned gain =
                    std::popcount(maximal_patterns[index] & uncovered);
                if (gain > best_gain ||
                    (gain == best_gain && best < maximal_patterns.size() &&
                     maximal_patterns[index] < maximal_patterns[best])) {
                    best_gain = gain;
                    best = index;
                }
            }
            require(best < maximal_patterns.size() && best_gain > 0,
                    "greedy cover stuck");
            greedy.push_back(best);
            uncovered &= ~maximal_patterns[best];
        }
        unsigned maximum_cover = 0;
        for (u64 pattern : maximal_patterns)
            maximum_cover = std::max(maximum_cover,
                                     static_cast<unsigned>(std::popcount(pattern)));
        const unsigned elementary_lower =
            (lost.size() + maximum_cover - 1) / maximum_cover;
        std::vector<std::tuple<unsigned, u64, u64, u64>> search_rows;
        std::vector<std::size_t> exact_choice;
        unsigned exact_minimum = 0;
        for (unsigned size = elementary_lower; size <= greedy.size(); ++size) {
            CoverSearch search{maximal_patterns};
            const bool feasible = search.solve(full_universe, size);
            search_rows.push_back({size, search.nodes, search.bound_prunes,
                                   search.memo_prunes});
            if (feasible) {
                exact_minimum = size;
                exact_choice = search.chosen;
                break;
            }
        }
        require(exact_minimum > 0, "exact cover search found no witness");
        u64 exact_union = 0;
        std::vector<u32> replacements;
        for (std::size_t index : exact_choice) {
            exact_union |= maximal_patterns[index];
            replacements.push_back(class_map[maximal_patterns[index]].least_mask);
        }
        require(exact_union == full_universe,
                "exact cover witness misses obligation");
        std::sort(replacements.begin(), replacements.end());
        require(std::adjacent_find(replacements.begin(), replacements.end()) ==
                    replacements.end(),
                "duplicate exact replacement");

        const int dual_numerator_sum = std::accumulate(
            kDualNumerators.begin(), kDualNumerators.end(), 0);
        int maximum_dual_class_sum = 0;
        for (const auto& [pattern, response_class] : class_map) {
            int class_sum = 0;
            u64 bits = pattern;
            while (bits) {
                const unsigned obligation = std::countr_zero(bits);
                bits &= bits - 1;
                class_sum += kDualNumerators[obligation];
            }
            maximum_dual_class_sum = std::max(maximum_dual_class_sum,
                                               class_sum);
        }
        const unsigned dual_lower =
            (dual_numerator_sum + kDualDenominator - 1) / kDualDenominator;
        require(dual_numerator_sum == 57 &&
                    maximum_dual_class_sum <= kDualDenominator &&
                    dual_lower == exact_minimum,
                "exact dual lower-bound certificate failed");

        std::array<u64, 64> co_covered{};
        for (u64 pattern : maximal_patterns) {
            u64 bits = pattern;
            while (bits) {
                const unsigned vertex = std::countr_zero(bits);
                bits &= bits - 1;
                co_covered[vertex] |= pattern;
            }
        }
        std::array<u64, 64> incompatible{};
        for (std::size_t vertex = 0; vertex < lost.size(); ++vertex)
            incompatible[vertex] =
                full_universe & ~co_covered[vertex] &
                ~(UINT64_C(1) << vertex);
        std::vector<unsigned> packing;
        unsigned packing_size = exact_minimum;
        while (packing_size > 0 &&
               !surgery_find_packing(incompatible, full_universe,
                                     packing_size, packing)) {
            --packing_size;
            packing.clear();
        }

        std::vector<u32> rebuilt;
        rebuilt.reserve(deck.size() - deleted_indices.size() + replacements.size());
        std::set<u32> rebuilt_distinct;
        for (std::size_t index = 0; index < deck.size(); ++index) {
            if (deleted[index]) continue;
            require(rebuilt_distinct.insert(deck[index]).second,
                    "retained deck duplicate");
            rebuilt.push_back(deck[index]);
        }
        for (u32 replacement : replacements) {
            const u64 rank = colex_rank8_local(replacement);
            require(active.active[rank] && rebuilt_distinct.insert(replacement).second,
                    "replacement inactive or already retained");
            rebuilt.push_back(replacement);
        }
        FnvLocal rebuilt_ledger;
        for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
        surgery_write_deck(argv[6], rebuilt);
        const CoverageAudit coverage = surgery_rescan_body_coverage(rebuilt);
        require(coverage.failures == 0,
                "rebuilt deck failed full body coverage");

        std::cout << "THM4281_GENERIC_SIGNATURE_SURGERY_ATLAS_V1\n"
                  << "TARGET " << kTargetQ << ',' << kTargetR
                  << " SIGNATURE_ROWS " << kSurgerySignatureRows
                  << " SIGNATURE_INACTIVE " << signature.inactive_count
                  << " DELETED_INDICES 57,107,222,275,345\n"
                  << "DELETED_MASKS";
        for (std::size_t index : deleted_indices)
            std::cout << ' ' << std::hex << std::setw(8) << std::setfill('0')
                      << deck[index] << std::dec << std::setfill(' ');
        std::cout << "\nBODIES " << kSurgeryBodyCount
                  << " RETAINED_CHECKS " << retained_checks
                  << " LOST_OBLIGATIONS " << lost.size()
                  << " OBLIGATION_FNV " << std::hex << obligation_fnv
                  << std::dec << '\n';
        for (unsigned response = 1; response < response_histogram.size(); ++response)
            if (response_histogram[response])
                std::cout << "DELETED_RESPONSE " << std::hex << response
                          << std::dec << " COUNT "
                          << response_histogram[response] << '\n';
        std::cout << "RANK8_UNIVERSE " << EXPECTED_REPAIRS
                  << " ACTIVE " << active.count << " ACTIVE_FNV "
                  << std::hex << active.fnv << std::dec
                  << " ZETA_OPERATIONS " << active.zeta_operations << '\n'
                  << "COMPLEMENT_CHECKS " << complement_checks
                  << " ACTIVE_INCIDENCES " << active_incidences
                  << " RESPONSE_CLASSES_ALL " << class_map.size()
                  << " NONEMPTY " << class_map.size() - class_map.contains(0)
                  << " MAXIMAL " << maximal_patterns.size()
                  << " CLASS_FNV " << std::hex << class_ledger.state
                  << " ACTIVE_RESPONSE_FNV " << active_response_ledger.state
                  << std::dec << '\n'
                  << "RESPONSE_UNION " << std::hex << response_union
                  << std::dec << " COVERED " << std::popcount(response_union)
                  << " OF " << lost.size() << " MAX_COVER " << maximum_cover
                  << " ELEMENTARY_LOWER " << elementary_lower
                  << " GREEDY_UPPER " << greedy.size() << '\n';
        for (const auto& [size, nodes, bound_prunes, memo_prunes] : search_rows)
            std::cout << "EXACT_SEARCH K " << size << " NODES " << nodes
                      << " BOUND_PRUNES " << bound_prunes
                      << " MEMO_PRUNES " << memo_prunes << '\n';
        std::cout << "MINIMUM_REPLACEMENTS " << exact_minimum
                  << " WITNESS_MASKS";
        for (u32 mask : replacements)
            std::cout << ' ' << std::hex << std::setw(8) << std::setfill('0')
                      << mask << std::dec << std::setfill(' ');
        std::cout << "\nWITNESS_RESPONSES";
        for (u32 mask : replacements) {
            const u64 pattern = response_by_rank[colex_rank8_local(mask)];
            std::cout << ' ' << std::hex << std::setw(8) << std::setfill('0')
                      << mask << ':' << std::setw(14) << pattern << std::dec
                      << std::setfill(' ');
        }
        std::cout << "\nDUAL_LOWER_BOUND NUMERATOR_SUM "
                  << dual_numerator_sum << " DENOMINATOR " << kDualDenominator
                  << " MAX_CLASS_NUMERATOR " << maximum_dual_class_sum
                  << " CEILING " << dual_lower << " WEIGHTS";
        for (std::size_t obligation = 0;
             obligation < kDualNumerators.size(); ++obligation)
            if (kDualNumerators[obligation] != 0)
                std::cout << ' ' << obligation << ':'
                          << kDualNumerators[obligation];
        std::cout << "\nPACKING_LOWER_BOUND " << packing_size
                  << " OBLIGATIONS";
        for (unsigned obligation : packing)
            std::cout << ' ' << obligation << ':' << std::hex << std::setw(8)
                      << std::setfill('0') << lost[obligation].body << std::dec
                      << std::setfill(' ');
        std::cout << "\nREBUILT_DECK " << rebuilt.size() << " FNV "
                  << std::hex << rebuilt_ledger.state << std::dec
                  << " TARGET_ALL_ACTIVE 1\n"
                  << "REBUILT_BODY_SCAN " << kSurgeryBodyCount
                  << " CHECKS " << coverage.checks
                  << " MAX_PREFIX " << coverage.maximum_prefix
                  << " FAILURES " << coverage.failures << " ROW_FNV "
                  << std::hex << coverage.row_fnv << std::dec << '\n'
                  << "VERDICT PASS EXACT_SIGNATURE_SURGERY_AND_MINIMUM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SIGNATURE_SURGERY_ERROR " << error.what() << '\n';
        return 1;
    }
}
