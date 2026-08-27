// Exact multiword active-response cover for the 358 bodies left at
// (256,650) by the frozen THM-4281 carrier plus the first eleven THM-4282
// continuation masks.  The greedy cover is an explicit upper bound only;
// no minimum-cardinality claim is made.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace {

constexpr int kQ = 256;
constexpr int kR = 650;
constexpr std::size_t kFailures = 358;
constexpr std::size_t kWords = (kFailures + 63) / 64;
constexpr u64 kFailureFnv = UINT64_C(0x609ee065cee1f551);
constexpr std::size_t kBaseCarrierCount = 8951;
constexpr u64 kBaseCarrierFnv = UINT64_C(0x188f82ab9dd1695a);
constexpr std::size_t kPriorAdditionCount = 11;
constexpr u64 kPriorCarrierFnv = UINT64_C(0x5b9e9a81c1582d06);

using Pattern = std::array<u64, kWords>;

struct Candidate {
    u64 rank = 0;
    u32 mask = 0;
};

std::vector<u32> read_masks(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "mask rank/distinctness changed");
        masks.push_back(mask);
    }
    return masks;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure CSV");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure CSV header changed");
    std::vector<u32> failures;
    FnvLocal ledger;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        int q = 0;
        int r = 0;
        char comma1 = 0;
        char comma2 = 0;
        std::string body_hex;
        row >> q >> comma1 >> r >> comma2 >> body_hex;
        require(static_cast<bool>(row) && q == kQ && r == kR &&
                    comma1 == ',' && comma2 == ',',
                "failure row coordinates changed");
        std::size_t used = 0;
        const u64 wide = std::stoull(body_hex, &used, 16);
        require(used == body_hex.size() && wide < (UINT64_C(1) << 30),
                "failure body token changed");
        const u32 body = static_cast<u32>(wide);
        require(std::popcount(body) == 9 &&
                    (failures.empty() || failures.back() < body),
                "failure body rank/order changed");
        failures.push_back(body);
        ledger.add(body);
    }
    require(failures.size() == kFailures && ledger.state == kFailureFnv,
            "failure ledger identity changed");
    return failures;
}

unsigned gain(const Pattern& response, const Pattern& uncovered) {
    unsigned answer = 0;
    for (std::size_t word = 0; word < kWords; ++word)
        answer += std::popcount(response[word] & uncovered[word]);
    return answer;
}

bool empty(const Pattern& pattern) {
    for (u64 word : pattern)
        if (word != 0) return false;
    return true;
}

void unite(Pattern& target, const Pattern& source) {
    for (std::size_t word = 0; word < kWords; ++word)
        target[word] |= source[word];
}

void subtract(Pattern& target, const Pattern& source) {
    for (std::size_t word = 0; word < kWords; ++word)
        target[word] &= ~source[word];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 5,
                "usage: endpoint650-large FAILURE_CSV BASE_CARRIER "
                "PRIOR_ADDITIONS NEW_ADDITIONS_OUT");
        init_choose8_local();
        const std::vector<u32> failures = read_failures(argv[1]);
        std::vector<u32> carrier = read_masks(argv[2]);
        require(carrier.size() == kBaseCarrierCount &&
                    mask_fnv(carrier) == kBaseCarrierFnv,
                "base carrier identity changed");
        const std::vector<u32> prior = read_masks(argv[3]);
        require(prior.size() == kPriorAdditionCount,
                "prior addition count changed");
        std::set<u32> carrier_set(carrier.begin(), carrier.end());
        for (u32 mask : prior) {
            require(carrier_set.insert(mask).second,
                    "prior addition overlaps carrier");
            carrier.push_back(mask);
        }
        require(mask_fnv(carrier) == kPriorCarrierFnv,
                "prior augmented carrier identity changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        const ActiveUniverse active = build_active_universe(cells, kQ, kR);
        for (u32 body : failures)
            for (u32 mask : carrier)
                require(!active.active[colex_rank8_local(mask)] ||
                            (body & mask) != 0,
                        "failure body already covered by active carrier");

        std::vector<Pattern> response(EXPECTED_REPAIRS);
        u64 complement_checks = 0;
        u64 active_incidences = 0;
        for (std::size_t obligation = 0; obligation < failures.size();
             ++obligation) {
            enumerate_disjoint_repairs(failures[obligation],
                [&](u32, u64 rank) {
                    ++complement_checks;
                    if (!active.active[rank]) return;
                    response[rank][obligation / 64] |=
                        UINT64_C(1) << (obligation % 64);
                    ++active_incidences;
                });
        }
        require(complement_checks ==
                    failures.size() * DISJOINT_REPAIRS_PER_BODY,
                "complement enumeration changed");

        std::vector<Candidate> candidates;
        candidates.reserve(active.count);
        FnvLocal response_ledger;
        unsigned maximum_response = 0;
        u32 least_maximum = 0;
        Pattern response_union{};
        std::vector<Pattern> co_covered(kFailures);
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if (!active.active[rank] || empty(response[rank])) continue;
            const u32 mask = unrank_colex8(rank);
            candidates.push_back({rank, mask});
            response_ledger.add(mask);
            unsigned cover = 0;
            for (u64 word : response[rank]) {
                response_ledger.add(word);
                cover += std::popcount(word);
            }
            if (cover > maximum_response ||
                (cover == maximum_response &&
                 (least_maximum == 0 || mask < least_maximum))) {
                maximum_response = cover;
                least_maximum = mask;
            }
            unite(response_union, response[rank]);
            for (std::size_t word = 0; word < kWords; ++word) {
                u64 bits = response[rank][word];
                while (bits) {
                    const unsigned bit = std::countr_zero(bits);
                    bits &= bits - 1;
                    unite(co_covered[64 * word + bit], response[rank]);
                }
            }
        }
        Pattern full{};
        full.fill(UINT64_MAX);
        full.back() = (UINT64_C(1) << (kFailures % 64)) - 1;
        require(response_union == full,
                "some failure has no active rank-eight response");

        Pattern uncovered = full;
        std::vector<Candidate> greedy;
        std::vector<unsigned> greedy_gains;
        while (!empty(uncovered)) {
            unsigned best_gain = 0;
            Candidate best;
            for (const Candidate candidate : candidates) {
                const unsigned candidate_gain =
                    gain(response[candidate.rank], uncovered);
                if (candidate_gain > best_gain ||
                    (candidate_gain == best_gain && candidate_gain != 0 &&
                     candidate.mask < best.mask)) {
                    best_gain = candidate_gain;
                    best = candidate;
                }
            }
            require(best_gain != 0, "greedy response cover stuck");
            greedy.push_back(best);
            greedy_gains.push_back(best_gain);
            subtract(uncovered, response[best.rank]);
        }

        std::array<unsigned, kFailures> cover_counts{};
        for (const Candidate candidate : greedy)
            for (std::size_t word = 0; word < kWords; ++word) {
                u64 bits = response[candidate.rank][word];
                while (bits) {
                    const unsigned bit = std::countr_zero(bits);
                    bits &= bits - 1;
                    ++cover_counts[64 * word + bit];
                }
            }
        for (std::size_t index = greedy.size(); index-- > 0;) {
            bool removable = true;
            for (std::size_t word = 0; word < kWords && removable; ++word) {
                u64 bits = response[greedy[index].rank][word];
                while (bits) {
                    const unsigned bit = std::countr_zero(bits);
                    bits &= bits - 1;
                    if (cover_counts[64 * word + bit] == 1) {
                        removable = false;
                        break;
                    }
                }
            }
            if (!removable) continue;
            for (std::size_t word = 0; word < kWords; ++word) {
                u64 bits = response[greedy[index].rank][word];
                while (bits) {
                    const unsigned bit = std::countr_zero(bits);
                    bits &= bits - 1;
                    --cover_counts[64 * word + bit];
                }
            }
            greedy.erase(greedy.begin() + index);
            greedy_gains.erase(greedy_gains.begin() + index);
        }
        require(std::ranges::all_of(cover_counts,
                                    [](unsigned count) { return count > 0; }),
                "pruned witness lost an obligation");

        Pattern packing_candidates = full;
        std::vector<std::size_t> packing;
        while (!empty(packing_candidates)) {
            std::size_t best = kFailures;
            unsigned best_neighbours = UINT_MAX;
            for (std::size_t word = 0; word < kWords; ++word) {
                u64 bits = packing_candidates[word];
                while (bits) {
                    const unsigned bit = std::countr_zero(bits);
                    bits &= bits - 1;
                    const std::size_t vertex = 64 * word + bit;
                    const unsigned neighbours =
                        gain(co_covered[vertex], packing_candidates);
                    if (neighbours < best_neighbours) {
                        best_neighbours = neighbours;
                        best = vertex;
                    }
                }
            }
            require(best < kFailures, "packing selection stuck");
            packing.push_back(best);
            subtract(packing_candidates, co_covered[best]);
        }

        std::ofstream output(argv[4]);
        require(static_cast<bool>(output), "cannot create new additions");
        std::vector<u32> new_masks;
        for (const Candidate candidate : greedy) {
            require(carrier_set.insert(candidate.mask).second,
                    "greedy response overlaps prior carrier");
            new_masks.push_back(candidate.mask);
        }
        std::sort(new_masks.begin(), new_masks.end());
        for (u32 mask : new_masks)
            output << std::hex << std::setw(8) << std::setfill('0') << mask
                   << '\n';
        require(output.good(), "failed writing new additions");
        carrier.insert(carrier.end(), new_masks.begin(), new_masks.end());

        std::cout << "THM4282_ENDPOINT650_LARGE_RESPONSE_GREEDY_V1\n"
                  << "PAIR 256,650 FAILURES " << failures.size()
                  << " FAILURE_FNV " << std::hex << kFailureFnv << std::dec
                  << " ACTIVE " << active.count << " ACTIVE_FNV " << std::hex
                  << active.fnv << std::dec << '\n'
                  << "COMPLEMENT_CHECKS " << complement_checks
                  << " ACTIVE_INCIDENCES " << active_incidences
                  << " NONEMPTY_CANDIDATES " << candidates.size()
                  << " RESPONSE_FNV " << std::hex << response_ledger.state
                  << std::dec << " MAX_RESPONSE " << maximum_response
                  << " LEAST_MAXIMUM " << std::hex << std::setw(8)
                  << std::setfill('0') << least_maximum << std::dec
                  << std::setfill(' ') << '\n'
                  << "PACKING_LOWER_BOUND " << packing.size()
                  << " OBLIGATIONS";
        for (std::size_t obligation : packing)
            std::cout << ' ' << obligation << ':' << std::hex
                      << std::setw(8) << std::setfill('0')
                      << failures[obligation] << std::dec << std::setfill(' ');
        std::cout << "\nGREEDY_PRUNED_UPPER " << new_masks.size()
                  << " MASKS";
        for (u32 mask : new_masks)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << mask << std::dec
                      << std::setfill(' ');
        std::cout << "\nREBUILT_CARRIER " << carrier.size() << " FNV "
                  << std::hex << mask_fnv(carrier) << std::dec << '\n'
                  << "VERDICT PASS EXPLICIT_ACTIVE_RESPONSE_COVER_NO_MINIMUM_CLAIM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT650_LARGE_RESPONSE_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
