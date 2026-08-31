// Exact full response quotient and cardinality-minimum cover for the 101
// pair-body failures inherited from THM-4287.
#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <charconv>
#include <fstream>
#include <map>
#include <unordered_map>

namespace {

#ifndef FAILURE_COUNT
#define FAILURE_COUNT 101
#endif
#ifndef FAILURE_ENDPOINT
#define FAILURE_ENDPOINT 636
#endif
constexpr std::size_t kFailureCount = FAILURE_COUNT;
constexpr int kFailureEndpoint = FAILURE_ENDPOINT;
constexpr std::size_t kWords = 2;
static_assert(kFailureCount > 64 && kFailureCount < 128);

struct Failure {
    int q = 0;
    int r = 0;
    u32 body = 0;
};

std::vector<Failure> read_failure_csv(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure csv");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure csv header changed");
    std::vector<Failure> rows;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank failure row");
        const std::size_t a = line.find(',');
        const std::size_t b = line.find(',', a + 1);
        require(a != std::string::npos && b != std::string::npos &&
                    line.find(',', b + 1) == std::string::npos,
                "bad failure row");
        Failure row;
        row.q = std::stoi(line.substr(0, a));
        row.r = std::stoi(line.substr(a + 1, b - a - 1));
        row.body = static_cast<u32>(std::stoul(line.substr(b + 1), nullptr, 16));
        require(std::popcount(row.body) == 9 &&
                    ((row.q == 100 && row.r == kFailureEndpoint) ||
                     (row.q == 256 && row.r == kFailureEndpoint)),
                "failure row outside expected universe");
        rows.push_back(row);
    }
    require(rows.size() == kFailureCount, "failure count changed");
    return rows;
}

struct PatternHash {
    std::size_t operator()(const std::array<u64, kWords>& x) const noexcept {
        u64 h = x[0] ^ (x[1] + UINT64_C(0x9e3779b97f4a7c15) +
                        (x[0] << 6) + (x[0] >> 2));
        return static_cast<std::size_t>(h);
    }
};

struct PatternInfo {
    u64 multiplicity = 0;
    u32 least = 0;
};

struct MaximalPattern {
    std::array<u64, kWords> bits{};
    u32 least = 0;
    u64 multiplicity = 0;
};

unsigned cover(const std::array<u64, kWords>& p) {
    return std::popcount(p[0]) + std::popcount(p[1]);
}

bool subset(const std::array<u64, kWords>& a,
            const std::array<u64, kWords>& b) {
    return (a[0] & ~b[0]) == 0 && (a[1] & ~b[1]) == 0;
}

std::array<u64, kWords> join(const std::array<u64, kWords>& a,
                             const std::array<u64, kWords>& b) {
    return {a[0] | b[0], a[1] | b[1]};
}

unsigned gain(const std::array<u64, kWords>& pattern,
              const std::array<u64, kWords>& covered) {
    return std::popcount(pattern[0] & ~covered[0]) +
           std::popcount(pattern[1] & ~covered[1]);
}

struct ExactCoverSearch {
    const std::vector<MaximalPattern>& patterns;
    std::array<u64, kWords> full{};
    std::array<std::vector<std::size_t>, kFailureCount> by_obligation;
    std::unordered_map<std::array<u64, kWords>, int, PatternHash> seen;
    std::vector<std::size_t> path;
    u64 calls = 0;

    explicit ExactCoverSearch(const std::vector<MaximalPattern>& input)
        : patterns(input) {
        full = {UINT64_MAX,
                (UINT64_C(1) << (kFailureCount - 64)) - 1};
        for (std::size_t p = 0; p < patterns.size(); ++p) {
            for (std::size_t i = 0; i < kFailureCount; ++i) {
                if ((patterns[p].bits[i >> 6] >> (i & 63)) & 1)
                    by_obligation[i].push_back(p);
            }
        }
    }

    bool dfs(const std::array<u64, kWords>& covered, int remaining) {
        ++calls;
        if (covered == full) return true;
        if (remaining == 0) return false;
        const auto prior = seen.find(covered);
        if (prior != seen.end() && prior->second >= remaining) return false;
        seen[covered] = remaining;
        const unsigned debt = kFailureCount - cover(covered);
        unsigned maximum_gain = 0;
        for (const auto& pattern : patterns)
            maximum_gain = std::max(maximum_gain, gain(pattern.bits, covered));
        if (maximum_gain == 0 ||
            (debt + maximum_gain - 1) / maximum_gain >
                static_cast<unsigned>(remaining))
            return false;

        std::size_t chosen = kFailureCount;
        std::size_t choices = patterns.size() + 1;
        for (std::size_t i = 0; i < kFailureCount; ++i) {
            if ((covered[i >> 6] >> (i & 63)) & 1) continue;
            std::size_t viable = 0;
            for (std::size_t p : by_obligation[i])
                viable += gain(patterns[p].bits, covered) != 0;
            if (viable < choices) {
                choices = viable;
                chosen = i;
            }
        }
        require(chosen < kFailureCount && choices > 0,
                "uncovered obligation has no maximal responder");
        std::vector<std::size_t> candidates = by_obligation[chosen];
        std::sort(candidates.begin(), candidates.end(), [&](std::size_t a,
                                                            std::size_t b) {
            const unsigned ga = gain(patterns[a].bits, covered);
            const unsigned gb = gain(patterns[b].bits, covered);
            if (ga != gb) return ga > gb;
            return patterns[a].least < patterns[b].least;
        });
        for (std::size_t p : candidates) {
            const auto next = join(covered, patterns[p].bits);
            if (next == covered) continue;
            path.push_back(p);
            if (dfs(next, remaining - 1)) return true;
            path.pop_back();
        }
        return false;
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2 || argc == 3,
                "usage: failure-activity FAILURES.csv [MAXIMAL.tsv]");
        init_choose8_local();
        const auto failures = read_failure_csv(argv[1]);
        const auto cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const ActiveUniverse active100 =
            build_active_universe(cells, 100, kFailureEndpoint);
        const ActiveUniverse active256 =
            build_active_universe(cells, 256, kFailureEndpoint);
        std::cout << "THM4295_ENDPOINT636_FAILURE_RESPONSE_V1\n"
                  << "ACTIVE 100," << kFailureEndpoint << " COUNT "
                  << active100.count << " FNV "
                  << std::hex << active100.fnv << std::dec << '\n'
                  << "ACTIVE 256," << kFailureEndpoint << " COUNT "
                  << active256.count << " FNV "
                  << std::hex << active256.fnv << std::dec << '\n';

        std::vector<std::array<u64, kWords>> patterns(EXPECTED_REPAIRS);
        std::array<u64, kFailureCount> own_counts{};
        std::array<u64, kFailureCount> common_counts{};
        std::array<u32, kFailureCount> least_own{};
        std::array<u32, kFailureCount> least_common{};
        FnvLocal incidence_ledger;
        u64 total_own = 0;
        u64 total_common = 0;
        for (std::size_t i = 0; i < failures.size(); ++i) {
            const Failure& failure = failures[i];
            const auto& own = failure.q == 100 ? active100.active
                                               : active256.active;
            enumerate_disjoint_repairs(failure.body, [&](u32 mask, u64 rank) {
                if (own[rank]) {
                    patterns[rank][i >> 6] |= UINT64_C(1) << (i & 63);
                    ++own_counts[i];
                    ++total_own;
                    if (least_own[i] == 0 || mask < least_own[i])
                        least_own[i] = mask;
                    incidence_ledger.add(i);
                    incidence_ledger.add(mask);
                }
                if (active100.active[rank] && active256.active[rank]) {
                    ++common_counts[i];
                    ++total_common;
                    if (least_common[i] == 0 || mask < least_common[i])
                        least_common[i] = mask;
                }
            });
        }

        std::unordered_map<std::array<u64, kWords>, PatternInfo, PatternHash>
            classes;
        classes.reserve(EXPECTED_REPAIRS / 2);
        u64 nonempty = 0;
        u64 max_cover = 0;
        u32 max_cover_least = 0;
        std::array<u64, kWords> union_pattern{};
        FnvLocal nonempty_ledger;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const auto pattern = patterns[rank];
            if (pattern[0] == 0 && pattern[1] == 0) continue;
            const u32 mask = unrank_colex8(rank);
            ++nonempty;
            nonempty_ledger.add(mask);
            union_pattern[0] |= pattern[0];
            union_pattern[1] |= pattern[1];
            PatternInfo& info = classes[pattern];
            ++info.multiplicity;
            if (info.least == 0 || mask < info.least) info.least = mask;
            const u64 c = cover(pattern);
            if (c > max_cover || (c == max_cover && mask < max_cover_least)) {
                max_cover = c;
                max_cover_least = mask;
            }
        }
        const std::array<u64, kWords> full = {
            UINT64_MAX,
            (UINT64_C(1) << (kFailureCount - 64)) - 1,
        };
        const bool all_responsive = union_pattern == full;

        std::vector<std::pair<std::array<u64, kWords>, PatternInfo>> ordered;
        ordered.reserve(classes.size());
        for (const auto& item : classes) ordered.push_back(item);
        std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
            if (cover(a.first) != cover(b.first))
                return cover(a.first) > cover(b.first);
            if (a.second.least != b.second.least)
                return a.second.least < b.second.least;
            return a.first < b.first;
        });
        std::vector<std::pair<std::array<u64, kWords>, PatternInfo>> maximal_raw;
        for (const auto& item : ordered) {
            bool dominated = false;
            for (const auto& prior : maximal_raw) {
                if (subset(item.first, prior.first)) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) maximal_raw.push_back(item);
        }
        std::vector<MaximalPattern> maximal;
        maximal.reserve(maximal_raw.size());
        for (const auto& [bits, info] : maximal_raw)
            maximal.push_back({bits, info.least, info.multiplicity});
        if (argc == 3) {
            std::ofstream output(argv[2]);
            require(static_cast<bool>(output), "cannot create maximal tsv");
            output << "mask_hex\tword_hi\tword_lo\tcover\tmultiplicity\n";
            for (const auto& item : maximal) {
                output << std::hex << std::setw(8) << std::setfill('0')
                       << item.least << '\t' << std::setw(16) << item.bits[1]
                       << '\t' << std::setw(16) << item.bits[0] << std::dec
                       << std::setfill(' ') << '\t' << cover(item.bits) << '\t'
                       << item.multiplicity << '\n';
            }
            require(output.good(), "failed writing maximal tsv");
        }

        std::array<u64, kWords> greedy_covered{};
        std::vector<std::size_t> greedy;
        while (all_responsive && greedy_covered != full) {
            std::size_t best = maximal.size();
            unsigned best_gain = 0;
            for (std::size_t i = 0; i < maximal.size(); ++i) {
                const unsigned g = gain(maximal[i].bits, greedy_covered);
                if (g > best_gain ||
                    (g == best_gain && best < maximal.size() &&
                     maximal[i].least < maximal[best].least)) {
                    best = i;
                    best_gain = g;
                }
            }
            require(best < maximal.size() && best_gain > 0,
                    "greedy response cover stalled");
            greedy.push_back(best);
            greedy_covered = join(greedy_covered, maximal[best].bits);
        }

        std::vector<std::size_t> exact;
        int exact_depth = 0;
        u64 exact_calls = 0;
        const int global_lower =
            (kFailureCount + static_cast<int>(max_cover) - 1) /
            static_cast<int>(max_cover);
        for (int depth = global_lower;
             all_responsive && depth <= static_cast<int>(greedy.size()); ++depth) {
            ExactCoverSearch search(maximal);
            if (search.dfs({}, depth)) {
                exact = search.path;
                exact_depth = depth;
                exact_calls += search.calls;
                break;
            }
            exact_calls += search.calls;
        }
        require(!all_responsive || exact_depth > 0,
                "exact cover search missed greedy witness");

        std::cout << "FAILURES " << failures.size()
                  << " COMPLEMENT_CHECKS "
                  << failures.size() * DISJOINT_REPAIRS_PER_BODY
                  << " OWN_ACTIVE_INCIDENCES " << total_own
                  << " COMMON_ACTIVE_INCIDENCES " << total_common << '\n'
                  << "NONEMPTY_MASKS " << nonempty << " FNV " << std::hex
                  << nonempty_ledger.state << std::dec
                  << " RESPONSE_CLASSES " << classes.size()
                  << " MAXIMAL_CLASSES " << maximal.size()
                  << " MAX_COVER " << max_cover << " LEAST_MAX " << std::hex
                  << std::setw(8) << std::setfill('0') << max_cover_least
                  << std::dec << std::setfill(' ') << '\n';
        std::cout << "RESPONSIVE " << cover(union_pattern)
                  << " UNRESPONSIVE " << kFailureCount - cover(union_pattern)
                  << " INDICES";
        for (std::size_t i = 0; i < kFailureCount; ++i)
            if (((union_pattern[i >> 6] >> (i & 63)) & 1) == 0)
                std::cout << ' ' << i;
        std::cout << '\n' << "GREEDY_COVER " << greedy.size() << " MASKS";
        for (std::size_t p : greedy)
            std::cout << ' ' << std::hex << std::setw(8) << std::setfill('0')
                      << maximal[p].least << std::dec << std::setfill(' ');
        std::cout << '\n' << "EXACT_MINIMUM " << exact_depth << " CALLS "
                  << exact_calls << " MASKS";
        for (std::size_t p : exact)
            std::cout << ' ' << std::hex << std::setw(8) << std::setfill('0')
                      << maximal[p].least << std::dec << std::setfill(' ');
        std::cout << '\n';

        u32 union100 = 0;
        for (const Failure& failure : failures)
            if (failure.q == 100) union100 |= failure.body;
        const u32 complement100 = ((u32{1} << 30) - 1) ^ union100;
        std::cout << "UNION_COMPLEMENT_100 " << std::hex << std::setw(8)
                  << std::setfill('0') << complement100 << std::dec
                  << std::setfill(' ') << " RANK "
                  << std::popcount(complement100);
        if (std::popcount(complement100) == 8) {
            const u64 universal_rank = colex_rank8_local(complement100);
            std::cout << " ACTIVE100 "
                      << static_cast<int>(active100.active[universal_rank])
                      << " ACTIVE256 "
                      << static_cast<int>(active256.active[universal_rank])
                      << " COVER " << cover(patterns[universal_rank]);
        }
        std::cout << '\n';
        const u32 thm4286 = UINT32_C(0x042022c9);
        const u64 thm4286_rank = colex_rank8_local(thm4286);
        std::cout << "THM4286_ORDINAL_MASK 042022c9 ACTIVE100 "
                  << static_cast<int>(active100.active[thm4286_rank])
                  << " ACTIVE256 "
                  << static_cast<int>(active256.active[thm4286_rank])
                  << " COVER " << cover(patterns[thm4286_rank]) << '\n';

        for (std::size_t i = 0; i < failures.size(); ++i) {
            std::cout << "BODY " << i << " PAIR " << failures[i].q << ','
                      << failures[i].r << " MASK " << std::hex << std::setw(8)
                      << std::setfill('0') << failures[i].body << std::dec
                      << std::setfill(' ') << " OWN " << own_counts[i]
                      << " LEAST_OWN " << std::hex << std::setw(8)
                      << std::setfill('0') << least_own[i] << std::dec
                      << std::setfill(' ') << " COMMON " << common_counts[i]
                      << " LEAST_COMMON " << std::hex << std::setw(8)
                      << std::setfill('0') << least_common[i] << std::dec
                      << std::setfill(' ') << '\n';
        }

        std::cout << "TOP_MAXIMAL_CLASSES\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(20, maximal.size()); ++i) {
            const auto& item = maximal[i];
            std::cout << "  " << i << " COVER " << cover(item.bits)
                      << " MULT " << item.multiplicity << " LEAST "
                      << std::hex << std::setw(8) << std::setfill('0')
                      << item.least << " WORDS " << std::setw(16) << item.bits[1]
                      << ',' << std::setw(16) << item.bits[0] << std::dec
                      << std::setfill(' ') << '\n';
        }
        std::cout << "INCIDENCE_FNV " << std::hex << incidence_ledger.state
                  << std::dec << "\nVERDICT PASS FINITE_EXACT_AUDIT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4295_FAILURE_RESPONSE_ERROR " << error.what() << '\n';
        return 1;
    }
}
