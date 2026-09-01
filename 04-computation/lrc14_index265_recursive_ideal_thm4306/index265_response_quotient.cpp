// Exact audit of the largest failed singleton-signature ideal for THM-4306.
//
// The import supplies the maintained import-free literal-wall geometry and
// colex/zeta primitives.  Every index-265 object (ideal, private bodies,
// common-active intersection, response quotient, exact cover, and rebuilt
// deck) is derived afresh below.

#define main detached_j19_unused_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/detached_j19_audit.cpp"
#undef main

namespace {

constexpr std::size_t kTarget265 = 265;
constexpr u32 kOldMask265 = UINT32_C(0x20820649);
constexpr std::size_t kIdealCount265 = 367;
constexpr u64 kIdealFnv265 = UINT64_C(0xd422b161d94ebae4);

std::vector<u32> read_deck265(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "invalid or repeated deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof() && deck.size() == kDeckCount &&
                ledger.state == kDeckFnv,
            "joint deck identity changed");
    require(deck[kTarget265] == kOldMask265,
            "index-265 old mask changed");
    return deck;
}

struct Ideal265 {
    std::vector<Pair> rows;
    u64 fnv = 0;
    int maximum_endpoint = 0;
};

Ideal265 read_ideal265(const std::filesystem::path& path,
                       const std::set<Pair>& live,
                       const std::filesystem::path& output_path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    Ideal265 out;
    Fnv ledger;
    std::size_t atlas_rows = 0;
    Pair previous{-1, -1};
    while (std::getline(input, line)) {
        ++atlas_rows;
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 10, "malformed signature row");
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        require(previous < pair, "signature atlas order changed");
        previous = pair;
        std::array<u64, 7> words{};
        unsigned weight = 0;
        for (std::size_t word = 0; word < words.size(); ++word) {
            words[word] = std::stoull(fields[3 + word], nullptr, 16);
            weight += std::popcount(words[word]);
        }
        require(weight == std::stoul(fields[2]),
                "inactive count disagrees with signature words");
        if (!live.contains(pair) || weight != 1) continue;
        bool exact = true;
        for (std::size_t word = 0; word < words.size(); ++word) {
            const u64 wanted = word == kTarget265 / 64
                ? UINT64_C(1) << (kTarget265 % 64) : 0;
            exact &= words[word] == wanted;
        }
        if (!exact) continue;
        out.rows.push_back(pair);
        out.maximum_endpoint = std::max(out.maximum_endpoint, pair.r);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(input.eof() && atlas_rows == kAtlasRows,
            "signature atlas size changed");
    out.fnv = ledger.state;
    require(out.rows.size() == kIdealCount265 && out.fnv == kIdealFnv265,
            "index-265 ideal identity changed");

    std::ofstream output(output_path, std::ios::binary);
    require(static_cast<bool>(output), "cannot create ideal ledger");
    for (Pair pair : out.rows) output << pair.q << ',' << pair.r << '\n';
    require(output.good(), "ideal ledger write failed");
    return out;
}

struct Obligations265 {
    std::vector<u32> bodies;
    u64 target_disjoint = 0;
    u64 retained_checks = 0;
    u64 fnv = 0;
    u32 body_union = 0;
};

Obligations265 private_obligations265(const std::vector<u32>& deck) {
    Obligations265 out;
    Fnv ledger;
    u32 body = (u32{1} << 9) - 1;
    for (u64 ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        if ((body & deck[kTarget265]) == 0) {
            ++out.target_disjoint;
            bool retained = false;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == kTarget265) continue;
                ++out.retained_checks;
                if ((body & deck[index]) == 0) {
                    retained = true;
                    break;
                }
            }
            if (!retained) {
                out.bodies.push_back(body);
                out.body_union |= body;
                ledger.add(body);
            }
        }
        if (ordinal + 1 < kBodyCount) body = next_combination(body);
    }
    require(body == UINT32_C(0x3fe00000) &&
                out.target_disjoint == UINT64_C(497420) &&
                out.bodies.size() == 8,
            "index-265 private-body census changed");
    out.fnv = ledger.state;
    return out;
}

struct ResponseClass265 {
    u64 count = 0;
    u32 least = 0;
};

struct ExactCover265 {
    int minimum = 1000;
    std::vector<unsigned> patterns;
};

ExactCover265 exact_cover265(const std::vector<unsigned>& patterns,
                             unsigned full) {
    std::array<int, 256> distance{};
    std::array<int, 256> prior{};
    std::array<unsigned, 256> step{};
    distance.fill(1000);
    prior.fill(-1);
    distance[0] = 0;
    std::queue<unsigned> queue;
    queue.push(0);
    while (!queue.empty()) {
        const unsigned state = queue.front();
        queue.pop();
        for (unsigned pattern : patterns) {
            const unsigned next = state | pattern;
            if (distance[next] <= distance[state] + 1) continue;
            distance[next] = distance[state] + 1;
            prior[next] = static_cast<int>(state);
            step[next] = pattern;
            queue.push(next);
        }
    }
    require(distance[full] < 1000, "response family does not cover obligations");
    ExactCover265 out;
    out.minimum = distance[full];
    for (unsigned state = full; state != 0;
         state = static_cast<unsigned>(prior[state]))
        out.patterns.push_back(step[state]);
    std::reverse(out.patterns.begin(), out.patterns.end());
    return out;
}

struct Packing265 {
    unsigned support = 0;
    unsigned size = 0;
};

Packing265 largest_integral_packing265(
        const std::map<unsigned, ResponseClass265>& classes) {
    Packing265 best;
    for (unsigned support = 1; support < 256; ++support) {
        bool valid = true;
        for (const auto& [response, cls] : classes) {
            (void)cls;
            if (std::popcount(response & support) > 1) {
                valid = false;
                break;
            }
        }
        if (valid && std::popcount(support) > best.size)
            best = {support, static_cast<unsigned>(std::popcount(support))};
    }
    return best;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 7,
                "usage: index265 DECK SIGNATURES LIVE IDEAL_OUT QUOTIENT_OUT THREADS");
        init_choose();
        const std::vector<u32> deck = read_deck265(argv[1]);
        const std::set<Pair> live = read_live(argv[3]);
        const Ideal265 ideal = read_ideal265(argv[2], live, argv[4]);
        const Obligations265 obligations = private_obligations265(deck);
        require(obligations.bodies.size() <= 8,
                "response encoding exceeds one byte");

        const unsigned thread_count = static_cast<unsigned>(std::stoul(argv[6]));
        require(thread_count >= 1 && thread_count <= 8,
                "thread count outside 1..8");

        std::vector<u64> common((kRepairCount + 63) / 64, ~UINT64_C(0));
        if (kRepairCount % 64 != 0)
            common.back() = (UINT64_C(1) << (kRepairCount % 64)) - 1;
        Fnv activity_ledger;
        u64 operation_total = 0;
        u64 active_count_total = 0;
        const u64 target_rank = colex_rank(kOldMask265);

        // Process bounded batches so only THREADS full 5,852,925-bit activity
        // vectors coexist.  The intersection is order-independent, while the
        // audit ledger is consumed in canonical ideal order after every join.
        for (std::size_t begin = 0; begin < ideal.rows.size();
             begin += thread_count) {
            const std::size_t count = std::min<std::size_t>(
                thread_count, ideal.rows.size() - begin);
            std::vector<RowActivity> rows(count);
            std::atomic<std::size_t> next{0};
            std::vector<std::thread> workers;
            for (unsigned thread = 0; thread < thread_count; ++thread) {
                workers.emplace_back([&]() {
                    while (true) {
                        const std::size_t offset = next.fetch_add(1);
                        if (offset >= count) break;
                        rows[offset] = compute_activity(
                            ideal.rows[begin + offset], {});
                    }
                });
            }
            for (std::thread& worker : workers) worker.join();
            for (std::size_t offset = 0; offset < count; ++offset) {
                const RowActivity& row = rows[offset];
                const Pair expected = ideal.rows[begin + offset];
                require(row.pair == expected,
                        "parallel activity row order changed");
                require((row.active[target_rank / 64] &
                         (UINT64_C(1) << (target_rank % 64))) == 0,
                        "old index-265 mask is not strictly inactive");
                for (std::size_t word = 0; word < common.size(); ++word)
                    common[word] &= row.active[word];
                operation_total += row.operations;
                active_count_total += row.active_count;
                activity_ledger.add(row.pair.q);
                activity_ledger.add(row.pair.r);
                activity_ledger.add(row.grid);
                activity_ledger.add(row.cells);
                activity_ledger.add(row.failure_classes);
                activity_ledger.add(row.operations);
                activity_ledger.add(row.active_count);
            }
            std::cout << "PROGRESS " << (begin + count) << '/'
                      << ideal.rows.size() << '\n';
        }

        std::map<unsigned, ResponseClass265> classes;
        Fnv common_ledger;
        Fnv response_ledger;
        u64 common_count = 0;
        u64 nonempty_count = 0;
        u32 repair = (u32{1} << 8) - 1;
        for (u64 rank = 0; rank < kRepairCount; ++rank) {
            if ((common[rank / 64] & (UINT64_C(1) << (rank % 64))) != 0) {
                ++common_count;
                common_ledger.add(repair);
                unsigned response = 0;
                for (unsigned index = 0; index < obligations.bodies.size(); ++index)
                    if ((repair & obligations.bodies[index]) == 0)
                        response |= 1u << index;
                nonempty_count += response != 0;
                response_ledger.add(repair);
                response_ledger.add(response);
                ResponseClass265& cls = classes[response];
                ++cls.count;
                if (cls.count == 1 || repair < cls.least) cls.least = repair;
            }
            if (rank + 1 < kRepairCount) repair = next_combination(repair);
        }
        require(repair == UINT32_C(0x3fc00000),
                "rank-eight enumeration endpoint changed");

        std::vector<unsigned> maximal;
        for (const auto& [pattern, cls] : classes) {
            (void)cls;
            if (pattern == 0) continue;
            bool dominated = false;
            for (const auto& [other, other_cls] : classes) {
                (void)other_cls;
                if (other != pattern && (pattern & ~other) == 0) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) maximal.push_back(pattern);
        }
        const unsigned full = (1u << obligations.bodies.size()) - 1;
        const ExactCover265 cover = exact_cover265(maximal, full);
        const Packing265 packing = largest_integral_packing265(classes);

        std::vector<u32> witnesses;
        unsigned witness_union = 0;
        Fnv witness_ledger;
        for (unsigned pattern : cover.patterns) {
            const u32 witness = classes.at(pattern).least;
            witnesses.push_back(witness);
            witness_union |= pattern;
            witness_ledger.add(witness);
        }
        require(witness_union == full &&
                    std::set<u32>(witnesses.begin(), witnesses.end()).size() ==
                        witnesses.size(),
                "exact-cover reconstruction failed");

        std::vector<u32> rebuilt;
        rebuilt.reserve(420 + witnesses.size());
        for (std::size_t index = 0; index < deck.size(); ++index)
            if (index != kTarget265) rebuilt.push_back(deck[index]);
        rebuilt.insert(rebuilt.end(), witnesses.begin(), witnesses.end());
        require(std::set<u32>(rebuilt.begin(), rebuilt.end()).size() ==
                    rebuilt.size(),
                "rebuilt deck has a duplicate");
        Fnv rebuilt_ledger;
        for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
        const BodyAudit body_audit = scan_body_cover(rebuilt);
        require(body_audit.failures == 0,
                "rebuilt index-265 deck fails body cover");

        std::ofstream quotient(argv[5], std::ios::binary);
        require(static_cast<bool>(quotient), "cannot create quotient ledger");
        quotient << "response_hex,count,least_mask,maximal\n";
        for (const auto& [pattern, cls] : classes) {
            quotient << std::hex << pattern << std::dec << ',' << cls.count
                     << ',' << std::hex << std::setw(8) << std::setfill('0')
                     << cls.least << std::dec << std::setfill(' ') << ','
                     << (std::find(maximal.begin(), maximal.end(), pattern) !=
                         maximal.end() ? 1 : 0) << '\n';
        }
        require(quotient.good(), "quotient ledger write failed");

        std::cout << "INDEX265_EXACT_COMMON_RESPONSE_AUDIT_V1\n"
                  << "IDEAL INDEX 265 ROWS " << ideal.rows.size()
                  << " FNV " << std::hex << ideal.fnv << std::dec
                  << " MAX_ENDPOINT " << ideal.maximum_endpoint << '\n'
                  << "OLD_MASK " << std::hex << std::setw(8)
                  << std::setfill('0') << deck[kTarget265] << std::dec
                  << std::setfill(' ') << " OBLIGATIONS "
                  << obligations.bodies.size() << " FNV " << std::hex
                  << obligations.fnv << " UNION " << obligations.body_union
                  << std::dec << " TARGET_DISJOINT "
                  << obligations.target_disjoint << " RETAINED_CHECKS "
                  << obligations.retained_checks << " BODIES";
        for (u32 body : obligations.bodies)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << body << std::dec
                      << std::setfill(' ');
        std::cout << '\n'
                  << "ACTIVITY ROWS " << ideal.rows.size()
                  << " REPAIRS_PER_ROW " << kRepairCount
                  << " ZETA_OPERATIONS " << operation_total
                  << " ACTIVE_COUNT_SUM " << active_count_total
                  << " ACTIVITY_FNV " << std::hex << activity_ledger.state
                  << std::dec << '\n'
                  << "COMMON_ACTIVE " << common_count << " FNV " << std::hex
                  << common_ledger.state << " RESPONSE_FNV "
                  << response_ledger.state << std::dec << " CLASSES "
                  << classes.size() << " NONEMPTY_MASKS " << nonempty_count
                  << " FULL_RESPONDERS "
                  << (classes.contains(full) ? classes.at(full).count : 0)
                  << '\n';
        for (const auto& [pattern, cls] : classes)
            std::cout << "CLASS " << std::hex << std::setw(2)
                      << std::setfill('0') << pattern << " LEAST "
                      << std::setw(8) << cls.least << std::dec
                      << std::setfill(' ') << " COVER "
                      << std::popcount(pattern) << " COUNT " << cls.count
                      << " MAXIMAL "
                      << (std::find(maximal.begin(), maximal.end(), pattern) !=
                          maximal.end()) << '\n';
        std::cout << "MAXIMAL_ANTICHAIN";
        for (unsigned pattern : maximal)
            std::cout << ' ' << std::hex << pattern << std::dec;
        std::cout << "\nINTEGRAL_PACKING SUPPORT " << std::hex
                  << packing.support << std::dec << " VALUE " << packing.size
                  << " MAX_LOAD 1\n"
                  << "EXACT_COVER MINIMUM " << cover.minimum
                  << " DP_STATES 256 PATTERNS";
        for (unsigned pattern : cover.patterns)
            std::cout << ' ' << std::hex << pattern << std::dec;
        std::cout << " MASKS";
        for (u32 witness : witnesses)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << witness << std::dec
                      << std::setfill(' ');
        std::cout << " FNV " << std::hex << witness_ledger.state << std::dec
                  << '\n'
                  << "REBUILT_DECK " << rebuilt.size() << " FNV "
                  << std::hex << rebuilt_ledger.state << std::dec
                  << " BODY_SCAN " << kBodyCount << " CHECKS "
                  << body_audit.checks << " FAILURES " << body_audit.failures
                  << " BODY_FNV " << std::hex << body_audit.fnv << std::dec
                  << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_POOL_SEPARATE_COMMON_DECK_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS COMPLETE_INDEX265_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEX265_EXACT_ERROR " << error.what() << '\n';
        return 1;
    }
}
