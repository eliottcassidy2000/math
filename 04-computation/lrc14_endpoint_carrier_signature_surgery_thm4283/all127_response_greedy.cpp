// Exact 127-scenario signature-surgery response atlas for THM-4283.
//
// The seven endpoint-644 inactive signatures generate 127 nonempty unions.
// For each union this program derives the exact exposed-body family, intersects
// activity across the selected endpoint rows, and constructs a deterministic
// greedy active-response cover.  Greedy witnesses are explicit upper bounds;
// no minimum claim is made.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "../lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <atomic>
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <thread>

namespace {

constexpr std::size_t kDeckCount = 421;
constexpr u64 kDeckFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kResidualCount = 23373;
constexpr u64 kResidualFnv = UINT64_C(0xc6ab0ae49ee32273);
constexpr std::size_t kGlobalObligationCount = 401;
constexpr u64 kGlobalObligationFnv = UINT64_C(0xa149cb077a90ef39);
constexpr std::size_t kResponseWords = 7;

struct Pair {
    int q = 0;
    int r = 0;
};

constexpr std::array<Pair, 7> kTopPairs = {{
    {220, 644}, {256, 644}, {258, 644}, {294, 644},
    {366, 644}, {416, 644}, {512, 644},
}};

constexpr std::array<u64, 7> kTopActiveCounts = {{
    UINT64_C(3096611), UINT64_C(1465388), UINT64_C(3957432),
    UINT64_C(3290605), UINT64_C(2383426), UINT64_C(3098275),
    UINT64_C(2271751),
}};

constexpr std::array<u64, 7> kTopActiveFnvs = {{
    UINT64_C(0xe25cd7111fdbe71d), UINT64_C(0x7ee496dd5e3b67c6),
    UINT64_C(0x17fe914235c41d6a), UINT64_C(0x73c1d8e05b9d6db8),
    UINT64_C(0x3275a21c0f85fc9d), UINT64_C(0x790be8ea06091783),
    UINT64_C(0x9f37456f9594fb36),
}};

using Signature = std::array<u64, 7>;
using Response = std::array<u64, kResponseWords>;

u64 pair_key(int q, int r) {
    return (static_cast<u64>(static_cast<u32>(q)) << 32) |
           static_cast<u32>(r);
}

std::vector<u32> read_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    FnvLocal ledger;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "bad deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "deck rank or distinctness changed");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(deck.size() == kDeckCount && ledger.state == kDeckFnv,
            "joint deck identity changed");
    return deck;
}

std::vector<Pair> read_residual(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open final residual");
    std::vector<Pair> rows;
    FnvLocal ledger;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "bad residual row");
        std::size_t used_q = 0;
        std::size_t used_r = 0;
        Pair pair{std::stoi(line.substr(0, comma), &used_q),
                  std::stoi(line.substr(comma + 1), &used_r)};
        require(used_q == comma && used_r == line.size() - comma - 1 &&
                    pair.q > 0 && pair.q < pair.r,
                "invalid residual pair");
        if (!rows.empty())
            require(pair_key(rows.back().q, rows.back().r) <
                        pair_key(pair.q, pair.r),
                    "residual order changed");
        rows.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(rows.size() == kResidualCount && ledger.state == kResidualFnv,
            "final residual identity changed");
    return rows;
}

std::map<u64, Signature> read_signatures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::map<u64, Signature> signatures;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::stringstream stream(line);
        std::array<std::string, 10> fields{};
        for (std::size_t index = 0; index < fields.size(); ++index) {
            require(static_cast<bool>(std::getline(stream, fields[index],
                                                   index + 1 == fields.size()
                                                       ? '\n' : ',')),
                    "short signature row");
        }
        require(stream.peek() == std::char_traits<char>::eof(),
                "extra signature field");
        const int q = std::stoi(fields[0]);
        const int r = std::stoi(fields[1]);
        const unsigned stated = static_cast<unsigned>(std::stoul(fields[2]));
        Signature signature{};
        unsigned actual = 0;
        for (std::size_t word = 0; word < signature.size(); ++word) {
            std::size_t used = 0;
            signature[word] = std::stoull(fields[word + 3], &used, 16);
            require(used == fields[word + 3].size(), "bad signature word");
            actual += std::popcount(signature[word]);
        }
        require((signature[6] >> 37) == 0 && actual == stated,
                "signature weight or padding changed");
        require(signatures.emplace(pair_key(q, r), signature).second,
                "duplicate signature row");
    }
    require(signatures.size() == 24223, "signature row count changed");
    return signatures;
}

std::vector<unsigned> signature_indices(const Signature& signature) {
    std::vector<unsigned> indices;
    for (unsigned index = 0; index < kDeckCount; ++index)
        if ((signature[index / 64] >> (index % 64)) & 1)
            indices.push_back(index);
    return indices;
}

bool signature_subset(const Signature& left, const Signature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        if ((left[word] & ~right[word]) != 0) return false;
    return true;
}

Signature signature_union(Signature left, const Signature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        left[word] |= right[word];
    return left;
}

u64 response_fnv(const std::vector<u32>& witnesses) {
    FnvLocal ledger;
    for (u32 mask : witnesses) ledger.add(mask);
    return ledger.state;
}

unsigned response_gain(const Response& response, const Response& uncovered) {
    unsigned gain = 0;
    for (std::size_t word = 0; word < response.size(); ++word)
        gain += std::popcount(response[word] & uncovered[word]);
    return gain;
}

bool response_empty(const Response& response) {
    for (u64 word : response)
        if (word != 0) return false;
    return true;
}

struct ScenarioResult {
    Signature deleted{};
    u32 deleted_local = 0;
    Response obligations{};
    std::size_t obligation_count = 0;
    u64 obligation_fnv = 0;
    std::vector<u32> witnesses;
    std::vector<unsigned> gains;
    u64 witness_fnv = 0;
    u64 deck_fnv = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 6,
                "usage: all127-response JOINT421 SIGNATURES FINAL23373 "
                "OUTDIR SUMMARY_CSV");
        init_choose8_local();
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::map<u64, Signature> signatures = read_signatures(argv[2]);
        const std::vector<Pair> residual = read_residual(argv[3]);
        const std::filesystem::path output_root(argv[4]);
        std::filesystem::create_directories(output_root);

        std::array<Signature, 7> top_signatures{};
        Signature global_signature{};
        for (std::size_t row = 0; row < kTopPairs.size(); ++row) {
            const auto found = signatures.find(
                pair_key(kTopPairs[row].q, kTopPairs[row].r));
            require(found != signatures.end(), "top signature missing");
            top_signatures[row] = found->second;
            global_signature = signature_union(global_signature, found->second);
        }
        const std::array<std::vector<unsigned>, 7> expected_top = {{
            {25},
            {9,29,32,75,137,139,150,159,174,205,218,309,333,347,358,
             374,394,399,405,416,417},
            {412}, {173}, {396}, {236}, {107,374},
        }};
        for (std::size_t row = 0; row < kTopPairs.size(); ++row)
            require(signature_indices(top_signatures[row]) == expected_top[row],
                    "top inactive signature changed");
        const std::vector<unsigned> global_indices =
            signature_indices(global_signature);
        require(global_indices.size() == 27, "global deletion size changed");
        std::array<int, kDeckCount> global_position{};
        global_position.fill(-1);
        for (std::size_t local = 0; local < global_indices.size(); ++local)
            global_position[global_indices[local]] = static_cast<int>(local);

        std::array<Signature, 128> deleted_signatures{};
        std::array<u32, 128> deleted_local{};
        for (unsigned scenario = 1; scenario < 128; ++scenario) {
            Signature signature{};
            for (unsigned row = 0; row < 7; ++row)
                if ((scenario >> row) & 1)
                    signature = signature_union(signature, top_signatures[row]);
            deleted_signatures[scenario] = signature;
            u32 local_mask = 0;
            for (unsigned index : signature_indices(signature)) {
                require(global_position[index] >= 0,
                        "scenario index outside global union");
                local_mask |= u32{1} << global_position[index];
            }
            deleted_local[scenario] = local_mask;
        }

        std::vector<u32> global_bodies;
        std::vector<u32> global_responses;
        FnvLocal global_body_ledger;
        u32 body = (u32{1} << 9) - 1;
        for (u64 ordinal = 0; ordinal < EXPECTED_BODIES; ++ordinal) {
            bool retained_hit = false;
            for (unsigned index = 0; index < deck.size(); ++index) {
                if (global_position[index] >= 0) continue;
                if ((body & deck[index]) == 0) {
                    retained_hit = true;
                    break;
                }
            }
            if (!retained_hit) {
                u32 response = 0;
                for (std::size_t local = 0; local < global_indices.size();
                     ++local)
                    if ((body & deck[global_indices[local]]) == 0)
                        response |= u32{1} << local;
                require(response != 0, "original joint deck body failure");
                global_bodies.push_back(body);
                global_responses.push_back(response);
                global_body_ledger.add(body);
            }
            if (ordinal + 1 < EXPECTED_BODIES)
                body = next_combination(body);
        }
        require(body == UINT32_C(0x3fe00000), "body endpoint changed");
        require(global_bodies.size() == kGlobalObligationCount &&
                    global_body_ledger.state == kGlobalObligationFnv,
                "global obligation family changed");

        std::array<ScenarioResult, 128> results{};
        for (unsigned scenario = 1; scenario < 128; ++scenario) {
            ScenarioResult& result = results[scenario];
            result.deleted = deleted_signatures[scenario];
            result.deleted_local = deleted_local[scenario];
            FnvLocal ledger;
            for (std::size_t index = 0; index < global_bodies.size(); ++index) {
                if ((global_responses[index] & ~result.deleted_local) != 0)
                    continue;
                result.obligations[index / 64] |=
                    UINT64_C(1) << (index % 64);
                ++result.obligation_count;
                ledger.add(global_bodies[index]);
            }
            result.obligation_fnv = ledger.state;
            require(result.obligation_count != 0,
                    "empty nontrivial scenario obligation family");
        }

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        std::vector<unsigned char> activity_pattern(EXPECTED_REPAIRS, 0);
        for (unsigned row = 0; row < kTopPairs.size(); ++row) {
            ActiveUniverse active = build_active_universe(
                cells, kTopPairs[row].q, kTopPairs[row].r);
            require(active.count == kTopActiveCounts[row] &&
                        active.fnv == kTopActiveFnvs[row],
                    "top active universe changed");
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank)
                if (active.active[rank])
                    activity_pattern[rank] |=
                        static_cast<unsigned char>(1u << row);
        }

        std::vector<Response> responses(EXPECTED_REPAIRS);
        for (std::size_t index = 0; index < global_bodies.size(); ++index) {
            enumerate_disjoint_repairs(global_bodies[index],
                [&](u32, u64 rank) {
                    responses[rank][index / 64] |=
                        UINT64_C(1) << (index % 64);
                });
        }
        std::vector<u32> mask_values(EXPECTED_REPAIRS);
        std::array<std::vector<u32>, 128> pattern_ranks;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const u32 mask = unrank_colex8(rank);
            mask_values[rank] = mask;
            if (!response_empty(responses[rank]) &&
                activity_pattern[rank] != 0)
                pattern_ranks[activity_pattern[rank]].push_back(
                    static_cast<u32>(rank));
        }

        std::atomic<unsigned> next_scenario{1};
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned worker_count =
            std::max(1u, std::min(4u, hardware ? hardware : 1u));
        std::vector<std::thread> workers;
        for (unsigned worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const unsigned scenario = next_scenario.fetch_add(1);
                    if (scenario >= 128) break;
                    ScenarioResult& result = results[scenario];
                    Response uncovered = result.obligations;
                    while (!response_empty(uncovered)) {
                        unsigned best_gain = 0;
                        u32 best_mask = 0;
                        u32 best_rank = 0;
                        for (unsigned pattern = 1; pattern < 128; ++pattern) {
                            if ((pattern & scenario) != scenario) continue;
                            for (u32 rank : pattern_ranks[pattern]) {
                                const unsigned gain =
                                    response_gain(responses[rank], uncovered);
                                const u32 mask = mask_values[rank];
                                if (gain > best_gain ||
                                    (gain == best_gain && gain != 0 &&
                                     (best_mask == 0 || mask < best_mask))) {
                                    best_gain = gain;
                                    best_mask = mask;
                                    best_rank = rank;
                                }
                            }
                        }
                        require(best_gain != 0,
                                "greedy active cover became stuck");
                        result.witnesses.push_back(best_mask);
                        result.gains.push_back(best_gain);
                        for (std::size_t word = 0; word < uncovered.size();
                             ++word)
                            uncovered[word] &= ~responses[best_rank][word];
                    }
                    result.witness_fnv = response_fnv(result.witnesses);
                    FnvLocal deck_ledger;
                    for (unsigned index = 0; index < deck.size(); ++index)
                        if (((result.deleted[index / 64] >> (index % 64)) & 1)
                            == 0)
                            deck_ledger.add(deck[index]);
                    for (u32 mask : result.witnesses)
                        deck_ledger.add(mask);
                    result.deck_fnv = deck_ledger.state;
                }
            });
        }
        for (std::thread& worker : workers) worker.join();

        const std::array<std::pair<unsigned, u64>, 8> recovered_controls = {{
            {1, UINT64_C(0xbba0d34db7d3b89a)},
            {2, UINT64_C(0x3257b4f3c1de01bb)},
            {4, UINT64_C(0x7da398510a109775)},
            {8, UINT64_C(0x6e355f8a61f632c6)},
            {16, UINT64_C(0x980d14a50f5cce6)},
            {32, UINT64_C(0x66efb34937743a7b)},
            {64, UINT64_C(0xab70d4857912a764)},
            {127, UINT64_C(0x2872403c9d48f141)},
        }};
        for (const auto [scenario, expected] : recovered_controls)
            require(results[scenario].witness_fnv == expected,
                    "recovered scenario witness control changed");

        std::ofstream summary(argv[5]);
        require(static_cast<bool>(summary), "cannot create scenario summary");
        summary << "scenario,deleted,obligations,obligation_fnv,witnesses,"
                   "witness_fnv,deck_size,deck_fnv\n";
        std::cout << "THM4283_ALL127_RESPONSE_GREEDY_V1\n"
                  << "DECK " << deck.size() << " FNV " << std::hex
                  << kDeckFnv << std::dec << " GLOBAL_DELETED "
                  << global_indices.size() << " GLOBAL_OBLIGATIONS "
                  << global_bodies.size() << " BODY_FNV " << std::hex
                  << global_body_ledger.state << std::dec << '\n';
        FnvLocal scenario_ledger;
        for (unsigned scenario = 1; scenario < 128; ++scenario) {
            const ScenarioResult& result = results[scenario];
            const std::vector<unsigned> deleted =
                signature_indices(result.deleted);
            const std::size_t deck_size =
                deck.size() - deleted.size() + result.witnesses.size();
            summary << scenario << ',' << deleted.size() << ','
                    << result.obligation_count << ',' << std::hex
                    << result.obligation_fnv << std::dec << ','
                    << result.witnesses.size() << ',' << std::hex
                    << result.witness_fnv << std::dec << ',' << deck_size
                    << ',' << std::hex << result.deck_fnv << std::dec << '\n';
            std::ostringstream name;
            name << "scenario_" << std::setw(3) << std::setfill('0')
                 << scenario << "_witnesses.txt";
            std::ofstream witness_file(output_root / name.str());
            require(static_cast<bool>(witness_file),
                    "cannot create witness file");
            for (u32 mask : result.witnesses)
                witness_file << std::hex << std::setw(8) << std::setfill('0')
                             << mask << std::dec << std::setfill(' ') << '\n';
            require(witness_file.good(), "witness write failed");
            scenario_ledger.add(scenario);
            scenario_ledger.add(deleted.size());
            scenario_ledger.add(result.obligation_count);
            scenario_ledger.add(result.obligation_fnv);
            scenario_ledger.add(result.witnesses.size());
            scenario_ledger.add(result.witness_fnv);
            scenario_ledger.add(deck_size);
            scenario_ledger.add(result.deck_fnv);
            std::cout << "SCENARIO " << scenario << " DELETED "
                      << deleted.size() << " OBLIGATIONS "
                      << result.obligation_count << " OBLIGATION_FNV "
                      << std::hex << result.obligation_fnv << std::dec
                      << " WITNESSES " << result.witnesses.size()
                      << " WITNESS_FNV " << std::hex << result.witness_fnv
                      << " DECK " << std::dec << deck_size << " DECK_FNV "
                      << std::hex << result.deck_fnv << std::dec << '\n';
        }
        require(summary.good(), "summary write failed");
        std::cout << "SCENARIO_LEDGER_FNV " << std::hex
                  << scenario_ledger.state << std::dec << '\n'
                  << "VERDICT PASS EXACT_127_ACTIVE_GREEDY_COVERS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_ALL127_ERROR " << error.what() << '\n';
        return 1;
    }
}
