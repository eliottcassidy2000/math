// Target-aware response covers for the four global-signature candidates left
// outside the first 127 greedy common-row union.
//
// For each target, use the inclusion-minimal top-signature scenario containing
// its old inactive signature, intersect exact activity at the selected top
// row(s) with exact activity at the target, and greedily cover the scenario's
// exposed bodies.  The output is an explicit upper bound only.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "../lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace {

constexpr std::size_t kDeckCount = 421;
constexpr u64 kDeckFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kGlobalObligationCount = 401;
constexpr u64 kGlobalObligationFnv = UINT64_C(0xa149cb077a90ef39);
constexpr std::size_t kResponseWords = 7;

struct LiftPair {
    int q = 0;
    int r = 0;
};

struct Target {
    LiftPair pair;
    unsigned scenario = 0;
};

constexpr std::array<LiftPair, 7> kTopPairs = {{
    {220, 644}, {256, 644}, {258, 644}, {294, 644},
    {366, 644}, {416, 644}, {512, 644},
}};
constexpr std::array<Target, 4> kTargets = {{
    {{206, 263}, 2}, {{250, 256}, 10},
    {{256, 394}, 2}, {{256, 400}, 2},
}};

using Signature = std::array<u64, 7>;
using Response = std::array<u64, kResponseWords>;

u64 lift_pair_key(int q, int r) {
    return (static_cast<u64>(static_cast<u32>(q)) << 32) |
           static_cast<u32>(r);
}

std::vector<std::string> lift_split(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    return fields;
}

std::vector<u32> lift_read_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    FnvLocal ledger;
    std::string token;
    while (input >> token) {
        const u64 wide = std::stoull(token, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "bad deck token");
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

std::map<u64, Signature> lift_read_signatures(
        const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::map<u64, Signature> signatures;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = lift_split(line);
        require(fields.size() == 10, "malformed signature row");
        Signature signature{};
        unsigned weight = 0;
        for (std::size_t word = 0; word < signature.size(); ++word) {
            signature[word] = std::stoull(fields[word + 3], nullptr, 16);
            weight += std::popcount(signature[word]);
        }
        require(weight == std::stoul(fields[2]),
                "signature weight changed");
        require(signatures.emplace(
                    lift_pair_key(std::stoi(fields[0]), std::stoi(fields[1])),
                    signature).second,
                "duplicate signature row");
    }
    require(signatures.size() == 24223, "signature row count changed");
    return signatures;
}

Signature lift_union(Signature left, const Signature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        left[word] |= right[word];
    return left;
}

bool lift_subset(const Signature& left, const Signature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        if ((left[word] & ~right[word]) != 0) return false;
    return true;
}

std::vector<unsigned> lift_indices(const Signature& signature) {
    std::vector<unsigned> indices;
    for (unsigned index = 0; index < kDeckCount; ++index)
        if ((signature[index / 64] >> (index % 64)) & 1)
            indices.push_back(index);
    return indices;
}

unsigned lift_gain(const Response& response, const Response& uncovered) {
    unsigned answer = 0;
    for (std::size_t word = 0; word < response.size(); ++word)
        answer += std::popcount(response[word] & uncovered[word]);
    return answer;
}

bool lift_empty(const Response& response) {
    for (u64 word : response)
        if (word != 0) return false;
    return true;
}

u64 lift_mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 4,
                "usage: targeted-lift JOINT421 SIGNATURES OUTDIR");
        init_choose8_local();
        const std::vector<u32> deck = lift_read_deck(argv[1]);
        const std::map<u64, Signature> signatures =
            lift_read_signatures(argv[2]);
        const std::filesystem::path output_root(argv[3]);
        std::filesystem::create_directories(output_root);

        std::array<Signature, 7> top_signatures{};
        Signature global_signature{};
        for (std::size_t row = 0; row < kTopPairs.size(); ++row) {
            const auto found = signatures.find(
                lift_pair_key(kTopPairs[row].q, kTopPairs[row].r));
            require(found != signatures.end(), "top signature missing");
            top_signatures[row] = found->second;
            global_signature = lift_union(global_signature, found->second);
        }
        const std::vector<unsigned> global_indices =
            lift_indices(global_signature);
        require(global_indices.size() == 27,
                "global deletion union changed");
        std::array<int, kDeckCount> global_position{};
        global_position.fill(-1);
        for (std::size_t local = 0; local < global_indices.size(); ++local)
            global_position[global_indices[local]] = static_cast<int>(local);

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
        require(global_bodies.size() == kGlobalObligationCount &&
                    global_body_ledger.state == kGlobalObligationFnv,
                "global obligation family changed");

        std::vector<Response> responses(EXPECTED_REPAIRS);
        for (std::size_t index = 0; index < global_bodies.size(); ++index)
            enumerate_disjoint_repairs(global_bodies[index],
                [&](u32, u64 rank) {
                    responses[rank][index / 64] |=
                        UINT64_C(1) << (index % 64);
                });
        std::vector<u32> mask_values(EXPECTED_REPAIRS);
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank)
            mask_values[rank] = unrank_colex8(rank);

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        const ActiveUniverse active256 =
            build_active_universe(cells, 256, 644);
        const ActiveUniverse active294 =
            build_active_universe(cells, 294, 644);
        require(active256.count == UINT64_C(1465388) &&
                    active256.fnv == UINT64_C(0x7ee496dd5e3b67c6) &&
                    active294.count == UINT64_C(3290605) &&
                    active294.fnv == UINT64_C(0x73c1d8e05b9d6db8),
                "top activity controls changed");

        std::cout << "THM4283_TARGETED_RESPONSE_LIFT_V1\n"
                  << "GLOBAL_OBLIGATIONS " << global_bodies.size()
                  << " FNV " << std::hex << global_body_ledger.state
                  << std::dec << '\n';
        FnvLocal target_ledger;
        for (const Target target : kTargets) {
            const auto target_signature_found = signatures.find(
                lift_pair_key(target.pair.q, target.pair.r));
            require(target_signature_found != signatures.end(),
                    "target signature missing");
            Signature deleted{};
            for (unsigned row = 0; row < 7; ++row)
                if ((target.scenario >> row) & 1)
                    deleted = lift_union(deleted, top_signatures[row]);
            require(lift_subset(target_signature_found->second, deleted),
                    "target signature not contained in scenario");
            for (unsigned proper = 1; proper < target.scenario; ++proper) {
                if ((proper & target.scenario) != proper) continue;
                Signature proper_deleted{};
                for (unsigned row = 0; row < 7; ++row)
                    if ((proper >> row) & 1)
                        proper_deleted = lift_union(
                            proper_deleted, top_signatures[row]);
                require(!lift_subset(target_signature_found->second,
                                     proper_deleted),
                        "target scenario is not inclusion-minimal");
            }
            u32 deleted_local = 0;
            for (unsigned index : lift_indices(deleted)) {
                require(global_position[index] >= 0,
                        "target scenario escaped global union");
                deleted_local |= u32{1} << global_position[index];
            }
            Response obligations{};
            std::size_t obligation_count = 0;
            FnvLocal obligation_ledger;
            for (std::size_t index = 0; index < global_bodies.size(); ++index)
                if ((global_responses[index] & ~deleted_local) == 0) {
                    obligations[index / 64] |=
                        UINT64_C(1) << (index % 64);
                    ++obligation_count;
                    obligation_ledger.add(global_bodies[index]);
                }

            const ActiveUniverse target_active = build_active_universe(
                cells, target.pair.q, target.pair.r);
            std::vector<u32> allowed_ranks;
            allowed_ranks.reserve(EXPECTED_REPAIRS / 4);
            Response response_union{};
            FnvLocal allowed_ledger;
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
                if (!active256.active[rank] ||
                    (target.scenario == 10 && !active294.active[rank]) ||
                    !target_active.active[rank] ||
                    lift_gain(responses[rank], obligations) == 0)
                    continue;
                allowed_ranks.push_back(static_cast<u32>(rank));
                allowed_ledger.add(mask_values[rank]);
                for (std::size_t word = 0; word < response_union.size(); ++word)
                    response_union[word] |=
                        responses[rank][word] & obligations[word];
            }
            require(response_union == obligations,
                    "target-aware active response family is infeasible");

            Response uncovered = obligations;
            std::vector<u32> witnesses;
            std::vector<unsigned> gains;
            while (!lift_empty(uncovered)) {
                unsigned best_gain = 0;
                u32 best_mask = 0;
                u32 best_rank = 0;
                for (u32 rank : allowed_ranks) {
                    const unsigned gain = lift_gain(responses[rank], uncovered);
                    const u32 mask = mask_values[rank];
                    if (gain > best_gain ||
                        (gain == best_gain && gain != 0 &&
                         (best_mask == 0 || mask < best_mask))) {
                        best_gain = gain;
                        best_mask = mask;
                        best_rank = rank;
                    }
                }
                require(best_gain != 0, "target-aware greedy cover stuck");
                witnesses.push_back(best_mask);
                gains.push_back(best_gain);
                for (std::size_t word = 0; word < uncovered.size(); ++word)
                    uncovered[word] &= ~responses[best_rank][word];
            }
            // Reverse deletion makes the explicit greedy cover irredundant.
            for (std::size_t position = witnesses.size(); position-- > 0;) {
                Response covered{};
                for (std::size_t index = 0; index < witnesses.size(); ++index) {
                    if (index == position) continue;
                    const u64 rank = colex_rank8_local(witnesses[index]);
                    for (std::size_t word = 0; word < covered.size(); ++word)
                        covered[word] |= responses[rank][word] & obligations[word];
                }
                if (covered == obligations) witnesses.erase(
                    witnesses.begin() + position);
            }
            require(!witnesses.empty(), "target witness cover vanished");
            Response covered{};
            for (u32 witness : witnesses) {
                const u64 rank = colex_rank8_local(witness);
                require(active256.active[rank] && target_active.active[rank] &&
                            (target.scenario != 10 || active294.active[rank]),
                        "target witness activity changed");
                for (std::size_t word = 0; word < covered.size(); ++word)
                    covered[word] |= responses[rank][word] & obligations[word];
            }
            require(covered == obligations, "target witness cover changed");

            FnvLocal deck_ledger;
            for (unsigned index = 0; index < deck.size(); ++index)
                if (((deleted[index / 64] >> (index % 64)) & 1) == 0)
                    deck_ledger.add(deck[index]);
            for (u32 witness : witnesses) deck_ledger.add(witness);
            const std::size_t deck_size = deck.size() -
                lift_indices(deleted).size() + witnesses.size();
            std::ostringstream name;
            name << "target_" << target.pair.q << '_' << target.pair.r
                 << "_scenario_" << std::setw(3) << std::setfill('0')
                 << target.scenario << "_witnesses.txt";
            std::ofstream output(output_root / name.str());
            require(static_cast<bool>(output), "cannot create target witnesses");
            for (u32 witness : witnesses)
                output << std::hex << std::setw(8) << std::setfill('0')
                       << witness << std::dec << std::setfill(' ') << '\n';
            require(output.good(), "target witness write failed");

            target_ledger.add(target.pair.q);
            target_ledger.add(target.pair.r);
            target_ledger.add(target.scenario);
            target_ledger.add(obligation_count);
            target_ledger.add(obligation_ledger.state);
            target_ledger.add(allowed_ranks.size());
            target_ledger.add(allowed_ledger.state);
            target_ledger.add(witnesses.size());
            target_ledger.add(lift_mask_fnv(witnesses));
            target_ledger.add(deck_size);
            target_ledger.add(deck_ledger.state);
            std::cout << "TARGET " << target.pair.q << ',' << target.pair.r
                      << " SCENARIO " << target.scenario << " DELETED "
                      << lift_indices(deleted).size() << " OBLIGATIONS "
                      << obligation_count << " OBLIGATION_FNV " << std::hex
                      << obligation_ledger.state << std::dec << " ALLOWED "
                      << allowed_ranks.size() << " ALLOWED_FNV " << std::hex
                      << allowed_ledger.state << std::dec << " WITNESSES "
                      << witnesses.size() << " WITNESS_FNV " << std::hex
                      << lift_mask_fnv(witnesses) << std::dec << " DECK "
                      << deck_size << " DECK_FNV " << std::hex
                      << deck_ledger.state << std::dec << '\n';
        }
        std::cout << "TARGET_LEDGER_FNV " << std::hex << target_ledger.state
                  << std::dec << '\n'
                  << "SCOPE EXPLICIT_TARGET_AWARE_UPPER_BOUNDS_NO_MINIMALITY\n"
                  << "VERDICT PASS FOUR_TARGET_RESPONSE_LIFTS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_TARGETED_LIFT_ERROR " << error.what() << '\n';
        return 1;
    }
}
