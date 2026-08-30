// Exact common-row synthesis for the 127 endpoint-644 signature surgeries.
//
// The companion greedy program proves body coverage for every nonempty union
// of the seven top inactive signatures.  This audit rederives the exposed
// bodies, checks every serialized witness deck, evaluates its common rows by
// exact cocycle arithmetic, and extracts a compact scenario subfamily with
// the same common-row union.  Two explicitly named singleton variants are
// included because they trade lexicographic witness choice for a larger
// common-row family.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "../lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

constexpr std::size_t kDeckCountAudit = 421;
constexpr u64 kDeckFnvAudit = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kResidualCountAudit = 23373;
constexpr u64 kResidualFnvAudit = UINT64_C(0xc6ab0ae49ee32273);
constexpr std::size_t kSignatureCountAudit = 24223;
constexpr std::size_t kGlobalObligationCountAudit = 401;
constexpr u64 kGlobalObligationFnvAudit = UINT64_C(0xa149cb077a90ef39);
constexpr u64 kScenarioSummaryFnvAudit = UINT64_C(0xc3518d83d0178da0);

struct AuditPair {
    int q = 0;
    int r = 0;
};

constexpr std::array<AuditPair, 7> kTopPairsAudit = {{
    {220, 644}, {256, 644}, {258, 644}, {294, 644},
    {366, 644}, {416, 644}, {512, 644},
}};

using AuditSignature = std::array<u64, 7>;

struct SummaryRow {
    unsigned scenario = 0;
    std::size_t deleted = 0;
    std::size_t obligations = 0;
    u64 obligation_fnv = 0;
    std::size_t witnesses = 0;
    u64 witness_fnv = 0;
    std::size_t deck_size = 0;
    u64 deck_fnv = 0;
};

struct ScenarioAudit {
    unsigned id = 0;
    std::string label;
    AuditSignature deleted{};
    u32 deleted_local = 0;
    std::vector<u32> witnesses;
    std::size_t obligation_count = 0;
    u64 obligation_fnv = 0;
    std::size_t deck_size = 0;
    u64 deck_fnv = 0;
    std::size_t candidate_count = 0;
    std::vector<std::size_t> common_rows;
    u64 common_fnv = 0;
};

struct TargetVariantDefinition {
    unsigned id = 0;
    const char* label = nullptr;
    unsigned scenario = 0;
    AuditPair target;
    const char* filename = nullptr;
};

constexpr std::array<TargetVariantDefinition, 4> kTargetVariants = {{
    {130, "target_206_263_s002", 2, {206, 263},
     "target_206_263_scenario_002_witnesses.txt"},
    {131, "target_250_256_s010", 10, {250, 256},
     "target_250_256_scenario_010_witnesses.txt"},
    {132, "target_256_394_s002", 2, {256, 394},
     "target_256_394_scenario_002_witnesses.txt"},
    {133, "target_256_400_s002", 2, {256, 400},
     "target_256_400_scenario_002_witnesses.txt"},
}};

u64 audit_pair_key(int q, int r) {
    return (static_cast<u64>(static_cast<u32>(q)) << 32) |
           static_cast<u32>(r);
}

std::vector<std::string> audit_split(const std::string& line,
                                     char delimiter) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, delimiter)) fields.push_back(field);
    return fields;
}

std::vector<u32> audit_read_deck(const std::filesystem::path& path) {
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
    require(deck.size() == kDeckCountAudit &&
                ledger.state == kDeckFnvAudit,
            "joint deck identity changed");
    return deck;
}

std::vector<AuditPair> audit_read_residual(
        const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open final residual");
    std::vector<AuditPair> rows;
    FnvLocal ledger;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = audit_split(line, ',');
        require(fields.size() == 2, "bad residual row");
        AuditPair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        require(pair.q > 0 && pair.q < pair.r, "invalid residual pair");
        if (!rows.empty())
            require(audit_pair_key(rows.back().q, rows.back().r) <
                        audit_pair_key(pair.q, pair.r),
                    "residual order changed");
        rows.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(rows.size() == kResidualCountAudit &&
                ledger.state == kResidualFnvAudit,
            "final residual identity changed");
    return rows;
}

std::map<u64, AuditSignature> audit_read_signatures(
        const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::map<u64, AuditSignature> signatures;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = audit_split(line, ',');
        require(fields.size() == 10, "bad signature row");
        const int q = std::stoi(fields[0]);
        const int r = std::stoi(fields[1]);
        const unsigned stated = static_cast<unsigned>(std::stoul(fields[2]));
        AuditSignature signature{};
        unsigned actual = 0;
        for (std::size_t word = 0; word < signature.size(); ++word) {
            std::size_t used = 0;
            signature[word] = std::stoull(fields[word + 3], &used, 16);
            require(used == fields[word + 3].size(),
                    "bad signature word");
            actual += std::popcount(signature[word]);
        }
        require((signature[6] >> 37) == 0 && actual == stated,
                "signature weight or padding changed");
        require(signatures.emplace(audit_pair_key(q, r), signature).second,
                "duplicate signature row");
    }
    require(signatures.size() == kSignatureCountAudit,
            "signature row count changed");
    return signatures;
}

std::array<SummaryRow, 128> audit_read_summary(
        const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open scenario summary");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "scenario,deleted,obligations,obligation_fnv,"
                        "witnesses,witness_fnv,deck_size,deck_fnv",
            "scenario summary header changed");
    std::array<SummaryRow, 128> rows{};
    FnvLocal ledger;
    unsigned expected_scenario = 1;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = audit_split(line, ',');
        require(fields.size() == 8, "bad scenario summary row");
        SummaryRow row;
        row.scenario = static_cast<unsigned>(std::stoul(fields[0]));
        row.deleted = std::stoull(fields[1]);
        row.obligations = std::stoull(fields[2]);
        row.obligation_fnv = std::stoull(fields[3], nullptr, 16);
        row.witnesses = std::stoull(fields[4]);
        row.witness_fnv = std::stoull(fields[5], nullptr, 16);
        row.deck_size = std::stoull(fields[6]);
        row.deck_fnv = std::stoull(fields[7], nullptr, 16);
        require(row.scenario == expected_scenario,
                "scenario summary order changed");
        rows[row.scenario] = row;
        ledger.add(row.scenario);
        ledger.add(row.deleted);
        ledger.add(row.obligations);
        ledger.add(row.obligation_fnv);
        ledger.add(row.witnesses);
        ledger.add(row.witness_fnv);
        ledger.add(row.deck_size);
        ledger.add(row.deck_fnv);
        ++expected_scenario;
    }
    require(expected_scenario == 128 &&
                ledger.state == kScenarioSummaryFnvAudit,
            "scenario summary identity changed");
    return rows;
}

std::vector<u32> audit_read_witnesses(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open scenario witness file");
    std::vector<u32> witnesses;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "bad witness token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "witness rank or distinctness changed");
        witnesses.push_back(mask);
    }
    require(!witnesses.empty(), "empty scenario witness file");
    return witnesses;
}

std::vector<unsigned> audit_signature_indices(
        const AuditSignature& signature) {
    std::vector<unsigned> indices;
    for (unsigned index = 0; index < kDeckCountAudit; ++index)
        if ((signature[index / 64] >> (index % 64)) & 1)
            indices.push_back(index);
    return indices;
}

bool audit_signature_subset(const AuditSignature& left,
                            const AuditSignature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        if ((left[word] & ~right[word]) != 0) return false;
    return true;
}

AuditSignature audit_signature_union(AuditSignature left,
                                     const AuditSignature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        left[word] |= right[word];
    return left;
}

u64 audit_mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

u64 audit_family_fnv(const std::vector<std::size_t>& indices,
                     const std::vector<AuditPair>& residual) {
    FnvLocal ledger;
    for (std::size_t index : indices) {
        ledger.add(residual[index].q);
        ledger.add(residual[index].r);
    }
    return ledger.state;
}

std::vector<std::size_t> audit_union_indices(
        const std::vector<ScenarioAudit>& scenarios,
        const std::vector<std::size_t>& selected,
        std::size_t residual_size) {
    std::vector<unsigned char> present(residual_size, 0);
    for (std::size_t scenario_index : selected)
        for (std::size_t row : scenarios[scenario_index].common_rows)
            present[row] = 1;
    std::vector<std::size_t> answer;
    for (std::size_t row = 0; row < present.size(); ++row)
        if (present[row]) answer.push_back(row);
    return answer;
}

void audit_write_family(const std::filesystem::path& path,
                        const std::vector<std::size_t>& indices,
                        const std::vector<AuditPair>& residual) {
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create common-family ledger");
    for (std::size_t index : indices)
        output << residual[index].q << ',' << residual[index].r << '\n';
    require(output.good(), "common-family ledger write failed");
}

i128 audit_mask_margin(const AtomData& atoms, const PrimitivePair& primitive,
                       i64 scale, u32 mask) {
    i128 mass = 0;
    for (const auto& [atom, value] : atoms.mass)
        if ((atom & ~mask) == 0) mass += value;
    const i128 denominator =
        static_cast<i128>(primitive.grid) * scale * COMMON;
    return static_cast<i128>(63) * mass -
           static_cast<i128>(4) * denominator;
}

std::string audit_scenario_filename(unsigned scenario) {
    std::ostringstream name;
    name << "scenario_" << std::setw(3) << std::setfill('0') << scenario
         << "_witnesses.txt";
    return name.str();
}

std::string audit_common_filename(const ScenarioAudit& scenario) {
    std::ostringstream name;
    name << "scenario_" << std::setw(3) << std::setfill('0') << scenario.id;
    if (scenario.id > 127) name << '_' << scenario.label;
    name << "_common.csv";
    return name.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 9,
                "usage: common-audit JOINT421 SIGNATURES FINAL23373 "
                "WITNESS_DIR SUMMARY_CSV TARGET_WITNESS_DIR OUTDIR "
                "BODY_LEDGER");
        init_choose8_local();
        const std::vector<u32> deck = audit_read_deck(argv[1]);
        const std::map<u64, AuditSignature> signatures =
            audit_read_signatures(argv[2]);
        const std::vector<AuditPair> residual = audit_read_residual(argv[3]);
        const std::filesystem::path witness_root(argv[4]);
        const std::array<SummaryRow, 128> summary =
            audit_read_summary(argv[5]);
        const std::filesystem::path target_witness_root(argv[6]);
        const std::filesystem::path output_root(argv[7]);
        std::filesystem::create_directories(output_root);

        std::array<AuditSignature, 7> top_signatures{};
        AuditSignature global_signature{};
        for (std::size_t row = 0; row < kTopPairsAudit.size(); ++row) {
            const auto found = signatures.find(audit_pair_key(
                kTopPairsAudit[row].q, kTopPairsAudit[row].r));
            require(found != signatures.end(), "top signature missing");
            top_signatures[row] = found->second;
            global_signature = audit_signature_union(global_signature,
                                                     found->second);
        }
        const std::array<std::vector<unsigned>, 7> expected_top = {{
            {25},
            {9,29,32,75,137,139,150,159,174,205,218,309,333,347,358,
             374,394,399,405,416,417},
            {412}, {173}, {396}, {236}, {107,374},
        }};
        for (std::size_t row = 0; row < expected_top.size(); ++row)
            require(audit_signature_indices(top_signatures[row]) ==
                        expected_top[row],
                    "top inactive signature changed");
        const std::vector<unsigned> global_indices =
            audit_signature_indices(global_signature);
        require(global_indices.size() == 27,
                "global inactive-signature union changed");
        std::array<int, kDeckCountAudit> global_position{};
        global_position.fill(-1);
        for (std::size_t local = 0; local < global_indices.size(); ++local)
            global_position[global_indices[local]] = static_cast<int>(local);

        std::vector<ScenarioAudit> scenarios;
        scenarios.reserve(133);
        for (unsigned scenario = 1; scenario < 128; ++scenario) {
            ScenarioAudit audit;
            audit.id = scenario;
            audit.label = "greedy";
            for (unsigned row = 0; row < 7; ++row)
                if ((scenario >> row) & 1)
                    audit.deleted = audit_signature_union(
                        audit.deleted, top_signatures[row]);
            for (unsigned index : audit_signature_indices(audit.deleted)) {
                require(global_position[index] >= 0,
                        "scenario deletion escaped global union");
                audit.deleted_local |=
                    u32{1} << global_position[index];
            }
            audit.witnesses = audit_read_witnesses(
                witness_root / audit_scenario_filename(scenario));
            scenarios.push_back(std::move(audit));
        }
        ScenarioAudit optimized_258;
        optimized_258.id = 128;
        optimized_258.label = "optimized_q258";
        optimized_258.deleted = top_signatures[2];
        optimized_258.deleted_local =
            u32{1} << global_position[expected_top[2][0]];
        optimized_258.witnesses = {UINT32_C(0x02511601)};
        scenarios.push_back(optimized_258);

        ScenarioAudit optimized_366;
        optimized_366.id = 129;
        optimized_366.label = "optimized_q366";
        optimized_366.deleted = top_signatures[4];
        optimized_366.deleted_local =
            u32{1} << global_position[expected_top[4][0]];
        optimized_366.witnesses = {UINT32_C(0x042022c9)};
        scenarios.push_back(optimized_366);

        for (const TargetVariantDefinition definition : kTargetVariants) {
            ScenarioAudit targeted;
            targeted.id = definition.id;
            targeted.label = definition.label;
            for (unsigned row = 0; row < 7; ++row)
                if ((definition.scenario >> row) & 1)
                    targeted.deleted = audit_signature_union(
                        targeted.deleted, top_signatures[row]);
            for (unsigned index : audit_signature_indices(targeted.deleted)) {
                require(global_position[index] >= 0,
                        "target variant deletion escaped global union");
                targeted.deleted_local |=
                    u32{1} << global_position[index];
            }
            targeted.witnesses = audit_read_witnesses(
                target_witness_root / definition.filename);
            scenarios.push_back(std::move(targeted));
        }

        std::vector<u32> global_bodies;
        std::vector<u32> global_responses;
        FnvLocal global_body_ledger;
        std::ofstream body_output(argv[8]);
        require(static_cast<bool>(body_output),
                "cannot create global-obligation ledger");
        body_output << "body,response_local\n";
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
                require(response != 0,
                        "original joint deck has an uncovered body");
                global_bodies.push_back(body);
                global_responses.push_back(response);
                global_body_ledger.add(body);
                body_output << std::hex << std::setw(8) << std::setfill('0')
                            << body << ',' << std::setw(7) << response
                            << std::dec << std::setfill(' ') << '\n';
            }
            if (ordinal + 1 < EXPECTED_BODIES)
                body = next_combination(body);
        }
        require(body == UINT32_C(0x3fe00000), "body endpoint changed");
        require(global_bodies.size() == kGlobalObligationCountAudit &&
                    global_body_ledger.state == kGlobalObligationFnvAudit,
                "global obligation family changed");
        require(body_output.good(), "global-obligation write failed");

        const std::set<u32> original_masks(deck.begin(), deck.end());
        for (ScenarioAudit& scenario : scenarios) {
            const std::vector<unsigned> deleted_indices =
                audit_signature_indices(scenario.deleted);
            std::set<u32> appended;
            for (u32 witness : scenario.witnesses) {
                require(original_masks.count(witness) == 0 &&
                            appended.insert(witness).second,
                        "scenario witness overlaps deck or repeats");
            }
            FnvLocal obligation_ledger;
            for (std::size_t index = 0; index < global_bodies.size(); ++index) {
                if ((global_responses[index] & ~scenario.deleted_local) != 0)
                    continue;
                ++scenario.obligation_count;
                obligation_ledger.add(global_bodies[index]);
                bool covered = false;
                for (u32 witness : scenario.witnesses)
                    if ((witness & global_bodies[index]) == 0) {
                        covered = true;
                        break;
                    }
                require(covered, "serialized scenario witness misses body");
            }
            scenario.obligation_fnv = obligation_ledger.state;
            FnvLocal deck_ledger;
            for (unsigned index = 0; index < deck.size(); ++index)
                if (((scenario.deleted[index / 64] >> (index % 64)) & 1)
                    == 0)
                    deck_ledger.add(deck[index]);
            for (u32 witness : scenario.witnesses)
                deck_ledger.add(witness);
            scenario.deck_size = deck.size() - deleted_indices.size() +
                                 scenario.witnesses.size();
            scenario.deck_fnv = deck_ledger.state;
            if (scenario.id <= 127) {
                const SummaryRow& expected = summary[scenario.id];
                require(deleted_indices.size() == expected.deleted &&
                            scenario.obligation_count == expected.obligations &&
                            scenario.obligation_fnv == expected.obligation_fnv &&
                            scenario.witnesses.size() == expected.witnesses &&
                            audit_mask_fnv(scenario.witnesses) ==
                                expected.witness_fnv &&
                            scenario.deck_size == expected.deck_size &&
                            scenario.deck_fnv == expected.deck_fnv,
                        "serialized scenario disagrees with summary");
            }
        }

        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "pool-cell count changed");
        u64 evaluated_cells = 0;
        u64 equality_cells = 0;
        i128 minimum_active_margin = 0;
        i128 maximum_inactive_margin = 0;
        bool have_active = false;
        bool have_inactive = false;
        std::size_t global_candidate_rows = 0;
        std::vector<std::size_t> global_candidate_indices;

        for (std::size_t row_index = 0; row_index < residual.size();
             ++row_index) {
            const AuditPair& pair = residual[row_index];
            const auto signature_found = signatures.find(
                audit_pair_key(pair.q, pair.r));
            require(signature_found != signatures.end(),
                    "residual signature missing");
            const AuditSignature& row_signature = signature_found->second;
            if (!audit_signature_subset(row_signature, global_signature))
                continue;
            ++global_candidate_rows;
            global_candidate_indices.push_back(row_index);
            std::vector<std::size_t> applicable;
            for (std::size_t scenario_index = 0;
                 scenario_index < scenarios.size(); ++scenario_index) {
                if (audit_signature_subset(
                        row_signature, scenarios[scenario_index].deleted)) {
                    applicable.push_back(scenario_index);
                    ++scenarios[scenario_index].candidate_count;
                }
            }
            require(!applicable.empty(),
                    "global candidate has no containing scenario");

            const i64 scale = std::gcd(pair.q, pair.r);
            const PrimitivePair primitive =
                build_primitive(pair.q / scale, pair.r / scale);
            const AtomData atoms =
                build_cocycle_atoms(pool_cells, primitive, scale);
            std::unordered_map<u32, i128> margin_cache;
            auto active = [&](u32 mask) {
                const auto found = margin_cache.find(mask);
                i128 margin = 0;
                if (found == margin_cache.end()) {
                    margin = audit_mask_margin(atoms, primitive, scale, mask);
                    margin_cache.emplace(mask, margin);
                    ++evaluated_cells;
                    if (margin == 0) ++equality_cells;
                    if (margin >= 0) {
                        if (!have_active || margin < minimum_active_margin)
                            minimum_active_margin = margin;
                        have_active = true;
                    } else {
                        if (!have_inactive || margin > maximum_inactive_margin)
                            maximum_inactive_margin = margin;
                        have_inactive = true;
                    }
                } else {
                    margin = found->second;
                }
                return margin >= 0;
            };
            for (std::size_t scenario_index : applicable) {
                bool common = true;
                for (u32 witness : scenarios[scenario_index].witnesses)
                    if (!active(witness)) {
                        common = false;
                        break;
                    }
                if (common)
                    scenarios[scenario_index].common_rows.push_back(row_index);
            }
        }
        require(equality_cells == 0,
                "an evaluated common-family activity cell is equality");

        std::ofstream scenario_summary(output_root /
                                       "scenario_common_summary.csv");
        require(static_cast<bool>(scenario_summary),
                "cannot create common scenario summary");
        scenario_summary <<
            "id,label,deleted,witnesses,obligations,obligation_fnv,deck_size,"
            "deck_fnv,candidates,common,common_fnv\n";
        for (ScenarioAudit& scenario : scenarios) {
            scenario.common_fnv =
                audit_family_fnv(scenario.common_rows, residual);
            scenario_summary << scenario.id << ',' << scenario.label << ','
                << audit_signature_indices(scenario.deleted).size() << ','
                << scenario.witnesses.size() << ','
                << scenario.obligation_count << ',' << std::hex
                << scenario.obligation_fnv << std::dec << ','
                << scenario.deck_size << ',' << std::hex << scenario.deck_fnv
                << std::dec << ',' << scenario.candidate_count << ','
                << scenario.common_rows.size() << ',' << std::hex
                << scenario.common_fnv << std::dec << '\n';
            audit_write_family(output_root / audit_common_filename(scenario),
                               scenario.common_rows, residual);
        }
        require(scenario_summary.good(),
                "common scenario summary write failed");

        const std::array<std::pair<unsigned, std::size_t>, 8>
            recovered_common_controls = {{
                {1, 41}, {2, 245}, {4, 6}, {8, 29},
                {16, 35}, {32, 48}, {64, 13}, {127, 460},
            }};
        for (const auto [scenario, expected] : recovered_common_controls)
            require(scenarios[scenario - 1].common_rows.size() == expected,
                    "recovered singleton/all-seven common count changed");
        require(scenarios[127].common_rows.size() == 33,
                "optimized q258 common count changed");
        require(scenarios[128].common_rows.size() == 36,
                "optimized q366 common count changed");
        for (std::size_t offset = 0; offset < kTargetVariants.size(); ++offset) {
            const TargetVariantDefinition& definition =
                kTargetVariants[offset];
            const ScenarioAudit& scenario = scenarios[129 + offset];
            const auto target_row = std::find_if(
                residual.begin(), residual.end(), [&](const AuditPair& pair) {
                    return pair.q == definition.target.q &&
                           pair.r == definition.target.r;
                });
            require(target_row != residual.end(), "target variant row missing");
            const std::size_t target_index =
                static_cast<std::size_t>(target_row - residual.begin());
            require(std::binary_search(scenario.common_rows.begin(),
                                       scenario.common_rows.end(),
                                       target_index),
                    "target-aware deck is not common at its target");
        }

        const std::vector<std::size_t> eight_indices =
            {0, 1, 3, 7, 15, 31, 63, 126};
        const std::vector<std::size_t> eight_union = audit_union_indices(
            scenarios, eight_indices, residual.size());
        require(eight_union.size() == 565,
                "recovered eight-deck common union changed");
        std::vector<std::size_t> optimized_eight_indices = eight_indices;
        optimized_eight_indices.push_back(127);
        optimized_eight_indices.push_back(128);
        const std::vector<std::size_t> optimized_eight_union =
            audit_union_indices(scenarios, optimized_eight_indices,
                                residual.size());
        require(optimized_eight_union.size() == 566,
                "recovered optimized eight-deck union changed");

        std::vector<std::size_t> base_indices(127);
        std::iota(base_indices.begin(), base_indices.end(), 0);
        const std::vector<std::size_t> union127 = audit_union_indices(
            scenarios, base_indices, residual.size());
        std::vector<std::size_t> extended_indices(scenarios.size());
        std::iota(extended_indices.begin(), extended_indices.end(), 0);
        const std::vector<std::size_t> extended_union = audit_union_indices(
            scenarios, extended_indices, residual.size());
        std::vector<std::size_t> uncovered_global_candidates;
        std::set_difference(global_candidate_indices.begin(),
                            global_candidate_indices.end(),
                            extended_union.begin(), extended_union.end(),
                            std::back_inserter(uncovered_global_candidates));
        require(uncovered_global_candidates.empty() &&
                    extended_union == global_candidate_indices,
                "target variants do not close the exact global candidate fibre");
        audit_write_family(output_root / "global_signature_candidates.csv",
                           global_candidate_indices, residual);
        audit_write_family(output_root /
                               "uncovered_global_signature_candidates.csv",
                           uncovered_global_candidates, residual);
        audit_write_family(output_root / "all127_common_union.csv",
                           union127, residual);
        audit_write_family(output_root /
                               "full_signature_fibre_common_union.csv",
                           extended_union, residual);

        const std::size_t bit_words = (residual.size() + 63) / 64;
        std::vector<std::vector<u64>> family_bits(
            scenarios.size(), std::vector<u64>(bit_words, 0));
        for (std::size_t scenario = 0; scenario < scenarios.size(); ++scenario)
            for (std::size_t row : scenarios[scenario].common_rows)
                family_bits[scenario][row / 64] |=
                    UINT64_C(1) << (row % 64);
        std::vector<u64> uncovered(bit_words, 0);
        for (std::size_t row : extended_union)
            uncovered[row / 64] |= UINT64_C(1) << (row % 64);
        auto bit_empty = [](const std::vector<u64>& bits) {
            for (u64 word : bits)
                if (word != 0) return false;
            return true;
        };
        std::vector<std::size_t> selected;
        std::vector<unsigned char> already_selected(scenarios.size(), 0);
        while (!bit_empty(uncovered)) {
            std::size_t best = scenarios.size();
            std::size_t best_gain = 0;
            for (std::size_t scenario = 0; scenario < scenarios.size();
                 ++scenario) {
                if (already_selected[scenario]) continue;
                std::size_t gain = 0;
                for (std::size_t word = 0; word < bit_words; ++word)
                    gain += std::popcount(
                        family_bits[scenario][word] & uncovered[word]);
                if (gain > best_gain) {
                    best = scenario;
                    best_gain = gain;
                }
            }
            require(best != scenarios.size() && best_gain != 0,
                    "scenario union cover became stuck");
            selected.push_back(best);
            already_selected[best] = 1;
            for (std::size_t word = 0; word < bit_words; ++word)
                uncovered[word] &= ~family_bits[best][word];
        }
        for (std::size_t position = selected.size(); position-- > 0;) {
            std::vector<std::size_t> trial;
            for (std::size_t index = 0; index < selected.size(); ++index)
                if (index != position) trial.push_back(selected[index]);
            if (audit_union_indices(scenarios, trial, residual.size()) ==
                extended_union)
                selected.erase(selected.begin() + position);
        }
        require(audit_union_indices(scenarios, selected, residual.size()) ==
                    extended_union,
                "selected scenario family does not reproduce union");

        std::ofstream selected_output(output_root /
                                      "selected_scenario_cover.txt");
        require(static_cast<bool>(selected_output),
                "cannot create selected scenario cover");
        for (std::size_t index : selected) {
            const ScenarioAudit& scenario = scenarios[index];
            selected_output << scenario.id << ' ' << scenario.label
                << " COMMON " << scenario.common_rows.size()
                << " COMMON_FNV " << std::hex << scenario.common_fnv
                << " DECK_FNV " << scenario.deck_fnv << std::dec << '\n';
        }
        require(selected_output.good(), "selected scenario write failed");

        std::cout << "THM4283_ALL127_COMMON_FAMILY_AUDIT_V1\n"
                  << "GLOBAL_OBLIGATIONS " << global_bodies.size()
                  << " FNV " << std::hex << global_body_ledger.state
                  << std::dec << " SCENARIOS " << scenarios.size() << '\n'
                  << "GLOBAL_CANDIDATES " << global_candidate_rows
                  << " EVALUATED_WITNESS_CELLS " << evaluated_cells
                  << " EQUALITIES " << equality_cells
                  << " MIN_ACTIVE_MARGIN "
                  << decimal(minimum_active_margin)
                  << " MAX_INACTIVE_MARGIN "
                  << decimal(maximum_inactive_margin) << '\n';
        for (const ScenarioAudit& scenario : scenarios) {
            std::cout << "SCENARIO " << scenario.id << " LABEL "
                      << scenario.label << " DELETED "
                      << audit_signature_indices(scenario.deleted).size()
                      << " WITNESSES " << scenario.witnesses.size()
                      << " CANDIDATES " << scenario.candidate_count
                      << " COMMON " << scenario.common_rows.size()
                      << " COMMON_FNV " << std::hex << scenario.common_fnv
                      << " DECK_FNV " << scenario.deck_fnv << std::dec
                      << '\n';
        }
        std::cout << "EIGHT_UNION " << eight_union.size() << " FNV "
                  << std::hex << audit_family_fnv(eight_union, residual)
                  << std::dec << " OPTIMIZED_EIGHT_UNION "
                  << optimized_eight_union.size() << " FNV " << std::hex
                  << audit_family_fnv(optimized_eight_union, residual)
                  << std::dec << '\n'
                  << "ALL127_UNION " << union127.size() << " FNV "
                  << std::hex << audit_family_fnv(union127, residual)
                  << std::dec << " FULL_SIGNATURE_FIBRE "
                  << extended_union.size()
                  << " FNV " << std::hex
                  << audit_family_fnv(extended_union, residual) << std::dec
                  << " UNCOVERED_GLOBAL_CANDIDATES "
                  << uncovered_global_candidates.size() << " FNV "
                  << std::hex
                  << audit_family_fnv(uncovered_global_candidates, residual)
                  << std::dec << '\n'
                  << "SELECTED_SCENARIOS " << selected.size() << " IDS";
        for (std::size_t index : selected)
            std::cout << ' ' << scenarios[index].id;
        std::cout << "\nVERDICT PASS EXACT_COMMON_FAMILY_SYNTHESIS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_COMMON_FAMILY_ERROR " << error.what() << '\n';
        return 1;
    }
}
