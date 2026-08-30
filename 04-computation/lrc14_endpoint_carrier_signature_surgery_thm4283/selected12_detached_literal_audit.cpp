// Detached audit of the twelve explicit surgery decks selected by the THM-4283
// common-family synthesis.  This translation unit imports no project source
// and uses neither the primitive/cocycle nor atom/zeta activity engine.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr std::array<std::pair<int, int>, 7> kTopPairs = {{
    {220, 644}, {256, 644}, {258, 644}, {294, 644},
    {366, 644}, {416, 644}, {512, 644},
}};
constexpr std::array<unsigned, 12> kSelectedIds =
    {63, 87, 91, 110, 17, 102, 131, 130, 32, 83, 132, 133};
constexpr std::size_t kJointCount = 421;
constexpr u64 kJointFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr i64 kFixedGrid = INT64_C(18241159416480);

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 word) {
        for (unsigned shift = 0; shift < 64; shift += 8) {
            state ^= (word >> shift) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

void add_i128(Fnv& ledger, i128 value) {
    const u128 bits = static_cast<u128>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    u128 magnitude = static_cast<u128>(value);
    if (negative) magnitude = u128{0} - magnitude;
    std::string answer;
    while (magnitude != 0) {
        answer.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

std::string hex8(u32 value) {
    std::ostringstream output;
    output << std::hex << std::setw(8) << std::setfill('0') << value;
    return output.str();
}

u32 parse_hex32(const std::string& token) {
    std::size_t used = 0;
    const u64 wide = std::stoull(token, &used, 16);
    require(used == token.size() && wide < (UINT64_C(1) << 30),
            "bad hexadecimal mask token");
    const u32 mask = static_cast<u32>(wide);
    require(std::popcount(mask) == 8, "mask does not have rank eight");
    return mask;
}

std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, delimiter)) fields.push_back(field);
    return fields;
}

struct Pair {
    int q = 0;
    int r = 0;
    friend bool operator<(const Pair& left, const Pair& right) {
        return std::tie(left.q, left.r) < std::tie(right.q, right.r);
    }
    friend bool operator==(const Pair& left, const Pair& right) {
        return left.q == right.q && left.r == right.r;
    }
};

std::vector<u32> read_masks(const std::filesystem::path& path,
                            bool require_nonempty = true) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask file " + path.string());
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_hex32(token);
        require(distinct.insert(mask).second, "duplicate mask in input");
        masks.push_back(mask);
    }
    require(!require_nonempty || !masks.empty(), "empty mask file");
    return masks;
}

std::vector<u32> read_joint(const std::filesystem::path& path) {
    const std::vector<u32> deck = read_masks(path);
    Fnv ledger;
    for (u32 mask : deck) ledger.add(mask);
    require(deck.size() == kJointCount && ledger.state == kJointFnv,
            "joint deck identity changed");
    return deck;
}

using Signature = std::array<u64, 7>;

std::array<Signature, 7> read_top_signatures(
        const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::array<Signature, 7> answer{};
    std::array<unsigned char, 7> found{};
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 10, "malformed signature row");
        const int q = std::stoi(fields[0]);
        const int r = std::stoi(fields[1]);
        for (std::size_t top = 0; top < kTopPairs.size(); ++top) {
            if (kTopPairs[top] != std::pair<int, int>{q, r}) continue;
            unsigned weight = 0;
            for (std::size_t word = 0; word < 7; ++word) {
                answer[top][word] =
                    std::stoull(fields[word + 3], nullptr, 16);
                weight += std::popcount(answer[top][word]);
            }
            require(weight == std::stoul(fields[2]),
                    "top signature weight changed");
            found[top] = 1;
        }
    }
    for (unsigned char value : found)
        require(value == 1, "top signature missing");
    return answer;
}

Signature signature_union(Signature left, const Signature& right) {
    for (std::size_t word = 0; word < left.size(); ++word)
        left[word] |= right[word];
    return left;
}

std::vector<unsigned> signature_indices(const Signature& signature) {
    std::vector<unsigned> answer;
    for (unsigned index = 0; index < kJointCount; ++index)
        if ((signature[index / 64] >> (index % 64)) & 1)
            answer.push_back(index);
    return answer;
}

struct SelectedLine {
    unsigned id = 0;
    std::string label;
    std::size_t common_count = 0;
    u64 common_fnv = 0;
    u64 deck_fnv = 0;
};

std::vector<SelectedLine> read_selected(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open selected scenario file");
    std::vector<SelectedLine> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::stringstream stream(line);
        SelectedLine row;
        std::string common_tag;
        std::string common_fnv_tag;
        std::string common_fnv;
        std::string deck_fnv_tag;
        std::string deck_fnv;
        stream >> row.id >> row.label >> common_tag >> row.common_count
               >> common_fnv_tag >> common_fnv >> deck_fnv_tag >> deck_fnv;
        require(stream && common_tag == "COMMON" &&
                    common_fnv_tag == "COMMON_FNV" &&
                    deck_fnv_tag == "DECK_FNV",
                "malformed selected scenario row");
        row.common_fnv = std::stoull(common_fnv, nullptr, 16);
        row.deck_fnv = std::stoull(deck_fnv, nullptr, 16);
        rows.push_back(row);
    }
    require(rows.size() == kSelectedIds.size(),
            "selected scenario count changed");
    const std::array<std::string, 12> expected_labels = {
        "greedy", "greedy", "greedy", "greedy", "greedy", "greedy",
        "target_250_256_s010", "target_206_263_s002", "greedy", "greedy",
        "target_256_394_s002", "target_256_400_s002"};
    for (std::size_t index = 0; index < rows.size(); ++index)
        require(rows[index].id == kSelectedIds[index] &&
                    rows[index].label == expected_labels[index],
                "selected scenario order changed");
    return rows;
}

unsigned deletion_scenario(unsigned id) {
    switch (id) {
        case 130: return 2;
        case 131: return 10;
        case 132: return 2;
        case 133: return 2;
        default:
            require(id >= 1 && id <= 127, "unknown selected scenario id");
            return id;
    }
}

std::string witness_filename(const SelectedLine& line) {
    if (line.id <= 127) {
        std::ostringstream name;
        name << "scenario_" << std::setw(3) << std::setfill('0') << line.id
             << "_witnesses.txt";
        return name.str();
    }
    switch (line.id) {
        case 130: return "target_206_263_scenario_002_witnesses.txt";
        case 131: return "target_250_256_scenario_010_witnesses.txt";
        case 132: return "target_256_394_scenario_002_witnesses.txt";
        case 133: return "target_256_400_scenario_002_witnesses.txt";
        default: fail("unknown target witness id");
    }
}

std::vector<Pair> read_common(const std::filesystem::path& path,
                              std::size_t expected_count,
                              u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open common-row ledger");
    std::vector<Pair> rows;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty() && line.back() != '\r',
                "bad common-row line ending");
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 2, "malformed common-row ledger");
        Pair row{std::stoi(fields[0]), std::stoi(fields[1])};
        require(row.q > 0 && row.q < row.r &&
                    (rows.empty() || rows.back() < row),
                "common-row order changed");
        rows.push_back(row);
        ledger.add(row.q);
        ledger.add(row.r);
    }
    require(rows.size() == expected_count && ledger.state == expected_fnv,
            "common-row identity changed");
    return rows;
}

struct Scenario {
    SelectedLine expected;
    Signature deleted{};
    std::vector<u32> witnesses;
    std::vector<u32> deck;
    std::vector<Pair> common;
    u64 body_checks = 0;
    u64 maximum_prefix = 0;
};

u32 next_combination(u32 value) {
    const u32 low = value & (u32{0} - value);
    require(low != 0, "zero low bit in combination successor");
    const u32 ripple = value + low;
    return ripple | (((ripple ^ value) >> 2) / low);
}

void audit_body_cover(Scenario& scenario) {
    const u32 limit = u32{1} << 30;
    u64 bodies = 0;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++bodies;
        bool covered = false;
        u64 prefix = 0;
        for (u32 mask : scenario.deck) {
            ++prefix;
            ++scenario.body_checks;
            if ((mask & body) == 0) {
                covered = true;
                break;
            }
        }
        scenario.maximum_prefix = std::max(scenario.maximum_prefix, prefix);
        require(covered, "selected scenario deck misses body " + hex8(body));
    }
    require(bodies == kBodyCount, "nine-body universe count changed");
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(),
            "literal grid overflows i64");
    return static_cast<i64>(wide);
}

i64 fixed_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == kFixedGrid, "fixed-pool literal grid changed");
    return grid;
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

struct Geometry {
    i64 grid = 0;
    u64 cells = 0;
    i64 pair_ticks = 0;
    std::map<u32, i64> safe_width_by_failure;
};

Geometry build_geometry(Pair pair) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair.q);
    grid = checked_lcm(grid, 14LL * pair.r);
    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(pair.q);
    add_walls(pair.r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.front() == 0 && walls.back() == grid,
            "literal wall boundary changed");
    Geometry geometry;
    geometry.grid = grid;
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        require(left < right, "nonpositive literal cell");
        if (!safe_midpoint(pair.q, grid, left, right) ||
            !safe_midpoint(pair.r, grid, left, right))
            continue;
        u32 failed = 0;
        for (unsigned vertex = 0; vertex < kPool.size(); ++vertex)
            if (!safe_midpoint(kPool[vertex], grid, left, right))
                failed |= u32{1} << vertex;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        geometry.safe_width_by_failure[failed] += width;
    }
    return geometry;
}

struct DeckIndex {
    std::vector<std::vector<u64>> contains;
    std::vector<u64> all;
};

DeckIndex build_deck_index(const std::vector<u32>& deck) {
    const std::size_t words = (deck.size() + 63) / 64;
    DeckIndex index;
    index.contains.assign(30, std::vector<u64>(words, 0));
    index.all.assign(words, 0);
    for (std::size_t mask_index = 0; mask_index < deck.size(); ++mask_index) {
        index.all[mask_index / 64] |= UINT64_C(1) << (mask_index % 64);
        for (unsigned vertex = 0; vertex < 30; ++vertex)
            if ((deck[mask_index] >> vertex) & 1)
                index.contains[vertex][mask_index / 64] |=
                    UINT64_C(1) << (mask_index % 64);
    }
    return index;
}

std::vector<i64> classify_all(const Geometry& geometry,
                              const std::vector<u32>& deck,
                              const DeckIndex& index) {
    std::vector<i64> mass(deck.size(), 0);
    for (const auto& [failure, width] : geometry.safe_width_by_failure) {
        if (std::popcount(failure) > 8) continue;
        std::vector<u64> candidates = index.all;
        u32 remaining = failure;
        while (remaining != 0) {
            const unsigned vertex = std::countr_zero(remaining);
            remaining &= remaining - 1;
            for (std::size_t word = 0; word < candidates.size(); ++word)
                candidates[word] &= index.contains[vertex][word];
        }
        for (std::size_t word = 0; word < candidates.size(); ++word) {
            u64 bits = candidates[word];
            while (bits != 0) {
                const unsigned offset = std::countr_zero(bits);
                const std::size_t mask_index = 64 * word + offset;
                require(mask_index < mass.size(), "deck bitset tail escaped");
                mass[mask_index] += width;
                bits &= bits - 1;
            }
        }
    }
    return mass;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 8,
                "usage: detached-audit JOINT421 SIGNATURES WITNESS_DIR "
                "TARGET_WITNESS_DIR COMMON_DIR SELECTED_TXT OUTPUT_DIR");
        const std::filesystem::path witness_root(argv[3]);
        const std::filesystem::path target_witness_root(argv[4]);
        const std::filesystem::path common_root(argv[5]);
        const std::filesystem::path output_root(argv[7]);
        require(std::filesystem::is_directory(output_root),
                "output directory does not exist");
        const std::vector<u32> joint = read_joint(argv[1]);
        const std::array<Signature, 7> top_signatures =
            read_top_signatures(argv[2]);
        const std::vector<SelectedLine> selected = read_selected(argv[6]);

        std::vector<Scenario> scenarios;
        for (const SelectedLine& line : selected) {
            Scenario scenario;
            scenario.expected = line;
            const unsigned deletion = deletion_scenario(line.id);
            for (unsigned row = 0; row < 7; ++row)
                if ((deletion >> row) & 1)
                    scenario.deleted = signature_union(
                        scenario.deleted, top_signatures[row]);
            const std::filesystem::path selected_witness_root =
                line.id <= 127 ? witness_root : target_witness_root;
            scenario.witnesses = read_masks(
                selected_witness_root / witness_filename(line));
            const std::vector<unsigned> deleted_vector =
                signature_indices(scenario.deleted);
            const std::set<unsigned> deleted(deleted_vector.begin(),
                                             deleted_vector.end());
            for (unsigned index = 0; index < joint.size(); ++index)
                if (deleted.count(index) == 0)
                    scenario.deck.push_back(joint[index]);
            scenario.deck.insert(scenario.deck.end(),
                                 scenario.witnesses.begin(),
                                 scenario.witnesses.end());
            std::set<u32> distinct(scenario.deck.begin(), scenario.deck.end());
            require(distinct.size() == scenario.deck.size(),
                    "selected deck contains duplicate mask");
            Fnv deck_ledger;
            for (u32 mask : scenario.deck) deck_ledger.add(mask);
            require(deck_ledger.state == line.deck_fnv,
                    "selected deck FNV changed");
            std::ostringstream common_name;
            common_name << "scenario_" << std::setw(3) << std::setfill('0')
                        << line.id;
            if (line.id > 127) common_name << '_' << line.label;
            common_name << "_common.csv";
            scenario.common = read_common(common_root / common_name.str(),
                                          line.common_count,
                                          line.common_fnv);
            for (unsigned row = 0; row < 7; ++row)
                if ((deletion >> row) & 1)
                    require(std::binary_search(
                                scenario.common.begin(), scenario.common.end(),
                                Pair{kTopPairs[row].first,
                                     kTopPairs[row].second}),
                            "selected endpoint absent from common family");
            scenarios.push_back(std::move(scenario));
        }

        std::atomic<std::size_t> next_scenario{0};
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned worker_count =
            std::max(1u, std::min(4u, hardware ? hardware : 1u));
        std::vector<std::thread> workers;
        for (unsigned worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next_scenario.fetch_add(1);
                    if (index >= scenarios.size()) break;
                    audit_body_cover(scenarios[index]);
                }
            });
        }
        for (std::thread& worker : workers) worker.join();

        std::map<Pair, std::vector<std::size_t>> row_scenarios;
        for (std::size_t scenario = 0; scenario < scenarios.size(); ++scenario)
            for (Pair pair : scenarios[scenario].common)
                row_scenarios[pair].push_back(scenario);
        require(row_scenarios.size() == 640,
                "selected common union count changed");
        std::vector<DeckIndex> deck_indices;
        for (const Scenario& scenario : scenarios)
            deck_indices.push_back(build_deck_index(scenario.deck));

        std::ofstream literal_output(output_root /
                                     "selected12_literal_rows.csv");
        require(static_cast<bool>(literal_output),
                "cannot create detached literal ledger");
        literal_output << "scenario,q,r,deck_size,grid,cells,pair_ticks,"
                          "failure_patterns,min_mask,min_mass,min_margin\n";
        u64 literal_checks = 0;
        u64 equalities = 0;
        i128 global_minimum = 0;
        Pair global_minimum_pair;
        unsigned global_minimum_scenario = 0;
        u32 global_minimum_mask = 0;
        bool minimum_set = false;
        Fnv geometry_ledger;
        Fnv activity_ledger;
        for (const auto& [pair, containing_scenarios] : row_scenarios) {
            const Geometry geometry = build_geometry(pair);
            geometry_ledger.add(pair.q);
            geometry_ledger.add(pair.r);
            geometry_ledger.add(geometry.grid);
            geometry_ledger.add(geometry.cells);
            geometry_ledger.add(geometry.pair_ticks);
            geometry_ledger.add(geometry.safe_width_by_failure.size());
            for (std::size_t scenario_index : containing_scenarios) {
                const Scenario& scenario = scenarios[scenario_index];
                const std::vector<i64> mass = classify_all(
                    geometry, scenario.deck, deck_indices[scenario_index]);
                i128 row_minimum = 0;
                std::size_t row_minimum_index = 0;
                bool row_minimum_set = false;
                for (std::size_t mask_index = 0;
                     mask_index < scenario.deck.size(); ++mask_index) {
                    const i128 margin =
                        static_cast<i128>(63) * mass[mask_index] -
                        static_cast<i128>(4) * geometry.grid;
                    ++literal_checks;
                    equalities += margin == 0;
                    require(margin >= 0,
                            "claimed common cell is literal-inactive");
                    activity_ledger.add(scenario.expected.id);
                    activity_ledger.add(pair.q);
                    activity_ledger.add(pair.r);
                    activity_ledger.add(scenario.deck[mask_index]);
                    activity_ledger.add(mass[mask_index]);
                    add_i128(activity_ledger, margin);
                    if (!row_minimum_set || margin < row_minimum ||
                        (margin == row_minimum &&
                         scenario.deck[mask_index] <
                            scenario.deck[row_minimum_index])) {
                        row_minimum_set = true;
                        row_minimum = margin;
                        row_minimum_index = mask_index;
                    }
                    if (!minimum_set || margin < global_minimum ||
                        (margin == global_minimum &&
                         std::tie(scenario.expected.id, pair.q, pair.r,
                                  scenario.deck[mask_index]) <
                         std::tie(global_minimum_scenario,
                                  global_minimum_pair.q,
                                  global_minimum_pair.r,
                                  global_minimum_mask))) {
                        minimum_set = true;
                        global_minimum = margin;
                        global_minimum_pair = pair;
                        global_minimum_scenario = scenario.expected.id;
                        global_minimum_mask = scenario.deck[mask_index];
                    }
                }
                literal_output << scenario.expected.id << ',' << pair.q
                    << ',' << pair.r << ',' << scenario.deck.size() << ','
                    << geometry.grid << ',' << geometry.cells << ','
                    << geometry.pair_ticks << ','
                    << geometry.safe_width_by_failure.size() << ','
                    << hex8(scenario.deck[row_minimum_index]) << ','
                    << mass[row_minimum_index] << ','
                    << decimal(row_minimum) << '\n';
            }
        }
        require(literal_output.good(), "detached literal write failed");
        require(equalities == 0, "detached literal equality encountered");

        Fnv body_ledger;
        u64 total_body_checks = 0;
        for (const Scenario& scenario : scenarios) {
            body_ledger.add(scenario.expected.id);
            body_ledger.add(scenario.deck.size());
            body_ledger.add(scenario.body_checks);
            body_ledger.add(scenario.maximum_prefix);
            total_body_checks += scenario.body_checks;
        }
        std::cout << "THM4283_SELECTED12_DETACHED_LITERAL_AUDIT_V1\n"
                  << "SCENARIOS " << scenarios.size() << " COMMON_UNION "
                  << row_scenarios.size() << " BODY_UNIVERSE_EACH "
                  << kBodyCount << '\n';
        for (const Scenario& scenario : scenarios)
            std::cout << "SCENARIO " << scenario.expected.id << " DECK "
                      << scenario.deck.size() << " DECK_FNV " << std::hex
                      << scenario.expected.deck_fnv << std::dec
                      << " COMMON " << scenario.common.size()
                      << " COMMON_FNV " << std::hex
                      << scenario.expected.common_fnv << std::dec
                      << " BODY_CHECKS " << scenario.body_checks
                      << " MAX_PREFIX " << scenario.maximum_prefix << '\n';
        std::cout << "BODY_TOTAL_CHECKS " << total_body_checks
                  << " BODY_LEDGER_FNV " << std::hex << body_ledger.state
                  << std::dec << '\n'
                  << "LITERAL_CHECKS " << literal_checks
                  << " EQUALITIES " << equalities
                  << " GEOMETRY_FNV " << std::hex << geometry_ledger.state
                  << " ACTIVITY_FNV " << activity_ledger.state << std::dec
                  << '\n'
                  << "GLOBAL_MINIMUM SCENARIO " << global_minimum_scenario
                  << " ROW " << global_minimum_pair.q << ','
                  << global_minimum_pair.r << " MASK "
                  << hex8(global_minimum_mask) << " MARGIN "
                  << decimal(global_minimum) << '\n'
                  << "SCOPE EXPLICIT_DECKS_AND_CLAIMED_POSITIVE_CELLS_ONLY\n"
                  << "VERDICT PASS EXHAUSTIVE_BODY_AND_DETACHED_LITERAL_AUDIT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_SELECTED10_DETACHED_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
