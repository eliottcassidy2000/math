// Scratch exact rank-nine and mixed rank-eight/rank-nine response quotient.
// The literal-wall engine is the import-free scratch certificate beside this
// file; no maintained project source is imported.

#define main r632_detached_base_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

#include <unordered_map>

namespace {

constexpr u64 kRank9Count = UINT64_C(14307150);
constexpr u64 kDisjoint9 = UINT64_C(293930);

struct Pattern66 {
    u64 low = 0;
    u64 high = 0;
    auto operator<=>(const Pattern66&) const = default;
};

struct Class89 {
    u64 count8 = 0;
    u64 count9 = 0;
    u32 least8 = 0;
    u32 least9 = 0;
};

std::array<std::array<u64, 10>, 31> choose_table{};

void init_choose() {
    for (unsigned n = 0; n <= 30; ++n) {
        choose_table[n][0] = 1;
        for (unsigned k = 1; k <= 9; ++k) {
            if (k > n) choose_table[n][k] = 0;
            else if (k == n) choose_table[n][k] = 1;
            else choose_table[n][k] =
                choose_table[n - 1][k - 1] + choose_table[n - 1][k];
        }
    }
    require(choose_table[30][9] == kRank9Count,
            "rank-nine universe changed");
}

u64 colex9(u32 mask) {
    require(std::popcount(mask) == 9, "colex9 rank changed");
    u64 rank = 0;
    unsigned ordinal = 1;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((mask >> bit) & 1) rank += choose_table[bit][ordinal++];
    require(ordinal == 10 && rank < kRank9Count, "colex9 escaped");
    return rank;
}

template <class Callback>
void choose9_rec(const std::vector<unsigned char>& positions,
                 std::size_t start, unsigned chosen, u32 mask,
                 Callback& callback) {
    if (chosen == 9) {
        callback(mask);
        return;
    }
    const std::size_t needed = 9 - chosen;
    for (std::size_t index = start;
         index + needed <= positions.size(); ++index)
        choose9_rec(positions, index + 1, chosen + 1,
                    mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
u64 enumerate_disjoint9(u32 body, Callback callback) {
    std::vector<unsigned char> positions;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((body & (u32{1} << bit)) == 0)
            positions.push_back(static_cast<unsigned char>(bit));
    require(positions.size() == 21, "rank-nine body complement changed");
    u64 count = 0;
    auto counted = [&](u32 mask) {
        require(std::popcount(mask) == 9 && (mask & body) == 0,
                "rank-nine complement escaped");
        callback(mask);
        ++count;
    };
    choose9_rec(positions, 0, 0, 0, counted);
    require(count == kDisjoint9, "rank-nine disjoint count changed");
    return count;
}

void set_bit(Pattern66& pattern, std::size_t index) {
    if (index < 64) pattern.low |= UINT64_C(1) << index;
    else pattern.high |= UINT64_C(1) << (index - 64);
}

u64 count_bits(Pattern66 pattern) {
    return std::popcount(pattern.low) + std::popcount(pattern.high);
}

Pattern66 shift_global_b(u64 global_low, u64 global_high) {
    return {(global_low >> 6) | (global_high << 58), global_high >> 6};
}

u64 mass_subset(const Geometry& geometry, u32 mask) {
    i64 mass = 0;
    for (const auto& [failure, width] : geometry.classes)
        if ((failure & ~mask) == 0) mass += width;
    return static_cast<u64>(mass);
}

std::map<Pattern66, Class89> read_rank8_classes(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open rank-eight atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "a_hex\tb_hex\tcover\tcount\tleast_mask",
            "rank-eight atlas header changed");
    std::map<Pattern66, Class89> classes;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), '\t', ' ');
        std::istringstream row(line);
        std::string low_hex;
        std::string high_hex;
        std::string mask_hex;
        u64 cover = 0;
        u64 count = 0;
        row >> low_hex >> high_hex >> cover >> count >> mask_hex;
        require(static_cast<bool>(row), "malformed rank-eight atlas row");
        const u64 global_low = std::stoull(low_hex, nullptr, 16);
        const u64 global_high = std::stoull(high_hex, nullptr, 16);
        const Pattern66 pattern = shift_global_b(global_low, global_high);
        require(count_bits(pattern) == cover &&
                    (pattern.low != 0 || pattern.high != 0),
                "rank-eight response shift changed");
        Class89& c = classes[pattern];
        require(c.count8 == 0, "duplicate rank-eight pattern");
        c.count8 = count;
        c.least8 = static_cast<u32>(std::stoul(mask_hex, nullptr, 16));
        require(std::popcount(c.least8) == 8, "rank-eight least mask changed");
    }
    require(classes.size() == 395, "rank-eight class count changed");
    return classes;
}

}  // namespace


int main(int argc, char** argv) {
    try {
        require(argc == 5,
                "usage: rank9 FAILURES RANK8_ATLAS RANK9_TSV MIXED_TSV");
        init_choose();
        const std::vector<u32> failures = read_failures(argv[1]);
        const Geometry geometry = build_geometry();
        std::vector<Pattern66> raw(kRank9Count);
        u64 raw_incidences = 0;
        for (std::size_t index = 0; index < failures.size(); ++index)
            enumerate_disjoint9(failures[index], [&](u32 mask) {
                set_bit(raw[colex9(mask)], index);
                ++raw_incidences;
            });
        require(raw_incidences == failures.size() * kDisjoint9,
                "rank-nine raw incidence count changed");

        std::map<Pattern66, Class89> rank9_classes;
        Fnv active_response_ledger;
        u64 raw_candidates = 0;
        u64 active_candidates = 0;
        u64 active_incidences = 0;
        u64 full_responders = 0;
        u32 least_full = 0;
        u64 rank = 0;
        const u32 limit = u32{1} << 30;
        for (u32 mask = (u32{1} << 9) - 1; mask < limit;
             mask = next_combination(mask), ++rank) {
            require(rank == colex9(mask), "numeric/colex rank-nine order changed");
            const Pattern66 response = raw[rank];
            if (response.low == 0 && response.high == 0) continue;
            ++raw_candidates;
            const i64 mass = static_cast<i64>(mass_subset(geometry, mask));
            const i128 ticks = static_cast<i128>(63) * mass -
                               static_cast<i128>(4) * geometry.grid;
            if (ticks < 0) continue;
            ++active_candidates;
            active_incidences += count_bits(response);
            active_response_ledger.add(mask);
            active_response_ledger.add(response.low);
            active_response_ledger.add(response.high);
            Class89& c = rank9_classes[response];
            ++c.count9;
            if (c.count9 == 1 || mask < c.least9) c.least9 = mask;
            if (count_bits(response) == failures.size()) {
                ++full_responders;
                if (full_responders == 1 || mask < least_full) least_full = mask;
            }
        }
        require(rank == kRank9Count, "rank-nine enumeration count changed");

        std::ofstream rank9_output(argv[3]);
        require(static_cast<bool>(rank9_output), "cannot create rank-nine atlas");
        rank9_output << "a_hex\tb_hex\tcover\tcount\tleast_mask\tleast_rank\n";
        Fnv rank9_class_ledger;
        u64 maximum_cover = 0;
        u32 least_maximum = 0;
        for (const auto& [pattern, c] : rank9_classes) {
            require(c.count8 == 0 && c.count9 > 0, "rank-nine class type changed");
            rank9_class_ledger.add(pattern.low);
            rank9_class_ledger.add(pattern.high);
            rank9_class_ledger.add(c.count9);
            rank9_class_ledger.add(c.least9);
            const u64 cover = count_bits(pattern);
            if (cover > maximum_cover ||
                (cover == maximum_cover && c.least9 < least_maximum)) {
                maximum_cover = cover;
                least_maximum = c.least9;
            }
            rank9_output << std::hex << std::setw(16) << std::setfill('0')
                         << pattern.low << '\t' << std::setw(16)
                         << pattern.high << std::dec << std::setfill(' ')
                         << '\t' << cover << '\t' << c.count9 << '\t'
                         << hex8(c.least9) << "\t9\n";
        }
        require(rank9_output.good(), "rank-nine atlas write failed");

        std::map<Pattern66, Class89> mixed = read_rank8_classes(argv[2]);
        for (const auto& [pattern, rank9_class] : rank9_classes) {
            Class89& c = mixed[pattern];
            c.count9 = rank9_class.count9;
            c.least9 = rank9_class.least9;
        }
        std::ofstream mixed_output(argv[4]);
        require(static_cast<bool>(mixed_output), "cannot create mixed atlas");
        mixed_output << "a_hex\tb_hex\tcover\tcount\tleast_mask\tleast_rank"
                        "\tcount8\tcount9\n";
        Fnv mixed_ledger;
        for (const auto& [pattern, c] : mixed) {
            const u64 count = c.count8 + c.count9;
            const bool choose8 = c.count8 > 0 &&
                (c.count9 == 0 || c.least8 < c.least9);
            const u32 least = choose8 ? c.least8 : c.least9;
            const unsigned least_rank = choose8 ? 8 : 9;
            require(count > 0, "empty mixed class");
            mixed_ledger.add(pattern.low);
            mixed_ledger.add(pattern.high);
            mixed_ledger.add(c.count8);
            mixed_ledger.add(c.count9);
            mixed_ledger.add(least);
            mixed_ledger.add(least_rank);
            mixed_output << std::hex << std::setw(16) << std::setfill('0')
                         << pattern.low << '\t' << std::setw(16)
                         << pattern.high << std::dec << std::setfill(' ')
                         << '\t' << count_bits(pattern) << '\t' << count
                         << '\t' << hex8(least) << '\t' << least_rank << '\t'
                         << c.count8 << '\t' << c.count9 << '\n';
        }
        require(mixed_output.good(), "mixed atlas write failed");
        std::cout << "R632_RANK9_MIXED_RESPONSE_ATLAS_V1\n"
                  << "FAILURES 66 FNV " << std::hex << kFailuresFnv << std::dec
                  << " RAW_INCIDENCES " << raw_incidences
                  << " RAW_CANDIDATES " << raw_candidates << '\n'
                  << "RANK9_ACTIVE_CANDIDATES " << active_candidates
                  << " ACTIVE_INCIDENCES " << active_incidences
                  << " CLASSES " << rank9_classes.size()
                  << " RESPONSE_FNV " << std::hex
                  << active_response_ledger.state << " CLASS_FNV "
                  << rank9_class_ledger.state << std::dec
                  << " MAX_COVER " << maximum_cover << " LEAST_MAX "
                  << hex8(least_maximum) << " FULL_RESPONDERS "
                  << full_responders << " LEAST_FULL " << hex8(least_full)
                  << '\n'
                  << "MIXED_CLASSES " << mixed.size() << " LEDGER_FNV "
                  << std::hex << mixed_ledger.state << std::dec << '\n'
                  << "SCOPE LITERAL_WALL_FIXED_PAIR_LABELLED_RANK9_AND_"
                     "PINNED_RANK8_ATLAS_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_RANK9_AND_MIXED_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R632_RANK9_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}
