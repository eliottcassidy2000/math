// Scratch-only exact probe: test whether an active fixed-rank repair family
// is closed under elementary left-compressions.  This includes the maintained
// import-free literal-wall engine but does not modify any theorem artifact.

#define main thm4296_r632_hostile_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

#include <unordered_set>

namespace {

u64 colex_any(u32 mask) {
    std::array<std::array<u64, 10>, 31> c{};
    for (unsigned n = 0; n <= 30; ++n) {
        c[n][0] = 1;
        for (unsigned k = 1; k <= 9; ++k)
            c[n][k] = (k > n) ? 0 : ((k == n) ? 1 : c[n - 1][k] + c[n - 1][k - 1]);
    }
    u64 answer = 0;
    unsigned ordinal = 1;
    for (unsigned bit = 0; bit < 30; ++bit)
        if (mask & (u32{1} << bit)) answer += c[bit][ordinal++];
    return answer;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream in(line);
    std::string field;
    while (std::getline(in, field, '\t')) out.push_back(field);
    return out;
}

std::vector<u32> read_representatives(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)), "empty atlas");
    const std::vector<std::string> header = split_tab(line);
    auto find_column = [&](const std::string& name) {
        const auto it = std::find(header.begin(), header.end(), name);
        require(it != header.end(), "missing atlas column " + name);
        return static_cast<std::size_t>(it - header.begin());
    };
    const std::size_t mask_column = find_column("least_mask");
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_tab(line);
        require(fields.size() == header.size(), "malformed atlas row");
        const u32 mask = static_cast<u32>(std::stoul(fields[mask_column], nullptr, 16));
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    mask < (u32{1} << 30),
                "invalid representative mask");
        distinct.insert(mask);
    }
    return {distinct.begin(), distinct.end()};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4, "usage: compressed_probe Q R ATLAS_TSV");
        const int q = std::stoi(argv[1]);
        const int r = std::stoi(argv[2]);
        const Geometry geometry = build_geometry(q, r);
        const std::vector<u32> representatives = read_representatives(argv[3]);

        u64 shifts = 0;
        u64 violations = 0;
        bool have = false;
        u32 best_source = 0;
        u32 best_target = 0;
        unsigned best_i = 0;
        unsigned best_j = 0;
        Margin best_source_margin{};
        Margin best_target_margin{};
        Fnv ledger;

        for (u32 source : representatives) {
            const Margin source_margin = margin(geometry, source);
            require(source_margin.ticks >= 0, "atlas representative inactive");
            for (unsigned j = 0; j < 30; ++j) {
                if ((source & (u32{1} << j)) == 0) continue;
                for (unsigned i = 0; i < j; ++i) {
                    if (source & (u32{1} << i)) continue;
                    const u32 target = (source ^ (u32{1} << j)) | (u32{1} << i);
                    const Margin target_margin = margin(geometry, target);
                    ++shifts;
                    ledger.add(source);
                    ledger.add(target);
                    add_i128(ledger, source_margin.ticks);
                    add_i128(ledger, target_margin.ticks);
                    if (target_margin.ticks >= 0) continue;
                    ++violations;
                    const auto key = std::tuple{colex_any(source), colex_any(target), i, j};
                    const auto best_key = std::tuple{colex_any(best_source), colex_any(best_target), best_i, best_j};
                    if (!have || key < best_key) {
                        have = true;
                        best_source = source;
                        best_target = target;
                        best_i = i;
                        best_j = j;
                        best_source_margin = source_margin;
                        best_target_margin = target_margin;
                    }
                }
            }
        }
        require(have, "no compression violation among atlas representatives");
        std::array<u64, 10> source_atom_count{};
        std::array<u64, 10> target_atom_count{};
        std::array<i64, 10> source_atom_mass{};
        std::array<i64, 10> target_atom_mass{};
        u64 source_only_count = 0;
        u64 target_only_count = 0;
        i64 source_only_mass = 0;
        i64 target_only_mass = 0;
        Fnv atom_ledger;
        for (const auto& [atom, width] : geometry.classes) {
            const bool in_source = (atom & ~best_source) == 0;
            const bool in_target = (atom & ~best_target) == 0;
            if (in_source == in_target) continue;
            const unsigned arity = std::popcount(atom);
            atom_ledger.add(atom);
            atom_ledger.add(static_cast<u64>(width));
            atom_ledger.add(in_source ? 1 : 2);
            if (in_source) {
                ++source_only_count;
                source_only_mass += width;
                ++source_atom_count[arity];
                source_atom_mass[arity] += width;
            } else {
                ++target_only_count;
                target_only_mass += width;
                ++target_atom_count[arity];
                target_atom_mass[arity] += width;
            }
        }
        require(best_source_margin.mass - best_target_margin.mass ==
                    source_only_mass - target_only_mass,
                "failure-atom decomposition mismatch");
        const i128 denominator = static_cast<i128>(63) * geometry.grid;
        const i128 source_divisor = gcd128(best_source_margin.ticks, denominator);
        const i128 target_divisor = gcd128(best_target_margin.ticks, denominator);
        std::cout << "THM4296_COMPRESSED_REPRESENTATIVE_PROBE_V1\n"
                  << "PAIR " << q << ',' << r << " GRID " << geometry.grid
                  << " CLASSES_RANK_LE9 " << geometry.classes.size() << '\n'
                  << "REPRESENTATIVES " << representatives.size()
                  << " ELEMENTARY_LEFT_SHIFTS " << shifts
                  << " VIOLATIONS " << violations
                  << " LEDGER_FNV " << std::hex << ledger.state << std::dec << '\n'
                  << "LEAST_SOURCE_COLEX " << colex_any(best_source)
                  << " SOURCE " << hex8(best_source)
                  << " LABELS {" << labels(best_source) << "} RANK "
                  << std::popcount(best_source) << " MASS " << best_source_margin.mass
                  << " MARGIN_TICKS63 " << decimal(best_source_margin.ticks)
                  << " SURPLUS " << decimal(best_source_margin.ticks / source_divisor)
                  << '/' << decimal(denominator / source_divisor) << '\n'
                  << "SHIFT " << best_j << "->" << best_i
                  << " TARGET_COLEX " << colex_any(best_target)
                  << " TARGET " << hex8(best_target)
                  << " LABELS {" << labels(best_target) << "} MASS "
                  << best_target_margin.mass << " MARGIN_TICKS63 "
                  << decimal(best_target_margin.ticks)
                  << " SURPLUS " << decimal(best_target_margin.ticks / target_divisor)
                  << '/' << decimal(denominator / target_divisor) << '\n'
                  << "FAILURE_ATOM_SWAP SOURCE_ONLY_COUNT " << source_only_count
                  << " SOURCE_ONLY_MASS " << source_only_mass
                  << " TARGET_ONLY_COUNT " << target_only_count
                  << " TARGET_ONLY_MASS " << target_only_mass
                  << " NET_MASS " << (source_only_mass - target_only_mass)
                  << " LEDGER_FNV " << std::hex << atom_ledger.state << std::dec << '\n';
        for (unsigned arity = 0; arity <= 9; ++arity) {
            if (source_atom_count[arity] == 0 && target_atom_count[arity] == 0)
                continue;
            std::cout << "ATOM_ARITY " << arity
                      << " SOURCE_COUNT " << source_atom_count[arity]
                      << " SOURCE_MASS " << source_atom_mass[arity]
                      << " TARGET_COUNT " << target_atom_count[arity]
                      << " TARGET_MASS " << target_atom_mass[arity]
                      << " NET_MASS "
                      << (source_atom_mass[arity] - target_atom_mass[arity]) << '\n';
        }
        std::cout
                  << "VERDICT REFUTED_LEFT_COMPRESSED_CLOSURE_ON_ACTIVE_ATLAS_REPRESENTATIVES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "COMPRESSED_PROBE_ERROR " << error.what() << '\n';
        return 1;
    }
}
