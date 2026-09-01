// Independent literal-wall audit of the five-mask arbitrary-rank cover for
// the 145 THM-4303 endpoint-595 failures.  This freezes an upper witness only;
// the separate CP-SAT exhaustive search supplies the no-four lower bound.

#define main lrc595_five_mask_hidden_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

#include <map>

namespace {

constexpr std::array<int, 3> kQ595 = {96, 100, 210};
constexpr std::array<u32, 5> kFiveMasks = {
    UINT32_C(0x00612a76), UINT32_C(0x00a183f6),
    UINT32_C(0x024d8b32), UINT32_C(0x0280a1ae),
    UINT32_C(0x10110bf6)};
constexpr u64 kTaggedFailureFnv = UINT64_C(0xf3d7f95fc38e7b49);
constexpr u64 kFiveMaskFnv = UINT64_C(0xebea2511eb7fa46f);
constexpr u64 kResponseFnv = UINT64_C(0x8e13cc00c5bca917);

Geometry build_full_geometry(int pair_q, int pair_r) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair_q);
    grid = checked_lcm(grid, 14LL * pair_r);
    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 unit = grid / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(pair_q);
    add_walls(pair_r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> by_failure;
    Geometry geometry;
    geometry.grid = grid;
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(pair_q, grid, left, right) ||
            !safe_midpoint(pair_r, grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], grid, left, right))
                failure |= u32{1} << bit;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        by_failure[failure] += width;
    }
    for (const auto& entry : by_failure) geometry.classes.push_back(entry);
    return geometry;
}

struct TaggedBody {
    int q = 0;
    int r = 0;
    u32 body = 0;
};

std::vector<TaggedBody> read_tagged(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    std::vector<TaggedBody> rows;
    std::array<u64, 3> counts{};
    Fnv ledger;
    while (std::getline(input, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        TaggedBody tagged;
        std::string token;
        row >> tagged.q >> tagged.r >> token;
        require(row && tagged.r == 595, "malformed tagged body");
        const auto where = std::find(kQ595.begin(), kQ595.end(), tagged.q);
        require(where != kQ595.end(), "unexpected pair");
        tagged.body = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(std::popcount(tagged.body) == 9, "body rank changed");
        ++counts[where - kQ595.begin()];
        ledger.add(tagged.q);
        ledger.add(tagged.r);
        ledger.add(tagged.body);
        rows.push_back(tagged);
    }
    require(rows.size() == 145 && counts == std::array<u64, 3>{116, 13, 16} &&
                ledger.state == kTaggedFailureFnv,
            "tagged failure identity changed");
    return rows;
}

std::string hex16(u64 value) {
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: five_mask_literal_audit FAILURE_CSV");
        const std::vector<TaggedBody> failures = read_tagged(argv[1]);
        const std::array<Geometry, 3> geometry = {
            build_full_geometry(96, 595), build_full_geometry(100, 595),
            build_full_geometry(210, 595)};
        require(geometry[0].grid == INT64_C(36482318832960) &&
                    geometry[0].cells == 8515 &&
                    geometry[0].classes.size() == 2453 &&
                    geometry[0].pair_ticks == INT64_C(26803336285440) &&
                    geometry[1].grid == INT64_C(91205797082400) &&
                    geometry[1].cells == 8523 &&
                    geometry[1].classes.size() == 2519 &&
                    geometry[1].pair_ticks == INT64_C(67008340713600) &&
                    geometry[2].grid == INT64_C(18241159416480) &&
                    geometry[2].cells == 8743 &&
                    geometry[2].classes.size() == 2543 &&
                    geometry[2].pair_ticks == INT64_C(13412617218000),
                "full geometry identity changed");

        Fnv geometry_ledger;
        for (std::size_t p = 0; p < geometry.size(); ++p) {
            geometry_ledger.add(kQ595[p]);
            geometry_ledger.add(595);
            geometry_ledger.add(geometry[p].grid);
            geometry_ledger.add(geometry[p].cells);
            geometry_ledger.add(geometry[p].pair_ticks);
            geometry_ledger.add(geometry[p].classes.size());
            for (const auto& [failure, width] : geometry[p].classes) {
                geometry_ledger.add(failure);
                geometry_ledger.add(width);
            }
        }

        Fnv witness_ledger;
        for (u32 mask : kFiveMasks) witness_ledger.add(mask);
        require(witness_ledger.state == kFiveMaskFnv,
                "five-mask identity changed");

        std::array<std::array<u64, 3>, 5> response_bits{};
        std::array<std::array<u64, 3>, 5> pair_hits{};
        std::array<u64, 3> joined{};
        u64 incidences = 0;
        Fnv response_ledger;
        for (std::size_t ordinal = 0; ordinal < failures.size(); ++ordinal) {
            const TaggedBody tagged = failures[ordinal];
            const std::size_t p =
                std::find(kQ595.begin(), kQ595.end(), tagged.q) - kQ595.begin();
            bool covered = false;
            for (std::size_t m = 0; m < kFiveMasks.size(); ++m) {
                const Margin value = margin(geometry[p], kFiveMasks[m]);
                const bool hit = value.ticks >= 0 &&
                                 (kFiveMasks[m] & tagged.body) == 0;
                if (!hit) continue;
                response_bits[m][ordinal / 64] |=
                    UINT64_C(1) << (ordinal % 64);
                ++pair_hits[m][p];
                ++incidences;
                covered = true;
                response_ledger.add(kFiveMasks[m]);
                response_ledger.add(tagged.q);
                response_ledger.add(tagged.r);
                response_ledger.add(ordinal);
                response_ledger.add(tagged.body);
            }
            joined[p] += covered;
        }
        require(joined == std::array<u64, 3>{116, 13, 16} &&
                    incidences == 244 &&
                    response_ledger.state == kResponseFnv,
                "five-mask response cover changed");

        std::cout << "LRC595_FIVE_MASK_LITERAL_AUDIT_V1\n"
                  << "FAILURES 145 COUNTS 116,13,16 FNV " << std::hex
                  << kTaggedFailureFnv << std::dec << '\n'
                  << "GEOMETRY CLASSES 2453,2519,2543 UNION_SUPPORT 2783 FNV "
                  << std::hex << geometry_ledger.state << std::dec << '\n';
        for (std::size_t m = 0; m < kFiveMasks.size(); ++m) {
            std::cout << "MASK " << hex8(kFiveMasks[m]) << " RANK "
                      << std::popcount(kFiveMasks[m]) << " LABELS {"
                      << labels(kFiveMasks[m]) << "} SIG ";
            for (std::size_t p = 0; p < geometry.size(); ++p)
                std::cout << (margin(geometry[p], kFiveMasks[m]).ticks >= 0);
            std::cout << " COVER " << pair_hits[m][0] << ','
                      << pair_hits[m][1] << ',' << pair_hits[m][2]
                      << " BITS " << hex16(response_bits[m][2]) << ':'
                      << hex16(response_bits[m][1]) << ':'
                      << hex16(response_bits[m][0]);
            for (std::size_t p = 0; p < geometry.size(); ++p) {
                const Margin value = margin(geometry[p], kFiveMasks[m]);
                std::cout << " M" << kQ595[p] << ' ' << value.mass
                          << " T" << kQ595[p] << ' '
                          << decimal(value.ticks) << '/'
                          << decimal(static_cast<i128>(63) * geometry[p].grid);
            }
            std::cout << '\n';
        }
        std::cout << "UNION JOINED 116,13,16 TOTAL 145 UNCOVERED 0 "
                     "UNCOVERED_FNV cbf29ce484222325 INCIDENCES 244\n"
                  << "WITNESS_FNV " << std::hex << witness_ledger.state
                  << " RESPONSE_FNV " << response_ledger.state << std::dec
                  << "\nSCOPE FIXED_POOL_FULL_LITERAL_GEOMETRY_"
                     "ARBITRARY_RANK_RESPONSE_COVER_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT PASS FIVE_MASK_UPPER_WITNESS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FIVE_MASK_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
