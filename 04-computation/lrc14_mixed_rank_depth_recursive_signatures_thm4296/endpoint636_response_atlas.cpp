// Scratch-only exact response atlas for the two THM-4287 endpoint-636
// carrier failures.  Imports the frozen exact wall/activity implementation,
// but keeps the new consumer outside maintained sources.

#include <filesystem>

#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

#include <algorithm>
#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr std::size_t kRepairs = 5852925;
constexpr u64 kAExpectedFnv = UINT64_C(0xd611500ea833ff83);
constexpr u64 kBExpectedFnv = UINT64_C(0xee7792a8a2fd51c9);

struct Pattern {
    u64 a = 0;
    u64 b = 0;
    auto operator<=>(const Pattern&) const = default;
};

struct PatternClass {
    u64 count = 0;
    u32 least = 0;
};

struct Failures {
    std::vector<u32> a;
    std::vector<u32> b;
};

std::string hex8(u32 value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

std::string agent_labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (unsigned bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!first) out << ',';
        out << bit;
        first = false;
    }
    return out.str();
}

u64 body_fnv_agent(const std::vector<u32>& bodies) {
    FnvLocal ledger;
    for (u32 body : bodies) ledger.add(body);
    return ledger.state;
}

Failures read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure-ledger header changed");
    Failures failures;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int q = 0;
        int r = 0;
        std::string body_hex;
        row >> q >> r >> body_hex;
        require(row && r == 636, "malformed failure row");
        const u32 body = static_cast<u32>(std::stoul(body_hex, nullptr, 16));
        require(std::popcount(body) == 9, "failure body rank changed");
        if (q == 100) failures.a.push_back(body);
        else if (q == 256) failures.b.push_back(body);
        else require(false, "unexpected failure pair");
    }
    require(failures.a.size() == 64 &&
                body_fnv_agent(failures.a) == kAExpectedFnv,
            "(100,636) failures changed");
    require(failures.b.size() == 37 &&
                body_fnv_agent(failures.b) == kBExpectedFnv,
            "(256,636) failures changed");
    return failures;
}

std::set<u32> read_carrier(const std::filesystem::path& inherited,
                           const std::filesystem::path& additions,
                           const std::filesystem::path& witness) {
    const std::vector<u32> base =
        read_masks(inherited, kCarrierCount, kCarrierFnv);
    const std::vector<u32> add = read_additions(additions);
    const std::vector<u32> suffix =
        read_masks(witness, 9, UINT64_C(0x02b936529030e4bc));
    std::set<u32> out(base.begin(), base.end());
    for (u32 mask : add) require(out.insert(mask).second, "carrier overlap");
    require(out.insert(UINT32_C(0x014c9084)).second, "prior repair overlap");
    for (u32 mask : suffix)
        require(out.insert(mask).second, "witness overlap");
    require(out.size() == 9006, "carrier size changed");
    return out;
}

u64 popcount(Pattern value) {
    return std::popcount(value.a) + std::popcount(value.b);
}

Pattern select_pattern(Pattern raw, unsigned mode,
                       bool active_a, bool active_b) {
    if (mode == 0) return {active_a ? raw.a : 0, 0};
    if (mode == 1) return {0, active_b ? raw.b : 0};
    if (mode == 2)
        return {active_a && active_b ? raw.a : 0,
                active_a && active_b ? raw.b : 0};
    require(mode == 3, "unknown atlas mode");
    return {active_a ? raw.a : 0, active_b ? raw.b : 0};
}

const char* mode_name(unsigned mode) {
    constexpr std::array<const char*, 4> names = {
        "LOCAL_100", "LOCAL_256", "COMMON_ACTIVE", "CARRIER_UNION"};
    return names.at(mode);
}

void write_atlas(const std::filesystem::path& path, unsigned mode,
                 const std::vector<Pattern>& raw,
                 const ActiveUniverse& active_a,
                 const ActiveUniverse& active_b,
                 const std::set<u32>& carrier) {
    std::map<Pattern, PatternClass> classes;
    FnvLocal response_ledger;
    u64 candidates = 0;
    u64 incidences = 0;
    u64 carrier_overlap = 0;
    u64 maximum_cover = 0;
    u32 least_maximum = 0;
    for (u64 rank = 0; rank < kRepairs; ++rank) {
        const Pattern response = select_pattern(
            raw[rank], mode, active_a.active[rank], active_b.active[rank]);
        if (response.a == 0 && response.b == 0) continue;
        const u32 mask = unrank_colex8(rank);
        ++candidates;
        incidences += popcount(response);
        response_ledger.add(mask);
        response_ledger.add(response.a);
        response_ledger.add(response.b);
        PatternClass& c = classes[response];
        ++c.count;
        if (c.count == 1 || mask < c.least) c.least = mask;
        if (carrier.contains(mask)) ++carrier_overlap;
        const u64 covered = popcount(response);
        if (covered > maximum_cover ||
            (covered == maximum_cover && mask < least_maximum)) {
            maximum_cover = covered;
            least_maximum = mask;
        }
    }
    FnvLocal class_ledger;
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create atlas");
    output << "a_hex\tb_hex\tcover\tcount\tleast_mask\tlabels\n";
    for (const auto& [pattern, c] : classes) {
        class_ledger.add(pattern.a);
        class_ledger.add(pattern.b);
        class_ledger.add(c.count);
        class_ledger.add(c.least);
        output << std::hex << std::setw(16) << std::setfill('0') << pattern.a
               << '\t' << std::setw(16) << pattern.b << std::dec
               << std::setfill(' ') << '\t' << popcount(pattern) << '\t'
               << c.count << '\t' << hex8(c.least) << '\t'
               << agent_labels(c.least) << '\n';
    }
    require(output.good(), "atlas write failed");
    std::cout << "ATLAS " << mode_name(mode)
              << " CANDIDATES " << candidates
              << " INCIDENCES " << incidences
              << " CLASSES " << classes.size()
              << " RESPONSE_FNV " << std::hex << response_ledger.state
              << " CLASS_FNV " << class_ledger.state << std::dec
              << " MAX_COVER " << maximum_cover
              << " LEAST_MAX " << hex8(least_maximum)
              << " LEAST_MAX_LABELS {" << agent_labels(least_maximum) << "}"
              << " CARRIER_OVERLAP " << carrier_overlap << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 7,
                "usage: atlas FAILURES BASE8951 ADD45 WITNESS9 OUTDIR THREADS");
        init_choose8_local();
        const Failures failures = read_failures(argv[1]);
        const std::set<u32> carrier = read_carrier(argv[2], argv[3], argv[4]);
        const std::filesystem::path outdir = argv[5];
        std::filesystem::create_directories(outdir);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        std::cout << "ENDPOINT636_RESPONSE_ATLAS_V1\n"
                  << "FAILURES A " << failures.a.size() << " FNV "
                  << std::hex << body_fnv_agent(failures.a) << " B " << std::dec
                  << failures.b.size() << " FNV " << std::hex
                  << body_fnv_agent(failures.b) << std::dec << '\n';

        const ActiveUniverse active_a = build_active_universe(cells, 100, 636);
        const ActiveUniverse active_b = build_active_universe(cells, 256, 636);
        FnvLocal both_ledger;
        u64 both = 0;
        u64 either = 0;
        for (u64 rank = 0; rank < kRepairs; ++rank) {
            const bool aa = active_a.active[rank];
            const bool ab = active_b.active[rank];
            if (aa && ab) {
                ++both;
                both_ledger.add(unrank_colex8(rank));
            }
            if (aa || ab) ++either;
        }
        std::cout << "ACTIVITY A " << active_a.count << " FNV " << std::hex
                  << active_a.fnv << " B " << std::dec << active_b.count
                  << " FNV " << std::hex << active_b.fnv << std::dec
                  << " BOTH " << both << " BOTH_FNV " << std::hex
                  << both_ledger.state << std::dec << " EITHER " << either
                  << '\n';

        std::vector<Pattern> raw(kRepairs);
        u64 incidences_a = 0;
        for (std::size_t i = 0; i < failures.a.size(); ++i) {
            enumerate_disjoint_repairs(failures.a[i], [&](u32, u64 rank) {
                raw[rank].a |= UINT64_C(1) << i;
                ++incidences_a;
            });
        }
        u64 incidences_b = 0;
        for (std::size_t i = 0; i < failures.b.size(); ++i) {
            enumerate_disjoint_repairs(failures.b[i], [&](u32, u64 rank) {
                raw[rank].b |= UINT64_C(1) << i;
                ++incidences_b;
            });
        }
        require(incidences_a == 64 * UINT64_C(203490) &&
                    incidences_b == 37 * UINT64_C(203490),
                "raw disjoint-incidence count changed");
        std::cout << "RAW_INCIDENCES A " << incidences_a
                  << " B " << incidences_b << '\n';

        for (unsigned mode = 0; mode < 4; ++mode) {
            const std::filesystem::path path =
                outdir / (std::string(mode_name(mode)) + ".tsv");
            write_atlas(path, mode, raw, active_a, active_b, carrier);
        }
        std::cout << "VERDICT PASS EXACT_ACTIVITY_AND_RESPONSE_ATLASES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT636_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}
