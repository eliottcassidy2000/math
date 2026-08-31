// Exact rank-eight response atlases at the first post-exchange failure layer.

#include <filesystem>

#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

#include <algorithm>
#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kRepairs = 5852925;
constexpr std::array<u32, 14> kDelete = {
    0x00003e1a, 0x000132a3, 0x00017464, 0x00033388,
    0x000a16c2, 0x000f8118, 0x00142a1a, 0x00154348,
    0x00184ba0, 0x001aa260, 0x00202c2b, 0x002066a4,
    0x002b018a, 0x0030c2a2};
constexpr std::array<u32, 14> kAdd = {
    0x18468880, 0x080e8281, 0x22081017, 0x08422a82,
    0x004cac40, 0x19c04044, 0x00c08ec0, 0x10443016,
    0x01609124, 0x10413209, 0x01611640, 0x00606449,
    0x0128d084, 0x08806449};

struct Pattern72 {
    u64 low = 0;
    u64 high = 0;
    auto operator<=>(const Pattern72&) const = default;
};

struct Class72 {
    u64 count = 0;
    u32 least = 0;
};

struct Input72 {
    Pair first;
    Pair second;
    std::vector<u32> bodies;
    std::vector<unsigned char> row;
    Pattern72 first_mask;
    Pattern72 second_mask;
};

void set_bit(Pattern72& pattern, std::size_t index) {
    if (index < 64) pattern.low |= UINT64_C(1) << index;
    else pattern.high |= UINT64_C(1) << (index - 64);
}

Pattern72 bit_and(Pattern72 left, Pattern72 right) {
    return {left.low & right.low, left.high & right.high};
}

Pattern72 bit_or(Pattern72 left, Pattern72 right) {
    return {left.low | right.low, left.high | right.high};
}

u64 pattern_count(Pattern72 pattern) {
    return std::popcount(pattern.low) + std::popcount(pattern.high);
}

std::string hex8_agent(u32 value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

Input72 read_input(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure header changed");
    Input72 out;
    std::vector<Pair> order;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        Pair pair;
        std::string token;
        row >> pair.q >> pair.r >> token;
        require(row && pair.r == 632, "malformed r632 row");
        if (order.empty() ||
            !(order.back().q == pair.q && order.back().r == pair.r))
            order.push_back(pair);
        require(order.size() <= 2, "more than two failure rows");
        const std::size_t index = out.bodies.size();
        out.bodies.push_back(
            static_cast<u32>(std::stoul(token, nullptr, 16)));
        out.row.push_back(static_cast<unsigned char>(order.size() - 1));
        set_bit(order.size() == 1 ? out.first_mask : out.second_mask, index);
    }
    require(order.size() == 2 && out.bodies.size() == 72,
            "r632 failure shape changed");
    out.first = order[0];
    out.second = order[1];
    require(out.first.q == 100 && out.second.q == 256,
            "r632 pair identities changed");
    FnvLocal ledger;
    for (std::size_t i = 0; i < out.bodies.size(); ++i) {
        require(std::popcount(out.bodies[i]) == 9, "body rank changed");
        ledger.add(out.row[i]);
        ledger.add(out.bodies[i]);
    }
    std::cout << "FAILURES FIRST 6 SECOND 66 TOTAL 72 FNV "
              << std::hex << ledger.state << std::dec << '\n';
    return out;
}

std::set<u32> build_replaced(const std::filesystem::path& inherited,
                             const std::filesystem::path& additions,
                             const std::filesystem::path& witness) {
    std::vector<u32> carrier = read_masks(inherited, kCarrierCount, kCarrierFnv);
    const std::vector<u32> add45 = read_additions(additions);
    const std::vector<u32> suffix =
        read_masks(witness, 9, UINT64_C(0x02b936529030e4bc));
    std::set<u32> distinct(carrier.begin(), carrier.end());
    for (u32 mask : add45) distinct.insert(mask);
    distinct.insert(UINT32_C(0x014c9084));
    for (u32 mask : suffix) distinct.insert(mask);
    for (u32 mask : kDelete)
        require(distinct.erase(mask) == 1, "delete absent");
    for (u32 mask : kAdd)
        require(distinct.insert(mask).second, "add overlap");
    require(distinct.size() == 9006, "replaced carrier size changed");
    return distinct;
}

Pattern72 selected(Pattern72 raw, unsigned mode, bool active_first,
                   bool active_second, const Input72& input) {
    const Pattern72 first = bit_and(raw, input.first_mask);
    const Pattern72 second = bit_and(raw, input.second_mask);
    if (mode == 0) return active_first ? first : Pattern72{};
    if (mode == 1) return active_second ? second : Pattern72{};
    if (mode == 2)
        return active_first && active_second ? raw : Pattern72{};
    require(mode == 3, "mode changed");
    return bit_or(active_first ? first : Pattern72{},
                  active_second ? second : Pattern72{});
}

const char* mode_name72(unsigned mode) {
    constexpr std::array<const char*, 4> names = {
        "R632_LOCAL_100", "R632_LOCAL_256",
        "R632_COMMON_ACTIVE", "R632_CARRIER_UNION"};
    return names.at(mode);
}

void write_mode(const std::filesystem::path& path, unsigned mode,
                const std::vector<Pattern72>& raw, const Input72& input,
                const ActiveUniverse& first, const ActiveUniverse& second,
                const std::set<u32>& carrier) {
    std::map<Pattern72, Class72> classes;
    FnvLocal response_ledger;
    u64 candidates = 0;
    u64 incidences = 0;
    u64 overlap = 0;
    u64 max_cover = 0;
    u32 least_max = 0;
    for (u64 rank = 0; rank < kRepairs; ++rank) {
        const Pattern72 response = selected(raw[rank], mode,
            first.active[rank], second.active[rank], input);
        if (response.low == 0 && response.high == 0) continue;
        const u32 mask = unrank_colex8(rank);
        ++candidates;
        incidences += pattern_count(response);
        response_ledger.add(mask);
        response_ledger.add(response.low);
        response_ledger.add(response.high);
        Class72& c = classes[response];
        ++c.count;
        if (c.count == 1 || mask < c.least) c.least = mask;
        if (carrier.contains(mask)) ++overlap;
        if (pattern_count(response) > max_cover ||
            (pattern_count(response) == max_cover && mask < least_max)) {
            max_cover = pattern_count(response);
            least_max = mask;
        }
    }
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create atlas");
    output << "a_hex\tb_hex\tcover\tcount\tleast_mask\n";
    FnvLocal class_ledger;
    for (const auto& [pattern, c] : classes) {
        class_ledger.add(pattern.low);
        class_ledger.add(pattern.high);
        class_ledger.add(c.count);
        class_ledger.add(c.least);
        output << std::hex << std::setw(16) << std::setfill('0') << pattern.low
               << '\t' << std::setw(16) << pattern.high << std::dec
               << std::setfill(' ') << '\t' << pattern_count(pattern) << '\t'
               << c.count << '\t' << hex8_agent(c.least) << '\n';
    }
    std::cout << "ATLAS " << mode_name72(mode) << " CANDIDATES "
              << candidates << " INCIDENCES " << incidences << " CLASSES "
              << classes.size() << " RESPONSE_FNV " << std::hex
              << response_ledger.state << " CLASS_FNV " << class_ledger.state
              << std::dec << " MAX_COVER " << max_cover << " LEAST_MAX "
              << hex8_agent(least_max) << " CARRIER_OVERLAP " << overlap
              << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 6,
                "usage: atlas FAILURES BASE8951 ADD45 WITNESS9 OUTDIR");
        init_choose8_local();
        const Input72 input = read_input(argv[1]);
        const std::set<u32> carrier = build_replaced(argv[2], argv[3], argv[4]);
        const std::vector<Cell> cells = build_pool_cells();
        const ActiveUniverse first =
            build_active_universe(cells, input.first.q, input.first.r);
        const ActiveUniverse second =
            build_active_universe(cells, input.second.q, input.second.r);
        u64 both = 0;
        FnvLocal both_ledger;
        for (u64 rank = 0; rank < kRepairs; ++rank)
            if (first.active[rank] && second.active[rank]) {
                ++both;
                both_ledger.add(unrank_colex8(rank));
            }
        std::cout << "R632_RESPONSE_ATLAS_V1\n"
                  << "ACTIVITY FIRST " << first.count << " FNV " << std::hex
                  << first.fnv << " SECOND " << std::dec << second.count
                  << " FNV " << std::hex << second.fnv << std::dec
                  << " BOTH " << both << " BOTH_FNV " << std::hex
                  << both_ledger.state << std::dec << '\n';
        std::vector<Pattern72> raw(kRepairs);
        for (std::size_t index = 0; index < input.bodies.size(); ++index)
            enumerate_disjoint_repairs(input.bodies[index],
                [&](u32, u64 rank) { set_bit(raw[rank], index); });
        const std::filesystem::path outdir = argv[5];
        std::filesystem::create_directories(outdir);
        for (unsigned mode = 0; mode < 4; ++mode)
            write_mode(outdir / (std::string(mode_name72(mode)) + ".tsv"),
                       mode, raw, input, first, second, carrier);
        std::cout << "VERDICT PASS EXACT_R632_RESPONSE_ATLASES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R632_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}
