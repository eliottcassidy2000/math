// Independent integer-only replay of the endpoint-636 13-cover obstruction.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Row {
    u128 pattern = 0;
    u32 mask = 0;
    unsigned load = 0;
};

constexpr unsigned kScale = 528;
constexpr std::array<unsigned, 101> kWeight = [] {
    std::array<unsigned, 101> out{};
    out[2]=264; out[5]=454; out[7]=65; out[8]=417; out[12]=227;
    out[19]=74; out[22]=236; out[26]=264; out[30]=37; out[33]=181;
    out[36]=37; out[37]=528; out[44]=301; out[45]=74; out[51]=227;
    out[53]=37; out[56]=454; out[61]=46; out[64]=264; out[65]=12;
    out[66]=9; out[67]=213; out[71]=210; out[73]=84; out[74]=63;
    out[75]=42; out[76]=42; out[77]=117; out[79]=74; out[80]=60;
    out[82]=232; out[83]=139; out[84]=264; out[85]=105; out[86]=97;
    out[87]=75; out[89]=84; out[90]=21; out[91]=180; out[92]=273;
    out[94]=63; out[95]=96; out[98]=23; out[100]=65;
    return out;
}();

constexpr std::array<u32, 14> kWitness{{
    UINT32_C(0x18468880), UINT32_C(0x080e8281),
    UINT32_C(0x22081017), UINT32_C(0x08422a82),
    UINT32_C(0x004cac40), UINT32_C(0x19c04044),
    UINT32_C(0x00c08ec0), UINT32_C(0x10443016),
    UINT32_C(0x01609124), UINT32_C(0x10413209),
    UINT32_C(0x01611640), UINT32_C(0x00606449),
    UINT32_C(0x0128d084), UINT32_C(0x08806449)}};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}
void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

unsigned popcount(u128 value) {
    return std::popcount(static_cast<u64>(value)) +
           std::popcount(static_cast<u64>(value >> 64));
}

unsigned load(u128 pattern) {
    unsigned answer = 0;
    for (unsigned bit = 0; bit < kWeight.size(); ++bit)
        if (((pattern >> bit) & 1) != 0) answer += kWeight[bit];
    return answer;
}

u64 parse_hex64(const std::string& text) {
    std::size_t used = 0;
    const u64 value = std::stoull(text, &used, 16);
    require(used == text.size(), "malformed hex field");
    return value;
}

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::stringstream stream(text);
    std::vector<std::string> fields;
    std::string field;
    while (std::getline(stream, field, delimiter)) fields.push_back(field);
    return fields;
}

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "a_hex\tb_hex\tcover\tcount\tleast_mask\tlabels",
            "atlas header changed");
    std::vector<Row> rows;
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = split(line, '\t');
        require(fields.size() == 6, "malformed atlas row");
        const u64 a = parse_hex64(fields[0]);
        const u64 b = parse_hex64(fields[1]);
        const u32 mask = static_cast<u32>(parse_hex64(fields[4]));
        const u128 pattern = static_cast<u128>(a) |
                             (static_cast<u128>(b) << 64);
        require(popcount(pattern) == std::stoul(fields[2]),
                "declared cover changed");
        rows.push_back({pattern, mask, load(pattern)});
    }
    require(rows.size() == 1835, "atlas class count changed");
    return rows;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: dual-gap-independent CARRIER_UNION.tsv");
        const std::vector<Row> rows = read_rows(argv[1]);
        require(std::accumulate(kWeight.begin(), kWeight.end(), 0u) == 6830,
                "dual total changed");
        u128 full = 0;
        unsigned maximum_load = 0;
        std::map<u32, u128> by_mask;
        for (const Row& row : rows) {
            full |= row.pattern;
            maximum_load = std::max(maximum_load, row.load);
            require(by_mask.emplace(row.mask, row.pattern).second,
                    "least-mask labels stopped being unique");
        }
        require(popcount(full) == 101 && maximum_load == kScale,
                "universe or dual feasibility changed");

        std::vector<Row> survivors;
        for (const Row& row : rows)
            if (kScale - row.load <= 34) survivors.push_back(row);
        require(survivors.size() == 103,
                "dual-gap survivor count changed");
        std::vector<Row> maximal;
        for (const Row& row : survivors) {
            bool dominated = false;
            for (const Row& other : survivors)
                if (row.pattern != other.pattern &&
                    (row.pattern | other.pattern) == other.pattern) {
                    dominated = true;
                    break;
                }
            if (!dominated) maximal.push_back(row);
        }
        require(maximal.size() == 58,
                "maximal survivor count changed");

        std::array<std::vector<std::size_t>, 101> through;
        for (std::size_t index = 0; index < maximal.size(); ++index)
            for (unsigned bit = 0; bit < 101; ++bit)
                if (((maximal[index].pattern >> bit) & 1) != 0)
                    through[bit].push_back(index);

        using State = std::tuple<u64, u64, int, int>;
        std::set<State> memo;
        std::vector<std::size_t> path;
        std::vector<std::size_t> answer;
        u64 nodes = 0;
        std::function<bool(u128, int, int)> search =
            [&](u128 covered, int budget, int charged) -> bool {
            ++nodes;
            if (covered == full) {
                if (budget == 0) {
                    answer = path;
                    return true;
                }
                return false;
            }
            if (budget == 0 || charged > 34) return false;
            const State state{static_cast<u64>(covered),
                              static_cast<u64>(covered >> 64),
                              budget, charged};
            if (memo.contains(state)) return false;
            const u128 uncovered = full & ~covered;
            unsigned maximum_gain = 0;
            for (const Row& row : maximal)
                maximum_gain = std::max(maximum_gain,
                                        popcount(row.pattern & uncovered));
            if ((popcount(uncovered) + maximum_gain - 1) / maximum_gain >
                static_cast<unsigned>(budget)) {
                memo.insert(state);
                return false;
            }

            unsigned pivot = 0;
            std::size_t fewest = std::numeric_limits<std::size_t>::max();
            for (unsigned bit = 0; bit < 101; ++bit) {
                if (((uncovered >> bit) & 1) == 0) continue;
                std::size_t count = 0;
                for (std::size_t index : through[bit]) {
                    const int extra =
                        static_cast<int>(kScale - maximal[index].load +
                                         load(maximal[index].pattern & covered));
                    count += charged + extra <= 34;
                }
                if (count < fewest) {
                    fewest = count;
                    pivot = bit;
                }
            }
            struct Choice {
                unsigned gain;
                int extra;
                u32 mask;
                std::size_t index;
            };
            std::vector<Choice> choices;
            for (std::size_t index : through[pivot]) {
                const int extra =
                    static_cast<int>(kScale - maximal[index].load +
                                     load(maximal[index].pattern & covered));
                if (charged + extra > 34) continue;
                choices.push_back({popcount(maximal[index].pattern & uncovered),
                                   extra, maximal[index].mask, index});
            }
            std::sort(choices.begin(), choices.end(),
                      [](const Choice& left, const Choice& right) {
                if (left.gain != right.gain) return left.gain > right.gain;
                if (left.extra != right.extra) return left.extra < right.extra;
                return left.mask < right.mask;
            });
            for (const Choice& choice : choices) {
                path.push_back(choice.index);
                if (search(covered | maximal[choice.index].pattern,
                           budget - 1, charged + choice.extra))
                    return true;
                path.pop_back();
            }
            memo.insert(state);
            return false;
        };
        require(!search(0, 13, 0), "unexpected 13-cover");

        u128 witness_union = 0;
        for (u32 mask : kWitness) witness_union |= by_mask.at(mask);
        require(witness_union == full, "14-mask witness stopped covering");

        std::cout << "ENDPOINT636_DUAL_GAP_INDEPENDENT_V1\n"
                  << "CLASSES " << rows.size() << " UNIVERSE "
                  << popcount(full) << " DUAL_SCALE " << kScale
                  << " DUAL_TOTAL 6830 MAX_LOAD " << maximum_load << '\n'
                  << "THIRTEEN_CHARGE_CAP 34 SURVIVORS "
                  << survivors.size() << " MAXIMAL " << maximal.size()
                  << " SEARCH_NODES " << nodes << " MEMO " << memo.size()
                  << " FEASIBLE 0\nWITNESS14";
        for (u32 mask : kWitness)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << mask;
        std::cout << std::dec << std::setfill(' ') << " COVER "
                  << popcount(witness_union)
                  << "\nVERDICT PASS EXACT_INTEGER_DUAL_GAP_AND_WITNESS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DUAL_GAP_INDEPENDENT_ERROR " << error.what() << '\n';
        return 1;
    }
}
